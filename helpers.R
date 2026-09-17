##### global.R for SpliceR

##### Set options for environment
options(shiny.sanitize.errors = FALSE)

##### Load packages
library(shiny)
library(Biostrings)
library(magrittr)
library(stringi)
library(dplyr)
library(tidyr)
library(ggplot2)
library(grr)
library(printr)
library(plyr)
library(readr)
library(rmarkdown)
library(DT)
library(httr)
library(jsonlite)

# Global GraphQL & Refget endpoints
ENSEMBL_GRAPHQL_URL <- "https://beta.ensembl.org/data/graphql"
ENSEMBL_REFGET_URL  <- "https://beta.ensembl.org/data/refget/sequence/"

# Query Ensembl GraphQL to extract exon genomic coordinates and strand for a given transcript ID
getTranscriptExonCoordinates = function(id, upstream = 0, downstream = 0, species = "Homo_sapiens"){
  # GraphQL query to retrieve transcript coordinates, strand, and associated exon slices
  gql_query <- '
  query GetTranscriptExons($transcriptId: String!) {
    transcript(by_id: { stable_id: $transcriptId }) {
      stable_id
      gene {
        stable_id
      }
      slice {
        region {
          name
        }
        strand {
          code
          value
        }
      }
      exons {
        slice {
          location {
            start
            end
          }
        }
      }
    }
  }'
  
  response <- POST(
    url = ENSEMBL_GRAPHQL_URL,
    body = list(query = gql_query, variables = list(transcriptId = id)),
    encode = "json",
    content_type_json()
  )
  
  stop_for_status(response)
  res_content <- content(response, as = "parsed", simplifyVector = FALSE)
  
  tx_data <- res_content$data$transcript
  if (is.null(tx_data)) {
    stop(paste("Transcript ID not found via Ensembl GraphQL:", id))
  }
  
  chr_name  <- tx_data$slice$region$name
  strand_val <- as.numeric(tx_data$slice$strand$value)
  gene_id   <- tx_data$gene$stable_id
  
  exons_list <- tx_data$exons
  
  df_exons <- do.call(rbind, lapply(exons_list, function(ex) {
    data.frame(
      exon_chrom_start = ex$slice$location$start,
      exon_chrom_end   = ex$slice$location$end,
      strand           = strand_val,
      chromosome_name  = chr_name,
      gene_id          = gene_id,
      transcript_id    = id,
      stringsAsFactors = FALSE
    )
  }))
  
  # Apply upstream/downstream adjustments based on strand definition
  if (upstream != 0 || downstream != 0) {
    if (strand_val == 1) {
      df_exons$exon_chrom_start <- df_exons$exon_chrom_start - upstream
      df_exons$exon_chrom_end   <- df_exons$exon_chrom_end + downstream
    } else {
      df_exons$exon_chrom_start <- df_exons$exon_chrom_start - downstream
      df_exons$exon_chrom_end   <- df_exons$exon_chrom_end + upstream
    }
  }
  
  return(as_tibble(df_exons))
}

# Fetch genomic DNA sequence via GraphQL / Refget and return a Biostrings DNAString
coordinatesToDNAString = function(start, end, strand, chromosome, species = "Homo_sapiens", upstream = 0, downstream = 0){
  seq_str <- coordinatesToDNAChar(start, end, strand, chromosome, species, upstream, downstream)
  return(Biostrings::DNAString(seq_str))
}

# Fetch genomic DNA sequence via GraphQL / Refget and return a raw character string
coordinatesToDNAChar = function(start, end, strand, chromosome, species = "Homo_sapiens", upstream = 0, downstream = 0){
  # Sort coordinates
  c_start <- min(start, end) - upstream
  c_end   <- max(start, end) + downstream
  
  # GraphQL query to resolve chromosome sequence checksum ID for refget API
  gql_query <- '
  query GetRegionSequenceChecksum($regionName: String!) {
    region(by_name: { name: $regionName }) {
      sequence {
        checksum
      }
    }
  }'
  
  response <- POST(
    url = ENSEMBL_GRAPHQL_URL,
    body = list(query = gql_query, variables = list(regionName = as.character(chromosome))),
    encode = "json",
    content_type_json()
  )
  
  stop_for_status(response)
  res_content <- content(response, as = "parsed", simplifyVector = FALSE)
  checksum <- res_content$data$region$sequence$checksum
  
  # refget protocol uses 0-based start coordinates
  refget_start <- c_start - 1
  refget_end   <- c_end
  
  refget_url <- paste0(ENSEMBL_REFGET_URL, checksum, "?start=", refget_start, "&end=", refget_end)
  
  seq_resp <- GET(refget_url, add_headers(Accept = "text/plain"))
  stop_for_status(seq_resp)
  
  dna_seq <- content(seq_resp, as = "text", encoding = "UTF-8")
  
  # Reverse-complement if on the negative strand
  if (strand == -1) {
    dna_seq <- revcom(dna_seq)
  }
  
  return(dna_seq)
}

# Vectorized 'matchPatterns' function
matchPatterns = Vectorize(matchPattern, vectorize.args = "subject")

# Function to pull out start site of the guides with respect to the target sequence
extractGuideStart = Vectorize(
  FUN = function(alignments, exon){alignments[[exon]]@ranges@start},
  vectorize.args = "exon"
)

# Function to pull out the guide as a character
extractGuide = Vectorize(
  FUN = function(alignments, exon){alignments[[exon]] %>% as.character},
  vectorize.args = "exon"
)

revcom = function(x){
  x %>%
    gsub("A", 't', .) %>%
    gsub("C", 'g', .) %>%
    gsub("G", 'c', .) %>%
    gsub("T", 'a', .) %>%
    toupper(.) %>%
    stringi::stri_reverse()
}

probability = function(logit){mapply(FUN = function(l){exp(l)/(1 + exp(l))}, l = logit)}

logit = function(probability){mapply(FUN = function(p){log(p/(1-p))}, p = probability)}

# Add genetic coordinates to the protospacer
addProtospacerCoordinates = function(data, guide_length){
  data %>% 
    mutate(
      cbe_position_tmp = {
        ifelse(strand == 1,
               {
                 ifelse((enzyme == "CBE and ABE") & (splice_site == "donor"),
                        ifelse(cbe_position <= 1, cbe_position + 1, cbe_position),
                        {
                          ifelse((enzyme == "CBE") & (splice_site == "acceptor"),
                                 ifelse(cbe_position <= 2, cbe_position + 1, cbe_position),
                                 ifelse(cbe_position <= 2, cbe_position + 1, cbe_position))
                        }) 
               },
               {
                 ifelse((enzyme == "CBE and ABE") & (splice_site == "donor"),
                        ifelse(cbe_position <= 1, cbe_position + 1, cbe_position),
                        {
                          ifelse((enzyme == "CBE") & (splice_site == "acceptor"),
                                 ifelse(cbe_position <= 2, cbe_position + 1, cbe_position),
                                 ifelse(cbe_position <= 1, cbe_position + 1, cbe_position))
                        }) 
               })
      },
      abe_position_tmp = {
        ifelse(strand == 1,
               {
                 ifelse((enzyme == "CBE and ABE") & (splice_site == "donor"),
                        ifelse(abe_position <= 1, abe_position + 1, abe_position),
                        {
                          ifelse((enzyme == "CBE") & (splice_site == "acceptor"),
                                 ifelse(abe_position <= 1, abe_position + 1, abe_position),
                                 ifelse(abe_position <= 0, abe_position + 1, abe_position))
                        }) 
               },
               {
                 ifelse((enzyme == "CBE and ABE") & (splice_site == "donor"),
                        ifelse(abe_position <= 0, abe_position + 1, abe_position),
                        {
                          ifelse((enzyme == "CBE") & (splice_site == "acceptor"),
                                 ifelse(abe_position <= 0, abe_position + 1, abe_position),
                                 ifelse(abe_position <= 0, abe_position + 1, abe_position))
                        }) 
               })
      }
    ) %>%
    mutate(
      chromStart = {
        ifelse(strand == 1,
               {
                 ifelse((enzyme == "CBE and ABE") & (splice_site == "donor"),
                        (exon_chrom_end - (guide_length - cbe_position_tmp - 1)),
                        {
                          ifelse((enzyme == "CBE") & (splice_site == "acceptor"),
                                 (exon_chrom_start - (guide_length - cbe_position_tmp)) - 1,
                                 (exon_chrom_start - abe_position_tmp) - 1)
                        }) 
               },
               {
                 ifelse((enzyme == "CBE and ABE") & (splice_site == "donor"),
                        (exon_chrom_start - cbe_position_tmp),
                        {
                          ifelse((enzyme == "CBE") & (splice_site == "acceptor"),
                                 (exon_chrom_end - (cbe_position_tmp - 2)),
                                 (exon_chrom_end - (guide_length - abe_position_tmp - 1)) + 1)
                        }) 
               })
      },
      chromEnd = {
        ifelse(strand == 1,
               {
                 ifelse((enzyme == "CBE and ABE") & (splice_site == "donor"),
                        (exon_chrom_end + cbe_position_tmp),
                        {
                          ifelse((enzyme == "CBE") & (splice_site == "acceptor"),
                                 (exon_chrom_start + (cbe_position_tmp - 2)),
                                 (exon_chrom_start + (guide_length - abe_position_tmp - 1)) - 1)
                        }) 
               },
               {
                 ifelse((enzyme == "CBE and ABE") & (splice_site == "donor"),
                        (exon_chrom_start + (guide_length - cbe_position_tmp - 1)),
                        {
                          ifelse((enzyme == "CBE") & (splice_site == "acceptor"),
                                 (exon_chrom_end + (guide_length - cbe_position_tmp + 1)),
                                 (exon_chrom_end + abe_position_tmp) + 1)
                        }) 
               })
      },
      chromStrand = {
        ifelse(strand == 1,
               {
                 ifelse((enzyme == "CBE and ABE") & (splice_site == "donor"),
                        ("-"),
                        {
                          ifelse((enzyme == "CBE") & (splice_site == "acceptor"),
                                 ("-"),
                                 ("+"))
                        }) 
               },
               {
                 ifelse((enzyme == "CBE and ABE") & (splice_site == "donor"),
                        ("+"),
                        {
                          ifelse((enzyme == "CBE") & (splice_site == "acceptor"),
                                 ("+"),
                                 ("-"))
                        }) 
               })
      }
    ) %>%
    dplyr::select(-abe_position_tmp, -cbe_position_tmp, -exon_chrom_start, -exon_chrom_end, -strand) %>%
    dplyr::rename(chrom = chromosome_name, strand = chromStrand)
} 

# Load weight data to calculate guide scores
motif_weights = read_tsv("motif_weights.tsv")
position_weights = read_tsv("position_weights.tsv")

max_weight = max(motif_weights$motif_weight) + max(position_weights$position_weight)
min_weight = min(motif_weights$motif_weight) + min(position_weights$position_weight)

cbe_motif_weights = motif_weights %>%
  filter(enzyme == cbe) %>%
  dplyr::rename(cbe_motif = motif, cbe_motif_weight = motif_weight) %>%
  dplyr::select(cbe_motif, cbe_motif_weight)

abe_motif_weights = motif_weights %>%
  filter(enzyme == abe) %>%
  dplyr::rename(abe_motif = motif, abe_motif_weight = motif_weight) %>%
  dplyr::select(abe_motif, abe_motif_weight)

cbe_position_weights = position_weights %>%
  filter(enzyme == cbe) %>%
  dplyr::rename(cbe_position = position, cbe_position_weight = position_weight) %>%
  dplyr::select(cbe_position, cbe_position_weight)

abe_position_weights = position_weights %>%
  filter(enzyme == abe) %>%
  dplyr::rename(abe_position = position, abe_position_weight = position_weight) %>%
  dplyr::select(abe_position, abe_position_weight)

filterGuides = function(runSpliceR.React,
                        enzymeClass.React,
                        min_editing_window,
                        max_editing_window,
                        splice_site.React,
                        strictFilter){
  runSpliceR.React %>%
    {if(enzymeClass.React == "CBE and-or ABE") {
      filter(.,
             (abe_position >= min_editing_window | cbe_position >= min_editing_window) &
               (abe_position <= max_editing_window | cbe_position <= max_editing_window)
      )
    } else {
      if(enzymeClass.React == "CBE and ABE") {
        filter(., enzyme == "CBE and ABE") %>%
          filter(
            (abe_position >= min_editing_window| cbe_position >= min_editing_window) &
              (abe_position <= max_editing_window | cbe_position <= max_editing_window)
          )
      } else {
        if(enzymeClass.React == "CBE") {
          filter(., enzyme == "CBE" | enzyme == "CBE and ABE") %>%
            filter(
              (cbe_position >= min_editing_window) & (cbe_position <= max_editing_window)
            )
        } else {
          filter(., enzyme == "ABE" | enzyme == "CBE and ABE") %>%
            filter(
              (abe_position >= min_editing_window) & (abe_position <= max_editing_window)
            )
        }
      }
    }} %>%
    {
      if(splice_site.React == "splice-donors") {
        filter(., splice_site == "donor")
      } else {
        if(splice_site.React == "splice-acceptors") {
          filter(., splice_site == "acceptor")
        } else {
          .
        }
      }
    } %>%
    {
      if(strictFilter) {
        filter(.,
               (abe_position >= min_editing_window) &
                 (abe_position <= max_editing_window) &
                 (cbe_position >= min_editing_window) &
                 (cbe_position <= max_editing_window)
        )
      } else {
        .
      }
    } %>%
    dplyr::rename(Exon = 1, `Splice Site` = 2, `Protospacer` = 3, `PAM` = 4, `Enzyme` = 5, `cDNA Disruption` = 6,
                  `CBE Position` = 7, `CBE Score` = 8, `ABE Position` = 9, `ABE Score` = 10, `Transcript ID` = 11, `Gene ID` = 12)
}