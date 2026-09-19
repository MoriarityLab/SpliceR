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
library(tibble)
library(ggplot2)
library(printr)
library(readr)
library(rmarkdown)
library(DT)
library(httr)
library(jsonlite)

# Global Ensembl REST endpoint
ENSEMBL_REST_URL <- "https://rest.ensembl.org"

# Query Ensembl REST API to extract exon genomic coordinates and strand for a given transcript ID
getTranscriptExonCoordinates = function(id, upstream = 0, downstream = 0, species = "Homo_sapiens"){
  clean_id <- sub("\\..*", "", id)
  
  url <- paste0(ENSEMBL_REST_URL, "/lookup/id/", clean_id, "?expand=1")
  
  response <- GET(
    url = url,
    add_headers(`User-Agent` = "SpliceR-App/2.0 (R/httr)"),
    content_type_json()
  )
  
  if (status_code(response) >= 400) {
    err_body <- content(response, as = "text", encoding = "UTF-8")
    stop(paste("Ensembl REST Request failed [HTTP", status_code(response), "]:", err_body))
  }
  
  tx_data <- content(response, as = "parsed", simplifyVector = FALSE)
  
  if (is.null(tx_data) || is.null(tx_data$Exon)) {
    stop(paste("Transcript ID or exons not found via Ensembl REST:", id))
  }
  
  chr_name  <- as.character(tx_data$seq_region_name)
  strand_val <- as.numeric(tx_data$strand)
  gene_id   <- tx_data$Parent
  
  exons_list <- tx_data$Exon
  
  df_exons <- do.call(rbind, lapply(exons_list, function(ex) {
    data.frame(
      exon_chrom_start = as.numeric(ex$start),
      exon_chrom_end   = as.numeric(ex$end),
      strand           = strand_val,
      chromosome_name  = chr_name,
      gene_id          = gene_id,
      transcript_id    = id,
      stringsAsFactors = FALSE
    )
  }))
  
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

# Fetch genomic DNA sequence via Ensembl REST and return a Biostrings DNAString
coordinatesToDNAString = function(start, end, strand, chromosome, species = "Homo_sapiens", upstream = 0, downstream = 0){
  seq_str <- coordinatesToDNAChar(start, end, strand, chromosome, species, upstream, downstream)
  return(Biostrings::DNAString(seq_str))
}

# Fetch genomic DNA sequence via Ensembl REST and return a raw character string
coordinatesToDNAChar = function(start, end, strand, chromosome, species = "Homo_sapiens", upstream = 0, downstream = 0){
  c_start <- min(start, end) - upstream
  c_end   <- max(start, end) + downstream
  
  # Ensure strand is represented as 1 or -1
  strand_param <- ifelse(strand == -1, -1, 1)
  
  region_str <- paste0(chromosome, ":", c_start, "..", c_end, ":", strand_param)
  url <- paste0(ENSEMBL_REST_URL, "/sequence/region/", species, "/", region_str)
  
  response <- GET(
    url = url,
    add_headers(`User-Agent` = "SpliceR-App/2.0 (R/httr)"),
    content_type_json()
  )
  
  if (status_code(response) >= 400) {
    err_body <- content(response, as = "text", encoding = "UTF-8")
    stop(paste("Ensembl Sequence Request failed [HTTP", status_code(response), "]:", err_body))
  }
  
  res_content <- content(response, as = "parsed", simplifyVector = FALSE)
  dna_seq <- res_content$seq
  
  return(dna_seq)
}

matchPatterns = Vectorize(matchPattern, vectorize.args = "subject")

extractGuideStart = Vectorize(
  FUN = function(alignments, exon){alignments[[exon]]@ranges@start},
  vectorize.args = "exon"
)

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

motif_weights = read_tsv("motif_weights.tsv")
position_weights = read_tsv("position_weights.tsv")

max_weight = max(motif_weights$motif_weight) + max(position_weights$position_weight)
min_weight = min(motif_weights$motif_weight) + min(position_weights$position_weight)

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