##### global.R for SpliceR

##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
# Copyright (C) 2019-2020 Mitchell Kluesner
#  
# This file is part of SpliceR
# 
# Please only copy and/or distribute this script with proper citation of SpliceR publication
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 

##### Set options for environment
options(shiny.sanitize.errors = FALSE)
httr::set_config(httr::config(ssl_verifypeer = FALSE))

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
library(printr)
library(rmarkdown)
library(DT)
library(httr)
library(curl)


# Given a chromosome, species, and genomic coordinates this function extracts genomic sequences. 1.12.18
# formerly getTranscript
getTranscriptExonCoordinates = function(id, upstream = 0, downstream = 0, species = "Homo_sapiens"){
  require(magrittr)
  require(tidyverse)
  url = paste0(
    "http://useast.ensembl.org/",
    species,
    "/Export/Output/Gene?db=core",";",
    "flank3_display=", as.character(upstream),";",
    "flank5_display=",as.character(downstream),";",
    "output=tab",";",
    "strand=feature",";",
    "t=",id,";",
    "param=gene",";",
    "miscset_ABC=yes",";",
    "miscset_RPCI-11=yes",";",
    "miscset_CHORI-17=yes",";",
    "miscset_WIBR-2=yes",";",
    "miscset_encode=yes",";",
    "miscset_encode_excluded=yes",";",
    "miscset_tilepath=yes",";",
    "_format=Text"
  )
  
  ensembl_data = url %>% read_tsv()
  ensembl_data %>%
    dplyr::filter(., transcript_id == id) %>%
    dplyr::select(start, end, strand, seqname, gene_id) %>%
    dplyr::rename(exon_chrom_start = start, exon_chrom_end = end, chromosome_name = seqname) %>%
    dplyr::mutate(strand = as.numeric(paste0(strand, 1))) %>%
    as_tibble()
}

### Function for getting sequence from genomic coordinates
### Formerly get_seq
coordinatesToDNAString = function(start, end, strand, chromosome, species = "Homo_sapiens", upstream = 0, downstream = 0){
  require(magrittr)
  require(tidyverse)
  require(Biostrings)
  seq = if(strand == 1) {
    start = start - upstream
    end = end + downstream
    paste0("http://useast.ensembl.org/", species, "/Export/Output/Location?db=core;flank3_display=0;flank5_display=0;output=fasta;r=",
           chromosome, ":", start, "-", end, ";strand=", strand, 
           ";utr5=yes;cdna=yes;intron=yes;utr3=yes;peptide=yes;coding=yes;genomic=unmasked;exon=yes;_format=Text") %>%
      readDNAStringSet() %>% .[[1]]
  } else {
    anti_end = start + upstream
    anti_start = end - downstream
    paste0("http://useast.ensembl.org/", species, "/Export/Output/Location?db=core;flank3_display=0;flank5_display=0;output=fasta;r=",
           chromosome, ":", anti_start, "-", anti_end, ";strand=", strand, 
           ";utr5=yes;cdna=yes;intron=yes;utr3=yes;peptide=yes;coding=yes;genomic=unmasked;exon=yes;_format=Text") %>%
      readDNAStringSet() %>% .[[1]]
  } 
  return(seq)
}

coordinatesToDNAChar = function(start, end, strand, chromosome, species = "Homo_sapiens", upstream = 0, downstream = 0){
  require(magrittr)
  require(tidyverse)
  require(Biostrings)
  # seq = if(strand == 1) {
  start = min(c(start, end))
  end = max(c(start, end))
    start = start - upstream
    end = end + downstream
   seq = paste0("http://useast.ensembl.org/", species, "/Export/Output/Location?db=core;flank3_display=0;flank5_display=0;output=fasta;r=",
           chromosome, ":", start, "-", end, ";strand=", strand, 
           ";utr5=yes;cdna=yes;intron=yes;utr3=yes;peptide=yes;coding=yes;genomic=unmasked;exon=yes;_format=Text") %>%
      readDNAStringSet() %>% .[[1]] %>% as.character()
  # } else {
  #   anti_end = start + upstream
  #   anti_start = end - downstream
  #   paste0("http://useast.ensembl.org/", species, "/Export/Output/Location?db=core;flank3_display=0;flank5_display=0;output=fasta;r=",
  #          chromosome, ":", anti_start, "-", anti_end, ";strand=", strand, 
  #          ";utr5=yes;cdna=yes;intron=yes;utr3=yes;peptide=yes;coding=yes;genomic=unmasked;exon=yes;_format=Text") %>%
  #     readDNAStringSet() %>% .[[1]] %>% as.character()
  # } 
  return(seq)
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
        # IF gene is sense...
        ifelse(strand == 1,
               {
                 # THEN assign as sense...
                 ifelse(
                   
                   # IF CBE and ABE AND a splice donor...
                   (enzyme == "CBE and ABE") & (splice_site == "donor"),
                   
                   # THEN assign as CBE and ABE splice donor...
                   ifelse(cbe_position <= 1, cbe_position + 1, cbe_position),
                   
                   
                   {
                     # ELSE IF CBE AND splice acceptor...
                     ifelse(
                       (enzyme == "CBE") & (splice_site == "acceptor"),
                       
                       # THEN assign as CBE and splice acceptor
                       ifelse(cbe_position <= 2, cbe_position + 1, cbe_position),
                       
                       # ELSE assing as ABE and splice acceptor
                       ifelse(cbe_position <= 2, cbe_position+ 1, cbe_position)
                     )
                   }
                 ) 
               },
               
               {
                 # THEN assign as antisense...
                 ifelse(
                   
                   # IF CBE and ABE AND a splice donor...
                   (enzyme == "CBE and ABE") & (splice_site == "donor"),
                   
                   # THEN assign as CBE and ABE splice donor...
                   ifelse(cbe_position <= 1, cbe_position + 1, cbe_position),
                   
                   
                   {
                     # ELSE IF CBE AND splice acceptor...
                     ifelse(
                       (enzyme == "CBE") & (splice_site == "acceptor"),
                       
                       # THEN assign as CBE and splice acceptor
                       ifelse(cbe_position <= 2, cbe_position + 1, cbe_position),
                       
                       # ELSE assing as ABE and splice acceptor
                       ifelse(cbe_position <= 1, cbe_position + 1, cbe_position)
                     )
                   }
                 ) 
               }
        )
      },
      
      abe_position_tmp = {
        # IF gene is sense...
        ifelse(strand == 1,
               {
                 # THEN assign as sense...
                 ifelse(
                   
                   # IF CBE and ABE AND a splice donor...
                   (enzyme == "CBE and ABE") & (splice_site == "donor"),
                   
                   # THEN assign as CBE and ABE splice donor...
                   ifelse(abe_position <= 1, abe_position + 1, abe_position),
                   
                   
                   {
                     # ELSE IF CBE AND splice acceptor...
                     ifelse(
                       (enzyme == "CBE") & (splice_site == "acceptor"),
                       
                       # THEN assign as CBE and splice acceptor
                       ifelse(abe_position <= 1, abe_position + 1, abe_position),
                       
                       # ELSE assesing as ABE and splice acceptor
                       ifelse(abe_position <= 0, abe_position + 1, abe_position)
                     )
                   }
                 ) 
               },
               
               {
                 # THEN assign as antisense...
                 ifelse(
                   
                   # IF CBE and ABE AND a splice donor...
                   (enzyme == "CBE and ABE") & (splice_site == "donor"),
                   
                   # THEN assign as CBE and ABE splice donor...
                   ifelse(abe_position <= 0, abe_position + 1, abe_position),
                   
                   
                   {
                     # ELSE IF CBE AND splice acceptor...
                     ifelse(
                       (enzyme == "CBE") & (splice_site == "acceptor"),
                       
                       # THEN assign as CBE and splice acceptor
                       ifelse(abe_position <= 0, abe_position + 1, abe_position),
                       
                       # ELSE assing as ABE and splice acceptor
                       ifelse(abe_position <= 0, abe_position + 1, abe_position)
                     )
                   }
                 ) 
               }
        )
      },
    ) %>%
    mutate(
      chromStart = {
        # IF gene is sense...
        ifelse(strand == 1,
               {
                 # THEN assign as sense...
                 ifelse(
                   
                   # IF CBE and ABE AND a splice donor...
                   (enzyme == "CBE and ABE") & (splice_site == "donor"),
                   
                   # THEN assign as CBE and ABE splice donor...
                   (exon_chrom_end - (guide_length - cbe_position_tmp - 1)),
                   
                   
                   {
                     # ELSE IF CBE AND splice acceptor...
                     ifelse(
                       (enzyme == "CBE") & (splice_site == "acceptor"),
                       
                       # THEN assign as CBE and splice acceptor
                       (exon_chrom_start - (guide_length - cbe_position_tmp)) - 1,
                       
                       # ELSE assign as ABE and splice acceptor
                       (exon_chrom_start - abe_position_tmp) - 1
                     )
                   }
                 ) 
               },
               
               {
                 # THEN assign as antisense...
                 ifelse(
                   
                   # IF CBE and ABE AND a splice donor...
                   (enzyme == "CBE and ABE") & (splice_site == "donor"),
                   
                   # THEN assign as CBE and ABE splice donor...
                   (exon_chrom_start - cbe_position_tmp),
                   
                   
                   {
                     # ELSE IF CBE AND splice acceptor...
                     ifelse(
                       (enzyme == "CBE") & (splice_site == "acceptor"),
                       
                       # THEN assign as CBE and splice acceptor
                       (exon_chrom_end - (cbe_position_tmp - 2)),
                       
                       # ELSE assign as ABE and splice acceptor
                       (exon_chrom_end - (guide_length - abe_position_tmp - 1)) + 1
                     )
                   }
                 ) 
               }
        )
      },
      
      chromEnd = {
        # IF gene is sense...
        ifelse(strand == 1,
               {
                 # THEN assign as sense...
                 ifelse(
                   
                   # IF CBE and ABE AND a splice donor...
                   (enzyme == "CBE and ABE") & (splice_site == "donor"),
                   
                   # THEN assign as CBE and ABE splice donor...
                   (exon_chrom_end + cbe_position_tmp),
                   
                   
                   {
                     # ELSE IF CBE AND splice acceptor...
                     ifelse(
                       (enzyme == "CBE") & (splice_site == "acceptor"),
                       
                       # THEN assign as CBE and splice acceptor
                       (exon_chrom_start + (cbe_position_tmp - 2)),
                       
                       # ELSE assign as ABE and splice acceptor
                       (exon_chrom_start + (guide_length - abe_position_tmp - 1)) - 1
                     )
                   }
                 ) 
               },
               
               {
                 # THEN assign as antisense...
                 ifelse(
                   
                   # IF CBE and ABE AND a splice donor...
                   (enzyme == "CBE and ABE") & (splice_site == "donor"),
                   
                   # THEN assign as CBE and ABE splice donor...
                   (exon_chrom_start + (guide_length - cbe_position_tmp - 1)),
                   
                   
                   {
                     # ELSE IF CBE AND splice acceptor...
                     ifelse(
                       (enzyme == "CBE") & (splice_site == "acceptor"),
                       
                       # THEN assign as CBE and splice acceptor
                       (exon_chrom_end + (guide_length - cbe_position_tmp + 1)),
                       
                       # ELSE assign as ABE and splice acceptor
                       (exon_chrom_end + abe_position_tmp) + 1
                     )
                   }
                 ) 
               }
        )
      },
      
      chromStrand = {
        # IF gene is sense...
        ifelse(strand == 1,
               {
                 # THEN assign as sense...
                 ifelse(
                   
                   # IF CBE and ABE AND a splice donor...
                   (enzyme == "CBE and ABE") & (splice_site == "donor"),
                   
                   # THEN assign as CBE and ABE splice donor...
                   ("-"),
                   
                   
                   {
                     # ELSE IF CBE AND splice acceptor...
                     ifelse(
                       (enzyme == "CBE") & (splice_site == "acceptor"),
                       
                       # THEN assign as CBE and splice acceptor
                       ("-"),
                       
                       # ELSE assign as ABE and splice acceptor
                       ("+")
                     )
                   }
                 ) 
               },
               
               {
                 # THEN assign as antisense...
                 ifelse(
                   
                   # IF CBE and ABE AND a splice donor...
                   (enzyme == "CBE and ABE") & (splice_site == "donor"),
                   
                   # THEN assign as CBE and ABE splice donor...
                   ("+"),
                   
                   
                   {
                     # ELSE IF CBE AND splice acceptor...
                     ifelse(
                       (enzyme == "CBE") & (splice_site == "acceptor"),
                       
                       # THEN assign as CBE and splice acceptor
                       ("+"),
                       
                       # ELSE assign as ABE and splice acceptor
                       ("-")
                     )
                   }
                 ) 
               }
        )
      }
      
    ) %>%
    dplyr::select(-abe_position_tmp, -cbe_position_tmp,-exon_chrom_start, -exon_chrom_end, -strand) %>%
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


#### Shiny Server Functions
# IF CBE or ABE is selected, THEN do not filter, else

filterGuides = function(runSpliceR.React,
         enzymeClass.React,
         min_editing_window,
         max_editing_window,
         splice_site.React,
         strictFilter
         ){
runSpliceR.React %>%
{if(enzymeClass.React == "CBE and-or ABE") {
  {.} %>%
    filter(.,
           (abe_position >= min_editing_window | cbe_position >= min_editing_window) &
             (abe_position <= max_editing_window | cbe_position <= max_editing_window)
    )
} else {
  
  # IF CBE and ABE is selected, THEN include only CBE and ABE
  if(enzymeClass.React == "CBE and ABE") {
    filter(., enzyme == "CBE and ABE") %>%
      filter(
        (abe_position >= min_editing_window| cbe_position >= min_editing_window) &
          (abe_position <= max_editing_window | cbe_position <= max_editing_window)
      )
  } else {
    
    # IF CBE is only base editor selected, THEN include only CBE AND CBE and ABE
    if(enzymeClass.React == "CBE") {
      filter(., enzyme == "CBE" | enzyme == "CBE and ABE") %>%
        filter(
          (cbe_position >= min_editing_window) & (cbe_position <= max_editing_window)
        )
    } else {
      
      # IF ABE is only base editor selected, THEN include only ABE AND CBE and ABE
      filter(., enzyme == "ABE" | enzyme == "CBE and ABE") %>%
        filter(
          (abe_position >= min_editing_window) & (abe_position <= max_editing_window)
        )
    }
  }
}
  } %>%
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




