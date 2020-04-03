### v1.0.0

### Packages
library(shiny)
library(Biostrings)
library(magrittr)
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

### Parameters

### General parameters
flank = 20
downstream_gene_seq = 20
remove_end_guides = FALSE
max_pos = 13

### Scoring weights
general_filter = FALSE
normalization = TRUE
drop_low_guides = FALSE
low_cutoff = 1
efficiency_weight = 0.5
position_weight = 0.25
motif_weight = 0.1
splice_site_weight = 0.25
total_weight = if(normalization == TRUE) {efficiency_weight + position_weight + motif_weight + splice_site_weight} else {1}


# Given a chromosome, species, and genomic coordinates this function extracts genomic sequences. 1.12.18
get_transcript = function(id, upstream = 0, downstream = 0, species = "Homo_sapiens"){
  require(magrittr)
  require(tidyverse)
  paste0(
    "https://useast.ensembl.org/",
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
  ) %>%
    readr::read_delim(., delim = "\t") %>% 
    dplyr::filter(., transcript_id == id) %>%
    dplyr::select(start, end, strand, seqname, gene_id) %>%
    dplyr::rename(exon_chrom_start = start, exon_chrom_end = end, chromosome_name = seqname) %>%
    dplyr::mutate(strand = as.numeric(paste0(strand, 1))) %>%
    as.data.frame()
}

### Function for getting sequence from genomic coordinates
get_seq = function(start, end, strand, chromosome, species = "Homo_sapiens", upstream = 0, downstream = 0){
  require(magrittr)
  require(tidyverse)
  require(Biostrings)
  seq = if(strand == 1) {
    start = start - upstream
    end = end + downstream
    paste0("https://useast.ensembl.org/", species, "/Export/Output/Location?db=core;flank3_display=0;flank5_display=0;output=fasta;r=",
           chromosome, ":", start, "-", end, ";strand=", strand, 
           ";utr5=yes;cdna=yes;intron=yes;utr3=yes;peptide=yes;coding=yes;genomic=unmasked;exon=yes;_format=Text") %>%
      readDNAStringSet %>% .[[1]]
  } else {
    anti_end = start + upstream
    anti_start = end - downstream
    paste0("https://useast.ensembl.org/", species, "/Export/Output/Location?db=core;flank3_display=0;flank5_display=0;output=fasta;r=",
           chromosome, ":", anti_start, "-", anti_end, ";strand=", strand, 
           ";utr5=yes;cdna=yes;intron=yes;utr3=yes;peptide=yes;coding=yes;genomic=unmasked;exon=yes;_format=Text") %>%
      readDNAStringSet %>% .[[1]]
  }
  return(seq)
}

# Extract the bases from a DNAStringSet object (target), given the position in the string set (position in string), the start within the substring (start) and the end within the substring (end). If only a single base is required, it is signified by the start argument where end is default to "". This function is vectorized to pull out multiple bases from multiple strings in a DNAStringSet (Bioconductor). This function returns a list object. 

extract_bases = Vectorize(
  function(string_set, sub_string, start, end){
    string_set[[sub_string]][start:end]},
  vectorize.args = c("sub_string","start", "end")
)

# General function to convert from a DNAStringSet object to a character vector of the DNA sequences.
# Function modified 1.14.18 to incorporate apply statements instead of for loops
unlist_DNAStringSet = function(DNAStringSet_object){
  as.vector(
    ldply(mapply("[", DNAStringSet_object), .fun = as.character)[,1]
  )
}

# Function to subset the target sequence surrounding each intron-exon or exon-intron boundary given a single DNAString
subset_target_seq = Vectorize(
  function(sequence, upstream_index, downstream_index){
    sequence[upstream_index:downstream_index]
  }, vectorize.args = c("upstream_index", "downstream_index")
)

# General vectorized function to extract some subsequence given a subject, exon number and coordinates
# **Need to revise to use extract_bases()**
extract_subsequence = Vectorize(
  function(targets, exons, subsequence){targets[[exons]][subsequence]},
  vectorize.args = "exons"
)

# Vectorized 'matchPatterns' function
matchPatterns = Vectorize(matchPattern, vectorize.args = "subject")

# Function to pull out start site of the guides with respect to the target sequence
extract_start = Vectorize(
  FUN = function(alignments, exon){alignments[[exon]]@ranges@start},
  vectorize.args = "exon"
)

# Function to pull out the guide as a character
extract_guide = Vectorize(
  FUN = function(alignments, exon){alignments[[exon]] %>% as.character},
  vectorize.args = "exon"
)

# Function to pull out the exon targetted for a given guide
extract_exon = Vectorize(
  FUN = function(alignments, exon){alignments[[exon]] %>% as.character},
  vectorize.args = "exon"
)

# Takes in a character vector and puts out a character vector
reverseComplements = function(x){
  DNAStrings = Vectorize(DNAString, vectorize.args = c("x"))
  reverseComplement_s = Vectorize(reverseComplement, vectorize.args = c("x"))
  reverseComplement_s(DNAStrings(x)) %>% unlist_DNAStringSet()
}

countPatterns = countPattern %>% Vectorize(., vectorize.args = c("pattern", "subject"))

abe_dinucleotide.fn = function(x, splice_site){
  if(splice_site == "acceptor") {x %>% DNAStringSet %>% subseq(., 1, 2) %>% as.character} else {x %>% DNAStringSet %>% subseq(., 5, 6) %>% reverseComplement %>% as.character}
}

score_motif = function(splice_site, motif, guide_motifs){
  guide_number = length(splice_site)
  ifelse(splice_site == "acceptor", 
         ifelse(countPatterns(extract_bases(string_set = guide_motifs, sub_string = 1:guide_number, start = 1, end = 4)[1,] %>% DNAStringSet, DNAStringSet("NHNN"), fixed = FALSE) == 1, 0,
                ifelse(countPatterns(extract_bases(string_set = guide_motifs, sub_string = 1:guide_number, start = 1, end = 4)[1,] %>% DNAStringSet, DNAStringSet("NGHN"), fixed = FALSE) == 1, 1,
                       ifelse(countPatterns(extract_bases(string_set = guide_motifs, sub_string = 1:guide_number, start = 1, end = 4)[1,] %>% DNAStringSet, DNAStringSet("NGGN"), fixed = FALSE) == 1, 2,0))),
         ifelse(countPatterns(extract_bases(string_set = guide_motifs, sub_string = 1:guide_number, start = 1, end = 4)[1,] %>% DNAStringSet, DNAStringSet("NNHN"), fixed = FALSE) == 1, 0,
                ifelse(countPatterns(extract_bases(string_set = guide_motifs, sub_string = 1:guide_number, start = 1, end = 4)[1,] %>% DNAStringSet, DNAStringSet("NHGN"), fixed = FALSE) == 1, 1,
                       ifelse(countPatterns(extract_bases(string_set = guide_motifs, sub_string = 1:guide_number, start = 1, end = 4)[1,] %>% DNAStringSet, DNAStringSet("NGGN"), fixed = FALSE) == 1, 2,0)))
  )
}

revcom = function(x){
  as.character(Biostrings::reverseComplement(Biostrings::DNAString(x)))
}


## Scoring Efficiencies
# Scoring efficiencies based on position according to komor et al., 2016

be_efficiencies = data.frame(stringsAsFactors=FALSE,
                             be_dinucleotide = c("CC", "CC", "CC", "CC", "CC", "CC", "CC", "CC", "CC", "CC",
                                                 "CC", "CC", "CC", "CC", "CC", "CC", "CC", "CC", "CC", "CC",
                                                 "CC", "CC", "CC", "CC", "CC", "CC", "CC", "CC", "CC", "CC",
                                                 "TC", "TC", "TC", "TC", "TC", "TC", "TC", "TC", "TC", "TC",
                                                 "TC", "TC", "TC", "TC", "TC", "TC", "TC", "TC", "TC", "TC",
                                                 "TC", "TC", "TC", "TC", "TC", "TC", "TC", "TC", "TC", "TC",
                                                 "AC", "AC", "AC", "AC", "AC", "AC", "AC", "AC", "AC", "AC",
                                                 "AC", "AC", "AC", "AC", "AC", "AC", "AC", "AC", "AC", "AC", "AC",
                                                 "AC", "AC", "AC", "AC", "AC", "AC", "AC", "AC", "AC", "GC",
                                                 "GC", "GC", "GC", "GC", "GC", "GC", "GC", "GC", "GC", "GC",
                                                 "GC", "GC", "GC", "GC", "GC", "GC", "GC", "GC", "GC", "GC",
                                                 "GC", "GC", "GC", "GC", "GC", "GC", "GC", "GC", "GC", "NN"),
                             be_position = c(-9L, -8L, -7L, -6L, -5L, -4L, -3L, -2L, -1L, 0L, 1L, 2L,
                                             3L, 4L, 5L, 6L, 7L, 8L, 9L, 10L, 11L, 12L, 13L, 14L, 15L,
                                             16L, 17L, 18L, 19L, 20L, -9L, -8L, -7L, -6L, -5L, -4L, -3L,
                                             -2L, -1L, 0L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L, 10L, 11L,
                                             12L, 13L, 14L, 15L, 16L, 17L, 18L, 19L, 20L, -9L, -8L, -7L, -6L,
                                             -5L, -4L, -3L, -2L, -1L, 0L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L,
                                             9L, 10L, 11L, 12L, 13L, 14L, 15L, 16L, 17L, 18L, 19L, 20L,
                                             -9L, -7L, -8L, -6L, -5L, -4L, -3L, -2L, -1L, 0L, 1L, 2L, 3L,
                                             4L, 5L, 6L, 7L, 8L, 9L, 10L, 11L, 12L, 13L, 14L, 15L, 16L,
                                             17L, 18L, 19L, 20L, 0L),
                             be_efficiency = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.807387863, 21.93403694,
                                               11.16886544, 18.30079156, 41.176781, 43.3298153,
                                               45.21372032, 17.08970976, 8.073878628, 2.960422164, 1.07651715, 0, 0,
                                               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2.556728232, 0,
                                               7.939313984, 12.24538259, 32.02638522, 43.73350923,
                                               48.98153034, 53.69129288, 51, 51, 48.98153034, 32.02638522,
                                               12.24538259, 7.939313984, 2.556728232, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                               0, 0, 0, 0, 0, 0, 1.07651715, 4.575197889, 10.7651715,
                                               23.01055409, 44.13720317, 26.91292876, 38.21635884, 13.32189974,
                                               5.517150396, 2.691292876, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                               0, 0, 0, 0, 0, 0, 0, 0, 0, 0.134564644, 0.269129288,
                                               0.941952507, 8.612137203, 26.91292876, 34.17941953, 5.248021108,
                                               1.480211082, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
)

abe_efficiencies = data.frame(stringsAsFactors=FALSE,
                              abe_dinucleotide = c("CA", "CA", "CA", "CA", "CA", "CA", "CA", "CA", "CA",
                                                   "CA", "CA", "CA", "CA", "CA", "CA", "CA", "CA", "CA",
                                                   "CA", "CA", "CA", "CA", "CA", "CA", "CA", "CA", "CA",
                                                   "CA", "CA", "CA", "TA", "TA", "TA", "TA", "TA", "TA", "TA",
                                                   "TA", "TA", "TA", "TA", "TA", "TA", "TA", "TA", "TA",
                                                   "TA", "TA", "TA", "TA", "TA", "TA", "TA", "TA", "TA", "TA",
                                                   "TA", "TA", "TA", "TA", "AA", "AA", "AA", "AA", "AA",
                                                   "AA", "AA", "AA", "AA", "AA", "AA", "AA", "AA", "AA", "AA",
                                                   "AA", "AA", "AA", "AA", "AA", "AA", "AA", "AA", "AA",
                                                   "AA", "AA", "AA", "AA", "AA", "AA", "GA", "GA", "GA", "GA",
                                                   "GA", "GA", "GA", "GA", "GA", "GA", "GA", "GA", "GA",
                                                   "GA", "GA", "GA", "GA", "GA", "GA", "GA", "GA", "GA", "GA",
                                                   "GA", "GA", "GA", "GA", "GA", "GA", "GA", "NN"),
                              abe_position = c(-9L, -8L, -7L, -6L, -5L, -4L, -3L, -2L, -1L, 0L, 1L,
                                               2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L, 10L, 11L, 12L, 13L,
                                               14L, 15L, 16L, 17L, 18L, 19L, 20L, -9L, -8L, -7L, -6L, -5L,
                                               -4L, -3L, -2L, -1L, 0L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L,
                                               9L, 10L, 11L, 12L, 13L, 14L, 15L, 16L, 17L, 18L, 19L, 20L,
                                               -9L, -8L, -7L, -6L, -5L, -4L, -3L, -2L, -1L, 0L, 1L, 2L,
                                               3L, 4L, 5L, 6L, 7L, 8L, 9L, 10L, 11L, 12L, 13L, 14L,
                                               15L, 16L, 17L, 18L, 19L, 20L, -9L, -7L, -8L, -6L, -5L, -4L,
                                               -3L, -2L, -1L, 0L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L,
                                               10L, 11L, 12L, 13L, 14L, 15L, 16L, 17L, 18L, 19L, 20L, 0L),
                              abe_efficiency = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.807387863, 21.93403694,
                                                 11.16886544, 18.30079156, 41.176781, 43.3298153,
                                                 45.21372032, 17.08970976, 8.073878628, 2.960422164, 1.07651715,
                                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                 2.556728232, 0, 7.939313984, 12.24538259, 32.02638522,
                                                 43.73350923, 48.98153034, 53.69129288, 51, 51, 48.98153034,
                                                 32.02638522, 12.24538259, 7.939313984, 2.556728232, 0, 0, 0,
                                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.07651715,
                                                 4.575197889, 10.7651715, 23.01055409, 44.13720317,
                                                 26.91292876, 38.21635884, 13.32189974, 5.517150396, 2.691292876,
                                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                 0, 0, 0.134564644, 0.269129288, 0.941952507, 8.612137203,
                                                 26.91292876, 34.17941953, 5.248021108, 1.480211082, 0, 0,
                                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
)

pam_efficiency = data.frame(stringsAsFactors=FALSE,
                            pam = c("AGG", "CGG", "GGG", "TGG", "AGT", "CGT", "TGT", "GGT",
                                    "AGC", "CGC", "TGC", "GGC", "AGA", "CGA", "TGA", "GGA"),
                            pam_weight = c(1, 1, 1, 1, 0.75, 0.75, 0.75, 0.75, 0.22, 0.22, 0.22, 0.22,
                                           0.55, 0.55, 0.55, 0.55)
                            
                            
)

generate_guides = function(ensembl_id, pam, species){
  transcript = ensembl_id #%>% strsplit(., split = ".", fixed = TRUE) %>% unlist(.) %>% .[1]
  
  pam_char = pam
  pam %<>% DNAString()
  pam_length = as.numeric(length(pam))
  guide = pam %>% as.character() %>% paste("NNNNNNNNNNNNNNNNNNNN", ., sep = "") %>% DNAString()
  # This line was altered 8.5.19, Ensembl must have changed their querying to work with the version number instead of just the transcrip id 
  
  gene = get_transcript(transcript, species = species)
  sense = gene$strand[1]
  
  # The one base is to account for the fact that nucleotide counting is base 1
  # Odd facet that sometimes the exons are out of order, need to pay attention to this, not sure why?
  if(sense == 1){ 
    gene = gene %>% arrange(., exon_chrom_start)
    exon_number = length(gene$exon_chrom_start)
    exon_intron_length = (gene$exon_chrom_end[exon_number] - gene$exon_chrom_start[1])+ 1 + 20 + pam_length
    gene$exon_start = (gene$exon_chrom_start-gene$exon_chrom_start[1]) + 1 + 20 + pam_length
    gene$exon_end = (gene$exon_chrom_end-gene$exon_chrom_start[1]) + 1 + 20 + pam_length
    seq_chrom_start = min(gene$exon_chrom_start) - (flank+pam_length)
    gene$sa_seq_coord = gene$exon_chrom_start - (flank+pam_length)
    gene$sd_seq_coord = gene$exon_chrom_end - (flank+pam_length)
    gene %<>% mutate(exon_length = exon_end - exon_start) %<>%
      mutate(sd_position = 1-cumsum(exon_length)/sum(exon_length)) %<>%
      mutate(sa_position = 1-(cumsum(exon_length)-exon_length)/sum(exon_length))
    
    gene_seq = get_seq(start = gene$exon_chrom_start[1], 
                       end = gene$exon_chrom_end[exon_number],
                       strand = gene$strand[1],
                       chromosome = gene$chromosome_name[1],
                       upstream = flank + pam_length,
                       downstream = downstream_gene_seq + pam_length,
                       species = species
    )
  } else
  {# Keep in mind this is in reverse orientation and the ends really signify the starts of exons.
    gene = gene %>% arrange(., desc(exon_chrom_start))
    exon_number = length(gene$exon_chrom_start)
    exon_intron_length = (gene$exon_chrom_end[1] - gene$exon_chrom_start[exon_number]) + 1 + 20 + pam_length 
    # Sense to the transcript
    gene$exon_start = (gene$exon_chrom_end[1] - gene$exon_chrom_end) + 1 + 20 + pam_length
    gene$exon_end = (gene$exon_chrom_end[1] - gene$exon_chrom_start) + 1 + 20 + pam_length
    seq_chrom_start = max(gene$exon_chrom_end) + (flank+pam_length)
    gene$sa_seq_coord = gene$exon_chrom_end + (flank+pam_length)
    gene$sd_seq_coord = gene$exon_chrom_start + (flank+pam_length)
    gene %<>% mutate(exon_length = exon_end- exon_start) %<>%
      mutate(sd_position = 1-cumsum(exon_length)/sum(exon_length)) %<>%
      mutate(sa_position = 1-(cumsum(exon_length)-exon_length)/sum(exon_length))
    
    gene_seq = get_seq(gene$exon_chrom_end[1], 
                       gene$exon_chrom_start[exon_number],
                       gene$strand[1],
                       chromosome = gene$chromosome_name[1],
                       upstream = flank + pam_length,
                       downstream = downstream_gene_seq + pam_length,
                       species = species
    )
  }
  
  # NAs are used because there is no intron upstream and downstream of the first and last exon respectively. Update 11.26.17, altered to simply call the 3' UTR and exon as it is possible that there may be an overlapping exon with another transcript, thus it is important to identify guides in this region as well.
  gene$exon = 1:exon_number
  
  gene$intron_start = c((gene$exon_end[1:(exon_number)])+1)
  
  if(exon_number > 1) {gene$intron_end = c((gene$exon_start[2:exon_number])-1, length(gene_seq))} else
  {gene$intron_end = length(gene_seq)}
  
  # Isolate splice acceptor and splice donor target sequences, **Need to revise to use extract_bases()**
  ### Error initiates here for sample 6
  # Appears error is the result of not having enough bases upstream in the last exon
  # thus the upstream index of the is going into nothing-ness?
  
  # 1.10.18 error starts here? Changed from %>% DNAStringSet to %>% unlist_DNAStringSet
  gene$sa_seq = subset_target_seq(sequence = gene_seq, 
                                  upstream_index = gene$exon_start-flank-pam_length,
                                  downstream_index = gene$exon_start+flank+pam_length) %>% unlist_DNAStringSet
  
  # Modified on 11.26.17 to account for the 3' UTR.
  # For ABE, this code works fine for subsetting 20 bases downstream, but not 15 bases.
  # Maybe would be best to drop the upstream/downstream all together and just use a "flank" value.
  gene$sd_seq = subset_target_seq(sequence = gene_seq, 
                                  upstream_index = gene$exon_end-flank-pam_length,
                                  downstream_index = c(gene$exon_end[1:exon_number]+flank+pam_length)) %>%
    unlist_DNAStringSet
  
  
  # Isolate motifs of splice acceptor and splice donor.
  # Modified 11.26.17 to account for the 3' UTR guides across multiple transcript variants.
  # Modified 1.10.18 to account for reversion back to DNAStringSet. Added 2nd DNAStringSet and last unlist_DNAStringSet
  gene$acceptor = gene$sa_seq %>% DNAStringSet %>% extract_bases(., sub_string = 1:exon_number,
                                                                 start = ((flank - 2) + pam_length),
                                                                 end = ((flank + 3) + pam_length)) %>% unlist_DNAStringSet
  
  gene$donor = gene$sd_seq %>% DNAStringSet %>% extract_bases(., sub_string = 1:exon_number, 
                                                              start = (flank + pam_length-1), 
                                                              end = ((flank + 4) + pam_length)) %>% unlist_DNAStringSet
  
  ##################################################################
  ### Identify guides
  
  # Identify raw splice acceptor and splice donor guides
  # Modified 1.10.18 to account for reversion back to DNAStringSet. %>% DNAStringSet added.
  # To solve the issue of the flank over-extending the region for SAs, the subseq function was used to restrict the range
  splice_acceptors = gene$sa_seq %>%
    DNAStringSet %>%
    reverseComplement %>%
    matchPatterns(pattern = guide, subject = ., fixed = FALSE)
  
  splice_acceptors_abe = gene$sa_seq %>%
    DNAStringSet %>%
    matchPatterns(pattern = guide, subject = ., fixed = FALSE)
  
  splice_donors= gene$sd_seq %>%
    DNAStringSet %>%
    reverseComplement %>%
    matchPatterns(pattern = guide, subject = ., fixed = FALSE)
  
  # Pull out guide character, start position and number the guide for both acceptors and donors
  unscored_sa_guides = extract_guide(splice_acceptors, exon = 1:exon_number)%>% 
    unlist %>%
    DNAStringSet %>%
    reverseComplement %>%
    unlist_DNAStringSet # still 20 + pam_length long
  sa_starts = extract_start(splice_acceptors, exon = 1:exon_number)
  sa_guide_number = length(unlist(sa_starts))
  
  # still 20 + pam_length long
  # Extra code is to reverse complement the abe guides.
  unscored_sa_abe_guides = extract_guide(splice_acceptors_abe, exon = 1:exon_number) %>% 
    unlist %>%
    DNAStringSet %>%
    reverseComplement %>%
    unlist_DNAStringSet
  sa_abe_starts = extract_start(splice_acceptors_abe, exon = 1:exon_number)
  sa_abe_guide_number = length(unlist(sa_abe_starts))
  
  unscored_sd_guides = extract_guide(splice_donors, exon = 1:exon_number) %>% 
    unlist %>%
    DNAStringSet %>%
    reverseComplement %>%
    unlist_DNAStringSet
  sd_starts = extract_start(splice_donors, exon = 1:exon_number)
  sd_guide_number = length(unlist(sd_starts))
  
  # Establish a dataframe with the information about the guides
  guides = data.frame(
    exon = c(sapply(splice_acceptors, length) %>% rep(c(1:exon_number), times = .),
             sapply(splice_acceptors_abe, length) %>% rep(c(1:exon_number), times = .),
             sapply(splice_donors, length) %>% rep(c(1:exon_number), times = .)),
    guide = c(unlist(unscored_sa_guides),unlist(unscored_sa_abe_guides),unlist(unscored_sd_guides)),
    guide_position = c(unlist(sa_starts), unlist(sa_abe_starts), unlist(sd_starts)),
    splice_site = c("acceptor", "acceptor", "donor") %>% rep(., times = c(sa_guide_number, sa_abe_guide_number, sd_guide_number)),
    be = c("BE3", "ABE", "BE3") %>% rep(., times = c(sa_guide_number, sa_abe_guide_number, sd_guide_number)),
    enzyme = c("CBE", "ABE", "CBE or ABE") %>% rep(., times = c(sa_guide_number, sa_abe_guide_number, sd_guide_number))
  ) # Guide is still 20 + pam_length nt
  
  # Convert the sequences from DNAStringSet objects to characters
  gene$sa_seq = gene$sa_seq %>% unlist_DNAStringSet()
  gene$sd_seq = gene$sd_seq %>% unlist_DNAStringSet()
  gene$acceptor = gene$acceptor %>% unlist_DNAStringSet()
  gene$donor %<>% unlist_DNAStringSet() %<>% {ifelse(.== "NANA", NA, .)}
  
  #############################################################
  ### Conditional if statement regarding number of guides
  
  if(length(guides$exon) > 0) { # Edited 
    guides$guide_number = 1:length(guides$exon)
    
    guides = gene %>% 
      dplyr::select(exon, acceptor, donor) %>% 
      tidyr::gather(exon) %>%
      mutate(splice_site = exon, exon = rep(1:exon_number, times = 2)) %>%
      set_names(c("exon","splice_site","motif")) %>%
      bind_cols(., (gene %>% 
                      dplyr::select(exon, sa_position, sd_position) %>% 
                      tidyr::gather(exon) %>% 
                      mutate(exon = {ifelse(exon == "sa_position", "acceptor", "donor")}) %>%
                      mutate(splice_site = exon, exon = rep(1:exon_number, times = 2)) %>%
                      set_names(c("exon","position_score", "splice_site")) %>%
                      dplyr::select(position_score))
      ) %>%
      set_names(c("exon","motif","splice_site","position_score")) %>%
      inner_join(guides, .)
    guide_number = length(guides$exon)
    
    guides = guides %>% 
      mutate(pam = guide %>% DNAStringSet %>% reverseComplement %>% subseq(., start = 21) %>% unlist_DNAStringSet(),
             guide = guide %>% DNAStringSet %>% reverseComplement %>% subseq(., end = 20) %>% unlist_DNAStringSet())
    
    ####################################################### 
    ### Score Motif
    guides = guides %>% mutate(full_motif = motif, motif = motif %>% DNAStringSet %>% subseq(., 2,5) %>% unlist_DNAStringSet())
    guides = guides %>% mutate(be_motif_score = score_motif(splice_site, motif, guide_motifs = guides$motif))
    
    ### For ABE, this same code can be copied and then simply change the start and end respectively
    ### However this code should also be deprectated as there is the subseq function from Biostrings which appears to be much more efficient
    guides = guides %>% 
      dplyr::select(splice_site) %>% 
      {ifelse(.=="acceptor",
              extract_bases(string_set = guides$motif %>% DNAStringSet(),
                            sub_string = 1:length(guides$guide_number),
                            start = 2, end = 3),
              extract_bases(string_set = guides$motif %>% DNAStringSet(),
                            sub_string = 1:length(guides$guide_number),
                            start = 3, end = 4)
      )
      } %>% 
      reverseComplements(.) %>%
      mutate(guides, be_dinucleotide = .)
    
    guides %<>%
      mutate(., abe_dinucleotide = mapply(abe_dinucleotide.fn, x = full_motif, splice_site = splice_site))
    
    guides %<>%
      mutate(be_position = {ifelse(splice_site == "acceptor",
                                   (23+pam_length-guide_position),
                                   (21+pam_length-guide_position))}) %<>%
      mutate(abe_position = 20+pam_length-guide_position) %<>%
      mutate(abe_position = {ifelse(abe_position > 0, abe_position, abe_position-1)}, 
             be_position = {ifelse(abe_position > 0, be_position, be_position-1)})
    
    guides %<>% filter(., (exon > 1 | splice_site != "acceptor"), (exon < max(exon) | splice_site != "donor"))
    
    ############################################################
    #Calculate editing efficiencies
    
    data = guides %>% 
      inner_join(., be_efficiencies) %>%
      mutate(splice_site_score = {ifelse(.$splice_site == "donor", 1, 0)}) %>%
      mutate(., 
             total_be_score = 100*((
               (be_motif_score/2)*motif_weight + 
                 be_efficiency/max(be_efficiencies$efficiency)*efficiency_weight +
                 position_score*position_weight +
                 splice_site_score*splice_site_weight)
               / total_weight
             )%>% 
               round(., digits = 4)
      ) %>%
      inner_join(., abe_efficiencies) %>%
      mutate(., 
             total_abe_score = 100*((
               abe_efficiency/max(abe_efficiencies$efficiency)*efficiency_weight +
                 position_score*position_weight)
               / (efficiency_weight + position_weight)
             )%>% 
               round(., digits = 4)
      ) %>%
      mutate(gene_id = gene$gene_id[1], id = ensembl_id, pam_class = pam_char, chromosome_name = gene$chromosome_name[1], strand = gene$strand[1])
    
    data %<>% mutate(efficiency = {ifelse(be =="BE3", be_efficiency, abe_efficiency)})
    
    data = arrange(data, desc(efficiency)) %>% 
      filter(efficiency > 0) %>%
      dplyr::select(be, efficiency, splice_site, be_efficiency, abe_efficiency, guide, exon, position_score, 
                    full_motif, be_dinucleotide, be_position,
                    abe_dinucleotide, abe_position,
                    strand, chromosome_name, gene_id, id, pam, enzyme)
    ###
    data$last_exon = ifelse(data$exon == exon_number, TRUE, FALSE)
    data$one_exon = rep(FALSE, times = length(data$exon))
    data$upstream_utr_guide = ifelse(data$exon == 1 & data$splice_site == "acceptor", TRUE, FALSE)
    data$downstream_utr_guide = ifelse(data$exon == exon_number & data$splice_site == "donor", TRUE, FALSE)
    
    
    gene$last_exon = ifelse(gene$exon == exon_number, TRUE, FALSE)
    gene$one_exon = rep(FALSE, times = length(gene$exon))
    
    test_output = data
    
    test_output %<>% mutate(abe_dinucleotide = {ifelse(be == "BE3" & splice_site == "acceptor",
                                                       NA,
                                                       abe_dinucleotide)},
                            abe_efficiency = {ifelse(be == "BE3" & splice_site == "acceptor",
                                                     NA,
                                                     abe_efficiency)},
                            abe_position = {ifelse(be == "BE3" & splice_site == "acceptor",
                                                   NA,
                                                   abe_position)})
    
    test_output %<>% mutate(be_dinucleotide = {ifelse(be == "ABE",
                                                      NA,
                                                      be_dinucleotide)},
                            be_efficiency = {ifelse(be == "ABE",
                                                    NA,
                                                    be_efficiency)},
                            be_position = {ifelse(be == "ABE",
                                                  NA,
                                                  be_position)})
    test_output %<>% mutate(efficiency = pmax(abe_efficiency, be_efficiency, na.rm = TRUE))
    
    test_output %<>%
      inner_join(., pam_efficiency) %<>%
      mutate(weighted_efficiency = pam_weight*efficiency) %<>% 
      arrange(desc(weighted_efficiency)) %<>%
      as_tibble()
    
    ensembl = paste0("https://useast.ensembl.org/",
                     species,
                     "/Gene/Summary?db=core;g=",
                     gene$gene_id[1],
                     ";r=",
                     gene$chromosome_name[1],
                     ":",
                     min(c(gene$exon_chrom_start, gene$exon_chrom_end)),
                     "-",
                     max(c(gene$exon_chrom_start, gene$exon_chrom_end))
    )
    
    
    # ensembl = "https://useast.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000139618.18;r=13:32315086-32400266"
    
  } else
  {
    data = data.frame(
      total_score = NA, 
      guide = NA, 
      exon = NA, 
      position_score = NA,
      splice_site = NA, 
      motif = NA, 
      efficiency = NA, 
      dinucleotide = NA,
      guide_start_coord = NA, 
      guide_end_coord = NA, 
      guide_position = NA,
      strand = NA,
      chromosome_name = NA,
      gene_id = NA,
      last_exon = NA, 
      one_exon = TRUE, 
      upstream_utr_guide = NA, 
      downstream_utr_guide = NA, 
      id = ensembl_id,
      pam = NA,
      pam_class = pam_char);
    gene$last_exon = ifelse(gene$exon == exon_number, TRUE, FALSE);
    gene$one_exon = TRUE;
    gene$id = ensembl_id
  }
  
  return(list(test_output, ensembl))
  
}


