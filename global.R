# CBE variant
cbe = "BE4"
# ABE variant
abe = "ABE7.10"

source("helpers.R")

##### Main function
runSpliceR = function(
  logistic_adjust = 1,
  ensembl_transcript_id,
  species,
  pam,
  guide_length = 20,
  flank_length = 20,
  downstream_gene_seq = 20
){

##### Operations
# Load data
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

# Establish PAM information
pam.char = pam
pam.DNAString = pam %>% DNAString()
pam_length = nchar(pam.char)

# Establish guide information
guide.DNAString = rep("N", guide_length) %>%
  paste0(., collapse = "") %>%
  paste0(., pam.char) %>%
  DNAString()

# Grab gene information from Ensembl via GraphQL
gene_coordinates = getTranscriptExonCoordinates(ensembl_transcript_id, species = species)

# Establish the sense of the gene
sense = gene_coordinates$strand[1]
ensembl_gene_id = gene_coordinates$gene_id[1]

# Generate URL for loading the iframe (Updated Ensembl URL target)
ensembl = paste0("https://beta.ensembl.org/genome-browser/gene?id=",
                 ensembl_gene_id
)

if(sense == 1){ 
  # Positive sense strand
  gene_coordinates %<>% arrange(., exon_chrom_start)
  exon_number = nrow(gene_coordinates)
  
  gene_sequence_length = gene_coordinates %$%
    c(.$exon_chrom_start, .$exon_chrom_end) %>%
    range() %>%
    {.[2] - .[1]} %>%
    {. + 1 + flank_length + pam_length}
  
  gene_coordinates %<>%
    mutate(exon_start_index = exon_chrom_start - min(exon_chrom_start) + 1 + flank_length + pam_length) %>%
    mutate(exon_end_index = exon_chrom_end - min(exon_chrom_start) + 1 + flank_length + pam_length)
  
  gene_seq_chrom_start = min(gene_coordinates$exon_chrom_start) - (flank_length + pam_length)
  
  gene_coordinates %<>%
    mutate(sa_seq_coordinates = exon_chrom_start - (flank_length + pam_length)) %<>%
    mutate(sd_seq_coordinates = exon_chrom_end - (flank_length + pam_length))
  
  gene_coordinates %<>% mutate(exon_length = exon_end_index - exon_start_index) %<>%
    mutate(sd_position = 1-cumsum(exon_length)/sum(exon_length)) %<>%
    mutate(sa_position = 1-(cumsum(exon_length)-exon_length)/sum(exon_length))
  
  gene.DNAString = coordinatesToDNAString(start = gene_coordinates$exon_chrom_start[1], 
                                          end = gene_coordinates$exon_chrom_end[exon_number],
                                          strand = gene_coordinates$strand[1],
                                          chromosome = gene_coordinates$chromosome_name[1],
                                          upstream = flank_length + pam_length,
                                          downstream = downstream_gene_seq + pam_length,
                                          species = species
  )
} else {
  # Negative sense strand
  gene_coordinates %<>% arrange(., desc(exon_chrom_start))
  exon_number = nrow(gene_coordinates)
  
  gene_sequence_length = gene_coordinates %$%
    c(.$exon_chrom_start, .$exon_chrom_end) %>%
    range() %>%
    {.[2] - .[1]} %>%
    {. + 1 + flank_length + pam_length}
  
  gene_coordinates %<>%
    mutate(exon_start_index = (max(exon_chrom_end) - exon_chrom_end) + 1 + flank_length + pam_length) %>%
    mutate(exon_end_index = (max(exon_chrom_end) - exon_chrom_start) + 1 + flank_length + pam_length)
  
  gene_seq_chrom_start = max(gene_coordinates$exon_chrom_end) + (flank_length + pam_length)
  
  gene_coordinates %<>%
    mutate(sa_seq_coordinates = exon_chrom_end + (flank_length + pam_length)) %<>%
    mutate(sd_seq_coordinates = exon_chrom_start + (flank_length + pam_length))
  
  gene_coordinates %<>% mutate(exon_length = exon_end_index - exon_start_index) %<>%
    mutate(sd_position = cumsum(exon_length)/sum(exon_length)) %<>%
    mutate(sa_position = (cumsum(exon_length)-exon_length)/sum(exon_length))
  
  gene.DNAString = coordinatesToDNAString(start = gene_coordinates$exon_chrom_end[1], 
                                          end = gene_coordinates$exon_chrom_start[exon_number],
                                          strand = gene_coordinates$strand[1],
                                          chromosome = gene_coordinates$chromosome_name[1],
                                          upstream = flank_length + pam_length,
                                          downstream = downstream_gene_seq + pam_length,
                                          species = species
  )
}

gene.char = as.character(gene.DNAString)

gene_coordinates = gene_coordinates %>%
  rowid_to_column("exon") %>%
  mutate(intron_start_index = exon_end_index + 1) %>%
  mutate(intron_end_index = if(exon_number > 1) {
    c(
      tail(exon_start_index, -1) - 1,
      length(gene.DNAString)
    )
  } else {
    length(gene.DNAString)
  })

gene_coordinates = gene_coordinates %>%
  mutate(sa_seq = mapply(FUN = substr,
                         x = gene.char,
                         start = exon_start_index - flank_length - pam_length,
                         stop = exon_start_index + flank_length + pam_length)) %>%
  mutate(sd_seq = mapply(FUN = substr,
                         x = gene.char,
                         start = exon_end_index - flank_length - pam_length,
                         stop = exon_end_index + flank_length + pam_length)) %>%
  mutate(acceptor = mapply(FUN = substr,
                           x = sa_seq,
                           start = flank_length - 4 + pam_length,
                           stop = flank_length  + 5 + pam_length)) %>%
  mutate(donor = mapply(FUN = substr,
                        x = sd_seq,
                        start = flank_length  + pam_length - 3,
                        stop = flank_length  + 6 + pam_length))

### Identify guides
cbe_splice_acceptors = gene_coordinates %>%
  .$sa_seq %>%
  DNAStringSet %>%
  reverseComplement %>%
  matchPatterns(pattern = guide.DNAString, subject = ., fixed = FALSE)

abe_splice_acceptors = gene_coordinates %>%
  .$sa_seq %>%
  DNAStringSet %>%
  matchPatterns(pattern = guide.DNAString, subject = ., fixed = FALSE)

splice_donors = gene_coordinates %>%
  .$sd_seq %>%
  DNAStringSet %>%
  reverseComplement %>%
  matchPatterns(pattern = guide.DNAString, subject = ., fixed = FALSE)

unscored_cbe_splice_acceptors = cbe_splice_acceptors %>%
  extractGuide(., 1:length(.)) %>%
  unlist() %>%
  as.vector()

unscored_abe_splice_acceptors = abe_splice_acceptors %>%
  extractGuide(., 1:length(.)) %>%
  unlist() %>%
  as.vector()

unscored_splice_donors = splice_donors %>%
  extractGuide(., 1:length(.)) %>%
  unlist() %>%
  as.vector()

cbe_splice_acceptors_positions = cbe_splice_acceptors %>% extractGuideStart(., exon = 1:length(.)) %>% unlist() %>% as.vector
abe_splice_acceptors_positions = abe_splice_acceptors %>% extractGuideStart(., exon = 1:length(.)) %>% unlist() %>% as.vector
splice_donors_positions = splice_donors %>% extractGuideStart(., exon = 1:length(.)) %>% unlist() %>% as.vector

# Build guide data frames
cbe_splice_acceptor_guides = tibble(guide = unscored_cbe_splice_acceptors,
                                    guide_position = cbe_splice_acceptors_positions,
                                    splice_site = "acceptor",
                                    enzyme = "CBE",
                                    exon = sapply(cbe_splice_acceptors, length) %>% rep(c(1:exon_number), times = .))

abe_splice_acceptor_guides = tibble(guide = unscored_abe_splice_acceptors,
                                    guide_position = abe_splice_acceptors_positions,
                                    splice_site = "acceptor",
                                    enzyme = "ABE",
                                    exon = sapply(abe_splice_acceptors, length) %>% rep(c(1:exon_number), times = .))

cbe_and_abe_splice_acceptor_guides = tibble(guide = unscored_splice_donors,
                                            guide_position = splice_donors_positions,
                                            splice_site = "donor",
                                            enzyme = "CBE and ABE",
                                            exon = sapply(splice_donors, length) %>% rep(c(1:exon_number), times = .))

guides = bind_rows(cbe_splice_acceptor_guides,
                   abe_splice_acceptor_guides,
                   cbe_and_abe_splice_acceptor_guides)

if(nrow(guides) > 0) { 
  guides = gene_coordinates %>% 
    dplyr::select(exon, acceptor, donor) %>% 
    tidyr::gather(exon) %>%
    mutate(splice_site = exon, exon = rep(1:exon_number, times = 2)) %>%
    set_names(c("exon","splice_site","motif")) %>%
    bind_cols(.,
              (gene_coordinates %>% 
                 dplyr::select(exon, sa_position, sd_position) %>% 
                 tidyr::gather(exon) %>% 
                 mutate(exon = {ifelse(exon == "sa_position", "acceptor", "donor")}) %>%
                 mutate(splice_site = exon, exon = rep(1:exon_number, times = 2)) %>%
                 set_names(c("exon","disruption_position", "splice_site")) %>%
                 dplyr::select(disruption_position)
              )
    ) %>%
    set_names(c("exon","motif","splice_site","disruption_position")) %>%
    inner_join(guides, .) %>%
    filter(!(exon == 1 & splice_site == "acceptor")) %>%
    filter(!(exon == exon_number & splice_site == "donor")) %>%
    mutate(pam = stringr::str_sub(guide, -pam_length, -1)) %>%
    mutate(protospacer = stringr::str_sub(guide, 1, guide_length)) %>%
    mutate(cbe_position = {ifelse(splice_site == "acceptor",
                                  (guide_length + 3 + pam_length - guide_position),
                                  (guide_length + 1 + pam_length - guide_position)
    )}) %>%
    mutate(abe_position = guide_length + pam_length - guide_position) %>%
    mutate(abe_position = {ifelse(abe_position > 0, abe_position, abe_position-1)}, 
           cbe_position = {ifelse(abe_position > 0, cbe_position, cbe_position-1)}) %>%
    mutate(target_motif = {ifelse(enzyme == "ABE" & splice_site == "acceptor", motif, revcom(motif))}) %>%
    mutate(cbe_motif = {ifelse(splice_site ==  "donor",
                               substr(target_motif, 2, 8),
                               substr(target_motif, 3, 9))},
           abe_motif = substr(target_motif, 1, 7)
    ) %>%
    inner_join(., gene_coordinates %>% dplyr::select(exon, sa_seq, sd_seq)) %>%
    mutate(sa_seq = {ifelse(enzyme == "ABE" & splice_site == "acceptor", sa_seq, revcom(sa_seq))}) %>%
    mutate(sd_seq = revcom(sd_seq)) %>%
    
    # Swapped legacy plyr::join with dplyr::left_join to maintain tibble integrity
    dplyr::left_join(cbe_motif_weights, by = "cbe_motif") %>%
    dplyr::left_join(abe_motif_weights, by = "abe_motif") %>%
    dplyr::left_join(cbe_position_weights, by = "cbe_position") %>%
    dplyr::left_join(abe_position_weights, by = "abe_position") %>%
    tibble() %>%
    mutate(cbe_score = probability(cbe_motif_weight + cbe_position_weight - logistic_adjust)) %>%
    mutate(abe_score = probability(abe_motif_weight + abe_position_weight - logistic_adjust)) %>% 
    dplyr::select(exon, splice_site, protospacer, pam, enzyme,
                  disruption_position, cbe_position, cbe_score, abe_position, abe_score) %>% 
    mutate(transcript = ensembl_transcript_id, gene = ensembl_gene_id) %>%
    mutate(cbe_position = {ifelse(enzyme == "ABE", NA, cbe_position)}) %>%
    mutate(cbe_score = {ifelse(enzyme == "ABE", NA, cbe_score)}) %>%
    mutate(abe_position = {ifelse(enzyme == "CBE", NA, abe_position)}) %>%
    mutate(abe_score = {ifelse(enzyme == "CBE", NA, abe_score)}) %>%
    mutate_if(is.numeric, signif, 3) %>%
    inner_join(.,
               gene_coordinates %>%
                 dplyr::select(exon, chromosome_name, exon_chrom_start, exon_chrom_end, strand)
    ) %>%
    addProtospacerCoordinates(., guide_length = guide_length) %>%
    arrange(desc(cbe_score))
  
} else {
  guides = tibble(exon, splice_site, protospacer, pam, enzyme,
                  disruption_position, cbe_position, cbe_score, abe_position, abe_score,
                  transcript = ensembl_transcript_id, gene = ensembl_gene_id)
}

return(list(guides, gene_coordinates, ensembl))
}