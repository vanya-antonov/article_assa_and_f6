
source('lib.R')

# BiocManager::install("Biostrings")
library(Biostrings)

###

ali <- readDNAMultipleAlignment(paste0(DATA_DIR, 'chaserr_alignment.txt'), format="clustal")
class(ali)
?DNAMultipleAlignment
seqs <- readDNAStringSet(paste0(DATA_DIR, 'chaserr_seqs.txt'),format="fasta")
seqs

palign1 <- pairwiseAlignment(seqs[1], seqs[2])
palign1


?pid


?readDNAStringSet

s1 <- DNAString("AGTATAGATGATAGAT")
s2 <- DNAString("AGTAGATAGATGGATGATAGATA")

palign1 <- pairwiseAlignment(s1, s2)
class(palign1)
palign1
pid(palign1)


