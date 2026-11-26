library(Matrix)
library(pagoda2)
library(ggplot2)
library(parallel)
library(scde)
library(seqinr)



barlist <- read.csv('data/raw/barcode list.txt', sep = '\t')
tfs <- unique(barlist$Gene)

### process TF barcodes
upstream <- 'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNATCTTGTGGAAAGGACGA'
downstream <- 'ACGCGTTAAGTCgacaatcaacctctggattacaaaatttgtgaaagattgactggtattcttaactatgttgcNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN'
barcodes <- barlist$Sequence

#ACTGATCTTGTGGAAAGGACGAGCTGACTCCTTCTTGAGAATCTGCACGCGTTAAGTCTACAATCAACCTCTGGATTACAAAATTTGTGAAAGATTGACT
#ATCTTGTGGAAAGGACGA                        ACGCGTTAAGTCgacaatcaacctctggattacaaaatttgtgaaagattgactggtattcttaactatgttgc

#GCATCTTGTGGAAAGGACGAGCTGACTCCTTCTTGGGAATCTGCACGCGTTAAGTCGACAATCAACCTCTGGATTACAAAATTTGTGAAAGATTGACTGG
#ATCTTGTGGAAAGGACGA

# save fastq file with barcode sequences
faseq.list <- lapply(1:length(barcodes), function(i){
  paste(upstream, barcodes[i], downstream, sep = "")
})
barcodes1 <- barcodes
write.fasta(faseq.list, barcodes1, file.out = "data/processed/barcodes.fa")

# save gtf file with barcode coordinates
seq <- vector()
for (i in 1:nrow(barlist)){
  #nm <- barcodes1[i]
  nm <- barlist$Sequence[i]
  prop <- paste('gene_id "',nm,'"; transcript_id "',nm,'"; gene_name "',nm,'"; gene_biotype "protein_coding";', sep = '')

  seq[i] <- paste(nm, "unknown", "exon",
                  #nchar(sleft) + 1, nchar(sleft) + nchar(barcodes1[i]),
                  #nchar(sgfp) + 1, nchar(sleft) + nchar(barcodes1[i]) + len.sright.barcode,
                  1, nchar(upstream) + nchar(barcodes1[i]) + nchar(downstream),
                  '.', '+', '.', prop,
                  sep = '\t')

}
fileConn<-file("barcodes_full.gft")
writeLines(seq, fileConn)
close(fileConn)
