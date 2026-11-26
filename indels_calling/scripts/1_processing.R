library(GenomicRanges)
library(Rsamtools)
library(GenomicAlignments)
library(data.table)

hash.match <- readRDS(file = "../RNA_analysis/data/processed/hashtag_meta.rds")
hash.match

run <- 1
gene <- 'DNMT3A'
estimate_crispr_site_coverage('DNMT3A', 4)

estimate_crispr_site_coverage('TET2', 4)

cov_editl <- lapply(c('TET2', 'DNMT3A'), function(gene){
  y <- do.call(rbind, lapply(1:4, function(run){
    geno <- as.character(unique(hash.match[hash.match$RNA == run,]$genotype))
    geno <- as.character(hash.match[hash.match$RNA == run & hash.match$exp == ifelse(gene == 'DNMT3A', 123, 127),]$genotype)
    cov <- estimate_crispr_site_coverage(gene, run)
    if (length(unique(geno)) > 1){
      geno <- 'mix'
    }else{
      geno <- geno[1]
    }
    data.frame(cov[1], cov[2], geno)
  }))
  colnames(y) <- c('edited', 'total', 'genotype')
  y
})
names(cov_editl) <- c('TET2', 'DNMT3A')

i <- 1
cov_edit <- cov_editl[[i]]
maxsz <- 250
pl <- ggplot(cov_edit, aes(total, edited, fill = genotype)) + geom_abline(linetype = 2, size = 2) + 
  geom_point(size = 5, shape = 21, col = 'black') + 
  #lims(x = c(min(cov_edit[,1:2]), max(cov_edit[,1:2])), y = c(min(cov_edit[,1:2]), max(cov_edit[,1:2]))) + 
  lims(x = c(0, maxsz), y = c(0, maxsz)) + 
  labs(x = 'total reads', y = 'edited reads', title = paste('gene', names(cov_editl)[i])) + theme_bw() + 
  theme(text = element_text(size = 25), 
        legend.position = c(0.05, 0.95), legend.justification   = c("left", "top"), legend.background = element_rect(fill = "white", color = "gray50")) 
pl

fig_path <- paste('figures/scrna_editing_', names(cov_editl)[i], '.pdf', sep ='')
pdf(fig_path, width = 5, height = 5)
print(pl)
dev.off()

estimate_crispr_site_coverage <- function(gene, run){
  if (gene == 'DNMT3A'){
    bam_path <- paste("../../human_larry_exp123/RNA", run, "_augment/outs/possorted_genome_bam.bam", sep = "")
  }else{
    bam_path <- paste("../../human_larry_exp127/RNA", run, "_augment/outs/possorted_genome_bam.bam", sep = "")
  }

# Define complete genomic coordinates for both target loci
if (gene == 'DNMT3A'){
  targets <- GRanges(
    seqnames = c("chr2"),
    ranges = IRanges(
      start = c(25225874), 
      end   = c(25344483)
    ),
    strand = c("-"), # Both genes sit on the minus strand
    gene   = c("DNMT3A")
  )
}else{
  targets <- GRanges(
    seqnames = c("chr4"),
    ranges = IRanges(
      start = c(105146934), 
      end   = c(105276332)
    ),
    strand = c("-"), # Both genes sit on the minus strand
    gene   = c("TET2")
  )
}
names(targets) <- mcols(targets)$gene

# Quick validation check for index file visibility
if (!file.exists(paste0(bam_path, ".bai"))) {
  cat("Index file missing. Generating BAM index first...\n")
  indexBam(bam_path)
}

# Restrict data pull exclusively to target regions and the 10x UMI tag ('UB')
param <- ScanBamParam(
  which = targets,
  #tag = "UB", 
  tag   = c("UB", "CB"),
  what = c("cigar", "pos", "qname"),
  flag = scanBamFlag(isDuplicate = FALSE, isUnmappedQuery = FALSE)
)

# Extract into a structured alignment list object (grouped by read pairs)
alignments_list <- readGAlignmentsList(bam_path, param = param)

# Flatten all read fragments into a single vector instantly (bypasses list loops)
flat_alignments <- unlist(alignments_list, use.names = FALSE)

# Convert directly to a data.table using vectorized C-extractors
dt <- data.table(
  Gene     = gene,
  Chr      = as.character(seqnames(flat_alignments)),
  Strand   = as.character(strand(flat_alignments)),
  Cigar    = as.character(cigar(flat_alignments)),
  Pos      = start(flat_alignments),
  ReadName = mcols(flat_alignments)$qname,
  CellBC   = mcols(flat_alignments)$CB,
  UB       = mcols(flat_alignments)$UB
)

# Assign the UMI identifier directly
dt[, UMI := UB]
# Drop reads lacking a valid cell-ranger approved UMI tag
dt <- dt[!is.na(UMI) & UMI != ""]
dt <- dt[!is.na(CellBC) & CellBC != ""]
# Identify true exonic alignments vs background splicing/pre-mRNA traffic.
# Reads bridging introns contain an 'N' code in their CIGAR profile.
table(grepl("N", dt$Cigar))
dt[, IsSplicedBridge := grepl("N", Cigar)]

# Tag individual reads carrying genomic Insertions ('I') or Deletions ('D')
dt[, HasInsertion := grepl("I", Cigar)]
dt[, HasDeletion  := grepl("D", Cigar)]
dt[, HasIndel     := HasInsertion | HasDeletion]

# Filter down to true exon-mapping alignments for mutation summary
exonic_dt <- dt[IsSplicedBridge == FALSE]

# Calculate high-level molecular depth profiles
total_unique_umis  <- uniqueN(dt$UMI)
exonic_unique_umis <- uniqueN(exonic_dt$UMI)

# Identify unique mutated molecules (UMIs supporting an indel)
mutated_umis       <- unique(exonic_dt[HasIndel == TRUE, UMI])
total_mutated_umis <- length(mutated_umis)
capture_frequency  <- (total_mutated_umis / exonic_unique_umis) * 100
capture_frequency

if (total_mutated_umis > 0) {
  cat("\n  [+] Top Prevalent Mutation CIGAR Profiles:\n")
  mutation_profiles <- exonic_dt[HasIndel == TRUE, .(Reads = .N), by = Cigar]
  setorder(mutation_profiles, -Reads)
  print(head(mutation_profiles, 3), row.names = FALSE)
}

# ==============================================================================
# 6. EXTRACT EXACT INDEL GENOMIC POSITIONS & HOTSPOTS
# ==============================================================================
if (total_mutated_umis > 0) {
  cat("\nExtracting exact genomic positions for indels...\n")
  
  # Filter for exonic reads that specifically contain an Indel
  mut_reads <- exonic_dt[HasIndel == TRUE]
  
  # Extract explicit genomic coordinates of every CIGAR operation block
  # This returns a IRangesList matching the length of our mut_reads vector
  cigar_ranges <- cigarRangesAlongReferenceSpace(
    mut_reads$Cigar, 
    pos = mut_reads$Pos, 
    ops = c("I", "D") # Only pull out Insertions and Deletions
  )
  
  # Convert the ranges list into a flat, index-matched data.table
  indel_intervals <- as.data.table(cigar_ranges)
  
  # Clean up column names from the range coercion
  # 'group' matches the row index of our mut_reads data table
  setnames(indel_intervals, old = c("group", "start", "end", "width"), new = c("ReadIdx", "IndelStart", "IndelEnd", "IndelWidth"))
  
  # Map the original UMI and CIGAR string back onto each individual indel event
  indel_intervals[, UMI := mut_reads$UMI[ReadIdx]]
  indel_intervals[, OriginalCigar := mut_reads$Cigar[ReadIdx]]
  
  # Determine if the event is an Insertion or Deletion by inspecting the Cigar chunk
  # (Deletions take up reference space width > 0; Insertions have a reference space width of 0)
  indel_intervals[, MutationType := fifelse(IndelWidth > 0, "Deletion", "Insertion")]
  
  # Aggregate the events by exact genomic coordinate to find true molecular hotspots
  # We count raw reads AND collapse by unique UMIs to get clean biology
  position_summary <- indel_intervals[, .(
    RawReadCount = .N,
    UniqueUMIs   = uniqueN(UMI),
    IndelLength  = max(IndelWidth) # 0 for insertions, or N bases deleted
  ), by = .(IndelStart, MutationType)]
  
  # Sort positions by the highest number of unique biological molecules (UMIs)
  setorder(position_summary, -UniqueUMIs)
  
  # Print the Top Hotspot Coordinates
  cat("\n==================================================\n")
  cat(sprintf("          TOP MUTATIONAL HOTSPOTS: %s             \n", gene))
  cat("==================================================\n")
  
  # Dynamically fetch the top 10 most mutated positions
  top_hotspots <- head(position_summary, 10)
  
  for (r in 1:nrow(top_hotspots)) {
    chrom <- unique(mut_reads$Chr)[1]
    pos_type <- top_hotspots$MutationType[r]
    pos_coord <- top_hotspots$IndelStart[r]
    u_count <- top_hotspots$UniqueUMIs[r]
    r_count <- top_hotspots$RawReadCount[r]
    i_len <- top_hotspots$IndelLength[r]
    
    if (pos_type == "Deletion") {
      cat(sprintf("  [!] %s:%d | %s | Size: -%dbp | Supported by %d UMIs (%d raw reads)\n",
                  chrom, pos_coord, pos_type, i_len, u_count, r_count))
    } else {
      cat(sprintf("  [!] %s:%d | %s | Insertion Site | Supported by %d UMIs (%d raw reads)\n",
                  chrom, pos_coord, pos_type, u_count, r_count))
    }
  }
  cat("==================================================\n")
  
  # Optional: Save the full position table to a CSV file for plotting
  # fwrite(position_summary, sprintf("%s_exact_indel_positions.csv", gene))
}

# ==============================================================================
# 7. CALCULATE SITE-SPECIFIC CRISPR EDITING EFFICIENCY
# ==============================================================================
cat("\nCalculating exact target-site editing efficiency...\n")

# Define your exact CRISPR target window (or single anchor base)
# For a specific double-strand break site, we can use the exact cut coordinate:
# DNMT3A chr2:25247596-25247612 (Strand: +)
# TET2: chr4:105234023-105234039
if (gene == 'TET2'){
  target_chr   <- "chr4"
  target_start <-  105234020 
  target_end   <- 105234039 
}else if (gene == 'DNMT3A'){
  target_chr   <- "chr2"
  target_start <- 25247596
  target_end   <- 25247612
}

###
exonic_dt[, End := Pos + cigarWidthAlongReferenceSpace(Cigar) - 1]
robust_reads <- exonic_dt[Pos <= target_start & End >= target_end]
total_robust_reads <- nrow(robust_reads)
total_robust_umis  <- uniqueN(robust_reads$UMI) # Handles dual-barcode unique tracking
edited_robust_umis  <- uniqueN(robust_reads[HasIndel == TRUE, .(CellBC, UMI)])
cat("\n==================================================\n")
cat("          ROBUST EDGE-CLEANED METRICS             \n")
cat("==================================================\n")
cat(sprintf("  [+] Cleaned Reads spanning site       : %d\n", total_robust_reads))
cat(sprintf("  [+] Cleaned Molecules (UMIs) at site  : %d\n", total_robust_umis))
cat(sprintf("  [+] Edited Molecules (UMIs) at site   : %d\n", edited_robust_umis))
cat("  ------------------------------------------------\n")
cat(sprintf("  [=] True Edge-Corrected Frequency     : %.3f%%\n", 
            (edited_robust_umis / total_robust_umis) * 100))
cat("==================================================\n")
###
return(c(edited = edited_robust_umis, total = total_robust_umis))
}