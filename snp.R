#This is an attempt to implement the basic idea of k-mer generator for snp mutation in R

library(VariantAnnotation)
library(tidyverse)
library(BSgenome.Hsapiens.UCSC.hg19)

fl <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
vcf <- readVcf(fl, "hg19")
vcf
header(vcf)

head(rowRanges(vcf), 3)
df_vcf <- as.data.frame(rowRanges(vcf))
df_vcf_sample <-head(df_vcf) %>%
  rownames_to_column("snp_id")%>%
  mutate(chromosome = paste0("chr", seqnames),
         chr.ref_pos.snp_id = paste(chromosome, start, snp_id, sep = ":"))
row_rep <- function(df, n) {
  df[rep(1:nrow(df), each = n), ]
}

vcf_size <- nrow(df_vcf_sample)
kmer_relative_position <- c(0, -1, -2, -3, -4)
relative_ranges <- rep(kmer_relative_position, vcf_size)
df_vcf_sample_kmers <- row_rep(df_vcf_sample,5) %>%
    mutate(kmer_start = start + relative_ranges,
           kmer_end = end + relative_ranges)

kmer_info <- getSeq(Hsapiens, df_vcf_sample_kmers$chromosome,
                   start = df_vcf_sample_kmers$kmer_start,
                   width = rep(5, nrow(df_vcf_sample_kmers))
                   )
kmer_width <- width(kmer_info)
kmer_ref_seq <- as.character(kmer_info)
kmer_chr <- names(kmer_info)
kmer_start <- df_vcf_sample_kmers$start
kmer_end <- df_vcf_sample_kmers$end
kmer_mut_ref <- df_vcf_sample_kmers$REF
kmer_mut_alt <- sapply(df_vcf_sample_kmers$ALT, as.character)
kmer_mut_snpid <- df_vcf_sample_kmers$snp_id
kmer_table <- data.frame(kmer_chr,
                         kmer_mut_snpid,
                         kmer_start,
                         kmer_end,
#                         kmer_width,
                         kmer_mut_ref,
                         kmer_mut_alt,
                         kmer_ref_seq)

########
vcf_size <- nrow(df_vcf_sample)
kmer_relative_position <- c(0, -1, -2, -3, -4)
relative_ranges <- rep(kmer_relative_position, vcf_size)
df_vcf_sample_kmers <- row_rep(df_vcf_sample,5) %>%
  mutate(kmer_start = start + relative_ranges,
         kmer_end = end + relative_ranges)

kmer_info <- getSeq(Hsapiens, df_vcf_sample_kmers$chromosome,
                    start = df_vcf_sample_kmers$kmer_start,
                    width = rep(5, nrow(df_vcf_sample_kmers))
)
kmer_width <- width(kmer_info)
kmer_ref_seq <- as.character(kmer_info)
kmer_chr <- names(kmer_info)
kmer_start <- df_vcf_sample_kmers$start
kmer_end <- df_vcf_sample_kmers$end
kmer_mut_ref <- df_vcf_sample_kmers$REF
kmer_mut_alt <- sapply(df_vcf_sample_kmers$ALT, as.character)
kmer_mut_snpid <- df_vcf_sample_kmers$snp_id
kmer_table <- data.frame(kmer_chr,
                         kmer_mut_snpid,
                         kmer_start,
                         kmer_end,
                         #kmer_width,
                         kmer_mut_ref,
                         kmer_mut_alt,
                         kmer_ref_seq)

kmer_list <- split(kmer_table, factor(kmer_table$kmer_mut_snpid))
