# testing manhattanly funcitons with HapMap

library(manhattanly)
# manhattanly(HapMap)

set.seed(12345)
# HapMap.subset <- subset(HapMap, CHR %in% 4:7)
# for highlighting SNPs of interest
# significantSNP <- sample(HapMap.subset$SNP, 20)

# df <- subset(HapMap, CHR %in% 1:3)

# head(HapMap.subset)

# dim(HapMap.subset)

# manhattanly(HapMap.subset, snp = "SNP", gene = "GENE")

testr <- read.csv("C:/Users/arthu/Documents/Research/localGWASdata/GCST90681941.h.tsv/GCST90681941.h.tsv", sep = "\t")
testr.chr1to2 <- subset(testr, chromosome %in% 1:2)
testr.chr1 <- subset(testr, chromosome %in% 1)
head(testr.chr1)
dim(testr.chr1)
manhattanly(testr.chr1, chr = "chromosome", snp = "rsid", bp = "base_pair_location", p = "p_value")