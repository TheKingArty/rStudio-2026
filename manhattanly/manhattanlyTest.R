# testing manhattanly funcitons with HapMap

library(manhattanly)
manhattanly(HapMap)

set.seed(12345)
HapMap.subset <- subset(HapMap, CHR %in% 4:7)
# for highlighting SNPs of interest
significantSNP <- sample(HapMap.subset$SNP, 20)
head(HapMap.subset)

dim(HapMap.subset)

manhattanly(HapMap.subset, snp = "SNP", gene = "GENE")
