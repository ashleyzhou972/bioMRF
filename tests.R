chrm='1'
rnaseq_IMR = read.table('~/hic/IMR90/rnaseq/rnaseq_RIW.txt', header = T)
names_IMR = as.character(rnaseq_IMR[rnaseq_IMR$chromosome==chrm,1])

rnaseq_GM = read.table('~/hic/rao2014/GM12878_10kb/rnaseq/rnaseq_ENCFF680ZFZ.txt', header = T)
names_GM = as.character(rnaseq_GM[rnaseq_GM$chromosome==chrm,1])

sum(!names_IMR%in%names_GM)