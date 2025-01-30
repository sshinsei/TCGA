setwd("E:/LZ/24080/01_TCGA/files")  #modify work path

rt <- read.table("mRNAmatrix.txt",header = T,check.names = F,sep = "\t")

lncRNA <- rt[grep("lncRNA", rt$id), ]
mRNA <- rt[grep("protein_coding", rt$id), ]
miRNA <- rt[grep("miRNA", rt$id), ]

#strsplit(lncRNA$id, split = "\\|")

library(tidyr)
lncRNA <- separate(lncRNA, id, into = c("id","Encode","Type"), sep = "\\|")[ , -c(2, 3)]
mRNA <- separate(mRNA, id, into = c("id", "Encode", "Type"), sep = "\\|")[ , -c(2, 3)]
miRNA <- separate(miRNA, id, into = c("id", "Encode", "Type"), sep = "\\|")[ , -c(2, 3)]

write.table(mRNA, "mRNA_count.txt", row.names = F, col.names = T, quote = F, sep = "\t")
write.table(lncRNA, "lncRNA_count.txt", row.names = F, col.names = T, quote = F, sep = "\t")
write.table(miRNA, "miRNA_count.txt", row.names = F, col.names = T, quote = F, sep = "\t")
