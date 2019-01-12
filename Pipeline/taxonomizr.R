argv <- commandArgs(TRUE)
setwd(argv[1])
options(stringsAsFactors=F)
library(taxonomizr)
library(xlsx)
load(argv[2]) #taxaNodes,taxaNames

mapped_contig <- read.table(argv[3],header = F,sep = "\t",quote = NULL)
row.names(mapped_contig) <- 1 : dim(mapped_contig)[1]
colnames(mapped_contig) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "slen", "stitle","qcovs")
mapped_contig$taxaid <- accessionToTaxa(mapped_contig$sseqid,argv[4])
taxonomy <- getTaxonomy(mapped_contig$taxaid,taxaNodes,taxaNames,desiredTaxa = c("superkingdom","phylum", "class", "order", "family", "genus", "species"))
mapped_contig <- cbind(mapped_contig,taxonomy)
mapped_contig[is.na(mapped_contig)] <- "Unassigned"
mapped_contig[mapped_contig=="NA"] <- "Unassigned"

cdd <- read.table(argv[5],header = F,sep = "\t",quote = NULL)
colnames(cdd) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "slen", "stitle","qcovs")
orf <- read.table(argv[6],header = T,sep = "\t",quote = NULL)
write.table(mapped_contig,file = argv[7],sep = "\t",quote = F,row.names = F)
write.table(cdd,file = argv[8],sep = "\t",quote = F,row.names = F)
write.table(orf,file = argv[9],sep = "\t",quote = F,row.names = F)

write.xlsx(mapped_contig,file = argv[9],sheetName = "nr_nt",row.names = F)
write.xlsx(cdd,file = argv[10],sheetName = "cdd",append = T,row.names = F)
write.xlsx(orf,file = argv[11],sheetName = "orf",append = T,row.names = F)
