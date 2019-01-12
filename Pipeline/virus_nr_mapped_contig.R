options(stringsAsFactors=F)
argv <- commandArgs(TRUE)
setwd(argv[1])
library(taxonomizr)

load(argv[2]) #taxaNodes,taxaNames
contig_eval_info <- read.table(argv[3],sep = "\t",header = T,quote = NULL) 
row.names(contig_eval_info) <- contig_eval_info[ ,1]

mapped_contig <- read.table(argv[4],header = F,sep = "\t",quote = NULL)
mapped_contig <- mapped_contig[!duplicated(mapped_contig[ ,1]), ]
row.names(mapped_contig) <- mapped_contig[ ,1]
mapped_contig <- cbind(mapped_contig, contig_eval_info[mapped_contig[ ,1] ,2:4])
colnames(mapped_contig) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "slen", "stitle","qcovs","contig_length","covered_percent","reads_number")

mapped_contig$taxaid <- accessionToTaxa(mapped_contig$sseqid,argv[5])
taxonomy <- getTaxonomy(mapped_contig$taxaid,taxaNodes,taxaNames,desiredTaxa = c("superkingdom","phylum", "class", "order", "family", "genus", "species"))
mapped_contig <- cbind(mapped_contig,taxonomy)
mapped_contig[is.na(mapped_contig)] <- "Unassigned"
mapped_contig[mapped_contig=="NA"] <- "Unassigned"
write.table(mapped_contig,file = argv[6],sep = "\t",quote = F,row.names = F)

vi_mapped_contig <- mapped_contig[mapped_contig$superkingdom == "Viruses", ]
write.table(vi_mapped_contig,file = argv[7],sep = "\t",row.names = F,quote = F)

write.table(vi_mapped_contig[ ,1],file = argv[8],quote = F,row.names = F,col.names = F)

