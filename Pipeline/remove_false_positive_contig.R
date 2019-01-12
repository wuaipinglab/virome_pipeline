options(stringsAsFactors=F)
argv <- commandArgs(TRUE)
setwd(argv[1])
library(taxonomizr)
library(RColorBrewer)
library(ggplot2)
library(stringr)
nr <- read.table(argv[2],header = F,sep = "\t",quote = NULL)
colnames(nr) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "slen", "stitle","qcovs")
nr_unique <- nr[!duplicated(nr[ ,1]), ]

load(argv[3])
accession <- as.data.frame(str_split(nr_unique[ ,2],"\\|",simplify = T))[ ,4]
nr_unique$taxaid <- accessionToTaxa(accession,argv[4])
taxonomy <- getTaxonomy(nr_unique$taxaid,taxaNodes,taxaNames,desiredTaxa = c("superkingdom","phylum", "class", "order", "family", "genus", "species"))
nr_unique <- cbind(nr_unique,taxonomy)
nr_unique[is.na(nr_unique)] <- "Unassigned"
nr_unique[nr_unique=="NA"] <- "Unassigned"
write.table(nr_unique,file = argv[5],sep = "\t",row.names = F,quote = F)

#fp contig
fp_contig <- nr_unique[nr_unique$qcovs >= 10 & nr_unique$pident >= 60 & nr_unique$superkingdom != "Viruses" , ]
write.table(fp_contig,file = argv[6],sep = "\t",row.names = F,quote = F)
fp_contig_id <- fp_contig[ ,1]
write.table(fp_contig_id,file = argv[7],sep = "\t",row.names = F,col.names = F,quote = F)

#all virus-nr blast mapped contig
vi_mapped_contig <- read.table(argv[8],header = T,sep = "\t",quote = NULL)
row.names(vi_mapped_contig) <- vi_mapped_contig[ ,1]
mv_fp_contig <- vi_mapped_contig[setdiff(rownames(vi_mapped_contig),fp_contig_id), ]
write.table(mv_fp_contig,file = argv[9],sep = "\t",row.names = F,quote = F)

write.table(mv_fp_contig[ ,1],file = argv[10],quote = F,row.names = F,col.names = F)

