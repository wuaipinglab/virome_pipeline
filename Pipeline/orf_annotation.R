options(stringsAsFactors=F)
argv <- commandArgs(TRUE)
setwd(argv[1])
library(taxonomizr)
library(stringr)

nr <- read.table(argv[2],header = F,sep = "\t",quote = NULL)
colnames(nr) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "slen", "stitle","qcovs")
nr_unique <- nr[!duplicated(nr[ ,1]), ]

a <- str_split(nr_unique[ ,1],":",simplify = T)
b <- a[ ,1]
a <- as.data.frame(a)
a[,2] <- as.numeric(a[,2])
a[,3] <- as.numeric(a[,3])
colnames(a) <- c("id","start","stop")
strand <- rep("+",dim(a)[1])
strand[a[,2] > a[,3]] <- "-"
length_nt <- abs(a[,2] - a[,3]) + 1
length_aa <- (length_nt / 3) - 1
a <- cbind(a[,2],a[,3],strand,length_nt,length_aa)
colnames(a) <- c("start","stop","strand","length_nt","length_aa")
a <- as.data.frame(a)


b <- cbind(str_split(b,"_",simplify = T,n = 2)[ ,2], str_split(str_split(b,"_",simplify = T,n = 2)[ ,1],"\\|",simplify = T)[ ,2])
b <- as.data.frame(b)
colnames(b) <- c("contig_id","orf_label")

accession <- str_split(nr_unique[,2],"\\|",simplify = T)[,4]

nr_unique <- cbind(b,a,accession,nr_unique)
nr_unique$taxaid <- accessionToTaxa(nr_unique$accession,argv[3])
load(argv[4]) #taxaNodes,taxaNames
taxonomy <- getTaxonomy(nr_unique$taxaid,taxaNodes,taxaNames,desiredTaxa = c("superkingdom","phylum", "class", "order", "family", "genus", "species"))
nr_unique <- cbind(nr_unique,taxonomy)
nr_unique[is.na(nr_unique)] <- "Unassigned"
nr_unique[nr_unique=="NA"] <- "Unassigned"
nr_unique$start <- as.numeric(nr_unique$start)
nr_unique$stop <- as.numeric(nr_unique$stop)
nr_unique$length_nt <- as.numeric(nr_unique$length_nt)
nr_unique$length_aa <- as.numeric(nr_unique$length_aa)
write.table(nr_unique,file = argv[5],sep = "\t",quote = F,row.names = F)

new_matrix <- nr_unique[nr_unique$superkingdom == "Viruses" & nr_unique$length_aa >= 80, ]
tmp <- new_matrix[new_matrix$pident <= 40 , ]
rep_id <- new_matrix[duplicated(new_matrix$contig_id) , 1] 
for(i in rep_id){
	  x <- new_matrix[new_matrix$contig_id == i, ]
  if(length(unique(x)) >= 2){
	      tmp <- rbind(tmp,x)
    }
}
final_matrix <- unique(tmp)
write.table(final_matrix,file = argv[6],sep = "\t",row.names = F,quote = F)

