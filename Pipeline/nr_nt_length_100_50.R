argv <- commandArgs(TRUE)
setwd(argv[1])
options(stringsAsFactors=F)

nr <- read.table(argv[2],header = F,quote = NULL,sep = "\t")
nt <- read.table(argv[3],header = F,quote = NULL,sep = "\t")
colnames(nr) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "slen", "stitle","qcovs")
colnames(nt) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "slen", "stitle","qcovs")
nr_per <- nr[nr$length >= 50 , ]
nt_per <- nt[nt$length >= 100 , ]
write.table(nr_per,file = argv[4],sep = "\t",quote = F,row.names = F,col.names = F)
write.table(nt_per,file = argv[5],sep = "\t",quote = F,row.names = F,col.names = F)
