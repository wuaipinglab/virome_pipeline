argv <- commandArgs(TRUE)
setwd(argv[1])
options(stringsAsFactors=F)
cov <- read.csv(argv[2],sep = "\t",quote = NULL)
cov[,1] <- gsub("\\s\\S+","",cov[ ,1])
cov[ ,7] <- as.numeric(cov[ ,7])
cov[ ,8] <- as.numeric(cov[ ,8])
Reads_number <- as.numeric(apply(cov[ ,7:8],1,sum))
contig_eval_info <- cbind(cov[ ,c(1,3,5)],Reads_number)
row.names(contig_eval_info) <- contig_eval_info[ ,1]
write.table(contig_eval_info,file = argv[3],sep = "\t",quote = F,row.names = F,col.names = T)

