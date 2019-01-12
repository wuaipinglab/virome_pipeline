options(stringsAsFactors=F)
argv <- commandArgs(TRUE)
setwd(argv[1])

m <- read.table(argv[2],header = T,sep = "\t",quote = NULL)
family_uni <- unique(m$family)

for (i in 1:length(family_uni)) {
  write.table(m[m$family==family_uni[i] ,1],file = paste(family_uni[i],"_family_id.txt",sep = ""),row.names = F,col.names = F,quote = F)
}


