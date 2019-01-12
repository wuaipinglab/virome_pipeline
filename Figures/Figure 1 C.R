#Figure 1:#########################################      concept pictures of comparative virome pipeline     ################################################
#C:sequencing data preliminary statistics                                                                                                                   #                              
#############################################################################################################################################################
rm(list = ls())
options(stringsAsFactors=F)
library(ggplot2)
setwd("/Intermediate_data/culicoides/species")
number <- c()
name <- c()
fileNames <- dir(pattern = "*_species.xls") 
for (i in 1:length(fileNames)) {
  cullicoides <- read.table(file = fileNames[i],header = T,quote = NULL,sep = "\t")
  name <- c(name,gsub("_species.xls","",fileNames[i]))
  number <- c(number,sum(cullicoides$Reads_number))
}
culicoides_virus_number <- cbind(name,number)
culicoides_virus_number <- as.data.frame(culicoides_virus_number)
write.table(culicoides_virus_number,file = "culicoides_virus_number.xls",sep = "\t",row.names = F)

setwd("/Intermediate_data/mosquito/species")
number <- c()
name <- c()
fileNames <- dir(pattern = "*_species.xls") 
for (i in 1:length(fileNames)) {
  mosquito <- read.table(file = fileNames[i],header = T,quote = NULL,sep = "\t")
  name <- c(name,gsub("_species.xls","",fileNames[i]))
  number <- c(number,sum(mosquito$Reads_number))
}
mosquito_virus_number <- cbind(name,number)
mosquito_virus_number <- as.data.frame(mosquito_virus_number)
write.table(mosquito_virus_number,file = "mosquito_virus_number.xls",sep = "\t",row.names = F)

setwd("/Intermediate_data")
culicoides_mosquito <- read.table("culicoides_mosquito_virus_non_host_reads.txt",header = T,sep = "\t",quote = NULL)
culicoides_virus_number$host <- culicoides_mosquito[1:61,3]
mosquito_virus_number$host <- culicoides_mosquito[1:62,7]
culicoides_mosquito <- rbind(culicoides_virus_number,mosquito_virus_number)
culicoides_mosquito$number <- as.numeric(culicoides_mosquito$number)
culicoides_mosquito$number <- log2(culicoides_mosquito$number)
culicoides_mosquito$host <- log2(culicoides_mosquito$host)
culicoides_mosquito$color <- c(rep("C",61),rep("M",62))
culicoides_mosquito$color <- factor(culicoides_mosquito$color)

pdf("Culicoides_mosquito_virus_host_number.pdf",width = 9,height = 7)
p_bar <- ggplot(culicoides_mosquito,aes(x = number,fill=color) ) +
  labs(title="",fill="") +
  xlab("log2(reads counts)") +
  ylab("Number of samples") + 
  scale_x_continuous(limits = c(17,24)) + 
  geom_histogram(binwidth = 0.2,position = "stack",bins = 50,colour = "black") + 
  theme_bw() +
  scale_fill_manual(values = c("red","blue") ) +
  theme(panel.border=element_blank(), panel.grid=element_blank(),axis.line=element_line(colour = "black") ) + 
  theme(axis.text.x=element_text(face = "bold", colour="black", size=10,angle = 90) , axis.text.y=element_text(face = "bold", colour="black", size=10)) + 
  theme(axis.title = element_text(face = "bold", colour="black", size=12)) +
  guides(fill=guide_legend(reverse=F)) + 
  theme(legend.position="bottom",legend.text = element_text(colour = 'black', size = 8, hjust = 2, vjust = 2, face = 'bold'),legend.key.size = unit(0.5,'cm'))
p_bar
dev.off() 