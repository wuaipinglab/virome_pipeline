#Figure 4:#########################    Similarities and differences of viral groups in different space and time(jiangcheng)     #############################
#B:space difference                                                                                                                                         #                              
#############################################################################################################################################################
rm(list = ls())
options(stringsAsFactors=F)
library(plyr)
library(ggplot2)
library(stringr)
library(RColorBrewer)
#culicoides(pie)
setwd("/Intermediate_data/culicoides/strand")
load("culicoides_strand.RData")
culicoides_data <- culicoides_data[ , c("Genome_Composition_Family","position_name_out","Sample","Reads_number")]
a <- paste(culicoides_data[,1],culicoides_data[,2],sep = "_")
a1 <- cbind(a,culicoides_data$Reads_number)
a1 <- as.data.frame(a1)
colnames(a1) <- c("name","reads")
a1$name <- as.character(a1$name)
a1$reads <- as.character(a1$reads)
a1$reads <- as.numeric(a1$reads)
m <- aggregate(a1$reads,list(a1$name),sum)

b <- as.data.frame(str_split(m[,1],"_",simplify = T,n=2))
culicoides_data_new <- cbind(b,m[,2])
colnames(culicoides_data_new) <- c("Type","Position","Reads_number")
write.table(culicoides_data_new,file = "culicoides_position_strand.xls",sep = "\t",row.names = F,quote = F)

rm(list = ls())
options(stringsAsFactors=F)
setwd("E:/Intermediate_data/culicoides/strand")
culicoides_data_new <- read.table(file = "culicoides_position_strand.xls",sep = "\t",header = T,quote = NULL)
culicoides_data_new <- culicoides_data_new[culicoides_data_new$Type != "Unassigned" ,]
unique(culicoides_data_new$Position)
culicoides_data_new[culicoides_data_new$Position == "JiangCheng2013" ,2] <- "JiangCheng"
culicoides_data_new[culicoides_data_new$Position == "JiangCheng2014" ,2] <- "JiangCheng"
culicoides_data_new[culicoides_data_new$Position == "JiangCheng2015" ,2] <- "JiangCheng"
culicoides_data_new[culicoides_data_new$Position == "JiangCheng2016" ,2] <- "JiangCheng"

culicoides_data_new[culicoides_data_new$Position == "JiangCheng2017" ,2] <- "JiangCheng"
aaa <- data.frame(Type=rep("ssDNA(+)",length(unique(culicoides_data_new$Position))) , Position=unique(culicoides_data_new$Position) , Reads_number=rep(0,length(unique(culicoides_data_new$Position))))
culicoides_data_new <- join(culicoides_data_new,aaa,type="full")

temp <- data.frame(Type = rep(sort(unique(culicoides_data_new$Type)),10),Position = sort(rep(unique(culicoides_data_new$Position),9)),Reads_number = rep(0,90)) 
culicoides_data_new <- rbind(temp,culicoides_data_new)
culicoides_data_new <- culicoides_data_new[order(culicoides_data_new$Type), ]
a <- paste(culicoides_data_new[,1],culicoides_data_new[,2],sep = "_")
a1 <- cbind(a,culicoides_data_new$Reads_number)
a1 <- as.data.frame(a1)
colnames(a1) <- c("name","reads")
a1$name <- as.character(a1$name)
a1$reads <- as.character(a1$reads)
a1$reads <- as.numeric(a1$reads)
m <- aggregate(a1$reads,list(a1$name),sum)
b <- as.data.frame(str_split(m[,1],"_",simplify = T,n=2))
culicoides_data_new <- cbind(b,m[,2])
colnames(culicoides_data_new) <- c("Type","Position","Reads_number")

Set1 <- colorRampPalette(brewer.pal(9,"Set1"))(9)
col <- c(Set1[3:8],Set1[1:2],Set1[9]) 

for (i in unique(culicoides_data_new$Position)) {
  ggplot(culicoides_data_new[culicoides_data_new$Position == i, ], aes(x ="", y = Reads_number, fill = Type )) + 
    labs(title=paste("culicoides_",i,sep = ""),x = "", y = "",fill="Type") +
    geom_bar(stat = "identity",width = 1) + 
    coord_polar(theta = "y") + 
    theme(axis.ticks = element_blank(),axis.text = element_blank()) +
    theme(panel.grid=element_blank(),panel.border=element_blank()) +
    theme(panel.background = element_rect(fill = "transparent",colour = NA),plot.background = element_rect(fill = "transparent",colour = NA)) +
    guides(fill=guide_legend(reverse=F)) +
    scale_fill_manual(values = col ) +
    theme(legend.position="bottom",legend.text = element_text(colour = 'black', size = 7, hjust = 3, vjust = 3, face = 'bold')) 
  ggsave( file = paste("culicoides_",i,".pdf",sep = ""),width = 7,height = 7) 
}

#mosquito(ring)
rm(list = ls())
options(stringsAsFactors=F)
setwd("/Intermediate_data/mosquito/strand")
load("mosquito_strand.RData")
mosquito_data <- mosquito_data[ , c("Genome_Composition_Family","position_name_out","Sample","Reads_number")]
a <- paste(mosquito_data[,1],mosquito_data[,2],sep = "_")
a1 <- cbind(a,mosquito_data$Reads_number)
a1 <- as.data.frame(a1)
colnames(a1) <- c("name","reads")
a1$name <- as.character(a1$name)
a1$reads <- as.character(a1$reads)
a1$reads <- as.numeric(a1$reads)
m <- aggregate(a1$reads,list(a1$name),sum)

b <- as.data.frame(str_split(m[,1],"_",simplify = T,n=2))
mosquito_data_new <- cbind(b,m[,2])
colnames(mosquito_data_new) <- c("Type","Position","Reads_number")
write.table(mosquito_data_new,file = "mosquito_position_strand.xls",sep = "\t",row.names = F,quote = F)

rm(list = ls())
options(stringsAsFactors=F)
setwd("/Intermediate_data/mosquito/strand")
mosquito_data_new <- read.table(file = "mosquito_position_strand.xls",sep = "\t",header = T,quote = NULL)
mosquito_data_new <- mosquito_data_new[mosquito_data_new$Type != "Unassigned" ,]
unique(mosquito_data_new$Position)
mosquito_data_new[mosquito_data_new$Position == "JiangCheng2013" ,2] <- "JiangCheng"
mosquito_data_new[mosquito_data_new$Position == "JiangCheng2014" ,2] <- "JiangCheng"
mosquito_data_new[mosquito_data_new$Position == "JiangCheng2015" ,2] <- "JiangCheng"
mosquito_data_new[mosquito_data_new$Position == "JiangCheng2016" ,2] <- "JiangCheng"

mosquito_data_new[mosquito_data_new$Position == "JiangCheng2017" ,2] <- "JiangCheng"
aaa <- data.frame(Type=rep("ssDNA(+)",length(unique(mosquito_data_new$Position))) , Position=unique(mosquito_data_new$Position) , Reads_number=rep(0,length(unique(mosquito_data_new$Position))))
mosquito_data_new <- join(mosquito_data_new,aaa,type="full")

temp <- data.frame(Type = rep(sort(unique(mosquito_data_new$Type)),14),Position = sort(rep(unique(mosquito_data_new$Position),9)),Reads_number = rep(0,126)) 
mosquito_data_new <- rbind(temp,mosquito_data_new)
mosquito_data_new <- mosquito_data_new[order(mosquito_data_new$Type), ]
a <- paste(mosquito_data_new[,1],mosquito_data_new[,2],sep = "_")
a1 <- cbind(a,mosquito_data_new$Reads_number)
a1 <- as.data.frame(a1)
colnames(a1) <- c("name","reads")
a1$name <- as.character(a1$name)
a1$reads <- as.character(a1$reads)
a1$reads <- as.numeric(a1$reads)
m <- aggregate(a1$reads,list(a1$name),sum)
b <- as.data.frame(str_split(m[,1],"_",simplify = T,n=2))
mosquito_data_new <- cbind(b,m[,2])
colnames(mosquito_data_new) <- c("Type","Position","Reads_number")

Set1 <- colorRampPalette(brewer.pal(9,"Set1"))(9)
col <- c(Set1[3:8],Set1[1:2],Set1[9])  

for (i in unique(mosquito_data_new$Position)) {
  ggplot(mosquito_data_new[mosquito_data_new$Position == i, ], aes(x ="", y = Reads_number, fill = Type )) + 
    labs(title=paste("mosquito_",i,sep = ""),x = "", y = "",fill="Type") +
    geom_bar(stat = "identity",width = 0.15) + 
    coord_polar(theta = "y") + 
    theme(axis.ticks = element_blank(),axis.text = element_blank()) +
    theme(panel.grid=element_blank(),panel.border=element_blank()) +
    theme(panel.background = element_rect(fill = "transparent",colour = NA),plot.background = element_rect(fill = "transparent",colour = NA)) +
    guides(fill=guide_legend(reverse=F)) +
    scale_fill_manual(values = col ) +
    theme(legend.position="bottom",legend.text = element_text(colour = 'black', size = 7, hjust = 3, vjust = 3, face = 'bold')) 
  ggsave( file = paste("mosquito_",i,".pdf",sep = ""),width = 7,height = 7) 
}