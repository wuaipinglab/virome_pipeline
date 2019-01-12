#Figure 4:#########################    Similarities and differences of viral groups in different space and time(jiangcheng)     #############################
#A:time difference                                                                                                                                          #                              
#############################################################################################################################################################
rm(list = ls())
options(stringsAsFactors=F)
library(ggplot2)
library(RColorBrewer)
library(plyr)
#culicoides
setwd("/Intermediate_data/culicoides/strand")
culicoides_jiangcheng <- read.table(file = "culicoides_jiangcheng.txt",sep = "\t",quote = NULL,header = F)
colnames(culicoides_jiangcheng) <- c("type","position","reads")
culicoides_jiangcheng <- culicoides_jiangcheng[culicoides_jiangcheng$type != "Unassigned" , ]

a <- data.frame(type=rep("ssDNA(+)",length(unique(culicoides_jiangcheng$position))) , position=unique(culicoides_jiangcheng$position) , reads=rep(0,length(unique(culicoides_jiangcheng$position))))
culicoides_jiangcheng <- join(culicoides_jiangcheng,a,type="full")
culicoides_jiangcheng$percent <- 0

culicoides_jiangcheng[culicoides_jiangcheng$position == "JiangCheng2013" ,"percent"] <- 
  culicoides_jiangcheng[culicoides_jiangcheng$position == "JiangCheng2013" ,"reads"]/sum(culicoides_jiangcheng[culicoides_jiangcheng$position == "JiangCheng2013" ,"reads"])
culicoides_jiangcheng[culicoides_jiangcheng$position == "JiangCheng2014" ,"percent"] <- 
  culicoides_jiangcheng[culicoides_jiangcheng$position == "JiangCheng2014" ,"reads"]/sum(culicoides_jiangcheng[culicoides_jiangcheng$position == "JiangCheng2014" ,"reads"])
culicoides_jiangcheng[culicoides_jiangcheng$position == "JiangCheng2015" ,"percent"] <- 
  culicoides_jiangcheng[culicoides_jiangcheng$position == "JiangCheng2015" ,"reads"]/sum(culicoides_jiangcheng[culicoides_jiangcheng$position == "JiangCheng2015" ,"reads"])
culicoides_jiangcheng[culicoides_jiangcheng$position == "JiangCheng2016" ,"percent"] <- 
  culicoides_jiangcheng[culicoides_jiangcheng$position == "JiangCheng2016" ,"reads"]/sum(culicoides_jiangcheng[culicoides_jiangcheng$position == "JiangCheng2016" ,"reads"])
culicoides_jiangcheng[culicoides_jiangcheng$position == "JiangCheng2017" ,"percent"] <- 
  culicoides_jiangcheng[culicoides_jiangcheng$position == "JiangCheng2017" ,"reads"]/sum(culicoides_jiangcheng[culicoides_jiangcheng$position == "JiangCheng2017" ,"reads"])

Set1 <- colorRampPalette(brewer.pal(9,"Set1"))(9)
col <- c(Set1[3:8],Set1[1:2],Set1[9])

pdf("culicoides_jiangcheng_time.pdf",width = 3.5,height = 3.5)
ggplot(culicoides_jiangcheng,aes(x=position,y=percent)) + 
  geom_point(aes(colour=factor(type),size=percent)) + 
  theme_bw() +
  scale_colour_manual(values =  col) +
  theme(panel.border=element_blank(), panel.grid=element_blank(),axis.line=element_line(colour = "black") ) + 
  theme(axis.text.x=element_text(face = "bold", colour="black", size=10,angle = 90) , axis.text.y=element_text(face = "bold", colour="black", size=10)) + 
  theme(axis.title = element_text(face = "bold", colour="black", size=14)) +
  guides(fill=guide_legend(reverse=F)) + 
  theme(legend.position="bottom",legend.text = element_text(colour = 'black', size = 10, hjust = 2, vjust = 2, face = 'bold'),legend.key.size = unit(0.8,'cm')) 
dev.off()

#mosquito
rm(list = ls())
options(stringsAsFactors=F)
setwd("/Intermediate_data/mosquito/strand")
mosquito_jiangcheng <- read.table(file = "mosquito_jiangcheng.txt",sep = "\t",quote = NULL,header = F)
colnames(mosquito_jiangcheng) <- c("type","position","reads")
mosquito_jiangcheng <- mosquito_jiangcheng[mosquito_jiangcheng$type != "Unassigned" , ]

a <- data.frame(type=rep("ssDNA(+)",length(unique(mosquito_jiangcheng$position))) , position=unique(mosquito_jiangcheng$position) , reads=rep(0,length(unique(mosquito_jiangcheng$position))))
mosquito_jiangcheng <- join(mosquito_jiangcheng,a,type="full")

mosquito_jiangcheng$percent <- 0

mosquito_jiangcheng[mosquito_jiangcheng$position == "JiangCheng2013" ,"percent"] <- 
  mosquito_jiangcheng[mosquito_jiangcheng$position == "JiangCheng2013" ,"reads"]/sum(mosquito_jiangcheng[mosquito_jiangcheng$position == "JiangCheng2013" ,"reads"])
mosquito_jiangcheng[mosquito_jiangcheng$position == "JiangCheng2014" ,"percent"] <- 
  mosquito_jiangcheng[mosquito_jiangcheng$position == "JiangCheng2014" ,"reads"]/sum(mosquito_jiangcheng[mosquito_jiangcheng$position == "JiangCheng2014" ,"reads"])
mosquito_jiangcheng[mosquito_jiangcheng$position == "JiangCheng2015" ,"percent"] <- 
  mosquito_jiangcheng[mosquito_jiangcheng$position == "JiangCheng2015" ,"reads"]/sum(mosquito_jiangcheng[mosquito_jiangcheng$position == "JiangCheng2015" ,"reads"])
mosquito_jiangcheng[mosquito_jiangcheng$position == "JiangCheng2016" ,"percent"] <- 
  mosquito_jiangcheng[mosquito_jiangcheng$position == "JiangCheng2016" ,"reads"]/sum(mosquito_jiangcheng[mosquito_jiangcheng$position == "JiangCheng2016" ,"reads"])

b <- data.frame(type=mosquito_jiangcheng[1:9,1],position=rep("JiangCheng2017",9),reads=rep(0,9),percent=rep(0,9))
mosquito_jiangcheng <- rbind(mosquito_jiangcheng,b)
mosquito_jiangcheng$type <- factor(mosquito_jiangcheng$type)
mosquito_jiangcheng <- mosquito_jiangcheng[order(mosquito_jiangcheng$position) , ]

Set1 <- colorRampPalette(brewer.pal(9,"Set1"))(9)
col <- c(Set1[3:8],Set1[1:2],Set1[9]) 

pdf("mosquito_jiangcheng_time.pdf",width = 3.5,height = 3.5)
ggplot(mosquito_jiangcheng,aes(x=position,y=percent)) + 
  geom_point(aes(colour=factor(type),size=percent)) + 
  theme_bw() +
  scale_colour_manual(values =  col) +
  theme(panel.border=element_blank(), panel.grid=element_blank(),axis.line=element_line(colour = "black") ) + 
  theme(axis.text.x=element_text(face = "bold", colour="black", size=10,angle = 90) , axis.text.y=element_text(face = "bold", colour="black", size=10)) + 
  theme(axis.title = element_text(face = "bold", colour="black", size=14)) +
  guides(fill=guide_legend(reverse=F)) + 
  theme(legend.position="bottom",legend.text = element_text(colour = 'black', size = 10, hjust = 2, vjust = 2, face = 'bold'),legend.key.size = unit(0.8,'cm')) 
dev.off()