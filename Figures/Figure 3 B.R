#Figure 3:#######################    Vrial families reads distribution & common species numbers in culicoides and mosquito    ###############################
#B:common virus species number in different spaces                                                                                                          #                              
#############################################################################################################################################################
rm(list = ls())
options(stringsAsFactors=F)
setwd("/Intermediate_data/mosquito/species")
library(plyr)
library(stringr)
library (VennDiagram)
mosquito <- read.table("species_all.xls",header = T,sep = "\t",quote = NULL)
mosquito <- mosquito[mosquito$Type != "Type" , ]
sample_position <- read.table("county_sample.txt",header = T,sep = "\t",quote = NULL)
row.names(sample_position) <- sample_position[,1]
sample_position[sample_position$position_name_out == "JiangCheng2013" ,3] <- "JiangCheng"
sample_position[sample_position$position_name_out == "JiangCheng2014" ,3] <- "JiangCheng"
sample_position[sample_position$position_name_out == "JiangCheng2015" ,3] <- "JiangCheng"
sample_position[sample_position$position_name_out == "JiangCheng2016" ,3] <- "JiangCheng"
mosquito <- join(mosquito,sample_position,type="full")
mosquito <- mosquito[ ,c("Type","position_name_out")]
mosquito <- unique(mosquito)
mosquito1 <- aggregate(mosquito$Type,list(mosquito$position_name_out),length)
colnames(mosquito1) <- c("Position","Number")
row.names(mosquito1) <- mosquito1[,1]

setwd("/Intermediate_data/culicoides/species")
culicoides <- read.table("species_all.xls",header = T,sep = "\t",quote = NULL)
culicoides <- culicoides[culicoides$Type != "Type" , ]
sample_position <- read.table("county_sample.txt",header = T,sep = "\t",quote = NULL)
sample_position[sample_position$position_name_out == "JiangCheng2013" ,3] <- "JiangCheng"
sample_position[sample_position$position_name_out == "JiangCheng2014" ,3] <- "JiangCheng"
sample_position[sample_position$position_name_out == "JiangCheng2015" ,3] <- "JiangCheng"
sample_position[sample_position$position_name_out == "JiangCheng2016" ,3] <- "JiangCheng"
sample_position[sample_position$position_name_out == "JiangCheng2017" ,3] <- "JiangCheng"
culicoides <- join(culicoides,sample_position,type="full")
culicoides <- culicoides[ ,c("Type","position_name_out")]
culicoides <- unique(culicoides)
culicoides1 <- aggregate(culicoides$Type,list(culicoides$position_name_out),length)
colnames(culicoides1) <- c("Position","Number")
row.names(culicoides1) <- culicoides1[,1]

#"DeHong"      "DongChuan"   "FuNing"      "HuaNing"     "JiangCheng"  "MengLa"      "ShuangJiang" "SongMing"    "WenShan"
for (i in intersect(mosquito1[,1],culicoides1[,1]) ){
  pdf(paste(i,".pdf",sep = ""),height = 7,width = 7)
  p=venn.diagram(x = list(C.species=unique(culicoides[culicoides$position_name_out == i ,1]), M.species=unique(mosquito[mosquito$position_name_out == i ,1]) ),
                 filename=NULL,fill = c("dodgerblue", "goldenrod1"),col="transparent",main = i)
  grid.draw(p)
  dev.off()
}