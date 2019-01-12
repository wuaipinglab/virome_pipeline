#Figure 3:#######################    Vrial families reads distribution & common species numbers in culicoides and mosquito    ###############################
#A:viral families reads in C/M                                                                                                                              #                              
#############################################################################################################################################################
#script for all families, some were shown in result pictures
rm(list = ls())
options(stringsAsFactors=F)
library(RColorBrewer)
library(pheatmap)
library(reshape)
library(reshape2)
library(stringr)
library(plyr)
#input culicoides data
setwd("/Intermediate_data/culicoides/family")
culicoides <- read.table("family_all.xls",header = T,sep = "\t",quote = NULL)
culicoides <- culicoides[culicoides$Type != "Type" , ]
culicoides <- culicoides[order(culicoides$Type,decreasing = T), ]
culicoides$Type <- factor(culicoides$Type)   
culicoides$Reads_number <- as.numeric(culicoides$Reads_number)
culicoides$Percent <- as.numeric(culicoides$Percent)
length(unique(culicoides$Type)) 
length(unique(culicoides$Sample))
culicoides$Type <- as.character(culicoides$Type)

tmp <- data.frame(Type = rep(unique(culicoides$Type),61),Sample = sort(rep(unique(culicoides$Sample),79)),Reads_number = rep(0,4819))
tmp <- rbind(tmp,culicoides[ ,c(1,3,2)])
tmp <- tmp[order(tmp$Type), ]
a <- paste(tmp[,1],tmp[,2],sep = "_")
a1 <- cbind(a,tmp$Reads_number)
a1 <- as.data.frame(a1)
colnames(a1) <- c("name","reads")
a1$name <- as.character(a1$name)
a1$reads <- as.character(a1$reads)
a1$reads <- as.numeric(a1$reads)
m <- aggregate(a1$reads,list(a1$name),sum)
b <- as.data.frame(str_split(m[,1],"_",simplify = T,n=2))
culicoides_new <- cbind(b,m[,2])
colnames(culicoides_new) <- c("Type","Sample","Reads_number")
culicoides_new$Reads_number <- as.numeric(culicoides_new$Reads_number)
culicoides_matrix<-dcast(culicoides_new,culicoides_new$Sample~culicoides_new$Type,value.var = 'Reads_number')
row.names(culicoides_matrix) <- culicoides_matrix[ ,1]
culicoides_matrix <- culicoides_matrix[ ,-1]

position <- read.table("county_sample.txt",header = T,sep = "\t",quote = NULL)
position <- position[order(position$position_name_out) , ]
row.names(position) <- position$Sample

culicoides_matrix <- culicoides_matrix[rownames(position) , ]
culicoides_matrix <- cbind(culicoides_matrix,position)
write.table(culicoides_matrix,file = "culicoides_sample_family_reads_matrix.xls",sep = "\t",quote = F)

#input mosquito data
setwd("/Intermediate_data/mosquito/family")
mosquito <- read.table("family_all.xls",header = T,sep = "\t",quote = NULL) 
mosquito <- mosquito[mosquito$Type != "Type" , ]
mosquito <- mosquito[order(mosquito$Type,decreasing = T), ]
mosquito$Type <- factor(mosquito$Type)   
mosquito$Reads_number <- as.numeric(mosquito$Reads_number)
mosquito$Percent <- as.numeric(mosquito$Percent)
length(unique(mosquito$Type))   
length(unique(mosquito$Sample))
mosquito$Type <- as.character(mosquito$Type)

tmp <- data.frame(Type = rep(unique(mosquito$Type),62),Sample = sort(rep(unique(mosquito$Sample),72)),Reads_number = rep(0,4464)) 
tmp <- rbind(tmp,mosquito[ ,c(1,3,2)])
tmp <- tmp[order(tmp$Type), ]
a <- paste(tmp[,1],tmp[,2],sep = "_")
a1 <- cbind(a,tmp$Reads_number)
a1 <- as.data.frame(a1)
colnames(a1) <- c("name","reads")
a1$name <- as.character(a1$name)
a1$reads <- as.character(a1$reads)
a1$reads <- as.numeric(a1$reads)
m <- aggregate(a1$reads,list(a1$name),sum)
b <- as.data.frame(str_split(m[,1],"_",simplify = T,n=2))
mosquito_new <- cbind(b,m[,2])
colnames(mosquito_new) <- c("Type","Sample","Reads_number")
mosquito_new$Reads_number <- as.numeric(mosquito_new$Reads_number)
new_matrix<-dcast(mosquito_new,mosquito_new$Sample~mosquito_new$Type,value.var = 'Reads_number')
row.names(mosquito_matrix) <- mosquito_matrix[ ,1]
mosquito_matrix <- mosquito_matrix[ ,-1]

position <- read.table("county_sample.txt",header = T,sep = "\t",quote = NULL) 
position <- position[order(position$position_name_out) , ]
row.names(position) <- position$Sample

mosquito_matrix <- mosquito_matrix[rownames(position) , ]
mosquito_matrix <- cbind(mosquito_matrix,position)
write.table(mosquito_matrix,file = "mosquito_sample_family_reads_matrix.xls",sep = "\t",quote = F)

#input data over!start drawing!

#heatmap for culicoides
rm(list = ls())
setwd("/Intermediate_data/culicoides/family")
culicoides <- read.table("culicoides_sample_family_reads_matrix.xls",header = T,sep = "\t")
culicoides <- culicoides[ , 1:79] 
culicoides <- culicoides[ , colnames(culicoides)[-78]] 
culicoides <- as.matrix(culicoides)

ann <- data.frame(position = as.character(culicoides$position_name_out),row.names = rownames(culicoides))
ann$position <- as.character(ann$position)    
ann1 <- data.frame(position = matrix(unlist(ann), nrow=61, byrow=T),stringsAsFactors=FALSE,row.names = rownames(ann))  
ann1$position <- as.character(ann1$position)
ann1 <- data.frame(ann[order(ann1$position), 1],row.names = rownames(ann)[order(ann1$position)])
colnames(ann1) <- "position"

dark2 <- colorRampPalette(brewer.pal(8,"Dark2"))(8)
paired <- colorRampPalette(brewer.pal(12,"Paired"))(12)
accent <- colorRampPalette(brewer.pal(8,"Accent"))(8)
Set1 <- colorRampPalette(brewer.pal(9,"Set1"))(9)
color <- c(dark2[3:5],accent[4],dark2[6:8],paired[1:2],paired[10],paired[4:6],paired[8]) 

names(color) <- unique(ann1$position)
ann_color <- list(position = color)

culicoides <- culicoides[rownames(ann1) , ]
culicoides <- culicoides + 1
culicoides <- log2(culicoides)
culicoides_data <- culicoides
culicoides_ann <- ann1
save(culicoides_data,culicoides_ann,file = "culicoides_data.RData")

pdf("family_sample_reads_number_percent_culicoides.pdf",height = 12,width = 12)
pheatmap(t(culicoides),fontsize_row = 5,fontsize_col = 4 ,cluster_rows = T,cluster_cols = F,show_colnames = T,show_rownames = T ,
         border_color = NA,cellwidth = 6,cellheight = 6,
         color = colorRampPalette(c("white","steelblue"))(100) ,
         annotation_col = ann,annotation_colors = ann_color)
dev.off()

#heatmap for mosquito
rm(list = ls())
setwd("/Intermediate_data/mosquito/family")
mosquito <- read.table("mosquito_sample_family_reads_matrix.xls",header = T,sep = "\t")
mosquito <- mosquito[ , 1:72]  
mosquito <- mosquito[ , colnames(mosquito)[-71]] 
mosquito <- as.matrix(mosquito)

ann <- data.frame(position = as.character(mosquito$position_name_out),row.names = rownames(mosquito))
ann$position <- as.character(ann$position)       
ann1 <- data.frame(position = matrix(unlist(ann), nrow=62, byrow=T),stringsAsFactors=FALSE,row.names = rownames(ann))  
ann1$position <- as.character(ann1$position)
ann1 <- data.frame(ann[order(ann1$position), 1],row.names = rownames(ann)[order(ann1$position)])
colnames(ann1) <- "position"

dark2 <- colorRampPalette(brewer.pal(8,"Dark2"))(8)
paired <- colorRampPalette(brewer.pal(12,"Paired"))(12)
accent <- colorRampPalette(brewer.pal(8,"Accent"))(8)
Set1 <- colorRampPalette(brewer.pal(9,"Set1"))(9)
color <- c(dark2,paired[1:9])  

names(color) <- unique(ann1$position)
ann_color <- list(position = color)

mosquito <- mosquito[rownames(ann1) , ]
mosquito <- mosquito + 1
mosquito <- log2(mosquito)

mosquito_data <- mosquito
mosquito_ann <- ann1
save(mosquito_data,mosquito_ann,file = "mosquito_data.RData")

pdf("family_sample_reads_number_percent_mosquito.pdf",height = 12,width = 12)
pheatmap(t(mosquito),fontsize_row = 5,fontsize_col = 4 ,cluster_rows = T,cluster_cols = F,show_colnames = T,show_rownames = T ,
         border_color = NA,cellwidth = 6,cellheight = 6,
         color = colorRampPalette(c("white","steelblue"))(100) ,
         annotation_col = ann,annotation_colors = ann_color)
dev.off()

#heatmap for culicoides and mosquito(partial)
rm(list = ls())
setwd("/Intermediate_data/mosquito/family")
load("mosquito_data.RData")
setwd("/Intermediate_data/culicoides/family")
load("culicoides_data.RData")
mosquito_data1 <- matrix(data = 0,nrow = 62,ncol = 12,dimnames = list(rownames(mosquito_data),setdiff(colnames(culicoides_data),colnames(mosquito_data)) ) )
culicoides_data1 <- matrix(data = 0,nrow = 61,ncol = 5,dimnames = list(rownames(culicoides_data),setdiff(colnames(mosquito_data),colnames(culicoides_data)) ) )
mosquito_data <- cbind(mosquito_data,mosquito_data1)
culicoides_data <- cbind(culicoides_data,culicoides_data1)
mosquito_data <- as.data.frame(mosquito_data)
culicoides_data <- as.data.frame(culicoides_data)
culicoides_mosquito_data <- join(mosquito_data,culicoides_data,type = "full")
row.names(culicoides_mosquito_data) <- c(rownames(mosquito_data),rownames(culicoides_data))
culicoides_mosquito_ann <- rbind(wz_ann,km_ann)
culicoides_mosquito_ann$type <- c(rep("Mosquito",62),rep("Culicoides",61))

dark2 <- colorRampPalette(brewer.pal(8,"Dark2"))(8)
paired <- colorRampPalette(brewer.pal(12,"Paired"))(12)
accent <- colorRampPalette(brewer.pal(8,"Accent"))(8)
Set1 <- colorRampPalette(brewer.pal(9,"Set1"))(9)
color <- c(dark2[1:5],accent[4],dark2[6:8],paired[1:10])  #19
Set2 <- colorRampPalette(brewer.pal(8,"Set2"))(8)
color1 <- c(Set2[1:2])
names(color) <- sort(unique(culicoides_mosquito_ann$position))
names(color1) <- sort(unique(culicoides_mosquito_ann$type))
ann_color <- list(position = color,type = color1)

a <- culicoides_mosquito_ann[1:62 , ]
a <- a[order(a$position) , ]
b <- culicoides_mosquito_ann[63:123 , ]
b <- b[order(b$position), ]
culicoides_mosquito_ann <- rbind(a,b)
rm(a,b)
culicoides_mosquito_data <- culicoides_mosquito_data[rownames(culicoides_mosquito_ann) ,]

sort_family <- read.table("/Intermediate_data/sort_family.txt")
sort_family <- sort_family[,1]
new_data <- culicoides_mosquito_data[ ,sort_family]

pdf("family_sample_reads_number_percent_culicoides_mosquito.pdf",height = 14,width = 16)
pheatmap(t(new_data),fontsize_row = 10,fontsize_col = 4 ,cluster_rows = F,cluster_cols = F,show_colnames = T,show_rownames = T ,
         border_color = NA,cellwidth = 6,cellheight = 10,
         color = colorRampPalette(c("white","steelblue"))(100) ,
         annotation_col = culicoides_mosquito_ann,annotation_colors = ann_color)
dev.off()