#Figure 2:######################################     viral groups reads distribution for culicoides and mosquito     ########################################
#B:9 assigned groups reads in C/M                                                                                                                           #                              
#############################################################################################################################################################
rm(list = ls())
options(stringsAsFactors=F)
library(ggplot2)
library(stringr)
library(Rmisc)
setwd("/Intermediate_data/mosquito/strand")
load("mosquito_strand.Rdata")
mosquito_data <- mosquito_data[ , c("Genome_Composition_Family","Sample","Reads_number")]
mosquito_data <- mosquito_data[mosquito_data$Genome_Composition_Family != "Unassigned" , ]
a <- paste(mosquito_data[,1],mosquito_data[,2],sep = "_")
a1 <- cbind(a,mosquito_data$Reads_number)
a1 <- as.data.frame(a1)
colnames(a1) <- c("name","reads")
a1$name <- as.character(a1$name)
a1$reads <- as.numeric(a1$reads)
m <- aggregate(a1$reads,list(a1$name),sum)
b <- as.data.frame(str_split(m[,1],"_",simplify = T,n=2))
mosquito_data <- cbind(b,m[,2])
colnames(mosquito_data) <- c("Type","Sample","Reads_number")
mosquito_data$Reads_number <- log2(mosquito_data$Reads_number)
mosquito_data$host <- rep("M",length(mosquito_data$Type))

setwd("/Intermediate_data/culicoides/strand")
load("culicoides_strand.Rdata")
culicoides_data <- culicoides_data
culicoides_data <- culicoides_data[ , c("Genome_Composition_Family","Sample","Reads_number")]
culicoides_data <- culicoides_data[culicoides_data$Genome_Composition_Family != "Unassigned" , ]
a <- paste(culicoides_data[,1],culicoides_data[,2],sep = "_")
a1 <- cbind(a,culicoides_data$Reads_number)
a1 <- as.data.frame(a1)
colnames(a1) <- c("name","reads")
a1$name <- as.character(a1$name)
a1$reads <- as.numeric(a1$reads)
m <- aggregate(a1$reads,list(a1$name),sum)
b <- as.data.frame(str_split(m[,1],"_",simplify = T,n=2))
culicoides_data <- cbind(b,m[,2])
colnames(culicoides_data) <- c("Type","Sample","Reads_number")
culicoides_data$Reads_number <- log2(culicoides_data$Reads_number)
culicoides_data$host <- rep("C",length(culicoides_data$Type))

culicoides_mosquito <- rbind(mosquito_data,culicoides_data)
culicoides_mosquito$Type <- factor(culicoides_mosquito$Type)
culicoides_mosquito$host <- factor(culicoides_mosquito$host)
tgc <- summarySE(data = culicoides_mosquito,measurevar = "Reads_number",groupvars = c("host","Type"))

pdf("error_barplot_assigned_all_sample_identity_se.pdf",width = 7,height = 7)
p_bar <- ggplot(tgc,aes(x = Type, y = Reads_number,fill=host)) + 
  geom_bar(stat = "identity", position = "dodge",width = 0.5) + 
  geom_errorbar(aes(ymin=Reads_number-se,ymax=Reads_number+se),position = position_dodge(0.5), width = 0.3) + 
  theme_bw() +
  scale_fill_manual(values =  c("red","blue")) +
  theme(panel.border=element_blank(), panel.grid=element_blank(),axis.line=element_line(colour = "black") ) + 
  theme(axis.text.x=element_text(face = "bold", colour="black", size=10,angle = 90) , axis.text.y=element_text(face = "bold", colour="black", size=10)) + 
  theme(axis.title = element_text(face = "bold", colour="black", size=14)) +
  guides(fill=guide_legend(reverse=F)) + 
  theme(legend.position="top",legend.text = element_text(colour = 'black', size = 10,  hjust = 2, vjust = 2, face = 'bold'),
        legend.key.size = unit(0.8,'cm')) 
p_bar
dev.off()