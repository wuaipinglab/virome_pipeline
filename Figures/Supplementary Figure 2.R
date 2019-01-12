#Supplementary Figure 2:##############################     viral groups reads distribution for culicoides and mosquito     ##################################
#9 assigned groups reads in M                                                                                                                               # 
#############################################################################################################################################################
rm(list = ls())
options(stringsAsFactors=F)
library(ggplot2)
library(RColorBrewer)
library(stringr)
setwd("/Intermediate_data/mosquito/strand")
load("mosquito_strand.Rdata")
mosquito_data <- mosquito_data[ , c("Genome_Composition_Family","position_name_out","Sample","Reads_number")]
mosquito_data <- mosquito_data[mosquito_data$Genome_Composition_Family != "Unassigned" ,]
a <- paste(mosquito_data[,2],mosquito_data[,3],sep = "_")
a1 <- cbind(a,mosquito_data$Genome_Composition_Family,mosquito_data$Reads_number)
a1 <- as.data.frame(a1)
colnames(a1) <- c("name","type","reads")
a1$name <- as.character(a1$name)
a1$reads <- as.character(a1$reads)
a1$reads <- as.numeric(a1$reads)

a <- paste(a1[,1],a1[,2],sep = "_")
a1 <- cbind(a,a1$reads)
a1 <- as.data.frame(a1)
colnames(a1) <- c("name","reads")
a1$name <- as.character(a1$name)
a1$reads <- as.character(a1$reads)
a1$reads <- as.numeric(a1$reads)
m <- aggregate(a1$reads,list(a1$name),sum)

b <- as.data.frame(str_split(m[,1],"_",simplify = T,n=3))
mosquito_data_new <- cbind(b,m[,2])
colnames(mosquito_data_new) <- c("Position","Sample","Type","Reads_number")
mosquito_data_new$sum <- 0
for(i in unique(mosquito_data_new$Sample))
{
  mosquito_data_new[mosquito_data_new$Sample == i ,"sum"] <- sum(mosquito_data_new[mosquito_data_new$Sample == i ,4])
}
mosquito_data_new$percent <- mosquito_data_new$Reads_number/mosquito_data_new$sum

Set1 <- colorRampPalette(brewer.pal(9,"Set1"))(9)
col <- c(Set1[3:8],Set1[1:2],Set1[9])

pdf("mosquito_virome_strand_family.pdf",width = 10,height =11)
p_bar <- ggplot(mosquito_data_new,aes(x = Sample, y = percent  ) ) +
  geom_point(aes(colour=factor(Type),size=percent)) + 
  theme_bw() +
  scale_colour_manual(values =  col) +
  theme(panel.border=element_blank(), panel.grid=element_blank(),axis.line=element_line(colour = "black") ) + 
  theme(axis.text.x=element_text(face = "bold", colour="black", size=10,angle = 90) , axis.text.y=element_text(face = "bold", colour="black", size=10)) + 
  theme(axis.title = element_text(face = "bold", colour="black", size=14)) +
  guides(fill=guide_legend(reverse=F)) + 
  theme(legend.position="bottom",legend.text = element_text(colour = 'black', size = 10, hjust = 2, vjust = 2, face = 'bold'),legend.key.size = unit(0.8,'cm')) + 
  facet_wrap(~ Position,scales = 'free_x',ncol = 5)
p_bar
dev.off() 