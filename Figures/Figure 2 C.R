#Figure 2:######################################     viral groups reads distribution for culicoides and mosquito     ########################################
#C:9 assigned groups reads in C                                                                                                                             #                              
#############################################################################################################################################################
rm(list = ls())
options(stringsAsFactors=F)
library(ggplot2)
library(RColorBrewer)
library(stringr)
setwd("/Intermediate_data/culicoides/strand")
load("culicoides_strand.Rdata")
culicoides_data <- culicoides_data[ , c("Genome_Composition_Family","position_name_out","Sample","Reads_number")]
culicoides_data <- culicoides_data[culicoides_data$Genome_Composition_Family != "Unassigned" ,]
a <- paste(culicoides_data[,2],culicoides_data[,3],sep = "_")
a1 <- cbind(a,culicoides_data$Genome_Composition_Family,culicoides_data$Reads_number)
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
culicoides_data_new <- cbind(b,m[,2])
colnames(culicoides_data_new) <- c("Position","Sample","Type","Reads_number")
culicoides_data_new$sum <- 0
for(i in unique(culicoides_data_new$Sample))
{
  culicoides_data_new[culicoides_data_new$Sample == i ,"sum"] <- sum(culicoides_data_new[culicoides_data_new$Sample == i ,4])
}
culicoides_data_new$percent <- culicoides_data_new$Reads_number/culicoides_data_new$sum

Set1 <- colorRampPalette(brewer.pal(9,"Set1"))(9)
col <- c(Set1[3:5], Set1[7:8],Set1[1:2],Set1[9])

pdf("culicoides_virome_strand_family.pdf",width = 10,height = 9)
p_bar <- ggplot(culicoides_data_new,aes(x = Sample, y = percent  ) ) +
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