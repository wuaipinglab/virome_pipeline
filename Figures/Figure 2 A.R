#Figure 2:######################################     viral groups reads distribution for culicoides and mosquito     ########################################
#A:assigned reads in C/M                                                                                                                                    #                              
#############################################################################################################################################################
rm(list = ls())
options(stringsAsFactors=F)
library(ggplot2)
library(stringr)
library(plyr)
library(ggpubr)

#culicoides
setwd("/Intermediate_data/culicoides/strand")
#genome composition information table from ICTV
mydata <- read.table("strand_information.txt",header = T,sep = "\t")

family_strand <- mydata[ , c(2,7)]
family_strand <- family_strand[!duplicated(family_strand[,1]), ]
row.names(family_strand) <- family_strand[,1]

species_strand <- mydata[ , c(5,7,2)]
species_strand <- species_strand[!duplicated(species_strand[,1]), ]
row.names(species_strand) <- species_strand[,1]

rm(mydata)
#genome composition information of all mosquito samples identified by blast result
all_taxonomy <- read.table("all_taxonomy.xls",sep = "\t",header = T,quote = NULL) 
all_taxonomy <- unique(all_taxonomy)
row.names(all_taxonomy) <- all_taxonomy$species
#species information of all mosquito samples identified by blast result
mosquito_species <- read.table("species_all.xls",header = T,sep = "\t",quote = NULL)
mosquito_species <- mosquito_species[mosquito_species$Type != "Type" , ]
mosquito_species <- mosquito_species[order(mosquito_species$Type,decreasing = T), ]
mosquito_species$Type <- factor(mosquito_species$Type)   
mosquito_species$Reads_number <- as.numeric(mosquito_species$Reads_number)
mosquito_species$Percent <- as.numeric(mosquito_species$Percent)
length(unique(mosquito_species$Type))  
length(unique(mosquito_species$Sample)) 
mosquito_species$Type <- as.character(mosquito_species$Type)
mosquito_species$Family <- all_taxonomy[mosquito_species$Type , "family"]
rm(all_taxonomy)

mosquito_species_assigned <- mosquito_species[mosquito_species$Family != "Unassigned" , ]
mosquito_species_assigned$Genome_Composition_Family <- family_strand[mosquito_species_assigned$Family ,2]
mosquito_species_assigned[is.na(mosquito_species_assigned$Genome_Composition_Family) , 6] <- 
  species_strand[mosquito_species_assigned[is.na(mosquito_species_assigned$Genome_Composition_Family) , 1] ,2]
mosquito_species_assigned[is.na(mosquito_species_assigned$Genome_Composition_Family) , 6] <- "Unassigned"

mosquito_species_unassigned <- mosquito_species[mosquito_species$Family == "Unassigned" , ]
mosquito_species_unassigned$Genome_Composition_Family <- species_strand[mosquito_species_unassigned$Type ,2]
mosquito_species_unassigned[is.na(mosquito_species_unassigned)] <- "Unassigned"

mosquito_species_strand <- rbind(mosquito_species_assigned,mosquito_species_unassigned)
rm(mosquito_species_unassigned,mosquito_species_assigned)

#county information of each mosquito sample
position <- read.table("county_sample.txt",header = T,sep = "\t",quote = NULL) 
row.names(position) <- position$sample_name
mosquito_data <- join(mosquito_species_strand,position,type = "full")
mosquito_data$position_name <- iconv(mosquito_data$position_name, from="GBK")
mosquito_data$position_name <- as.factor(mosquito_data$position_name)
length(unique(mosquito_data$Genome_Composition_Family))
save(mosquito_data,file = "mosquito_strand.RData")

mosquito <- mosquito_data[ , c("Genome_Composition_Family","position_name_out","Sample","Reads_number")]
mosquito[mosquito$Genome_Composition_Family != "Unassigned" , 1] <- "Assigned"

mosquito <- mosquito[ ,-2]
a <- paste(mosquito[,1],mosquito[,2],sep = "_")
a1 <- cbind(a,mosquito$Reads_number)
a1 <- as.data.frame(a1)
colnames(a1) <- c("name","reads")
a1$name <- as.character(a1$name)
a1$reads <- as.numeric(a1$reads)
m <- aggregate(a1$reads,list(a1$name),sum)
b <- as.data.frame(str_split(m[,1],"_",simplify = T,n=2))
mosquito_new <- cbind(b,m[,2])
colnames(mosquito_new) <- c("Type","Sample","Reads_number")
mosquito_new$sum <- 0
for(i in unique(mosquito_new$Sample))
{
  mosquito_new[mosquito_new$Sample == i ,"sum"] <- sum(mosquito_new[mosquito_new$Sample == i ,3])
}
mosquito_new$percent <- mosquito_new$Reads_number/mosquito_new$sum

mosquito_assigned_frame <- mosquito_new[mosquito_new$Type == "Assigned" , ]
mosquito_assigned_frame$Host <- rep("M",length(mosquito_assigned_frame$Type))


#culicoides
setwd("/Intermediate_data/strand")
#genome composition information table from ICTV
mydata <- read.table("strand_information.txt",header = T,sep = "\t")

family_strand <- mydata[ , c(2,7)]
family_strand <- family_strand[!duplicated(family_strand[,1]), ]
row.names(family_strand) <- family_strand[,1]

species_strand <- mydata[ , c(5,7,2)]
species_strand <- species_strand[!duplicated(species_strand[,1]), ]
row.names(species_strand) <- species_strand[,1]

rm(mydata)
#genome composition information of all culicoides samples identified by blast result
all_taxonomy <- read.table("all_taxonomy.xls",sep = "\t",header = T,quote = NULL)  #all_taxonomy_km.xls  all_taxonomy.xls
all_taxonomy <- unique(all_taxonomy)
row.names(all_taxonomy) <- all_taxonomy$species
#species information of all culicoides samples identified by blast result
culicoides_species <- read.table("species_all.xls",header = T,sep = "\t",quote = NULL)  #61species_all_km.xls 62species_all.xls 
culicoides_species <- culicoides_species[culicoides_species$Type != "Type" , ]
culicoides_species <- culicoides_species[order(culicoides_species$Type,decreasing = T), ]
culicoides_species$Type <- factor(culicoides_species$Type)   
culicoides_species$Reads_number <- as.numeric(culicoides_species$Reads_number)
culicoides_species$Percent <- as.numeric(culicoides_species$Percent)
length(unique(culicoides_species$Type))  
length(unique(culicoides_species$Sample)) 
culicoides_species$Type <- as.character(culicoides_species$Type)
culicoides_species$Family <- all_taxonomy[culicoides_species$Type , "family"]
rm(all_taxonomy)

culicoides_species_assigned <- culicoides_species[culicoides_species$Family != "Unassigned" , ]
culicoides_species_assigned$Genome_Composition_Family <- family_strand[culicoides_species_assigned$Family ,2]
culicoides_species_assigned[is.na(culicoides_species_assigned$Genome_Composition_Family) , 6] <- 
  species_strand[culicoides_species_assigned[is.na(culicoides_species_assigned$Genome_Composition_Family) , 1] ,2]
culicoides_species_assigned[is.na(culicoides_species_assigned$Genome_Composition_Family) , 6] <- "Unassigned"

culicoides_species_unassigned <- culicoides_species[culicoides_species$Family == "Unassigned" , ]
culicoides_species_unassigned$Genome_Composition_Family <- species_strand[culicoides_species_unassigned$Type ,2]
culicoides_species_unassigned[is.na(culicoides_species_unassigned)] <- "Unassigned"

culicoides_species_strand <- rbind(culicoides_species_assigned,culicoides_species_unassigned)
rm(culicoides_species_unassigned,culicoides_species_assigned)
#county information of each culicoides sample
position <- read.table("county_sample.txt",header = T,sep = "\t",quote = NULL) #county_sample61_km.txt  county_sample62.txt
row.names(position) <- position$sample_name
culicoides_data <- join(culicoides_species_strand,position,type = "full")
culicoides_data$position_name <- iconv(culicoides_data$position_name, from="GBK")
culicoides_data$position_name <- as.factor(culicoides_data$position_name)
length(unique(culicoides_data$Genome_Composition_Family))
save(culicoides_data,file = "culicoides_strand.RData")

culicoides <- culicoides_data[ , c("Genome_Composition_Family","position_name_out","Sample","Reads_number")]
culicoides[culicoides$Genome_Composition_Family != "Unassigned" , 1] <- "Assigned"

culicoides <- culicoides[ ,-2]
a <- paste(culicoides[,1],culicoides[,2],sep = "_")
a1 <- cbind(a,culicoides$Reads_number)
a1 <- as.data.frame(a1)
colnames(a1) <- c("name","reads")
a1$name <- as.character(a1$name)
a1$reads <- as.numeric(a1$reads)
m <- aggregate(a1$reads,list(a1$name),sum)
b <- as.data.frame(str_split(m[,1],"_",simplify = T,n=2))
culicoides_new <- cbind(b,m[,2])
colnames(culicoides_new) <- c("Type","Sample","Reads_number")
culicoides_new$sum <- 0
for(i in unique(culicoides_new$Sample))
{
  culicoides_new[culicoides_new$Sample == i ,"sum"] <- sum(culicoides_new[culicoides_new$Sample == i ,3])
}
culicoides_new$percent <- culicoides_new$Reads_number/culicoides_new$sum

culicoides_assigned_frame <- culicoides_new[culicoides_new$Type == "Assigned" , ]
culicoides_assigned_frame$Host <- rep("C",length(culicoides_assigned_frame$Type))

culicoides_mosquito_data <- rbind(culicoides_assigned_frame,mosquito_assigned_frame)
col <- c("red","blue")

compare_means(percent~Host, data=culicoides_mosquito_data,method = "t.test")
# A tibble: 1 x 8
#.y. group1 group2           p  p.adj p.format p.signif method
#<chr>  <chr>  <chr>       <dbl>  <dbl>    <chr>    <chr>  <chr>
#1 percent      C      M 0.002536327 0.0025   0.0025       ** T-test

pdf("Culicoides_mosquito_percent_of_assigned_boxplot.pdf",width = 7,height = 5)
p_bar <- ggplot(culicoides_mosquito_data,aes(x = Host, y = percent)) + 
  geom_boxplot(width = 0.2,aes(colour=Host),outlier.colour = NA) + 
  geom_jitter(width = 0.1,aes(colour=Host)) +
  scale_y_continuous(limits = c(0,1)) + 
  theme_bw() +
  stat_compare_means(method = "t.test",label.x = 0.5,label.y = 1) + 
  stat_compare_means(label = "p.signif") + 
  scale_fill_manual(values =  col) +
  scale_colour_manual(values =  col) +
  theme(panel.border=element_blank(), panel.grid=element_blank(),axis.line=element_line(colour = "black") ) + 
  theme(axis.text.x=element_text(face = "bold", colour="black", size=10,angle = 90) , axis.text.y=element_text(face = "bold", colour="black", size=10)) + 
  theme(axis.title = element_text(face = "bold", colour="black", size=14)) +
  guides(fill=guide_legend(reverse=F)) + 
  theme(legend.position="top",legend.text = element_text(colour = 'black', size = 10, hjust = 2, vjust = 2, face = 'bold'),legend.key.size = unit(0.8,'cm')) 
p_bar
dev.off() 
