#Figure 5:###################################    Human,poultry,livestock-infected viruses in culicoides and mosquito    #####################################
#A:Viruses host distribution in C/M                                                                                                                         #                              
#############################################################################################################################################################
#culicoides
rm(list = ls())
options(stringsAsFactors=F)
library(ggplot2)
library(RColorBrewer)
setwd("/Intermediate_data/culicoides/host")
host <-read.table(file = "host.txt",header = F,sep = "\t",quote = NULL)
colnames(host) <- c("virus.name","refseq.id","Source_information","Number_of_segments","Genome_length","Numbe_of_proteins",
                    "Genome_Neighbors","host.name","Date_completed","Date_updated")
host <- host[host$refseq.id != "" , ]
host[duplicated(host$virus.name), 1] <- "Trichomonas vaginalis virus dsRNA satellite S1_"

hostdb <-read.table(file = "virushostdb.txt",header = T,sep = "\t",quote = NULL)
hostdb[hostdb$virus.name == "Trichomonas vaginalis virus dsRNA satellite S1'" , 2] <- "Trichomonas vaginalis virus dsRNA satellite S1_"
hostdb$virus.name <- gsub("\'" ,"",hostdb$virus.name)
a <- hostdb[ ,1:2]
a <- unique(a)
row.names(a) <- a[,2]

host_new <- host[ ,c("virus.name","refseq.id","host.name")]
row.names(host_new) <- host_new$virus.name
host_new$id <- a[host_new$virus.name,1]
b <- host_new[is.na(host_new$id),]
row.names(b) <- b$refseq.id
c <- hostdb[hostdb$refseq.id %in% b[,2] ,c(1,4)]
c <- unique(c)
row.names(c) <- c$refseq.id
b[rownames(c) ,4] <- c$virus.tax.id
b[b[,2] == "NC_030403" ,4] <- "1859135"
row.names(b) <- b$virus.name
host_new <- rbind(host_new,b)
host_new <- host_new [!is.na(host_new$id) ,]
row.names(host_new) <- host_new$id
rm(host,hostdb,a,b,c)

culicoides <- read.table(file = "culicoides_information.xls",header = T,sep = "\t",quote = NULL)
culicoides <- culicoides[ , c("taxaid","species","sample","reads_number")]
culicoides$host <- host_new[rownames(culicoides) ,3]
culicoides$host <- gsub("\"","",culicoides$host)
culicoides[is.na(culicoides$host) ,5] <- "Others"
culicoides[culicoides$host == "" ,5] <- "Others"

type <- c("plants", "vertebrates", "bacteria", "invertebrates", "fungi")
culicoides <- culicoides[culicoides$host %in% type , ]
culicoides1 <- aggregate(culicoides$reads_number,list(culicoides$host),sum)
colnames(culicoides1) <- c("Type","reads")
culicoides1$Type <- factor(culicoides1$Type,levels = c("vertebrates","invertebrates","plants","bacteria","fungi"))

color <- colorRampPalette(brewer.pal(8,"Set2"))(8)

pdf("culicoides_host.pdf",width = 17,height = 7)
ggplot(culicoides1, aes(x ="", y = reads, fill = Type )) + 
  labs(title="",x = "", y = "",fill="Family") +
  geom_bar(stat = "identity",width = 1) + 
  coord_polar(theta = "y") + 
  theme(axis.ticks = element_blank(),axis.text = element_blank()) +
  theme(panel.grid=element_blank(),panel.border=element_blank()) +
  theme(panel.background = element_rect(fill = "transparent",colour = NA),plot.background = element_rect(fill = "transparent",colour = NA)) +
  guides(fill=guide_legend(reverse=F)) +
  scale_fill_manual(values = color[1:6] ) +
  theme(legend.position="bottom",legend.text = element_text(colour = 'black', size = 7, hjust = 3, vjust = 3, face = 'bold')) 
dev.off()

#mosquito
rm(list = ls())
options(stringsAsFactors=F)
setwd("/Intermediate_data/mosquito/host")
host <-read.table(file = "host.txt",header = F,sep = "\t",quote = NULL)
colnames(host) <- c("virus.name","refseq.id","Source_information","Number_of_segments","Genome_length","Numbe_of_proteins",
                    "Genome_Neighbors","host.name","Date_completed","Date_updated")
host <- host[host$refseq.id != "" , ]
host[duplicated(host$virus.name), 1] <- "Trichomonas vaginalis virus dsRNA satellite S1_"

hostdb <-read.table(file = "virushostdb.txt",header = T,sep = "\t",quote = NULL)
hostdb[hostdb$virus.name == "Trichomonas vaginalis virus dsRNA satellite S1'" , 2] <- "Trichomonas vaginalis virus dsRNA satellite S1_"
hostdb$virus.name <- gsub("\'" ,"",hostdb$virus.name)
a <- hostdb[ ,1:2]
a <- unique(a)
row.names(a) <- a[,2]

host_new <- host[ ,c("virus.name","refseq.id","host.name")]
row.names(host_new) <- host_new$virus.name
host_new$id <- a[host_new$virus.name,1]
b <- host_new[is.na(host_new$id),]
row.names(b) <- b$refseq.id
c <- hostdb[hostdb$refseq.id %in% b[,2] ,c(1,4)]
c <- unique(c)
row.names(c) <- c$refseq.id
b[rownames(c) ,4] <- c$virus.tax.id
b[b[,2] == "NC_030403" ,4] <- "1859135"
row.names(b) <- b$virus.name
host_new <- rbind(host_new,b)
host_new <- host_new [!is.na(host_new$id) ,]
row.names(host_new) <- host_new$id
rm(host,hostdb,a,b,c)

mosquito <- read.table(file = "mosquito_information.xls",header = T,sep = "\t",quote = NULL)
mosquito <- mosquito[ , c("taxaid","species","sample","reads_number")]
mosquito$host <- host_new[rownames(mosquito) ,3]
mosquito$host <- gsub("\"","",mosquito$host)
mosquito[is.na(mosquito$host) ,5] <- "Others"
mosquito[mosquito$host == "" ,5] <- "Others"

type <- c("plants", "vertebrates", "bacteria", "invertebrates", "fungi")
mosquito <- mosquito[mosquito$host %in% type , ]
mosquito1 <- aggregate(mosquito$reads_number,list(mosquito$host),sum)
colnames(mosquito1) <- c("Type","reads")
mosquito1$Type <- factor(mosquito1$Type,levels = c("vertebrates","invertebrates","plants","bacteria","fungi"))

color <- colorRampPalette(brewer.pal(8,"Set2"))(8)

pdf("mosquito_host.pdf",width = 17,height = 7)
ggplot(mosquito1, aes(x ="", y = reads, fill = Type )) + 
  labs(title="",x = "", y = "",fill="Family") +
  geom_bar(stat = "identity",width = 1) + 
  coord_polar(theta = "y") + 
  theme(axis.ticks = element_blank(),axis.text = element_blank()) +
  theme(panel.grid=element_blank(),panel.border=element_blank()) +
  theme(panel.background = element_rect(fill = "transparent",colour = NA),plot.background = element_rect(fill = "transparent",colour = NA)) +
  guides(fill=guide_legend(reverse=F)) +
  scale_fill_manual(values = color[1:6] ) +
  theme(legend.position="bottom",legend.text = element_text(colour = 'black', size = 7, hjust = 3, vjust = 3, face = 'bold')) 
dev.off()
