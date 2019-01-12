#Figure 5:###################################    Human,poultry,livestock-infected viruses in culicoides and mosquito    #####################################
#B:human,poultry,livestock-infected viruses reads abundance                                                                                                 #                              
#############################################################################################################################################################
#culicoides
rm(list = ls())
options(stringsAsFactors=F)
library(data.table)
library(ggplot2)
setwd("/Intermediate_data/culicoides/host")
hostdb <-read.table(file = "virushostdb.txt",header = T,sep = "\t")
hostdb[is.na(hostdb)] <- "unknown"
culicoides <- read.table(file = "culicoides_identity_cutoff_60.txt",header = T,sep = "\t",quote = NULL)
colnames(culicoides)[19] <- "virus.tax.id"
all_host_virus <- merge(culicoides,hostdb,by="virus.tax.id")
unique(all_host_virus$host.name)
all_host_virus <- all_host_virus[all_host_virus$virus.name!="Dengue virus 1" ,]

#infect human viurs
human_host_virus <- all_host_virus[all_host_virus$host.name == "Homo sapiens" , ]
sort(unique(human_host_virus$virus.name))
sort(unique(human_host_virus$DISEASE))

#infect poultry virus
poultry_host_virus <- all_host_virus[all_host_virus$host.name == "Gallus gallus" | 
                                       all_host_virus$host.name == "Meleagris gallopavo" |
                                       all_host_virus$host.name == "Columba livia" |
                                       all_host_virus$host.name == "Anas platyrhynchos" |
                                       all_host_virus$host.name == "Anser sp." , ]
sort(unique(poultry_host_virus$virus.name))
sort(unique(poultry_host_virus$DISEASE))  

livestock_host_virus <- all_host_virus[all_host_virus$host.name == "Bos taurus" | 
                                         all_host_virus$host.name == "Equus caballus" |
                                         all_host_virus$host.name == "Sus scrofa" |
                                         all_host_virus$host.name == "Ovis aries", ]
sort(unique(livestock_host_virus$virus.name))
sort(unique(livestock_host_virus$DISEASE))  

culicoides_host_virus <- rbind(human_host_virus,poultry_host_virus,livestock_host_virus)
a <- aggregate(culicoides_host_virus$reads_number,list(culicoides_host_virus$virus.name),sum)    #total reads:31564(identity60) 
a$sample <- rep("C",length(a$Group.1))
unique(culicoides_host_virus$host.name)
#"Homo sapiens"   "Anser sp."      "Bos taurus"     "Equus caballus" "Sus scrofa"     "Ovis aries"

#mosquito
setwd("/Intermediate_data/mosquito/host")
hostdb <-read.table(file = "virushostdb.txt",header = T,sep = "\t")
hostdb[is.na(hostdb)] <- "unknown"
mosquito <- read.table(file = "mosquito_identity_cutoff_60.txt",header = T,sep = "\t",quote = NULL)
colnames(mosquito)[19] <- "virus.tax.id"
all_host_virus <- merge(mosquito,hostdb,by="virus.tax.id")
unique(all_host_virus$host.name)

#infect human viurs
human_host_virus <- all_host_virus[all_host_virus$host.name == "Homo sapiens" , ]
sort(unique(human_host_virus$virus.name))
sort(unique(human_host_virus$DISEASE))

#infect poultry virus
poultry_host_virus <- all_host_virus[all_host_virus$host.name == "Gallus gallus" | 
                                       all_host_virus$host.name == "Meleagris gallopavo" |
                                       all_host_virus$host.name == "Columba livia" |
                                       all_host_virus$host.name == "Anas platyrhynchos" |
                                       all_host_virus$host.name == "Anser sp." , ]
sort(unique(poultry_host_virus$virus.name))
sort(unique(poultry_host_virus$DISEASE))  

#infect livestock
livestock_host_virus <- all_host_virus[all_host_virus$host.name == "Bos taurus" | 
                                         all_host_virus$host.name == "Equus caballus" |
                                         all_host_virus$host.name == "Sus scrofa" |
                                         all_host_virus$host.name == "Ovis aries", ]
sort(unique(livestock_host_virus$virus.name))
sort(unique(livestock_host_virus$DISEASE))  

mosquito_host_virus <- rbind(human_host_virus,poultry_host_virus,livestock_host_virus)
b <- aggregate(mosquito_host_virus$reads_number,list(mosquito_host_virus$virus.name),sum)    #total reads:207361
b$sample <- rep("M",length(b$Group.1))
unique(mosquito_host_virus$host.name)
#"Homo sapiens"   "Gallus gallus"  "Anser sp."      "Sus scrofa"     "Bos taurus"     "Equus caballus"

culicoides_mosquito <- rbind(a,b)
colnames(culicoides_mosquito) <- c("virus_name","reads","sample")
culicoides_mosquito$reads <- log2(culicoides_mosquito$reads)

pdf("culicoides_mosquito_infect_host_virus.pdf",width = 7,height = 7)
p_bar <- ggplot(culicoides_mosquito,aes(x = virus_name, y = reads,fill=sample)) + 
  geom_bar(stat = "identity", position = "dodge",width = 0.5) + 
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