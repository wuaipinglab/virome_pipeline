#Figure 4:#########################    Similarities and differences of viral groups in different space and time(jiangcheng)     #############################
#CDE:interesting pattern in C/M                                                                                                                             #                              
#############################################################################################################################################################
rm(list = ls())
options(stringsAsFactors=F)
library(ggpubr)
library(Rmisc)
library(plyr)
library(stringr)
library(ggplot2)
setwd("/Intermediate_data/culicoides/strand")
load("km_strand.Rdata") #culicoides_data
setwd("/Intermediate_data/mosquito/strand")
load("wz_strand.Rdata")  #mosquito_data

#same position,different medium
#1.jiangcheng13-17 c/m 
mosquito_data <- mosquito_data[mosquito_data$position_name_out %in% c("JiangCheng2013" , "JiangCheng2014" , "JiangCheng2015" , "JiangCheng2016"), ]
mosquito_data <- mosquito_data[ , c("Reads_number","Sample","Genome_Composition_Family")]
mosquito_data <- mosquito_data[mosquito_data$Genome_Composition_Family != "Unassigned" , ]
a <- data.frame(Reads_number=rep(0,8*5),Sample=sort(rep(unique(mosquito_data$Sample),8)), Genome_Composition_Family=rep(unique(mosquito_data$Genome_Composition_Family),5))
mosquito_data <- rbind(mosquito_data,a)
a <- paste(mosquito_data$Sample,mosquito_data$Genome_Composition_Family, sep = "_")
a1 <- cbind(a,mosquito_data$Reads_number)
a1 <- as.data.frame(a1)
colnames(a1) <- c("name","reads")
a1$name <- as.character(a1$name)
a1$reads <- as.numeric(a1$reads)
m <- aggregate(a1$reads,list(a1$name),sum)
b <- as.data.frame(str_split(m[,1],"_",simplify = T,n=2))
mosquito_data <- cbind(b,m[,2])
colnames(mosquito_data) <- c("sample","type","reads")
rm(a,b,a1,m)

mosquito_data$sum <- 0
for(i in unique(mosquito_data$sample))
{
  mosquito_data[mosquito_data$sample == i ,"sum"] <- sum(mosquito_data[mosquito_data$sample == i ,3])
}
mosquito_data$percent <- mosquito_data$reads/mosquito_data$sum
mosquito_data <- mosquito_data[,c("type","percent","reads")]
mosquito_data$host <- rep("M",length(mosquito_data$type))

culicoides_data <- culicoides_data[culicoides_data$position_name_out %in% c("JiangCheng2013" , "JiangCheng2014" , "JiangCheng2015" , "JiangCheng2016","JiangCheng2017"), ]
culicoides_data <- culicoides_data[ , c("Reads_number","Sample","Genome_Composition_Family")]
culicoides_data <- culicoides_data[culicoides_data$Genome_Composition_Family != "Unassigned" , ]
a <- data.frame(Reads_number=rep(0,7*25),Sample=sort(rep(unique(culicoides_data$Sample),7)), Genome_Composition_Family=rep(unique(culicoides_data$Genome_Composition_Family),25))
culicoides_data <- rbind(culicoides_data,a)
a <- paste(culicoides_data$Sample,culicoides_data$Genome_Composition_Family, sep = "_")
a1 <- cbind(a,culicoides_data$Reads_number)
a1 <- as.data.frame(a1)
colnames(a1) <- c("name","reads")
a1$name <- as.character(a1$name)
a1$reads <- as.numeric(a1$reads)
m <- aggregate(a1$reads,list(a1$name),sum)
b <- as.data.frame(str_split(m[,1],"_",simplify = T,n=2))
culicoides_data <- cbind(b,m[,2])
colnames(culicoides_data) <- c("sample","type","reads")
rm(a,b,a1,m)

culicoides_data$sum <- 0
for(i in unique(culicoides_data$sample))
{
  culicoides_data[culicoides_data$sample == i ,"sum"] <- sum(culicoides_data[culicoides_data$sample == i ,3])
}
culicoides_data$percent <- culicoides_data$reads/culicoides_data$sum
culicoides_data <- culicoides_data[,c("type","percent","reads")]
culicoides_data$host <- rep("C",length(culicoides_data$type))

culicoides_mosquito_data <- rbind(culicoides_data,mosquito_data)
culicoides_mosquito_data$percent <- as.numeric(culicoides_mosquito_data$percent)
culicoides_mosquito_data$type <- factor(culicoides_mosquito_data$type)
culicoides_mosquito_data$host <- factor(culicoides_mosquito_data$host)

culicoides_mosquito_data$reads <- log2(culicoides_mosquito_data$reads+1)
tgc <- summarySE(data = culicoides_mosquito_data,measurevar = "percent",groupvars = c("host","type"))

culicoides_mosquito_data <- rbind(culicoides_mosquito_data,c("ssDNA(+)",0,0,"C"))
culicoides_mosquito_data$percent <- as.numeric(culicoides_mosquito_data$percent)
culicoides_mosquito_data$reads <- as.numeric(culicoides_mosquito_data$reads)

compare_means(percent~host,data=culicoides_mosquito_data, group.by = "type",method = "anova")
# A tibble: 8 x 7
#type     .y.           p p.adj p.format p.signif method
#<fctr>   <chr>       <dbl> <dbl>    <chr>    <chr>  <chr>
#1      dsDNA percent 0.200820593 0.800   0.2008       ns  Anova
#2   dsDNA-RT percent 0.263663486 0.800   0.2637       ns  Anova
#3      dsRNA percent 0.079266836 0.480   0.0793       ns  Anova
#4 ssDNA(+/-) percent 0.486076082 0.970   0.4861       ns  Anova
#5   ssRNA-RT percent 0.137148932 0.690   0.1371       ns  Anova
#6   ssRNA(-) percent 0.003481322 0.028   0.0035       **  Anova
#7   ssRNA(+) percent 0.018102334 0.130   0.0181        *  Anova
#8   ssDNA(+) percent 0.592978893 0.970   0.5930       ns  Anova

pdf("culicoides_mosquito_space_jiangcheng_boxplot.pdf",width = 7,height = 7)
p_bar <- ggplot(tgc,aes(x = type, y = percent,fill=host)) + 
  geom_bar(stat = "identity", position = "dodge",width = 0.5) + 
  geom_errorbar(aes(ymin=percent-se,ymax=percent+se),position = position_dodge(0.5), width = 0.3) + 
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

#same medium,different position
#close postion
#2.yuxi: JChuan2015/CJ2015/HN2015/TH2015
mosquito_data <- mosquito_data[mosquito_data$position_name_out %in% c("JiangChuan" , "ChengJiang" , "HuaNing" , "TongHai"), ]
mosquito_data <- mosquito_data[ , c("Reads_number","Sample","Genome_Composition_Family")]
mosquito_data <- mosquito_data[mosquito_data$Genome_Composition_Family != "Unassigned" , ]
a <- data.frame(Reads_number=rep(0,8*20),Sample=sort(rep(unique(mosquito_data$Sample),8)), Genome_Composition_Family=rep(unique(mosquito_data$Genome_Composition_Family),20))
mosquito_data <- rbind(mosquito_data,a)
a <- paste(mosquito_data$Sample,mosquito_data$Genome_Composition_Family, sep = "_")
a1 <- cbind(a,mosquito_data$Reads_number)
a1 <- as.data.frame(a1)
colnames(a1) <- c("name","reads")
a1$name <- as.character(a1$name)
a1$reads <- as.numeric(a1$reads)
m <- aggregate(a1$reads,list(a1$name),sum)
b <- as.data.frame(str_split(m[,1],"_",simplify = T,n=2))
mosquito_data <- cbind(b,m[,2])
colnames(mosquito_data) <- c("sample","type","reads")
rm(a,b,a1,m)

mosquito_data$sum <- 0
for(i in unique(mosquito_data$sample))
{
  mosquito_data[mosquito_data$sample == i ,"sum"] <- sum(mosquito_data[mosquito_data$sample == i ,3])
}
mosquito_data$percent <- mosquito_data$reads/mosquito_data$sum

a <- mosquito_data[ ,c("Sample","position_name_out")]
a <- unique(a)
rownames(a) <- a[,1]
b <- a[mosquito_data[,1] , 2]
mosquito_data$position <- b
rm(a,b,i)

mosquito_data$type <- factor(mosquito_data$type)
mosquito_data$position <- factor(mosquito_data$position)

mosquito_data$reads <- log2(mosquito_data$reads+1)
tgc <- summarySE(data = mosquito_data,measurevar = "percent",groupvars = c("position","type"))

compare_means(percent~position, data=mosquito_data,method = "anova",group.by = "type")
# A tibble: 8 x 7
#type     .y.         p p.adj p.format p.signif method
#<fctr>   <chr>     <dbl> <dbl>    <chr>    <chr>  <chr>
#1      dsDNA percent 0.1142049  0.91     0.11       ns  Anova
#2   dsDNA-RT percent 0.8452470  1.00     0.85       ns  Anova
#3      dsRNA percent 0.5925461  1.00     0.59       ns  Anova
#4 ssDNA(+/-) percent 0.4821585  1.00     0.48       ns  Anova
#5   ssRNA-RT percent 0.2120328  1.00     0.21       ns  Anova
#6   ssRNA(-) percent 0.1897710  1.00     0.19       ns  Anova
#7   ssRNA(+) percent 0.1252438  0.91     0.13       ns  Anova
#8 ssRNA(+/-) percent 0.4104270  1.00     0.41       ns  Anova

pdf("culicoides_mosquito_space_yuxi_boxplot.pdf",width = 7,height = 7)
p_bar <- ggplot(tgc,aes(x = type, y = percent,fill=position)) + 
  geom_bar(stat = "identity", position = "dodge",width = 1) + 
  geom_errorbar(aes(ymin=percent-se,ymax=percent+se),position = position_dodge(1), width = 0.7) + 
  theme_bw() +
  scale_fill_manual(values =  c("red","blue","yellow","brown")) +
  theme(panel.border=element_blank(), panel.grid=element_blank(),axis.line=element_line(colour = "black") ) + 
  theme(axis.text.x=element_text(face = "bold", colour="black", size=10,angle = 90) , axis.text.y=element_text(face = "bold", colour="black", size=10)) + 
  theme(axis.title = element_text(face = "bold", colour="black", size=14)) +
  guides(fill=guide_legend(reverse=F)) + 
  theme(legend.position="top",legend.text = element_text(colour = 'black', size = 10,  hjust = 2, vjust = 2, face = 'bold'),
        legend.key.size = unit(0.8,'cm')) 
p_bar
dev.off() 

#2.1
#JChuan2015/CJ2015/HN2015/TH2015
#jiangchuan&chengjiang > no significant
mosquito_data1 <- mosquito_data[mosquito_data$position %in% c("JiangChuan","ChengJiang"),]
compare_means(percent~position, data=mosquito_data1,method = "anova",group.by = "type")

#jiangchuan&huaning > no significant
mosquito_data1 <- mosquito_data[mosquito_data$position %in% c("JiangChuan","HuaNing"),]
compare_means(percent~position, data=mosquito_data1,method = "anova",group.by = "type")

#jiangchuan&tonghai > no significant
mosquito_data1 <- mosquito_data[mosquito_data$position %in% c("JiangChuan","TongHai"),]
compare_means(percent~position, data=mosquito_data1,method = "anova",group.by = "type")

#chengjiang&huaning > no significant
mosquito_data1 <- mosquito_data[mosquito_data$position %in% c("ChengJiang","HuaNing"),]
compare_means(percent~position, data=mosquito_data1,method = "anova",group.by = "type")

#chengjiang&tonghai > dsDNA *
mosquito_data1 <- mosquito_data[mosquito_data$position %in% c("ChengJiang","TongHai"),]
compare_means(percent~position, data=mosquito_data1,method = "anova",group.by = "type")

#huaning&tonghai > no significant
mosquito_data1 <- mosquito_data[mosquito_data$position %in% c("HuaNing","TongHai"),]
compare_means(percent~position, data=mosquito_data1,method = "anova",group.by = "type")

#same medium,different position
#close postion
#3.yongsheng & binchuan
mosquito_data <- mosquito_data[mosquito_data$position_name_out %in% c("YongSheng","BinChuan"), ]
mosquito_data <- mosquito_data[ , c("Reads_number","Sample","Genome_Composition_Family")]
mosquito_data <- mosquito_data[mosquito_data$Genome_Composition_Family != "Unassigned" , ]
a <- data.frame(Reads_number=rep(0,8*7),Sample=sort(rep(unique(mosquito_data$Sample),8)), Genome_Composition_Family=rep(unique(mosquito_data$Genome_Composition_Family),7))
mosquito_data <- rbind(mosquito_data,a)
a <- paste(mosquito_data$Sample,mosquito_data$Genome_Composition_Family, sep = "_")
a1 <- cbind(a,mosquito_data$Reads_number)
a1 <- as.data.frame(a1)
colnames(a1) <- c("name","reads")
a1$name <- as.character(a1$name)
a1$reads <- as.numeric(a1$reads)
m <- aggregate(a1$reads,list(a1$name),sum)
b <- as.data.frame(str_split(m[,1],"_",simplify = T,n=2))
mosquito_data <- cbind(b,m[,2])
colnames(mosquito_data) <- c("sample","type","reads")
rm(a,b,a1,m)

mosquito_data$sum <- 0
for(i in unique(mosquito_data$sample))
{
  mosquito_data[mosquito_data$sample == i ,"sum"] <- sum(mosquito_data[mosquito_data$sample == i ,3])
}
mosquito_data$percent <- mosquito_data$reads/mosquito_data$sum

a <- mosquito_data[ ,c("Sample","position_name_out")]
a <- unique(a)
rownames(a) <- a[,1]
b <- a[mosquito_data[,1] , 2]
mosquito_data$position <- b
rm(a,b,i)

mosquito_data$type <- factor(mosquito_data$type)
mosquito_data$position <- factor(mosquito_data$position)

mosquito_data$reads <- log2(mosquito_data$reads+1)
tgc <- summarySE(data = mosquito_data,measurevar = "percent",groupvars = c("position","type"))

compare_means(percent~position, data=mosquito_data,method = "anova",group.by = "type")
# A tibble: 8 x 7
#type     .y.            p   p.adj p.format p.signif method
#<fctr>   <chr>        <dbl>   <dbl>    <chr>    <chr>  <chr>
#1      dsDNA percent 4.743908e-01 1.00000   0.4744       ns  Anova
#2   dsDNA-RT percent 8.094635e-01 1.00000   0.8095       ns  Anova
#3      dsRNA percent 3.733208e-01 1.00000   0.3733       ns  Anova
#4 ssDNA(+/-) percent 5.545096e-01 1.00000   0.5545       ns  Anova
#5   ssRNA-RT percent 4.652685e-01 1.00000   0.4653       ns  Anova
#6   ssRNA(-) percent 1.533058e-03 0.01100   0.0015       **  Anova
#7   ssRNA(+) percent 6.725576e-05 0.00054  6.7e-05     ****  Anova
#8 ssRNA(+/-) percent 5.761317e-01 1.00000   0.5761       ns  Anova

pdf("culicoides_mosquito_space_yongsheng_binchuan_boxplot.pdf",width = 7,height = 7)
p_bar <- ggplot(tgc,aes(x = type, y = percent,fill=position)) + 
  geom_bar(stat = "identity", position = "dodge",width = 0.6) + 
  geom_errorbar(aes(ymin=percent-se,ymax=percent+se),position = position_dodge(0.6), width = 0.4) + 
  theme_bw() +
  scale_fill_manual(values =  c("red","blue","yellow","brown")) +
  theme(panel.border=element_blank(), panel.grid=element_blank(),axis.line=element_line(colour = "black") ) + 
  theme(axis.text.x=element_text(face = "bold", colour="black", size=10,angle = 90) , axis.text.y=element_text(face = "bold", colour="black", size=10)) + 
  theme(axis.title = element_text(face = "bold", colour="black", size=14)) +
  guides(fill=guide_legend(reverse=F)) + 
  theme(legend.position="top",legend.text = element_text(colour = 'black', size = 10,  hjust = 2, vjust = 2, face = 'bold'),
        legend.key.size = unit(0.8,'cm')) 
p_bar
dev.off() 

#same medium,different position
#close postion
#4.mengla & jiangcheng--culicoides
culicoides_data <- culicoides_data[culicoides_data$position_name_out %in% c("MengLa" , "JiangCheng2013","JiangCheng2014","JiangCheng2015","JiangCheng2016","JiangCheng2017"), ]
culicoides_data[culicoides_data$position_name_out !="MengLa" ,8] <- "JiangCheng"
culicoides_data <- culicoides_data[ , c("Reads_number","Sample","Genome_Composition_Family")]
culicoides_data <- culicoides_data[culicoides_data$Genome_Composition_Family != "Unassigned" , ]
a <- data.frame(Reads_number=rep(0,7*30),Sample=sort(rep(unique(culicoides_data$Sample),7)), Genome_Composition_Family=rep(unique(culicoides_data$Genome_Composition_Family),30))
culicoides_data <- rbind(culicoides_data,a)
a <- paste(culicoides_data$Sample,culicoides_data$Genome_Composition_Family, sep = "_")
a1 <- cbind(a,culicoides_data$Reads_number)
a1 <- as.data.frame(a1)
colnames(a1) <- c("name","reads")
a1$name <- as.character(a1$name)
a1$reads <- as.numeric(a1$reads)
m <- aggregate(a1$reads,list(a1$name),sum)
b <- as.data.frame(str_split(m[,1],"_",simplify = T,n=2))
culicoides_data <- cbind(b,m[,2])
colnames(culicoides_data) <- c("sample","type","reads")
rm(a,b,a1,m)

culicoides_data$sum <- 0
for(i in unique(culicoides_data$sample))
{
  culicoides_data[culicoides_data$sample == i ,"sum"] <- sum(culicoides_data[culicoides_data$sample == i ,3])
}
culicoides_data$percent <- culicoides_data$reads/culicoides_data$sum

a <- culicoides_data[ ,c("Sample","position_name_out")]
a <- unique(a)
rownames(a) <- a[,1]
b <- a[culicoides_data[,1] , 2]
culicoides_data$position <- b
rm(a,b,i)

culicoides_data[culicoides_data$position!="MengLa" ,6] <- "JiangCheng"
culicoides_data$type <- factor(culicoides_data$type)
culicoides_data$position <- factor(culicoides_data$position)

culicoides_data$reads <- log2(culicoides_data$reads+1)
tgc <- summarySE(data = culicoides_data,measurevar = "percent",groupvars = c("position","type"))

compare_means(percent~position, data=culicoides_data,method = "anova",group.by = "type")
# A tibble: 7 x 7
#type     .y.            p   p.adj p.format p.signif method
#<fctr>   <chr>        <dbl>   <dbl>    <chr>    <chr>  <chr>
#1      dsDNA percent 7.637390e-02 0.38000    0.076       ns  Anova
#2   dsDNA-RT percent 6.043859e-01 1.00000    0.604       ns  Anova
#3      dsRNA percent 5.790461e-01 1.00000    0.579       ns  Anova
#4 ssDNA(+/-) percent 5.387979e-01 1.00000    0.539       ns  Anova
#5   ssRNA-RT percent 2.141245e-05 0.00015  2.1e-05     ****  Anova
#6   ssRNA(-) percent 5.645564e-02 0.34000    0.056       ns  Anova
#7   ssRNA(+) percent 3.820600e-01 1.00000    0.382       ns  Anova

pdf("culicoides_mosquito_space_mengla_jiangcheng_boxplot.pdf",width = 7,height = 7)
p_bar <- ggplot(tgc,aes(x = type, y = percent,fill=position)) + 
  geom_bar(stat = "identity", position = "dodge",width = 0.5) + 
  geom_errorbar(aes(ymin=percent-se,ymax=percent+se),position = position_dodge(0.5), width = 0.3) + 
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

#same medium,different position
#close postion
#5.mengla & jiangcheng--mosquito
mosquito_data <- mosquito_data[mosquito_data$position_name_out %in% c("MengLa" , "JiangCheng2013","JiangCheng2014","JiangCheng2015","JiangCheng2016","JiangCheng2017"), ]
mosquito_data[mosquito_data$position_name_out !="MengLa" ,8] <- "JiangCheng"
mosquito_data <- mosquito_data[ , c("Reads_number","Sample","Genome_Composition_Family")]
mosquito_data <- mosquito_data[mosquito_data$Genome_Composition_Family != "Unassigned" , ]
a <- data.frame(Reads_number=rep(0,9*10),Sample=sort(rep(unique(mosquito_data$Sample),9)), Genome_Composition_Family=rep(unique(mosquito_data$Genome_Composition_Family),10))
mosquito_data <- rbind(mosquito_data,a)
a <- paste(mosquito_data$Sample,mosquito_data$Genome_Composition_Family, sep = "_")
a1 <- cbind(a,mosquito_data$Reads_number)
a1 <- as.data.frame(a1)
colnames(a1) <- c("name","reads")
a1$name <- as.character(a1$name)
a1$reads <- as.numeric(a1$reads)
m <- aggregate(a1$reads,list(a1$name),sum)
b <- as.data.frame(str_split(m[,1],"_",simplify = T,n=2))
mosquito_data <- cbind(b,m[,2])
colnames(mosquito_data) <- c("sample","type","reads")
rm(a,b,a1,m)

mosquito_data$sum <- 0
for(i in unique(mosquito_data$sample))
{
  mosquito_data[mosquito_data$sample == i ,"sum"] <- sum(mosquito_data[mosquito_data$sample == i ,3])
}
mosquito_data$percent <- mosquito_data$reads/mosquito_data$sum

a <- mosquito_data[ ,c("Sample","position_name_out")]
a <- unique(a)
rownames(a) <- a[,1]
b <- a[mosquito_data[,1] , 2]
mosquito_data$position <- b
rm(a,b,i)

mosquito_data[mosquito_data$position!="MengLa" ,6] <- "JiangCheng"
mosquito_data$type <- factor(mosquito_data$type)
mosquito_data$position <- factor(mosquito_data$position)

mosquito_data$reads <- log2(mosquito_data$reads+1)
tgc <- summarySE(data = mosquito_data,measurevar = "percent",groupvars = c("position","type"))
compare_means(percent~position, data=mosquito_data,method = "anova",group.by = "type")
# A tibble: 9 x 7
#type     .y.          p p.adj p.format p.signif method
#<fctr>   <chr>      <dbl> <dbl>    <chr>    <chr>  <chr>
#1      dsDNA percent 0.11255103  0.79    0.113       ns  Anova
#2   dsDNA-RT percent 0.07297616  0.66    0.073       ns  Anova
#3      dsRNA percent 0.29107640  1.00    0.291       ns  Anova
#4   ssDNA(+) percent 0.19314162  1.00    0.193       ns  Anova
#5 ssDNA(+/-) percent 0.98627944  1.00    0.986       ns  Anova
#6   ssRNA-RT percent 0.09796477  0.78    0.098       ns  Anova
#7   ssRNA(-) percent 0.22588561  1.00    0.226       ns  Anova
#8   ssRNA(+) percent 0.21801676  1.00    0.218       ns  Anova
#9 ssRNA(+/-) percent 0.34659351  1.00    0.347       ns  Anova

pdf("culicoides_mosquito_space_mengla_jiangcheng_boxplot.pdf",width = 7,height = 7)
p_bar <- ggplot(tgc,aes(x = type, y = percent,fill=position)) + 
  geom_bar(stat = "identity", position = "dodge",width = 0.5) + 
  geom_errorbar(aes(ymin=percent-se,ymax=percent+se),position = position_dodge(0.5), width = 0.3) + 
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

#far position
#6.JChuan2015/DH2013
mosquito_data <- mosquito_data[mosquito_data$position_name_out %in% c("JiangChuan" , "DeHong"), ]
mosquito_data <- mosquito_data[ , c("Reads_number","Sample","Genome_Composition_Family")]
mosquito_data <- mosquito_data[mosquito_data$Genome_Composition_Family != "Unassigned" , ]
a <- data.frame(Reads_number=rep(0,8*10),Sample=sort(rep(unique(mosquito_data$Sample),8)), Genome_Composition_Family=rep(unique(mosquito_data$Genome_Composition_Family),10))
mosquito_data <- rbind(mosquito_data,a)
a <- paste(mosquito_data$Sample,mosquito_data$Genome_Composition_Family, sep = "_")
a1 <- cbind(a,mosquito_data$Reads_number)
a1 <- as.data.frame(a1)
colnames(a1) <- c("name","reads")
a1$name <- as.character(a1$name)
a1$reads <- as.numeric(a1$reads)
m <- aggregate(a1$reads,list(a1$name),sum)
b <- as.data.frame(str_split(m[,1],"_",simplify = T,n=2))
mosquito_data <- cbind(b,m[,2])
colnames(mosquito_data) <- c("sample","type","reads")
rm(a,b,a1,m)

mosquito_data$sum <- 0
for(i in unique(mosquito_data$sample))
{
  mosquito_data[mosquito_data$sample == i ,"sum"] <- sum(mosquito_data[mosquito_data$sample == i ,3])
}
mosquito_data$percent <- mosquito_data$reads/mosquito_data$sum

a <- mosquito_data[ ,c("Sample","position_name_out")]
a <- unique(a)
rownames(a) <- a[,1]
b <- a[mosquito_data[,1] , 2]
mosquito_data$position <- b
rm(a,b,i)

mosquito_data$type <- factor(mosquito_data$type)
mosquito_data$position <- factor(mosquito_data$position)

mosquito_data$reads <- log2(mosquito_data$reads+1)
tgc <- summarySE(data = mosquito_data,measurevar = "percent",groupvars = c("position","type"))

compare_means(percent~position, data=mosquito_data,method = "anova",group.by = "type")
#identity60
# A tibble: 8 x 7
#type     .y.           p p.adj p.format p.signif method
#<fctr>   <chr>       <dbl> <dbl>    <chr>    <chr>  <chr>
#1      dsDNA percent 0.001837473 0.015   0.0018       **  Anova
#2   dsDNA-RT percent 0.182895940 0.910   0.1829       ns  Anova
#3      dsRNA percent 0.024400665 0.150   0.0244        *  Anova
#4 ssDNA(+/-) percent 0.010556154 0.074   0.0106        *  Anova
#5   ssRNA-RT percent 0.190522146 0.910   0.1905       ns  Anova
#6   ssRNA(-) percent 0.455476836 0.910   0.4555       ns  Anova
#7   ssRNA(+) percent 0.197600481 0.910   0.1976       ns  Anova
#8 ssRNA(+/-) percent 0.342398154 0.910   0.3424       ns  Anova

pdf("culicoides_mosquito_space_jiangchuan_dehong_boxplot.pdf",width = 7,height = 7)
p_bar <- ggplot(tgc,aes(x = type, y = percent,fill=position)) + 
  geom_bar(stat = "identity", position = "dodge",width = 0.5) + 
  geom_errorbar(aes(ymin=percent-se,ymax=percent+se),position = position_dodge(0.5), width = 0.3) + 
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

#far position
#7.ShuangJiang2015/DongChuan2016
culicoides_data <- km[km$position_name_out %in% c("ShuangJiang" , "DongChuan"), ]
culicoides_data <- culicoides_data[ , c("Reads_number","Sample","Genome_Composition_Family")]
culicoides_data <- culicoides_data[culicoides_data$Genome_Composition_Family != "Unassigned" , ]
a <- data.frame(Reads_number=rep(0,7*10),Sample=sort(rep(unique(culicoides_data$Sample),7)), Genome_Composition_Family=rep(unique(culicoides_data$Genome_Composition_Family),10))
culicoides_data <- rbind(culicoides_data,a)
a <- paste(culicoides_data$Sample,culicoides_data$Genome_Composition_Family, sep = "_")
a1 <- cbind(a,culicoides_data$Reads_number)
a1 <- as.data.frame(a1)
colnames(a1) <- c("name","reads")
a1$name <- as.character(a1$name)
a1$reads <- as.numeric(a1$reads)
m <- aggregate(a1$reads,list(a1$name),sum)
b <- as.data.frame(str_split(m[,1],"_",simplify = T,n=2))
culicoides_data <- cbind(b,m[,2])
colnames(culicoides_data) <- c("sample","type","reads")
rm(a,b,a1,m)

culicoides_data$sum <- 0
for(i in unique(culicoides_data$sample))
{
  culicoides_data[culicoides_data$sample == i ,"sum"] <- sum(culicoides_data[culicoides_data$sample == i ,3])
}
culicoides_data$percent <- culicoides_data$reads/culicoides_data$sum

a <- km[ ,c("Sample","position_name_out")]
a <- unique(a)
rownames(a) <- a[,1]
b <- a[culicoides_data[,1] , 2]
culicoides_data$position <- b
rm(a,b,i)

culicoides_data$type <- factor(culicoides_data$type)
culicoides_data$position <- factor(culicoides_data$position)

culicoides_data$reads <- log2(culicoides_data$reads+1)
tgc <- summarySE(data = culicoides_data,measurevar = "percent",groupvars = c("position","type"))
compare_means(percent~position, data=culicoides_data,method = "anova",group.by = "type")
# A tibble: 7 x 7
#type     .y.         p p.adj p.format p.signif method
#<fctr>   <chr>     <dbl> <dbl>    <chr>    <chr>  <chr>
#1      dsDNA percent 0.5088046  1.00     0.51       ns  Anova
#2   dsDNA-RT percent 0.6756228  1.00     0.68       ns  Anova
#3      dsRNA percent 0.6132724  1.00     0.61       ns  Anova
#4 ssDNA(+/-) percent 0.1416933  0.85     0.14       ns  Anova
#5   ssRNA-RT percent 0.8164036  1.00     0.82       ns  Anova
#6   ssRNA(-) percent 0.1132316  0.79     0.11       ns  Anova
#7   ssRNA(+) percent 0.2170000  1.00     0.22       ns  Anova

pdf("culicoides_mosquito_space_shuangjiang_dongchuan_boxplot.pdf",width = 7,height = 7)
p_bar <- ggplot(tgc,aes(x = type, y = percent,fill=position)) + 
  geom_bar(stat = "identity", position = "dodge",width = 0.5) + 
  geom_errorbar(aes(ymin=percent-se,ymax=percent+se),position = position_dodge(0.5), width = 0.3) + 
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