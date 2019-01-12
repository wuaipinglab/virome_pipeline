#plot virome within viral family level
argv <- commandArgs(TRUE)
options(stringsAsFactors=F)
setwd(argv[1])
library(RColorBrewer)
library(ggplot2)
mv_fp_contig <- read.table(argv[2],header = T,quote = NULL,sep = "\t")
family <- aggregate(mv_fp_contig$reads_number,list(mv_fp_contig$family),sum)
colnames(family) <- c("Type","Reads_number")
family$Sample <- argv[3]
family$Percent <- (family$Reads_number/sum(family$Reads_number))*100
write.table(family,file=argv[4],sep = "\t",row.names = F,quote = F)

pdf(argv[5],width = 7,height = 7)
p <- ggplot(family, aes(x = Sample, y = Percent, fill = Type)) +
	  labs(title=paste(argv[3]," virome composition by family")) +
	    geom_bar(stat = "identity",position = "stack",width = 0.2) + 
		  theme_bw() +
		    theme(panel.border=element_blank(), panel.grid=element_blank(),axis.line=element_line(colour = "black") ) + 
			  scale_fill_manual(values = colorRampPalette(brewer.pal(12,"Set3"))(length(family$Type)))
p
dev.off()

species <- aggregate(mv_fp_contig$reads_number,list(mv_fp_contig$species),sum)
colnames(species) <- c("Type","Reads_number")
species$Sample <- argv[3]
species$Percent <- (species$Reads_number/sum(species$Reads_number))*100
write.table(species,file=argv[6],sep = "\t",row.names = F,quote = F)

pdf(argv[7],width = 7,height = 7)
p <- ggplot(species, aes(x = Sample, y = Percent, fill = Type)) +
	  labs(title=paste(argv[3]," virome composition by species")) +
	    geom_bar(stat = "identity",position = "stack",width = 0.2) + 
		  theme_bw() +
		    theme(panel.border=element_blank(), panel.grid=element_blank(),axis.line=element_line(colour = "black") ) + 
			  scale_fill_manual(values = colorRampPalette(brewer.pal(12,"Set3"))(length(species$Type)))
p
dev.off()
