setwd("~/Google Drive/Archived/Hurwitz Lab/PATRIC-taxoner-species/now-with-fastq")

library(RColorBrewer)

raw_graph_data<-read.csv("source_data_for_R_based_on_new_normalization.csv")

graph_data <- raw_graph_data

sum_by_species <- rowsum(graph_data[,c("percDNA1","percDNA2","percDNA3","percDNA4")],graph_data$species,reorder = T)

first_graph <- sum_by_species[c("Mucispirillum schaedleri","Lactobacillus murinus","Lactobacillus plantarum","Lachnospiraceae bacterium A4","Helicobacter hepaticus","Parabacteroides distasonis"),]

first_graph$bacteria <- factor(row.names(first_graph),levels=c("Mucispirillum schaedleri","Lactobacillus murinus","Lactobacillus plantarum","Lachnospiraceae bacterium A4","Helicobacter hepaticus","Parabacteroides distasonis"),labels = c("M.schaedleri","L.murinus","L.plantarum","L.bacterium-A4","H.hepaticus","P.distasonis"))

row.names(first_graph)<-first_graph$bacteria

#trying ggplot just for fun
library(ggplot2)
library(reshape2)
melted<-melt(first_graph)
colnames(melted)<-c("Species","Sample","Fraction")
melted$Sample<-gsub("percDNA","S",melted$Sample)
melted$Sample <- factor(melted$Sample, c("S3","S4","S1","S2"))
melted <- with(melted, melted[order(Sample),])
colors<-brewer.pal(4,"RdBu")[4:1]
ggplot(melted, aes(x=Sample,y=Fraction,fill=Sample)) + geom_bar(stat="identity") + facet_grid(~Species) + scale_fill_brewer(palette = "RdBu",direction = -1)

first_graph<-first_graph[,c(1,2,3,4)]

first_graph<-t(first_graph)

first_graph_ord<-first_graph[c("percDNA3","percDNA4","percDNA1","percDNA2"),]

colors<-brewer.pal(4,"RdBu")[4:1]

barplot(as.matrix(first_graph_ord)
        ,ylab="Fraction of Total Bacterial Counts"
        ,xlab="Species"
        ,beside=TRUE
        ,col=colors
        ,ylim=c(0,0.25))

legend("topright"
       ,c("SMAD3(+/+),H.hep neg"
          ,"SMAD3(-/-),H.hep neg"
          ,"SMAD3(+/+),H.hep pos"
          ,"SMAD3(-/-),H.hep pos")
       ,cex=.75
       ,fill=colors)

sum_by_phylum <- rowsum(graph_data[,c("percDNA1","percDNA2","percDNA3","percDNA4")],graph_data$phylum,reorder = T)

second_graph <- sum_by_phylum[c("Bacteroidetes","Deferribacteres","Firmicutes","Proteobacteria"),]

second_graph<-t(second_graph)

second_graph_ord<-second_graph[c("percDNA3","percDNA4","percDNA1","percDNA2"),]

barplot(as.matrix(second_graph_ord)
        ,ylab="Fraction of Total Bacterial Counts"
        ,xlab="Phylum"
        ,beside=TRUE
        ,col=colors
        ,ylim=c(0,0.9))

legend("topleft"
       ,c("SMAD3(+/+),H.hep neg"
          ,"SMAD3(-/-),H.hep neg"
          ,"SMAD3(+/+),H.hep pos"
          ,"SMAD3(-/-),H.hep pos")
       ,cex=.75
       ,fill=colors)

sum_by_family <- rowsum(graph_data[,c("percDNA1","percDNA2","percDNA3","percDNA4")],graph_data$family,reorder = T)

third_graph <- sum_by_family[c("Lachnospiraceae","Lactobacillaceae","Helicobacteraceae","Bacteroidaceae","Deferribacteraceae"),]

third_graph<-t(third_graph)

third_graph_ord<-third_graph[c("percDNA3","percDNA4","percDNA1","percDNA2"),]

barplot(as.matrix(third_graph_ord)
        ,ylab="Fraction of Total Bacterial Counts"
        ,xlab="Family"
        ,beside=TRUE
        ,col=colors
        ,ylim=c(0,0.6))

legend("topright"
       ,c("SMAD3(+/+),H.hep neg"
          ,"SMAD3(-/-),H.hep neg"
          ,"SMAD3(+/+),H.hep pos"
          ,"SMAD3(-/-),H.hep pos")
       ,cex=.75
       ,fill=colors)

phylum_changes <- data.frame(phyla=colnames(second_graph),"smad3ko"=0,"inflammation"=0,"cancer"=0)
only_control_phyla<-second_graph_ord[1,]
phylum_changes$smad3ko=second_graph_ord[2,]/only_control_phyla
phylum_changes$inflammation=second_graph_ord[3,]/only_control_phyla
phylum_changes$cancer=second_graph_ord[4,]/only_control_phyla

#write.csv(phylum_changes,"phylum_changes_perc.csv")

family_changes <- data.frame(family=colnames(third_graph),"smad3ko"=0,"inflammation"=0,"cancer"=0)
family_changes$smad3ko<-third_graph_ord[2,]/third_graph_ord[1,]
family_changes$inflammation<-third_graph_ord[3,]/third_graph_ord[1,]
family_changes$cancer<-third_graph_ord[4,]/third_graph_ord[1,]

#write.csv(family_changes,"family_changes_perc.csv")

species_changes <- data.frame(species=colnames(first_graph),"smad3ko"=0,"inflammation"=0,"cancer"=0)
species_changes$smad3ko<-first_graph_ord[2,]/first_graph_ord[1,]
species_changes$inflammation<-first_graph_ord[3,]/first_graph_ord[1,]
species_changes$cancer<-first_graph_ord[4,]/first_graph_ord[1,]

#write.csv(species_changes,"species_changes_perc.csv")
