#for IDS sample validations


#data from scimap
#spatial kmeans neighborhoods
#cell cell distances

dist <- read.csv('/media/oncosys/Expansion/cellcycle/cellcycle_distances_RCNs.csv')
data <- read.csv('/media/oncosys/Expansion/cellcycle/cellcycle_phenotypes_RCNs.csv')

all_data <- cbind(data, dist)

#spatial kmeans categories:
#TSI 1, 14, 15
#IBA1CD163 17, 13
#IBA1CD11c 9
#CD11c 6
#mixed 12, 2
#T cell 10, 5

library(dplyr)

RCNs <- all_data %>% group_by(spatial_kmeans, phenotype) %>% summarise(n=n()) %>% mutate(freq=n/sum(n))

library(ggplot2)

RCNs$spatial_kmeans <- as.factor(RCNs$spatial_kmeans)
RCNs$phenotype <- factor(RCNs$phenotype, levels=c('CD8Tcells', 'CD4Tcells', 'CD11c.myeloid', 'IBA1.CD11.mac', 'IBA1.CD163.mac', 'CD163.mac', 'Stromal', 'Tumor', 'Unknown'))

#colors

group.colors <- c('Tumor' = "#627D87", 'CD8Tcells' = "#ff5a5f", 'CD4Tcells'="#57cc99", 'CD11c.myeloid'="#996fd6", 
                  'IBA1.CD11.mac'="#b59ce0", 'IBA1.CD163.mac'="#c2b3e0", 'CD163.mac'="#d0c9ea",
                  'Stromal'="#d6ccc2", 'Unknown'="#f5ebe0")


p<- ggplot(RCNs, aes(x=spatial_kmeans, y=freq, fill=phenotype)) + geom_bar(position="fill", stat="identity") + scale_fill_manual(values=group.colors) + theme_bw()

pdf('/media/oncosys/Expansion/cellcycle/plots/barplot_spatial_kmeans.pdf', width=7, height=4)
p
dev.off()
#TSI 1, 14, 15
#IBA1CD163 17, 13
#IBA1CD11c 9
#CD11c 16, 6
#mixed 12, 2
#T cell 10, 5
#tumor 3 11
#stromal 0 4 7 8 

#merge categories
all_data[which(all_data$spatial_kmeans == '1'|all_data$spatial_kmeans == '14'|all_data$spatial_kmeans == '15' ), 'spatial_kmeans'] <- 'RCN1'
all_data[which(all_data$spatial_kmeans == '17'| all_data$spatial_kmeans == '13'), 'spatial_kmeans'] <- 'RCN7'
all_data[which(all_data$spatial_kmeans == '9'), 'spatial_kmeans'] <- 'RCN9'
all_data[which(all_data$spatial_kmeans == '6'|all_data$spatial_kmeans == '16'), 'spatial_kmeans'] <- 'RCN6'
all_data[which(all_data$spatial_kmeans == '2' | all_data$spatial_kmeans == '12'), 'spatial_kmeans'] <- 'RCN2'
all_data[which(all_data$spatial_kmeans == '10'| all_data$spatial_kmeans == '5'), 'spatial_kmeans'] <- 'RCN5'
all_data[which(all_data$spatial_kmeans == '3'|all_data$spatial_kmeans == '11'), 'spatial_kmeans'] <- 'RCN3'
all_data[which(all_data$spatial_kmeans == '0'|all_data$spatial_kmeans == '4'|all_data$spatial_kmeans == '7'|all_data$spatial_kmeans == '8'), 'spatial_kmeans'] <- 'RCN4'

RCNs <- all_data %>% group_by(spatial_kmeans, phenotype) %>% summarise(n=n()) %>% mutate(freq=n/sum(n))

library(ggplot2)

RCNs$spatial_kmeans <- as.factor(RCNs$spatial_kmeans)
RCNs$phenotype <- factor(RCNs$phenotype, levels=c('CD8Tcells', 'CD4Tcells', 'CD11c.myeloid', 'IBA1.CD11.mac', 'IBA1.CD163.mac', 'CD163.mac', 'Stromal', 'Tumor', 'Unknown'))

p<- ggplot(RCNs, aes(x=spatial_kmeans, y=freq, fill=phenotype)) + geom_bar(position="fill", stat="identity") + scale_fill_manual(values=group.colors) + theme_bw()

pdf('/media/oncosys/Expansion/cellcycle/plots/barplot_merged_spatial_kmeans.pdf', width=7, height=4)
p
dev.off()

#
#cells belonging to each RCN+spatially interacting with the RCN

subset <- all_data[which(all_data$spatial_kmeans == 'RCN2' | all_data$X2 <= 45| all_data$X12 <= 45 ), c(1:8)]
write.csv(subset, '/media/oncosys/Expansion/cellcycle/RCN2.csv', row.names = F)

subset <- all_data[which(all_data$spatial_kmeans == 'RCN1' | all_data$X1 <= 45| all_data$X14 <= 45 | all_data$X15 <= 45 ), c(1:8)]
write.csv(subset, '/media/oncosys/Expansion/cellcycle/RCN1.csv', row.names = F)

subset <- all_data[which(all_data$spatial_kmeans == 'RCN5' | all_data$X5 <= 45| all_data$X10 <= 45 ),c(1:8) ]
write.csv(subset, '/media/oncosys/Expansion/cellcycle/RCN5.csv', row.names = F)

subset <- all_data[which(all_data$spatial_kmeans == 'RCN6' | all_data$X6 <= 45| all_data$X16 <= 45 ), c(1:8)]
write.csv(subset, '/media/oncosys/Expansion/cellcycle/RCN6.csv', row.names = F)

subset <- all_data[which(all_data$spatial_kmeans == 'RCN7' | all_data$X17 <= 45| all_data$X13 <= 45 ), c(1:8) ]
write.csv(subset, '/media/oncosys/Expansion/cellcycle/RCN7.csv', row.names = F)

subset <- all_data[which(all_data$spatial_kmeans == 'RCN9' | all_data$X9 <= 45 ), c(1:8)]
write.csv(subset, '/media/oncosys/Expansion/cellcycle/RCN9.csv', row.names = F)







