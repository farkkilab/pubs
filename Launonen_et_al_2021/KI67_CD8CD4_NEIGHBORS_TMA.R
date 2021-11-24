##########################################################################################################################################
#VIOLIN PLOTS OF KI67 EXPRESSION IN PROLIFERATING EPITHELIAL CELLS AND THE NUMBER OF CD4 AND CD8+T-CELL NEIGHBORS
##########################################################################################################################################

#For Figure 5 and Supplementary Figure 5

#as data you need all neighbor ids

#for each prof epi cell, calculate the count of cd4 or cd8 cells

#then plot expression of Ki67 in the prof epi cells vs count of neighbours

#who are the neighbours of prof epi cells

#cell ids for each cell type


setwd("neighbor_id_files")
temp= list.files(pattern=".csv")
myfiles = lapply(temp, read.csv)

for (i in c(1:112)){
  
  myfiles[[i]]$core <- paste(temp)[i]
}



library(dplyr)

all_neighbour_ids <- bind_rows(myfiles[[1]], myfiles[[2]], myfiles[[3]], myfiles[[4]],myfiles[[5]], myfiles[[6]], myfiles[[7]], myfiles[[8]],
                               myfiles[[9]], myfiles[[10]], myfiles[[11]], myfiles[[12]],myfiles[[13]], myfiles[[14]], myfiles[[15]], myfiles[[16]],
                               myfiles[[17]], myfiles[[18]], myfiles[[19]], myfiles[[20]],myfiles[[21]], myfiles[[22]], myfiles[[23]], myfiles[[24]],
                               myfiles[[25]], myfiles[[26]], myfiles[[27]], myfiles[[28]],myfiles[[29]], myfiles[[30]], myfiles[[31]], myfiles[[32]],
                               myfiles[[33]], myfiles[[34]], myfiles[[35]], myfiles[[36]],myfiles[[37]], myfiles[[38]], myfiles[[39]], myfiles[[40]],
                               myfiles[[41]], myfiles[[42]], myfiles[[43]], myfiles[[44]],myfiles[[45]], myfiles[[46]], myfiles[[47]], myfiles[[48]],
                               myfiles[[49]], myfiles[[50]], myfiles[[51]], myfiles[[52]],myfiles[[53]], myfiles[[54]], myfiles[[55]], myfiles[[56]],
                               myfiles[[57]], myfiles[[58]], myfiles[[59]],myfiles[[60]], myfiles[[61]], myfiles[[62]],myfiles[[63]], myfiles[[64]], myfiles[[65]], myfiles[[66]],
                               myfiles[[67]], myfiles[[68]], myfiles[[69]],myfiles[[70]], myfiles[[71]], myfiles[[72]],myfiles[[73]], myfiles[[74]], myfiles[[75]], myfiles[[76]],
                               myfiles[[77]], myfiles[[78]], myfiles[[79]],myfiles[[80]], myfiles[[81]], myfiles[[82]],myfiles[[83]], myfiles[[84]], myfiles[[85]], myfiles[[86]],
                               myfiles[[87]], myfiles[[88]], myfiles[[89]],myfiles[[90]], myfiles[[91]], myfiles[[92]],myfiles[[93]], myfiles[[94]], myfiles[[95]], myfiles[[96]],
                               myfiles[[97]], myfiles[[98]], myfiles[[99]],myfiles[[100]], myfiles[[101]], myfiles[[102]],myfiles[[103]], myfiles[[104]], myfiles[[105]], myfiles[[106]],
                               myfiles[[107]], myfiles[[108]], myfiles[[109]],myfiles[[110]],myfiles[[111]],myfiles[[112]])


all_neighbour_ids <- all_neighbour_ids[, c(1, 2, 14, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24)]

#all_neighbour_ids$core <- gsub("^.{0,4}", "", all_fractions$core)
#all_fractions$core <- gsub("^.{0,14}", "", all_fractions$core)
all_neighbour_ids$core <- substr(all_neighbour_ids$core, 1, nchar(all_neighbour_ids$core) -4)

all_neighbour_ids[, 1] <- paste0(all_neighbour_ids$core, "_", all_neighbour_ids[, 1])
all_neighbour_ids[, 4] <- paste0(all_neighbour_ids$core, "_", all_neighbour_ids[, 4])
all_neighbour_ids[, 5] <- paste0(all_neighbour_ids$core, "_", all_neighbour_ids[, 5])
all_neighbour_ids[, 6] <- paste0(all_neighbour_ids$core, "_", all_neighbour_ids[, 6])
all_neighbour_ids[, 7] <- paste0(all_neighbour_ids$core, "_", all_neighbour_ids[, 7])
all_neighbour_ids[, 8] <- paste0(all_neighbour_ids$core, "_", all_neighbour_ids[, 8])
all_neighbour_ids[, 9] <- paste0(all_neighbour_ids$core, "_", all_neighbour_ids[, 9])
all_neighbour_ids[, 10] <- paste0(all_neighbour_ids$core, "_", all_neighbour_ids[, 10])
all_neighbour_ids[, 11] <- paste0(all_neighbour_ids$core, "_", all_neighbour_ids[, 11])
all_neighbour_ids[, 12] <- paste0(all_neighbour_ids$core, "_", all_neighbour_ids[, 12])
all_neighbour_ids[, 13] <- paste0(all_neighbour_ids$core, "_", all_neighbour_ids[, 13])
all_neighbour_ids[, 14] <- paste0(all_neighbour_ids$core, "_", all_neighbour_ids[, 14])
all_neighbour_ids[, 15] <- paste0(all_neighbour_ids$core, "_", all_neighbour_ids[, 15])
all_neighbour_ids[, 16] <- paste0(all_neighbour_ids$core, "_", all_neighbour_ids[, 16])
all_neighbour_ids[, 17] <- paste0(all_neighbour_ids$core, "_", all_neighbour_ids[, 17])
all_neighbour_ids[, 18] <- paste0(all_neighbour_ids$core, "_", all_neighbour_ids[, 18])
all_neighbour_ids[, 19] <- paste0(all_neighbour_ids$core, "_", all_neighbour_ids[, 19])
all_neighbour_ids[, 20] <- paste0(all_neighbour_ids$core, "_", all_neighbour_ids[, 20])
all_neighbour_ids[, 21] <- paste0(all_neighbour_ids$core, "_", all_neighbour_ids[, 21])
all_neighbour_ids[, 22] <- paste0(all_neighbour_ids$core, "_", all_neighbour_ids[, 22])
all_neighbour_ids[, 23] <- paste0(all_neighbour_ids$core, "_", all_neighbour_ids[, 23])
all_neighbour_ids[, 24] <- paste0(all_neighbour_ids$core, "_", all_neighbour_ids[, 24])


prof_epi_center_cells <- all_neighbour_ids[which(all_neighbour_ids$cluster == "Proliferating epithelial"),]
cd4_center_cells <- all_neighbour_ids[which(all_neighbour_ids$cluster == "CD4"),]
cd8_center_cells <- all_neighbour_ids[which(all_neighbour_ids$cluster == "CD8"),]

#prof_epi_center_cells[, c(4:24)] <- as.character(prof_epi_center_cells[, c(4:24)])
#cd4_center_cells[, c(4:24)] <- as.character(cd4_center_cells[, c(4:24)])
#keep the neighbours which are CD4

#Keep all the rows but 
#look for CD4 T-cells as neighbors of proliferating epithelial cells
profepi_n1 <- prof_epi_center_cells[prof_epi_center_cells$neighbour_CellId1 %in% cd4_center_cells$CellId, c("CellId", "cluster", "core", "neighbour_CellId1")]
profepi_n2 <- prof_epi_center_cells[prof_epi_center_cells$neighbour_CellId2 %in% cd4_center_cells$CellId, c("CellId", "cluster", "core", "neighbour_CellId2")]
profepi_n3 <- prof_epi_center_cells[prof_epi_center_cells$neighbour_CellId3 %in% cd4_center_cells$CellId, c("CellId", "cluster", "core", "neighbour_CellId3")]
profepi_n4 <- prof_epi_center_cells[prof_epi_center_cells$neighbour_CellId4 %in% cd4_center_cells$CellId, c("CellId", "cluster", "core", "neighbour_CellId4")]
profepi_n5 <- prof_epi_center_cells[prof_epi_center_cells$neighbour_CellId5 %in% cd4_center_cells$CellId, c("CellId", "cluster", "core", "neighbour_CellId5")]
profepi_n6 <- prof_epi_center_cells[prof_epi_center_cells$neighbour_CellId6 %in% cd4_center_cells$CellId, c("CellId", "cluster", "core", "neighbour_CellId6")]
profepi_n7 <- prof_epi_center_cells[prof_epi_center_cells$neighbour_CellId7 %in% cd4_center_cells$CellId, c("CellId", "cluster", "core", "neighbour_CellId7")]
profepi_n8 <- prof_epi_center_cells[prof_epi_center_cells$neighbour_CellId8 %in% cd4_center_cells$CellId, c("CellId", "cluster", "core", "neighbour_CellId8")]
profepi_n9 <- prof_epi_center_cells[prof_epi_center_cells$neighbour_CellId9 %in% cd4_center_cells$CellId, c("CellId", "cluster", "core", "neighbour_CellId9")]
profepi_n10 <- prof_epi_center_cells[prof_epi_center_cells$neighbour_CellId10 %in% cd4_center_cells$CellId, c("CellId", "cluster", "core", "neighbour_CellId10")]
profepi_n11 <- prof_epi_center_cells[prof_epi_center_cells$neighbour_CellId11 %in% cd4_center_cells$CellId, c("CellId", "cluster", "core", "neighbour_CellId11")]
profepi_n13 <- prof_epi_center_cells[prof_epi_center_cells$neighbour_CellId12 %in% cd4_center_cells$CellId, c("CellId", "cluster", "core", "neighbour_CellId12")]
profepi_n14 <- prof_epi_center_cells[prof_epi_center_cells$neighbour_CellId13 %in% cd4_center_cells$CellId, c("CellId", "cluster", "core", "neighbour_CellId13")]
profepi_n15 <- prof_epi_center_cells[prof_epi_center_cells$neighbour_CellId14 %in% cd4_center_cells$CellId, c("CellId", "cluster", "core", "neighbour_CellId14")]
profepi_n16 <- prof_epi_center_cells[prof_epi_center_cells$neighbour_CellId15 %in% cd4_center_cells$CellId, c("CellId", "cluster", "core", "neighbour_CellId15")]
profepi_n17 <- prof_epi_center_cells[prof_epi_center_cells$neighbour_CellId16 %in% cd4_center_cells$CellId, c("CellId", "cluster", "core", "neighbour_CellId16")]
profepi_n18 <- prof_epi_center_cells[prof_epi_center_cells$neighbour_CellId17 %in% cd4_center_cells$CellId, c("CellId", "cluster", "core", "neighbour_CellId17")]
profepi_n19 <- prof_epi_center_cells[prof_epi_center_cells$neighbour_CellId18 %in% cd4_center_cells$CellId, c("CellId", "cluster", "core", "neighbour_CellId18")]
profepi_n20 <- prof_epi_center_cells[prof_epi_center_cells$neighbour_CellId20 %in% cd4_center_cells$CellId, c("CellId", "cluster", "core", "neighbour_CellId20")]
profepi_n21 <- prof_epi_center_cells[prof_epi_center_cells$neighbour_CellId21 %in% cd4_center_cells$CellId, c("CellId", "cluster", "core", "neighbour_CellId21")]
profepi_n12 <- prof_epi_center_cells[prof_epi_center_cells$neighbour_CellId19 %in% cd4_center_cells$CellId, c("CellId", "cluster", "core", "neighbour_CellId19")]



#change the column name of 4th column

colnames(profepi_n1)[4] <- "neighbour_CellId"
colnames(profepi_n2)[4] <- "neighbour_CellId"
colnames(profepi_n3)[4] <- "neighbour_CellId"
colnames(profepi_n4)[4] <- "neighbour_CellId"
colnames(profepi_n5)[4] <- "neighbour_CellId"
colnames(profepi_n6)[4] <- "neighbour_CellId"
colnames(profepi_n7)[4] <- "neighbour_CellId"
colnames(profepi_n8)[4] <- "neighbour_CellId"
colnames(profepi_n9)[4] <- "neighbour_CellId"
colnames(profepi_n10)[4] <- "neighbour_CellId"
colnames(profepi_n11)[4] <- "neighbour_CellId"
colnames(profepi_n12)[4] <- "neighbour_CellId"
colnames(profepi_n13)[4] <- "neighbour_CellId"
colnames(profepi_n14)[4] <- "neighbour_CellId"
colnames(profepi_n15)[4] <- "neighbour_CellId"
colnames(profepi_n16)[4] <- "neighbour_CellId"
colnames(profepi_n17)[4] <- "neighbour_CellId"
colnames(profepi_n18)[4] <- "neighbour_CellId"
colnames(profepi_n19)[4] <- "neighbour_CellId"
colnames(profepi_n20)[4] <- "neighbour_CellId"
colnames(profepi_n21)[4] <- "neighbour_CellId"


profepi_n1$neighbour_CellId <- as.character(profepi_n1$neighbour_CellId)
profepi_n2$neighbour_CellId <- as.character(profepi_n2$neighbour_CellId)
profepi_n3$neighbour_CellId <- as.character(profepi_n3$neighbour_CellId)
profepi_n4$neighbour_CellId <- as.character(profepi_n4$neighbour_CellId)
profepi_n5$neighbour_CellId <- as.character(profepi_n5$neighbour_CellId)
profepi_n6$neighbour_CellId <- as.character(profepi_n6$neighbour_CellId)
profepi_n7$neighbour_CellId <- as.character(profepi_n7$neighbour_CellId)
profepi_n8$neighbour_CellId <- as.character(profepi_n8$neighbour_CellId)
profepi_n9$neighbour_CellId <- as.character(profepi_n9$neighbour_CellId)
profepi_n10$neighbour_CellId <- as.character(profepi_n10$neighbour_CellId)
profepi_n11$neighbour_CellId <- as.character(profepi_n11$neighbour_CellId)
profepi_n12$neighbour_CellId <- as.character(profepi_n12$neighbour_CellId)
profepi_n13$neighbour_CellId <- as.character(profepi_n13$neighbour_CellId)
profepi_n14$neighbour_CellId <- as.character(profepi_n14$neighbour_CellId)
profepi_n15$neighbour_CellId <- as.character(profepi_n15$neighbour_CellId)
profepi_n16$neighbour_CellId <- as.character(profepi_n16$neighbour_CellId)
profepi_n17$neighbour_CellId <- as.character(profepi_n17$neighbour_CellId)
profepi_n18$neighbour_CellId <- as.character(profepi_n18$neighbour_CellId)
profepi_n19$neighbour_CellId <- as.character(profepi_n19$neighbour_CellId)
profepi_n20$neighbour_CellId <- as.character(profepi_n20$neighbour_CellId)
profepi_n21$neighbour_CellId <- as.character(profepi_n21$neighbour_CellId)



#rbind those cells


all_prof_epi_cd4_neighbours <- bind_rows(profepi_n1,profepi_n2,profepi_n3,profepi_n4,profepi_n5,profepi_n6,profepi_n7,profepi_n8,profepi_n9,profepi_n10,profepi_n11,profepi_n12,profepi_n13,profepi_n14,profepi_n15,profepi_n16,profepi_n17,profepi_n18,profepi_n19,profepi_n20,profepi_n21)

#these are the cd4 cells as neighbbours of profepi

#however i need also the prof epi cells as cell id and cluster 


#number of prof_epi_cells

n_of_neighbours_prof_epi <- all_prof_epi_cd4_neighbours %>% group_by(CellId) %>% summarise(number=n())

#922 cells with some neighbours

#prof_epi_center_cells has all cells, take those that don't have neighbours

no_neighbours_prof_epi <- prof_epi_center_cells[-which(prof_epi_center_cells$CellId %in% n_of_neighbours_prof_epi$CellId),]

no_neighbours_prof_epi$number <- "0"


number_of_neigbours_prof_epi <- rbind(n_of_neighbours_prof_epi, no_neighbours_prof_epi[, c("CellId", "number")])

all_celltypes_24092020 <- read.csv("TMA_annotated_single_cell_data.csv")


all_celltypes_24092020_for_neighbourhoods <- all_celltypes_24092020

all_celltypes_24092020_for_neighbourhoods$cores <- paste0("core", all_celltypes_24092020_for_neighbourhoods$cores,"_", all_celltypes_24092020_for_neighbourhoods$Cellid)


number_of_neigbours_prof_epi <- merge(number_of_neigbours_prof_epi, all_celltypes_24092020_for_neighbourhoods[, c("Sample", "cores", "Ki67", "HR_defect")], by.x="CellId", by.y = "cores")

#now scatter plot number of neighbours vs Ki67
number_of_neigbours_prof_epi$HR_defect <- as.character(number_of_neigbours_prof_epi$HR_defect)


library(ggplot2)
dodge <- position_dodge(width = 0.7)

number_of_neigbours_prof_epi[which(number_of_neigbours_prof_epi$HR_defect == "0"), "HR_defect"] <- "HRwt"
number_of_neigbours_prof_epi[which(number_of_neigbours_prof_epi$HR_defect == "1"), "HR_defect"] <- "BRCA1/2 mutated"

theme <- theme(panel.border = element_rect(colour = "black", size=1, fill=NA), 
               plot.title=element_text(hjust=0.5), panel.grid.major =element_blank(),
               panel.background = element_rect(fill = "white"), panel.grid.minor = element_line(size = 0.15, linetype = 'solid', colour = "grey")) 

dodge <- position_dodge(width = 0.5)

my_comparisons <- list(c("0BRCA1/2 mutated", "0HRwt"),c("1BRCA1/2 mutated", "1HRwt"),c("2BRCA1/2 mutated", "2HRwt"), c("3BRCA1/2 mutated", "3HRwt"))

number_of_neigbours_prof_epi$combined <- paste0(number_of_neigbours_prof_epi$number, number_of_neigbours_prof_epi$HR_defect)

p <- ggplot(number_of_neigbours_prof_epi, aes(x=combined, y=Ki67, fill=HR_defect)) + geom_violin(position=dodge, width=1.5) + geom_boxplot(position = dodge, width=0.2)+ theme+
  scale_fill_manual(values=c("#3366CC", "#CC0033")) + ggtitle("Number of CD4+ T-cells as neighbours") + stat_compare_means(method="wilcox", label = "p.signif", comparisons=my_comparisons)

p



# width=12, height=4.5)




dodge <- position_dodge(width = 0.7)
number_of_neigbours_prof_epi$number_of_neighbours <- number_of_neigbours_prof_epi$number

number_of_neigbours_prof_epi$number_of_neighbours[which(number_of_neigbours_prof_epi$number_of_neighbours >= 1)] <- "1 or more CD4+ T-cell as neighbour"

library(ggpubr)
number_of_neigbours_prof_epi$Sample <- as.factor(number_of_neigbours_prof_epi$Sample)

#do a separate variable with n of neighbours and HR status
#calculate values and add to plot

number_of_neigbours_prof_epi$combined_variables <- paste0(number_of_neigbours_prof_epi$HR_defect,"_", number_of_neigbours_prof_epi$number_of_neighbours)





my_comparisons <- list( c("BRCA1/2 mutated_0", "BRCA1/2 mutated_1 or more CD4+ T-cell as neighbour"), c("HRwt_0", "HRwt_1 or more CD4+ T-cell as neighbour"), c("BRCA1/2 mutated_0", "HRwt_0"), c("BRCA1/2 mutated_1 or more CD4+ T-cell as neighbour", "HRwt_1 or more CD4+ T-cell as neighbour") )

p1 <- ggplot(number_of_neigbours_prof_epi, aes(x=combined_variables, y=Ki67, fill=HR_defect)) + geom_violin(position=dodge) + geom_boxplot(position=dodge, width=0.2)+ theme+
  scale_fill_manual(values=c("#3366CC", "#CC0033")) + ggtitle("Number of CD4+ T-cells as neighbours\nof proliferating epithelial cells")  + xlab("HR status") + stat_compare_means(comparisons=my_comparisons, label = "p.signif", label.y = c(11.3,11.3, 11.7, 12.1)) + ylab("Ki67 expression in proliferating epithelial cells") + theme(text = element_text(size = 15))
p1

p1 <- ggplot(number_of_neigbours_prof_epi, aes(x=combined_variables, y=Ki67, fill=HR_defect)) + geom_violin(position=dodge) + geom_boxplot(position=dodge, width=0.2)+ theme+
  scale_fill_manual(values=c("#3366CC", "#CC0033")) + ggtitle("Number of CD4+ T-cells as neighbours\nof proliferating epithelial cells")  + xlab("HR status") + stat_compare_means(comparisons=my_comparisons, label.y = c(11.3,11.3, 11.7, 12.1)) + ylab("Ki67 expression in proliferating epithelial cells") + theme(text = element_text(size = 15))
p1

#width=9, height=5.5)


#p2 <- ggplot(number_of_neigbours_prof_epi, aes(x=Ki67, fill=number_of_neighbours)) + #geom_density(alpha=0.4)+ theme+ facet_grid(. ~ HR_defect) + ggtitle("Number of CD4+ T-cells as #neighbours")
#p2




#same with CD8+ T-cells


#look for CD8+T-cell neighbors
profepi_n1 <- prof_epi_center_cells[prof_epi_center_cells$neighbour_CellId1 %in% cd8_center_cells$CellId, c("CellId", "cluster", "core", "neighbour_CellId1")]
profepi_n2 <- prof_epi_center_cells[prof_epi_center_cells$neighbour_CellId2 %in% cd8_center_cells$CellId, c("CellId", "cluster", "core", "neighbour_CellId2")]
profepi_n3 <- prof_epi_center_cells[prof_epi_center_cells$neighbour_CellId3 %in% cd8_center_cells$CellId, c("CellId", "cluster", "core", "neighbour_CellId3")]
profepi_n4 <- prof_epi_center_cells[prof_epi_center_cells$neighbour_CellId4 %in% cd8_center_cells$CellId, c("CellId", "cluster", "core", "neighbour_CellId4")]
profepi_n5 <- prof_epi_center_cells[prof_epi_center_cells$neighbour_CellId5 %in% cd8_center_cells$CellId, c("CellId", "cluster", "core", "neighbour_CellId5")]
profepi_n6 <- prof_epi_center_cells[prof_epi_center_cells$neighbour_CellId6 %in% cd8_center_cells$CellId, c("CellId", "cluster", "core", "neighbour_CellId6")]
profepi_n7 <- prof_epi_center_cells[prof_epi_center_cells$neighbour_CellId7 %in% cd8_center_cells$CellId, c("CellId", "cluster", "core", "neighbour_CellId7")]
profepi_n8 <- prof_epi_center_cells[prof_epi_center_cells$neighbour_CellId8 %in% cd8_center_cells$CellId, c("CellId", "cluster", "core", "neighbour_CellId8")]
profepi_n9 <- prof_epi_center_cells[prof_epi_center_cells$neighbour_CellId9 %in% cd8_center_cells$CellId, c("CellId", "cluster", "core", "neighbour_CellId9")]
profepi_n10 <- prof_epi_center_cells[prof_epi_center_cells$neighbour_CellId10 %in% cd8_center_cells$CellId, c("CellId", "cluster", "core", "neighbour_CellId10")]
profepi_n11 <- prof_epi_center_cells[prof_epi_center_cells$neighbour_CellId11 %in% cd8_center_cells$CellId, c("CellId", "cluster", "core", "neighbour_CellId11")]
profepi_n13 <- prof_epi_center_cells[prof_epi_center_cells$neighbour_CellId12 %in% cd8_center_cells$CellId, c("CellId", "cluster", "core", "neighbour_CellId12")]
profepi_n14 <- prof_epi_center_cells[prof_epi_center_cells$neighbour_CellId13 %in% cd8_center_cells$CellId, c("CellId", "cluster", "core", "neighbour_CellId13")]
profepi_n15 <- prof_epi_center_cells[prof_epi_center_cells$neighbour_CellId14 %in% cd8_center_cells$CellId, c("CellId", "cluster", "core", "neighbour_CellId14")]
profepi_n16 <- prof_epi_center_cells[prof_epi_center_cells$neighbour_CellId15 %in% cd8_center_cells$CellId, c("CellId", "cluster", "core", "neighbour_CellId15")]
profepi_n17 <- prof_epi_center_cells[prof_epi_center_cells$neighbour_CellId16 %in% cd8_center_cells$CellId, c("CellId", "cluster", "core", "neighbour_CellId16")]
profepi_n18 <- prof_epi_center_cells[prof_epi_center_cells$neighbour_CellId17 %in% cd8_center_cells$CellId, c("CellId", "cluster", "core", "neighbour_CellId17")]
profepi_n19 <- prof_epi_center_cells[prof_epi_center_cells$neighbour_CellId18 %in% cd8_center_cells$CellId, c("CellId", "cluster", "core", "neighbour_CellId18")]
profepi_n20 <- prof_epi_center_cells[prof_epi_center_cells$neighbour_CellId20 %in% cd8_center_cells$CellId, c("CellId", "cluster", "core", "neighbour_CellId20")]
profepi_n21 <- prof_epi_center_cells[prof_epi_center_cells$neighbour_CellId21 %in% cd8_center_cells$CellId, c("CellId", "cluster", "core", "neighbour_CellId21")]
profepi_n12 <- prof_epi_center_cells[prof_epi_center_cells$neighbour_CellId19 %in% cd8_center_cells$CellId, c("CellId", "cluster", "core", "neighbour_CellId19")]



#change the column name of 4th column

colnames(profepi_n1)[4] <- "neighbour_CellId"
colnames(profepi_n2)[4] <- "neighbour_CellId"
colnames(profepi_n3)[4] <- "neighbour_CellId"
colnames(profepi_n4)[4] <- "neighbour_CellId"
colnames(profepi_n5)[4] <- "neighbour_CellId"
colnames(profepi_n6)[4] <- "neighbour_CellId"
colnames(profepi_n7)[4] <- "neighbour_CellId"
colnames(profepi_n8)[4] <- "neighbour_CellId"
colnames(profepi_n9)[4] <- "neighbour_CellId"
colnames(profepi_n10)[4] <- "neighbour_CellId"
colnames(profepi_n11)[4] <- "neighbour_CellId"
colnames(profepi_n12)[4] <- "neighbour_CellId"
colnames(profepi_n13)[4] <- "neighbour_CellId"
colnames(profepi_n14)[4] <- "neighbour_CellId"
colnames(profepi_n15)[4] <- "neighbour_CellId"
colnames(profepi_n16)[4] <- "neighbour_CellId"
colnames(profepi_n17)[4] <- "neighbour_CellId"
colnames(profepi_n18)[4] <- "neighbour_CellId"
colnames(profepi_n19)[4] <- "neighbour_CellId"
colnames(profepi_n20)[4] <- "neighbour_CellId"
colnames(profepi_n21)[4] <- "neighbour_CellId"


profepi_n1$neighbour_CellId <- as.character(profepi_n1$neighbour_CellId)
profepi_n2$neighbour_CellId <- as.character(profepi_n2$neighbour_CellId)
profepi_n3$neighbour_CellId <- as.character(profepi_n3$neighbour_CellId)
profepi_n4$neighbour_CellId <- as.character(profepi_n4$neighbour_CellId)
profepi_n5$neighbour_CellId <- as.character(profepi_n5$neighbour_CellId)
profepi_n6$neighbour_CellId <- as.character(profepi_n6$neighbour_CellId)
profepi_n7$neighbour_CellId <- as.character(profepi_n7$neighbour_CellId)
profepi_n8$neighbour_CellId <- as.character(profepi_n8$neighbour_CellId)
profepi_n9$neighbour_CellId <- as.character(profepi_n9$neighbour_CellId)
profepi_n10$neighbour_CellId <- as.character(profepi_n10$neighbour_CellId)
profepi_n11$neighbour_CellId <- as.character(profepi_n11$neighbour_CellId)
profepi_n12$neighbour_CellId <- as.character(profepi_n12$neighbour_CellId)
profepi_n13$neighbour_CellId <- as.character(profepi_n13$neighbour_CellId)
profepi_n14$neighbour_CellId <- as.character(profepi_n14$neighbour_CellId)
profepi_n15$neighbour_CellId <- as.character(profepi_n15$neighbour_CellId)
profepi_n16$neighbour_CellId <- as.character(profepi_n16$neighbour_CellId)
profepi_n17$neighbour_CellId <- as.character(profepi_n17$neighbour_CellId)
profepi_n18$neighbour_CellId <- as.character(profepi_n18$neighbour_CellId)
profepi_n19$neighbour_CellId <- as.character(profepi_n19$neighbour_CellId)
profepi_n20$neighbour_CellId <- as.character(profepi_n20$neighbour_CellId)
profepi_n21$neighbour_CellId <- as.character(profepi_n21$neighbour_CellId)


#rbind those cells

all_prof_epi_cd8_neighbours <- bind_rows(profepi_n1,profepi_n2,profepi_n3,profepi_n4,profepi_n5,profepi_n6,profepi_n7,profepi_n8,profepi_n9,profepi_n10,profepi_n11,profepi_n12,profepi_n13,profepi_n14,profepi_n15,profepi_n16,profepi_n17,profepi_n18,profepi_n19,profepi_n20,profepi_n21)

#these are the cd8 cells as neighbbours of profepi

#i need also the prof epi cells as cell id and cluster 


#number of prof_epi_cells

n_of_neighbours_prof_epi_cd8 <- all_prof_epi_cd8_neighbours %>% group_by(CellId) %>% summarise(number=n())


#prof_epi_center_cells has all cells, take those that don't have neighbours

no_neighbours_prof_epi_cd8 <- prof_epi_center_cells[-which(prof_epi_center_cells$CellId %in% n_of_neighbours_prof_epi_cd8$CellId),]

no_neighbours_prof_epi_cd8$number <- "0"


number_of_neigbours_prof_epi_cd8 <- rbind(n_of_neighbours_prof_epi_cd8, no_neighbours_prof_epi_cd8[, c("CellId", "number")])

all_celltypes_24092020_for_neighbourhoods <- all_celltypes_24092020

all_celltypes_24092020_for_neighbourhoods$cores <- paste0("core", all_celltypes_24092020_for_neighbourhoods$cores,"_", all_celltypes_24092020_for_neighbourhoods$Cellid)


number_of_neigbours_prof_epi_cd8 <- merge(number_of_neigbours_prof_epi_cd8, all_celltypes_24092020_for_neighbourhoods[, c("Sample", "cores", "Ki67","PDL1", "HR_defect")], by.x="CellId", by.y = "cores")


number_of_neigbours_prof_epi_cd8$HR_defect <- as.character(number_of_neigbours_prof_epi_cd8$HR_defect)

number_of_neigbours_prof_epi_cd8[which(number_of_neigbours_prof_epi_cd8$HR_defect == "0"), "HR_defect"] <- "HRwt"
number_of_neigbours_prof_epi_cd8[which(number_of_neigbours_prof_epi_cd8$HR_defect == "1"), "HR_defect"] <- "BRCA1/2 mutated"


number_of_neigbours_prof_epi_cd8$combined <- paste0(number_of_neigbours_prof_epi_cd8$number, number_of_neigbours_prof_epi_cd8$HR_defect)

p3 <- ggplot(number_of_neigbours_prof_epi_cd8, aes(x=combined, y=Ki67, fill=HR_defect)) + geom_violin(positio=dodge) + geom_boxplot(position=dodge, width=0.1)+ theme+
  scale_fill_manual(values=c("#3366CC", "#CC0033")) + ggtitle("Number of CD8+ T-cells as neighbours") + stat_compare_means(method="wilcox.test", label="p.signf", comparisons=my_comparisons)
p3


#width=12, height=4.5


number_of_neigbours_prof_epi_cd8$number_of_neighbours <- number_of_neigbours_prof_epi_cd8$number

number_of_neigbours_prof_epi_cd8$number_of_neighbours[which(number_of_neigbours_prof_epi_cd8$number_of_neighbours > 0)] <- "1 or more CD8+ T-cell as neighbour"


number_of_neigbours_prof_epi_cd8$combined_variables <- paste0(number_of_neigbours_prof_epi_cd8$HR_defect,"_", number_of_neigbours_prof_epi_cd8$number_of_neighbours)





my_comparisons <- list( c("BRCA1/2 mutated_0", "BRCA1/2 mutated_1 or more CD8+ T-cell as neighbour"), c("HRwt_0", "HRwt_1 or more CD8+ T-cell as neighbour"), c("BRCA1/2 mutated_0", "HRwt_0"), c("BRCA1/2 mutated_1 or more CD8+ T-cell as neighbour", "HRwt_1 or more CD8+ T-cell as neighbour") )

p4 <- ggplot(number_of_neigbours_prof_epi_cd8, aes(x=combined_variables, y=Ki67, fill=HR_defect)) + geom_violin(position=dodge) + geom_boxplot(position=dodge, width=0.2)+ theme+
  scale_fill_manual(values=c("#3366CC", "#CC0033")) + ggtitle("Number of CD8+ T-cells as neighbours\nof proliferating epithelial cells")  + xlab("HR status") + stat_compare_means(comparisons=my_comparisons, label = "p.signif", label.y = c(11.3,11.3, 11.7, 12.1)) + ylab("Ki67 expression in proliferating epithelial cells") + theme(text = element_text(size = 15))
p4

p4 <- ggplot(number_of_neigbours_prof_epi_cd8, aes(x=combined_variables, y=Ki67, fill=HR_defect)) + geom_violin(position=dodge) + geom_boxplot(position=dodge, width=0.2)+ theme+
  scale_fill_manual(values=c("#3366CC", "#CC0033")) + ggtitle("Number of CD8+ T-cells as neighbours\nof proliferating epithelial cells")  + xlab("HR status") + stat_compare_means(comparisons=my_comparisons, label.y = c(11.3,11.3, 11.7, 12.1)) + ylab("Ki67 expression in proliferating epithelial cells") + theme(text = element_text(size = 15))
p4

#width=9, height=5.5



