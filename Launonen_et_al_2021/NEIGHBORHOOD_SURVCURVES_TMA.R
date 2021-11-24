##########################################################################################################################################
################################ NEIGHBORHOOD FRACTION SURVIVAL CURVES ###################################################################
##########################################################################################################################################

#For Figure 5

library(readxl)
library(survival)
library(survminer)
library(dplyr)


ann <- read.csv("TMA_clinicaldata.csv")
ann <- all_data_TMA[, c("Identifier","Category", "PFI_time")]
ann$PFI_time[which(ann$`PFI_time`> 365)] <- "long"
ann$PFI_time[which(ann$`PFI_time`> 182 & ann$`PFI_time`< 365)] <- "short"
ann$PFI_time[which(ann$`PFI_time`< 182)] <- "short"

ann$Category <- as.character(ann$Category)
ann$Category[which(ann$`Category`== "HR")] <- "HRwt"
ann$Category[which(ann$`Category`== "BRCA1")] <- "BRCA1/2 mutated"
ann$Category[which(ann$`Category`== "BRCA2")] <- "BRCA1/2 mutated"

### GET NEIGHBORS

setwd("interaction_frequencies")
temp= list.files(pattern=".csv")
myfiles = lapply(temp, read.csv)

for (i in c(1:112)){
  
  myfiles[[i]]$core <- paste(temp)[i]
}



library(dplyr)

all_fractions <- bind_rows(myfiles[[1]], myfiles[[2]], myfiles[[3]], myfiles[[4]],myfiles[[5]], myfiles[[6]], myfiles[[7]], myfiles[[8]],
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





tma_mapping <- read_excel("TMA_mapping.xlsx")

clinical_data_TMA <- read_excel("TMA_clinicaldata.xlsx")

tma_mapping <- merge(tma_mapping, clinical_data_TMA[, c(1:2)], by.x= "Patient", by.y="Identifier")
tma_mapping <- tma_mapping[, c(1, 3, 6)]

tma_mapping <- na.omit(tma_mapping)


all_fractions$core <- gsub("^.{0,14}", "", all_fractions$core)
all_fractions$core <- substr(all_fractions$core, 1, nchar(all_fractions$core) -4)


all_fractions <- merge(all_fractions, tma_mapping, by.x = "core", by.y="TMA1", all=T)



all_fractions_survival <- merge(all_fractions, all_data_TMA, by.x="Patient", by.y="Identifier")

all_fractions_survival_BRCA <- all_fractions_survival[which(all_fractions_survival$Type == "BRCA1" | all_fractions_survival$Type == "BRCA2"),]

#proliferating epithelial center cells

all_fractions_survival_BRCA_profepi <- all_fractions_survival_BRCA[which(all_fractions_survival_BRCA$cluster == "Proliferating epithelial"),]


### CD4+T-CELL SURVIVAL CURVES

median_cd4_neighbors <- all_fractions_survival_BRCA_profepi %>% group_by(Patient) %>% summarise(mean_CD4_neighbors <- mean(CD4.x))


median_cd4_neighbors <- merge(median_cd4_neighbors, all_data_TMA, by.x="Patient", by.y="Identifier")
colnames(median_cd4_neighbors)[2] <- "CD4_neighbors"

median_all <- median_cd4_neighbors


for (i in seq_along(colnames(median_all))[c(2)]){
  
  all_data_TMA_cat <- median_all
  all_data_TMA_cat[, i] <- cut(all_data_TMA_cat[, i], breaks = c(-Inf, median(all_data_TMA_cat[, i]), Inf), labels=c("low", "high"))
  
  
  colnames(all_data_TMA_cat) <- make.names(colnames(all_data_TMA_cat))
  formula <- paste0("fit <- survfit(Surv(PFI_time,status_PFI) ~ ",colnames(all_data_TMA_cat)[i],", data = all_data_TMA_cat)")
  eval(parse(text=formula))
  
  p1 <- ggsurvplot(fit, size=1.5,data=all_data_TMA_cat, cencor=T, pval=T, pval.coord=c(1200, 0.80), 
                   title=paste("BRCA1/2mut\n",colnames(all_data_TMA_cat)[i],"and PFI"), risk.table=T, risk.table.height = 0.2, pval.size=6,
                   xlim = c(0, 2000), break.time.by=400, palette = c("deepskyblue", "sienna2"), risk.table.y.text=F,legend=c(0.15, 0.15),
                  legend.title="", xlab="PFI(days)", ylab="Survival",
                   legend.labs=c("low fraction of CD4+ T-cell neighbors", "high fraction of CD4+T-cell neighbors"), conf.int = F)+ guides(colour=guide_legend(nrow=2))
  p1 <- ggpar(p1, font.main = c(20),font.x = c(20),font.y = c(20),font.caption = c(16), font.legend = c(10),font.tickslab = c(16))
  
  print(p1)
  
}

all_patients_cox <- all_data_TMA_cat
more.than.2.years <- which(all_patients_cox$PFI_time > 800)
all_patients_cox$PFI_time[more.than.2.years] <- 800
all_patients_cox$status_PFI[more.than.2.years] <- 0


res.cox <- coxph(Surv(PFI_time, status_PFI) ~ CD4_neighbors, data =  all_patients_cox)

summary(res.cox)
cox.zph(res.cox)


### CD8+T-CELL SURVIVAL CURVES

median_cd8_neighbors <- all_fractions_survival_BRCA_profepi %>% group_by(Patient) %>% summarise(mean_CD8_neighbors <- mean(CD8.x))


median_cd8_neighbors <- merge(median_cd8_neighbors, all_data_TMA, by.x="Patient", by.y="Identifier")
colnames(median_cd8_neighbors)[2] <- "CD8_neighbors"

median_all <- median_cd8_neighbors

#0.012
median_all[which(median_all$CD8_neighbors > 0.011976), "CD8_neighbors"] <- "high"
median_all[which(median_all$CD8_neighbors < 0.011976), "CD8_neighbors"] <- "low"
all_data_TMA_cat <- median_all


for (i in seq_along(colnames(median_all))[c(2)]){
  
  #all_data_TMA_cat <- median_all
  #all_data_TMA_cat[, i] <- cut(all_data_TMA_cat[, i], breaks = c(-Inf, median(all_data_TMA_cat[, i]), Inf), labels=c("low", "high"))
  
  
  colnames(all_data_TMA_cat) <- make.names(colnames(all_data_TMA_cat))
  formula <- paste0("fit <- survfit(Surv(PFI_time,status_PFI) ~ ",colnames(all_data_TMA_cat)[i],", data = all_data_TMA_cat)")
  eval(parse(text=formula))
  
  p1 <- ggsurvplot(fit, size=1.5,data=all_data_TMA_cat, cencor=T, pval=T, pval.coord=c(1200, 0.80), 
                   title=paste("BRCA1/2mut\n",colnames(all_data_TMA_cat)[i],"and PFI"), risk.table=T, risk.table.height = 0.2, pval.size=6,
                   xlim = c(0, 2000), break.time.by=400, palette = c("sienna2", "deepskyblue"), risk.table.y.text=F,legend=c(0.15, 0.15),
                   legend.title="", xlab="PFI(days)", ylab="Survival",
                   legend.labs=c("high fraction of CD8+ T-cell neighbors", "low fraction of CD8+T-cell neighbors"), conf.int = F)+ guides(colour=guide_legend(nrow=2))
  p1 <- ggpar(p1, font.main = c(20),font.x = c(20),font.y = c(20),font.caption = c(16), font.legend = c(10),font.tickslab = c(16))
  
  print(p1)
  
}

all_patients_cox <- all_data_TMA_cat
more.than.2.years <- which(all_patients_cox$PFI_time > 800)
all_patients_cox$PFI_time[more.than.2.years] <- 800
all_patients_cox$status_PFI[more.than.2.years] <- 0


res.cox <- coxph(Surv(PFI_time, status_PFI) ~ CD4_neighbors, data =  all_patients_cox)

summary(res.cox)
cox.zph(res.cox)



