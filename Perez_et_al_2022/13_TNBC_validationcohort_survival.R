setwd("C:/Users/fernpere/HRD/TCGA_analysis/TNBC/JStaaf/")
library(ggplot2)
library(gridExtra)

#Importing other function inside the repository
extra.functions <- list.files("R/", full.names = TRUE)
sapply(extra.functions, source)
load("sysdata.rda")


################## Reading input files and editing info ######################################
patient.data <- read.table(file="PatientDataTable.csv", sep=",", header=TRUE)

Segments <- read.table(file="ascat.summary.csv", sep=",")

ploy = rep(2, nrow(Segments)) #Not really used this value
tnbc.segs <- data.frame(SampleID = Segments[,1],
                            Chromosome = paste("chr",Segments[,2],sep=""),
                            Start_position = Segments[,3],
                            End_position = Segments[,4],
                            total_cn = (Segments[,5] + Segments[,6]),
                            A_cn  = Segments[,5],
                            B_cn = Segments[,6],
                            ploidy = ploy)

#Calculating scars using new criteria
tnbc_scars <- get_HRDs(tnbc.segs, chrominfo = chrominfo_grch37, LOH_windos=c(10,30), LST_segSize=5e6, LST_mindistance=2e6)
df.tnbc_scars <- as.data.frame(tnbc_scars)
HRDscar <- df.tnbc_scars$HRDsum

#Calculating scars using previous criteria
tnbc_scars_prevdef <- get_HRDs(tnbc.segs, chrominfo = chrominfo_grch37, LOH_windos=c(15,1000), LST_segSize=10e6, LST_mindistance=3e6)
df.tnbc_scars_prevdef <- as.data.frame(tnbc_scars_prevdef)
Telli2016 <- df.tnbc_scars_prevdef$HRDsum

##Merging patient data with scar info
scars.clinical.info <- cbind(patient.data, HRDscar, Telli2016)

#Binarizing scar and BRCAness status
scars.clinical.info$BRCAness <- c("1","0")[match(scars.clinical.info$BRCA1Germline == "TRUE" |
                                  scars.clinical.info$BRCA2Germline == "TRUE" |
                                  scars.clinical.info$BRCA1_RearrangementsBiAllelic == "TRUE" |
                                  scars.clinical.info$BRCA2_RearrangementsBiAllelic == "TRUE", c('TRUE', 'FALSE'))]

scars.clinical.info$HRDscarBin <- ifelse(scars.clinical.info$HRDscar >= 53, 1, 0)
scars.clinical.info$Telli2016Bin <- ifelse(scars.clinical.info$Telli2016 >= 42, 1, 0)

#Selecting the patients for DRFI survival analysis
scars.clinical.info.selected <- scars.clinical.info[scars.clinical.info$Include_ACT_DRFI_Analysis == 1,]

surv_objectPFI <- Surv(time = scars.clinical.info.selected$Clinical.DRFI, event = scars.clinical.info.selected$Clinical.DRFIbin)
#Plotting Figure 5h, 5i, 5j
plot <- make.merge.survplots(surv_objectPFI, scars.clinical.info.selected, variables = c("BRCAness","Telli2016Bin", "HRDscarBin"),
                             break_time=2, xlabplot="DRFI time (years)", ylabplot="DRFI probability", palette = c("#34A0D3" ,"#FA5A41"))
print(plot)
ggsave(plot, filename = paste0("Kaplan-Meir_DRFI.svg"), width = 21, height = 12, units = "cm")


############################Comparing vs HRDDetect ####################################################################
scars.clinical.info.selected$HRDetectBin <- c(1,0)[match(scars.clinical.info.selected$HRDetect.prob >= 0.7, c('TRUE', 'FALSE'))]
scars.clinical.info.selected$HRDetectBin[which(scars.clinical.info.selected$HRDetect.prob > 0.2 &
                                                 scars.clinical.info.selected$HRDetect.prob < 0.7)] <- NA
scars.clinical.info.selected[which(scars.clinical.info.selected$HRDetect.groups..3.groups. == "[0.0,0.2)"),
                    c("HRDetect.groups..3.groups.", "HRDetectBin")]

#Ignoring all those patients with intermediate values
scars.clinical.info.selected2  <- scars.clinical.info.selected[!is.na(scars.clinical.info.selected$HRDetectBin),]

surv_objectPFI <- Surv(time = scars.clinical.info.selected2$Clinical.DRFI, event = scars.clinical.info.selected2$Clinical.DRFIbin)
#Plotting Supplementary Figure 5i, 5j
plot <- make.merge.survplots(surv_objectPFI, scars.clinical.info.selected2, variables = c("HRDetectBin", "HRDscarBin"),
                             break_time=1, xlabplot="DRFI time (years)", ylabplot="DRFI probability", palette = c("#34A0D3" ,"#FA5A41"))
print(plot)
ggsave(plot, filename = paste0("Kaplan-Meir_DRFI_HRDetect.svg"), width = 15, height = 12, units = "cm")


### Evaluate different cutoffs for HRDetect #######
HRDetect.max.cutoffs <- rev(seq(0.3, 0.8, by=0.1))
input.data <- scars.clinical.info.selected

p.values.cutoffs <- NULL
for (min in c(0.1, 0.15, 0.2, 0.25)){
  pvals <- sapply(HRDetect.max.cutoffs, function(x){
          input.data$HRDetectBin <- c(1,0)[match(input.data$HRDetect.prob >= x, c('TRUE', 'FALSE'))]
          #Intermediate values ignored
          input.data$HRDetectBin[which(input.data$HRDetect.prob > min &
                                         input.data$HRDetect.prob < x)] <- NA
          fits <- surv_fit(surv_objectPFI ~ HRDetectBin, data = input.data)
          x <- surv_pvalue(fits, input.data, method ="survdiff")
          return(x$pval)
  })
  p.values.df <- data.frame(minval=rep(min, length(HRDetect.max.cutoffs)),
                          max.val=HRDetect.max.cutoffs, p.values=pvals)
  p.values.cutoffs <- rbind(p.values.cutoffs, p.values.df)
}

p.values.cutoffs$minval <- as.factor(p.values.cutoffs$minval)

#Plotting Supplementary Figure 5k
p <- ggplot(data=p.values.cutoffs, aes(x=max.val, y=p.values, group=minval)) +  geom_line(aes(linetype=minval))+  geom_point(aes(shape=minval))
p <- p + geom_hline(yintercept = 0.0022, colour="blue", size=1)
p <- p + theme(axis.text=element_text(size=rel(1.2)), strip.placement = "outside",strip.background = element_blank(),
               axis.text.x=element_text(size=rel(1.2)), axis.text.y=element_text(size=rel(1.2)),
               legend.text=element_text(size=rel(1.2)), legend.title=element_text(size=rel(1.2)),
               axis.title=element_text(size=rel(1.2)),
               panel.spacing = unit(1, "lines"),
               panel.border = element_rect(linetype = "solid", fill = NA),
               panel.background = element_rect(fill = "white"),
               panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "lightgrey"),
               panel.grid.minor = element_line(size = 0.5, linetype = 'solid', colour = "lightgrey"),
               panel.grid.major.x = element_blank(),
               panel.grid.minor.x = element_blank())
p <- p + ylab("P-value") + xlab("Minimal value for HRDetect positive")
p <- p + ylim(0,0.03)
print(p)
ggsave(p, filename = paste0("HRDetect_positive.svg"), width = 12, height = 12, units = "cm")


###### Cox regressions for the different algorithms
df <- scars.clinical.info.selected2
res.cox <- Cox_regresion_variables(df, c("BRCAness", "Telli2016Bin" , "HRDetectBin", "HRDscarBin"), formula ='Surv(Clinical.DRFI, Clinical.DRFIbin)~ ')
res.cox$Hazzard_ratio <- paste0(res.cox[,c(2)], "(CI:", res.cox[,c(3)], "-", res.cox[,c(4)], ")")
res.cox$N <- c(sum(df$BRCAness == 1),
               sum(df$Telli2016Bin == 1),
               sum(df$HRDetectBin == 1),
               sum(df$HRDscarBin == 1))
res.cox$Prop <- c(sum(df$BRCAness == 1)/nrow(df),
                  sum(df$Telli2016Bin == 1)/nrow(df),
                  sum(df$HRDetectBin == 1)/nrow(df),
                  sum(df$HRDscarBin == 1)/nrow(df))
res.cox$Prop <- paste0(round(res.cox$Prop,2) * 100, "%")
res.cox$TestName <- rownames(res.cox)
res.cox.t <- res.cox[,c(10,8,9,7,6)]
colnames(res.cox.t)[c(4,5)] <- c("Hazard ratio","Pval")
res.cox.t[,1] <- c("BRCAmut/del","Telli2016", "HRDetect","tnbcHRDscar")

#Plotting table for Supplementary Figure  5h
tt2 <- ttheme_minimal(base_size = 8, padding = unit(c(1.3, 2.9), "mm"))
pdf("HazardR_scores_DRFI_HRDetect.pdf", height=10, width=15)
grid.table(res.cox.t, theme=tt2, rows=NULL)
dev.off()

coxvalues <- res.cox
coxvalues$p.value <- (-1 * log10(coxvalues$p.value))
coxvalues$attribute <- factor(rownames(coxvalues), levels = rev(rownames(coxvalues)))
coxvalues$xval <- rep(1,nrow(coxvalues))

#Plotting dot-plots for Supplementary Figure  5h
p <- ggplot(coxvalues, aes(x=xval, y=attribute))
p <- p + geom_point(aes(color=p.value, size=HR)) + theme_classic()
p <- p + scale_color_gradientn(name="-log10(p.val)", colours = c("blue4", "deepskyblue","red"),
                               values=c(0,0.90,1), guide = guide_colourbar(direction = "horizontal"))
p <- p + scale_radius(range = c(7,1), limits = c(0.25, 1.1), breaks = c(0.40, 0.50, 0.60, 0.7))
p <- p + theme(axis.text=element_text(size=rel(1)), strip.placement = "outside",strip.background = element_blank(),
               axis.text.x=element_blank(), axis.text.y=element_text(size=rel(1)),
               legend.text=element_text(size=rel(1)), legend.title=element_text(size=rel(1)),
               axis.title=element_text(size=rel(1)))
p <- p + ylab("") + xlab("")
ggsave(p, filename = "Dotpvalue_DRFI_HRDetect.svg", width = 12, height = 5.5, units = "cm")
