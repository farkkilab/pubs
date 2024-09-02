#bulkRNAseq primary vs interval
#boxplots for S.figure 4g


#

data_bulk_prim <- read.csv("D:/sciset/bulk/immune_tmm_primary.tsv/immune_tmm_primary.tsv", sep="\t")
data_bulk_prim <- data_bulk_prim[which(data_bulk_prim$geneName == "LAG3" | data_bulk_prim$geneName == "CXCL12" | data_bulk_prim$geneName == "CXCR4" | data_bulk_prim$geneName == "CD226"| data_bulk_prim$geneName == "NECTIN2"| data_bulk_prim$geneName == "TIGIT"| data_bulk_prim$geneName == "CD96"),]

data_bulk_int <- read.csv("D:/sciset/bulk/immune_tmm_interval.tsv/immune_tmm_interval.tsv", sep="\t")
data_bulk_int <- data_bulk_int[which(data_bulk_int$geneName == "LAG3" | data_bulk_int$geneName == "CXCL12" | data_bulk_int$geneName == "CXCR4" | data_bulk_int$geneName == "CD226" | data_bulk_int$geneName == "NECTIN2" | data_bulk_int$geneName == "TIGIT" | data_bulk_int$geneName == "CD96"),]


data_bulk_prim <- t(data_bulk_prim)
data_bulk_int <- t(data_bulk_int)

data_bulk_prim <- as.data.frame(data_bulk_prim)
data_bulk_int <- as.data.frame(data_bulk_int)

colnames(data_bulk_prim) <- data_bulk_prim[1,]
data_bulk_prim <- data_bulk_prim[-1,]

colnames(data_bulk_int) <- data_bulk_int[1,]
data_bulk_int <- data_bulk_int[-1,]

data_bulk_prim[,1] <- as.numeric(data_bulk_prim[,1])
data_bulk_prim[,2] <- as.numeric(data_bulk_prim[,2])
data_bulk_prim[,3] <- as.numeric(data_bulk_prim[,3])
data_bulk_prim[,4] <- as.numeric(data_bulk_prim[,4])
data_bulk_prim[,5] <- as.numeric(data_bulk_prim[,5])
data_bulk_prim[,6] <- as.numeric(data_bulk_prim[,6])
data_bulk_prim[,7] <- as.numeric(data_bulk_prim[,7])

data_bulk_int[,1] <- as.numeric(data_bulk_int[,1])
data_bulk_int[,2] <- as.numeric(data_bulk_int[,2])
data_bulk_int[,3] <- as.numeric(data_bulk_int[,3])
data_bulk_int[,4] <- as.numeric(data_bulk_int[,4])
data_bulk_int[,5] <- as.numeric(data_bulk_int[,5])
data_bulk_int[,6] <- as.numeric(data_bulk_int[,6])
data_bulk_int[,7] <- as.numeric(data_bulk_int[,7])

data_bulk_prim$Stage <- "primary"
data_bulk_int$Stage <- "interval"

data_bulk <- rbind(data_bulk_prim, data_bulk_int)

data_bulk$patient <- rownames(data_bulk)
data_bulk$patient <- substr(data_bulk$patient, start = 1, stop = 4)
pdf("D:/Sciset/scimap/plots/boxplots/boxplots_bulk_selected_LR_updated_v2.pdf", width = 3.5, height = 4)

ggplot(data_bulk, aes(x=factor(Stage, levels = c("primary", "interval")), y=LAG3, fill=Stage)) +
  geom_boxplot(width=0.3) + geom_point() + geom_line(aes(group=patient)) + xlab("Stage") + theme + 
  scale_fill_manual(values=sample(moma.colors("OKeeffe", 2))) + stat_compare_means(paired=T)


ggplot(data_bulk, aes(x=factor(Stage, levels = c("primary", "interval")), y=CXCL12, fill=Stage)) +
  geom_boxplot(width=0.3) + geom_point() + geom_line(aes(group=patient)) + theme + xlab("Stage")+
  scale_fill_manual(values=sample(moma.colors("OKeeffe", 2))) + stat_compare_means(paired=T)

ggplot(data_bulk, aes(x=factor(Stage, levels = c("primary", "interval")), y=CXCR4, fill=Stage)) +
  geom_boxplot(width=0.3) + geom_point() + geom_line(aes(group=patient)) + theme + xlab("Stage")+
  scale_fill_manual(values=sample(moma.colors("OKeeffe", 2))) + stat_compare_means(paired=T)


ggplot(data_bulk, aes(x=factor(Stage, levels = c("primary", "interval")), y=CD226, fill=Stage)) +
  geom_boxplot(width=0.3) + geom_point() + geom_line(aes(group=patient)) + theme + xlab("Stage")+
  scale_fill_manual(values=sample(moma.colors("OKeeffe", 2))) + stat_compare_means(paired=T)

ggplot(data_bulk, aes(x=factor(Stage, levels = c("primary", "interval")), y=TIGIT, fill=Stage)) +
  geom_boxplot(width=0.3) + geom_point() + geom_line(aes(group=patient)) + theme + xlab("Stage")+
  scale_fill_manual(values=sample(moma.colors("OKeeffe", 2))) + stat_compare_means(paired=T)

ggplot(data_bulk, aes(x=factor(Stage, levels = c("primary", "interval")), y=CD96, fill=Stage)) +
  geom_boxplot(width=0.3) + geom_point() + geom_line(aes(group=patient)) + theme + xlab("Stage")+
  scale_fill_manual(values=sample(moma.colors("OKeeffe", 2))) + stat_compare_means(paired=T)

ggplot(data_bulk, aes(x=factor(Stage, levels = c("primary", "interval")), y=NECTIN2, fill=Stage)) +
  geom_boxplot(width=0.3) + geom_point() + geom_line(aes(group=patient)) + theme + xlab("Stage")+
  scale_fill_manual(values=sample(moma.colors("OKeeffe", 2))) + stat_compare_means(paired=T)

dev.off()

#
#then dotplot

data_bulk_long <- pivot_longer(data_bulk, cols=CXCR4:NECTIN2, values_to='expression', names_to='gene')
metacluster_percentages <- data_bulk_long

library(dplyr)
test <- metacluster_percentages %>%
  group_by(gene) %>%
  summarise(
    Stage_pval = wilcox.test(expression[Stage == 'interval'], 
                             expression[Stage == 'primary'], paired=TRUE)$p.value)


df2 <- metacluster_percentages %>%
  group_by(gene, Stage) %>%
  summarize(expression = mean(expression))

df2

df3 <- df2 %>%
  pivot_wider(names_from = 'Stage', values_from = 'expression')
df3

df3 <- df3 %>%
  mutate(foldChange = log2(`interval`/`primary`))

df2_bulk <- merge(test, df3)


df2_bulk <- df2_bulk[,-c(3,4)]

colnames(df2_bulk)[2] <- "pvalue"
colnames(df2_bulk)[3] <- "log_foldchange"

truncate.df <- function(df, na.cutoff=0.1, na.var="pvalue", na.var.boundary=1e-50, range.lims=c(-2, 2), range.var="log_foldchange"){
  df[which(df[,na.var] > na.cutoff), na.var] <- NA
  df[which(df[,na.var] < na.var.boundary), na.var] <- na.var.boundary
  df[,range.var] <- pmax( range.lims[1], pmin( df[,range.var], range.lims[2]))
  return(df)
}
to.plot <- truncate.df(df2_bulk) 
to.plot$comparison <- "bulk"
#row_order <- rev(c("tcycif", "scRNAseq", "bulkRNAseq"))
p <- ggplot(to.plot,  aes(x=factor(gene, levels=c("TIGIT", "CD96","CD226","NECTIN2", "LAG3", "CXCL12", "CXCR4")), y=factor(comparison))) +
  geom_point(aes(color=log_foldchange, size=pvalue)) + theme_classic() + xlab(NULL) + ylab(NULL) + labs(color="fold change (log2)") +
  scale_color_gradient2(low = "#0b71b0", high="#cc2127", mid = "gray90",  breaks=waiver()) +
  scale_size_area("p-values", trans="log10",max_size = 5.5, breaks=c(1e-6, 1e-2, 0.1), limits=c(1e-10,0.1)) +
  theme(axis.text=element_text(size=rel(1.3)), axis.text.x = element_text(angle=90, hjust=1), strip.placement = "outside",strip.background = element_blank())+
  theme(plot.title = element_text(size = 12, face = "bold"),  panel.border = element_rect(colour = "black", size=1.5, fill=NA))

print(p)

pdf("D:/sciset/revision/dotplot_bulk_LR.pdf")
p
dev.off()


