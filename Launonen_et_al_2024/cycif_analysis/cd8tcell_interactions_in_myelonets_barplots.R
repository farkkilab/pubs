#download myelonets

data = read.csv('/media/oncosys/T7/sciset/data_with_myelonets.csv')

#download interactions
iba1cd163 = read.csv('/media/oncosys/T7/sciset/IBA1.CD163_interactions.csv')
iba1cd11c = read.csv('/media/oncosys/T7/sciset/IBA1.CD11c_interactions.csv')
cd163 = read.csv('/media/oncosys/T7/sciset/CD163_interactions.csv')
cd11c = read.csv('/media/oncosys/T7/sciset/CD11c.myeloid_interactions.csv')

#combine interactions with myelonets
data = cbind(data, iba1cd163, iba1cd11c, cd163, cd11c)

#remove duplicate columns
test = data
test = test[,-c(53, 55, 57, 59)]
#then: where do interactions happen out of all interactions

#x akseli interaktoivat solutyypit
#y akseli missa interaktoi - myelonet vai muualla

test$myelonet_all <- NA
test[which(test$myelonet_macs == "myelonet_mac"), "myelonet_all"] <- "myelonet_RCN14"
test[which(test$myelonet_IBA1CD163 == "myelonet_IBA1CD163"), "myelonet_all"] <- "myelonet_RCN15"
test[which(test$myelonet_CD11c == "myelonet_CD11c"), "myelonet_all"] <- "myelonet_RCN16"
test[which(test$myelonet_mixed_immune == "myelonet_mixed_immune"), "myelonet_all"] <- "myelonet_RCN18"

library(ggplot2)
library(cowplot)

cd8t <- test[which(test$GlobalCellType2 == "CD8.T.cells"),]
p <- ggplot(cd8t[-which(cd8t$spatial_pscore_IBA1.CD163.Macrophages == "other"),],
            aes(x=factor(spatial_pscore_IBA1.CD163.Macrophages), fill=myelonet_all)) + geom_bar(position="fill", stat="count") + scale_fill_manual(values=c(  '#3e92ccff' ,  '#d8315bff', '#0a2463ff', '#fffaffff', 'seashell3')) +theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position='None', axis.title.x=element_blank()) + xlab("")

p

p1 <- ggplot(cd8t[-which(cd8t$spatial_pscore_IBA1.CD11c.Macrophages == "other"),],
            aes(x=factor(spatial_pscore_IBA1.CD11c.Macrophages), fill=myelonet_all)) + geom_bar(position="fill", stat="count") + scale_fill_manual(values=c(  '#3e92ccff' ,  '#d8315bff', '#0a2463ff', '#fffaffff', 'seashell3')) +theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position='None', axis.title.x=element_blank()) + xlab("")

p1

p2 <- ggplot(cd8t[-which(cd8t$spatial_pscore_CD163.Macrophages == "other"),],
            aes(x=factor(spatial_pscore_CD163.Macrophages), fill=myelonet_all))+ geom_bar(position="fill", stat="count") + scale_fill_manual(values=c(  '#3e92ccff' ,  '#d8315bff', '#0a2463ff', '#fffaffff', 'seashell3')) +theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position='None', axis.title.x=element_blank()) + xlab("")

p2

p3 <- ggplot(cd8t[-which(cd8t$spatial_pscore_CD11c.myeloid == "other"),],
            aes(x=factor(spatial_pscore_CD11c.myeloid), fill=myelonet_all)) + geom_bar(position="fill", stat="count") + scale_fill_manual(values=c(  '#3e92ccff' ,  '#d8315bff', '#0a2463ff', '#fffaffff', 'seashell3')) +theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position='None', axis.title.x=element_blank()) + xlab("")

p3

plot_grid(p, p1, p2, p3, ncol=4, rel_heights = 1)


pdf("/media/oncosys/T7/sciset/interactions_inside_myelonets.pdf", width=7, height=4)
plot_grid(p, p1, p2, p3, ncol=4, rel_heights = 1)
dev.off()
#3e92ccff , #d8315bff,#0a2463ff, #fffaffff, seashell3

#then stacked

p <- ggplot(cd8t[-which(cd8t$spatial_pscore_IBA1.CD163.Macrophages == "other"),],
            aes(x=factor(spatial_pscore_IBA1.CD163.Macrophages), fill=myelonet_all)) + geom_bar() + scale_fill_manual(values=c(  '#3e92ccff' ,  '#d8315bff', '#0a2463ff', '#fffaffff', 'seashell3')) +theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position='None', axis.title.x=element_blank()) + xlab("") + ylim(0, 180000)

p

p1 <- ggplot(cd8t[-which(cd8t$spatial_pscore_IBA1.CD11c.Macrophages == "other"),],
             aes(x=factor(spatial_pscore_IBA1.CD11c.Macrophages), fill=myelonet_all)) + geom_bar() + scale_fill_manual(values=c(  '#3e92ccff' ,  '#d8315bff', '#0a2463ff', '#fffaffff', 'seashell3')) +theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position='None', axis.title.x=element_blank()) + xlab("")+ ylim(0, 180000)

p1

p2 <- ggplot(cd8t[-which(cd8t$spatial_pscore_CD163.Macrophages == "other"),],
             aes(x=factor(spatial_pscore_CD163.Macrophages), fill=myelonet_all))+ geom_bar() + scale_fill_manual(values=c(  '#3e92ccff' ,  '#d8315bff', '#0a2463ff', '#fffaffff', 'seashell3')) +theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position='None', axis.title.x=element_blank()) + xlab("")+ ylim(0, 180000)

p2

p3 <- ggplot(cd8t[-which(cd8t$spatial_pscore_CD11c.myeloid == "other"),],
             aes(x=factor(spatial_pscore_CD11c.myeloid), fill=myelonet_all)) + geom_bar() + scale_fill_manual(values=c(  '#3e92ccff' ,  '#d8315bff', '#0a2463ff', '#fffaffff', 'seashell3')) +theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position='None', axis.title.x=element_blank()) + xlab("")+ ylim(0, 180000)

p3

plot_grid(p, p3, p1, p2, ncol=4, rel_heights = 1)
#iba1cd613, cd11c, iba1cd11c, cd163












