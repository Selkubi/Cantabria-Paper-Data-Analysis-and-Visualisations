setwd("C:/Users/kubilay/Nextcloud/Results/Result_1")

formula_crosstab<-read.csv("Formula_crosstab.csv", header=T, sep = ",")
dup<-formula_crosstab[duplicated(formula_crosstab$mz),1]
rows<-formula_crosstab[,1] %in% dup
formula_noRep<-formula_crosstab[!rows,]
sum(duplicated(formula_noRep$mz))

formula<-formula_noRep[,79:248]
rownames(formula)<-formula_noRep$mz
formula[is.na(formula)]<-0 #turns the NA's to 0's
formula_norm<-formula
sum_vector<-colSums(formula_norm)

for(i in 1:ncol(formula)){
  formula_norm[,i]<-(formula)[,i]/sum_vector[i]*100
}
colSums(formula_norm)

#Analyses of a season

#February PCA for temperate

info<-read.delim("Info-2.txt", header=T)
formula_norm2<-data.frame(t(formula_norm))
formula_norm2 <- cbind(rownames(formula_norm2), data.frame(formula_norm2, row.names=NULL))
names(formula_norm2)[names(formula_norm2) == "rownames(formula_norm2)"] <- "Sample.Name"
formula_norm2<-merge.data.frame(x= formula_norm2, y=info, by.x = "Sample.Name", by.y = "Sample.Name")

formula_feb<-subset(x=formula_norm2, subset=formula_norm2$Sampling_time=="18-Feb"& formula_norm2$Climate=="Temperate")
rownames(formula_feb)<-(formula_feb$Sample.Name)
formula_feb<-formula_feb[,-1]
formula_feb<-formula_feb[,-(4970:4979)]
formula_feb<-formula_feb[,colSums(formula_feb)>0]

for_PCA<-(formula_feb)
PCA_formula<-prcomp(for_PCA, scale=T)
plot(PCA_formula$x[,1:2])

PCA_formula$x <- cbind(rownames(PCA_formula$x), data.frame(PCA_formula$x, row.names=NULL))
names(PCA_formula$x)[names(PCA_formula$x) == "rownames(PCA_formula$x)"] <- "Sample.Name"
PCA_all<-merge.data.frame(x= PCA_formula$x, y=info, by.x = "Sample.Name", by.y = "Sample.Name")

summary_dsh<-data.frame(PCA_all[,1:3])
summary_dsh<-cbind(summary_dsh, PCA_all)
summary_dsh[,2:4]<-NULL

theme_nogrid <- function (base_size = 12, base_family = "") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace% 
    theme(
      panel.grid = element_blank()
    )   
}

pcpc<-ggplot(summary_dsh, aes(x=PC1, y=PC2, group=Site, col=Alteration))+ ggtitle("Formula PCA February Alteration")
pcpc<-pcpc+geom_point(aes(shape=Alteration, color=Class_alteration, size=3)) +theme_nogrid()+xlab("PC1")+ylab("PC2")+scale_colour_manual("Site",values=c("darkslategray1", "orangered", "dodgerblue3", "goldenrod1","deepskyblue2", "darkorange1"))
pcpc+theme(axis.text=element_text(size=12),axis.title=element_text(size=10,face="bold"), title=element_text(size=16), legend.text=element_text(size=14,face="bold"))

pcpc<-ggplot(summary_dsh, aes(x=PC3, y=PC4, group=Site, col=Alteration))+ ggtitle("Formula PCA February Alteration")
pcpc<-pcpc+geom_point(aes(shape=Alteration, color=Class_alteration, size=3)) +theme_nogrid()+xlab("PC3")+ylab("PC4")+scale_colour_manual("Site",values=c("darkslategray1", "orangered", "dodgerblue3", "goldenrod1","deepskyblue2", "darkorange1"))
pcpc+theme(axis.text=element_text(size=12),axis.title=element_text(size=10,face="bold"), title=element_text(size=16), legend.text=element_text(size=14,face="bold"))

summary(PCA_formula)

#Feb PCA for medditeranean
formula_feb<-subset(x=formula_norm2, subset=formula_norm2$Sampling_time=="18-Feb"& formula_norm2$Climate=="Mediterranean")
rownames(formula_feb)<-(formula_feb$Sample.Name)
formula_feb<-formula_feb[,-1]
formula_feb<-formula_feb[,-(4970:4979)]
formula_feb<-formula_feb[-19,]
formula_feb<-formula_feb[,colSums(formula_feb)>0]

for_PCA<-(formula_feb)
PCA_formula<-prcomp(for_PCA, scale=T)
plot(PCA_formula$x[,1:2])

PCA_formula$x <- cbind(rownames(PCA_formula$x), data.frame(PCA_formula$x, row.names=NULL))
names(PCA_formula$x)[names(PCA_formula$x) == "rownames(PCA_formula$x)"] <- "Sample.Name"
PCA_all<-merge.data.frame(x= PCA_formula$x, y=info, by.x = "Sample.Name", by.y = "Sample.Name")

summary_dsh<-data.frame(PCA_all[,1:3])
summary_dsh<-cbind(summary_dsh, PCA_all)
summary_dsh[,2:4]<-NULL

theme_nogrid <- function (base_size = 12, base_family = "") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace% 
    theme(
      panel.grid = element_blank()
    )   
}

pcpc<-ggplot(summary_dsh, aes(x=PC1, y=PC2, group=Site, col=Alteration))+ ggtitle("Formula PCA February Alteration")
pcpc<-pcpc+geom_point(aes(shape=Alteration, color=Class_alteration, size=3)) +theme_nogrid()+xlab("PC1")+ylab("PC2")+scale_colour_manual("Site",values=c("darkslategray1", "orangered", "dodgerblue3", "goldenrod1","deepskyblue2", "darkorange1"))
pcpc+theme(axis.text=element_text(size=12),axis.title=element_text(size=10,face="bold"), title=element_text(size=16), legend.text=element_text(size=14,face="bold"))


cor.test(summary_dsh$PC4, summary_dsh$DOC_orig.µmol.,method = "spearman")
 
#Feb total
formula_feb<-subset(x=formula_norm2, subset=formula_norm2$Sampling_time=="18-Feb")
rownames(formula_feb)<-(formula_feb$Sample.Name)
formula_feb<-formula_feb[,-1]
formula_feb<-formula_feb[,-(4970:4981)]
formula_feb<-formula_feb[-19,]
formula_feb<-formula_feb[,colSums(formula_feb)>0]


for_PCA<-(formula_feb)
PCA_formula<-prcomp(for_PCA, scale=T)
plot(PCA_formula$x[,1:2])

PCA_formula$x <- cbind(rownames(PCA_formula$x), data.frame(PCA_formula$x, row.names=NULL))
names(PCA_formula$x)[names(PCA_formula$x) == "rownames(PCA_formula$x)"] <- "Sample.Name"
PCA_all<-merge.data.frame(x= PCA_formula$x, y=info, by.x = "Sample.Name", by.y = "Sample.Name")

summary_dsh<-data.frame(PCA_all[,1:3])
summary_dsh<-cbind(summary_dsh, PCA_all)
summary_dsh[,2:4]<-NULL

theme_nogrid <- function (base_size = 12, base_family = "") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace% 
    theme(
      panel.grid = element_blank()
    )   
}

require(ggplot2)
pcpc<-ggplot(summary_dsh, aes(x=PC1, y=PC2, group=Site, col=Alteration))+ ggtitle("Formula PCA February Alteration")
pcpc<-pcpc+geom_point(aes(shape=Alteration, color=Class_alteration, size=3)) +theme_nogrid()+xlab("PC1")+ylab("PC2")+scale_colour_manual("Site",values=c("darkslategray1", "orangered", "dodgerblue3", "goldenrod1","deepskyblue2", "darkorange1"))
pcpc+theme(axis.text=element_text(size=12),axis.title=element_text(size=10,face="bold"), title=element_text(size=16), legend.text=element_text(size=14,face="bold"))

#adonis test for the climate (hydrological class) - February
formula_feb<-subset(x=formula_norm2, subset=formula_norm2$Sampling_time=="18-Feb")
rownames(formula_feb)<-(formula_feb$Sample.Name)
formula_feb<-formula_feb[,-1]

data<-formula_feb
data<-data[-19,]
groups <-(data$Climate)
data<-data[,-(4970:4981)]
data<-data[,colSums(data)>0]

require(vegan)
dis <- vegdist(data, method="bray")
mod <- betadisper(dis, groups)
anova(mod) #p value not significant. so the groups are homogenous
TukeyHSD(mod)
plot(mod)
adonis(data ~ groups, perm=999)

#adonis test for the climate (alteration)- February

formula_feb<-subset(x=formula_norm2, subset=formula_norm2$Sampling_time=="18-Feb")
rownames(formula_feb)<-(formula_feb$Sample.Name)
formula_feb<-formula_feb[,-1]

data<-formula_feb
data<-data[-19,]
groups <-(data$Class_alteration)
data<-data[,-(4970:4981)]
data<-data[,colSums(data)>0]

require(vegan)
dis <- vegdist(data, method="bray")
mod <- betadisper(dis, groups)
anova(mod) #p value not significant. so the groups are homogenous
TukeyHSD(mod)
plot(mod, main="Bray-Curtis distances between samples", sub="February")
adonis(data ~ groups, perm=999)

#April
info<-read.delim("Info-2.txt", header=T)
formula_norm2<-data.frame(t(formula_norm))
formula_norm2 <- cbind(rownames(formula_norm2), data.frame(formula_norm2, row.names=NULL))
names(formula_norm2)[names(formula_norm2) == "rownames(formula_norm2)"] <- "Sample.Name"
formula_norm2<-merge.data.frame(x= formula_norm2, y=info, by.x = "Sample.Name", by.y = "Sample.Name")

formula_apr<-subset(x=formula_norm2, subset=formula_norm2$Sampling_time=="18-Apr")
rownames(formula_apr)<-(formula_apr$Sample.Name)
formula_apr<-formula_apr[,-1]
formula_apr<-formula_apr[,-(4970:4981)]
formula_apr<-formula_apr[-c(4,26,27),]
formula_apr<-formula_apr[,colSums(formula_apr)>0]



for_PCA<-(formula_apr)
PCA_formula<-prcomp(for_PCA, scale=T)
plot(PCA_formula$x[,1:2])

PCA_formula$x <- cbind(rownames(PCA_formula$x), data.frame(PCA_formula$x, row.names=NULL))
names(PCA_formula$x)[names(PCA_formula$x) == "rownames(PCA_formula$x)"] <- "Sample.Name"
PCA_all<-merge.data.frame(x= PCA_formula$x, y=info, by.x = "Sample.Name", by.y = "Sample.Name")

summary_dsh<-data.frame(PCA_all[,1:3])
summary_dsh<-cbind(summary_dsh, PCA_all)
summary_dsh[,2:4]<-NULL

require(ggplot2)
theme_nogrid <- function (base_size = 12, base_family = "") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace% 
    theme(
      panel.grid = element_blank()
    )   
}

pcpc<-ggplot(summary_dsh, aes(x=PC1, y=PC2, group=Site, col=Alteration))+ ggtitle("Formula PCA February Alteration")
pcpc<-pcpc+geom_point(aes(shape=Alteration, color=Class_alteration, size=3)) +theme_nogrid()+xlab("PC1")+ylab("PC2")+scale_colour_manual("Site",values=c("darkslategray1", "orangered", "dodgerblue3", "goldenrod1","deepskyblue2", "darkorange1"))
pcpc+theme(axis.text=element_text(size=12),axis.title=element_text(size=10,face="bold"), title=element_text(size=16), legend.text=element_text(size=14,face="bold"))

formula_apr<-subset(x=formula_norm2, subset=formula_norm2$Sampling_time=="18-Apr")
rownames(formula_apr)<-(formula_apr$Sample.Name)
formula_apr<-formula_apr[,-1]

data<-formula_apr
data<-data[-c(4,26,27),]
groups <-(data$Climate)
data<-data[,-(4970:4981)]
data<-data[,colSums(data)>0]

require(vegan)
dis <- vegdist(data, method="bray")
mod <- betadisper(dis, groups)
anova(mod) #p value not significant. so the groups are homogenous
TukeyHSD(mod)
plot(mod)
adonis(data ~ groups, perm=999)


#October
formula_oct<-subset(x=formula_norm2, subset=formula_norm2$Sampling_time=="17-Oct")
rownames(formula_oct)<-(formula_oct$Sample.Name)
formula_oct<-formula_oct[,-1]
formula_oct<-formula_oct[,-(4970:4979)]
formula_oct<-formula_oct[,colSums(formula_oct)>0]


for_PCA<-(formula_oct)
PCA_formula<-prcomp(for_PCA, scale=T)
plot(PCA_formula$x[,1:2])

PCA_formula$x <- cbind(rownames(PCA_formula$x), data.frame(PCA_formula$x, row.names=NULL))
names(PCA_formula$x)[names(PCA_formula$x) == "rownames(PCA_formula$x)"] <- "Sample.Name"
PCA_all<-merge.data.frame(x= PCA_formula$x, y=info, by.x = "Sample.Name", by.y = "Sample.Name")

summary_dsh<-data.frame(PCA_all[,1:3])
summary_dsh<-cbind(summary_dsh, PCA_all)
summary_dsh[,2:4]<-NULL

theme_nogrid <- function (base_size = 12, base_family = "") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace% 
    theme(
      panel.grid = element_blank()
    )   
}

pcpc<-ggplot(summary_dsh, aes(x=PC1, y=PC2, group=Site, col=Climate))+ ggtitle("Formula PCA October Climate")
pcpc<-pcpc+geom_point(data=summary_dsh,aes(x=PC1, y=PC2), size=5)+
  theme_nogrid()+xlab("PC1")+ylab("PC2")+scale_colour_manual("Site",values=c("deepskyblue2", "darkorange1"))
pcpc+theme(axis.text=element_text(size=12),axis.title=element_text(size=10,face="bold"), title=element_text(size=16), legend.text=element_text(size=14,face="bold"))


pcpc<-ggplot(summary_dsh, aes(x=PC1, y=PC2, group=Site, col=Alteration))+ ggtitle("Formula PCA October Alteration")
pcpc<-pcpc+geom_text(data=summary_dsh,aes(x=PC1, y=PC2, label=Site), size=5, position=position_jitter(width=0, height=3))+theme_nogrid()+xlab("PC1")+ylab("PC2")+scale_colour_manual("Site",values=c("deepskyblue2", "darkorange1"))
pcpc+theme(axis.text=element_text(size=12),axis.title=element_text(size=10,face="bold"), title=element_text(size=16), legend.text=element_text(size=14,face="bold"))

#December
formula_dec<-subset(x=formula_norm2, subset=formula_norm2$Sampling_time=="17-Dec")
rownames(formula_dec)<-(formula_dec$Sample.Name)
formula_dec<-formula_dec[,-1]
formula_dec<-formula_dec[,-(4970:4979)]
formula_dec<-formula_dec[,colSums(formula_dec)>0]


for_PCA<-(formula_dec)
PCA_formula<-prcomp(for_PCA, scale=T)
plot(PCA_formula$x[,1:2])

PCA_formula$x <- cbind(rownames(PCA_formula$x), data.frame(PCA_formula$x, row.names=NULL))
names(PCA_formula$x)[names(PCA_formula$x) == "rownames(PCA_formula$x)"] <- "Sample.Name"
PCA_all<-merge.data.frame(x= PCA_formula$x, y=info, by.x = "Sample.Name", by.y = "Sample.Name")

summary_dsh<-data.frame(PCA_all[,1:3])
summary_dsh<-cbind(summary_dsh, PCA_all)
summary_dsh[,2:4]<-NULL

theme_nogrid <- function (base_size = 12, base_family = "") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace% 
    theme(
      panel.grid = element_blank()
    )   
}

pcpc<-ggplot(summary_dsh, aes(x=PC1, y=PC2, group=Site, col=Climate))+ ggtitle("Formula PCA December Climate")
pcpc<-pcpc+geom_text(data=summary_dsh,aes(x=PC1, y=PC2, label=Site), size=5, position=position_jitter(width=0, height=3))+theme_nogrid()+xlab("PC1")+ylab("PC2")+scale_colour_manual("Site",values=c("deepskyblue2", "darkorange1"))
pcpc+theme(axis.text=element_text(size=12),axis.title=element_text(size=10,face="bold"), title=element_text(size=16), legend.text=element_text(size=14,face="bold"))


pcpc<-ggplot(summary_dsh, aes(x=PC1, y=PC2, group=Site, col=Alteration))+ ggtitle("Formula PCA December Alteration")
pcpc<-pcpc+geom_text(data=summary_dsh,aes(x=PC1, y=PC2, label=Site), size=5, position=position_jitter(width=0, height=3))+theme_nogrid()+xlab("PC1")+ylab("PC2")+scale_colour_manual("Site",values=c("deepskyblue2", "darkorange1"))
pcpc+theme(axis.text=element_text(size=12),axis.title=element_text(size=10,face="bold"), title=element_text(size=16), legend.text=element_text(size=14,face="bold"))

#February
formula_feb<-subset(x=formula_norm2, subset=formula_norm2$Sampling_time=="18-Feb"& formula_norm2$Climate=="Temperate")
rownames(formula_feb)<-(formula_feb$Sample.Name)
formula_feb<-formula_feb[,-1]
formula_feb<-formula_feb[,-(4970:4979)]
formula_feb<-formula_feb[-19,]
formula_feb<-formula_feb[,colSums(formula_feb)>0]

for_PCA<-(formula_feb)
PCA_formula<-prcomp(for_PCA, scale=T)
plot(PCA_formula$x[,1:2])

PCA_formula$x <- cbind(rownames(PCA_formula$x), data.frame(PCA_formula$x, row.names=NULL))
names(PCA_formula$x)[names(PCA_formula$x) == "rownames(PCA_formula$x)"] <- "Sample.Name"
PCA_all<-merge.data.frame(x= PCA_formula$x, y=info, by.x = "Sample.Name", by.y = "Sample.Name")

summary_dsh<-data.frame(PCA_all[,1:3])
summary_dsh<-cbind(summary_dsh, PCA_all)
summary_dsh[,2:4]<-NULL

theme_nogrid <- function (base_size = 12, base_family = "") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace% 
    theme(
      panel.grid = element_blank()
    )   
}
summary(PCA_formula)


cor.test(summary_dsh$PC4, summary_dsh$DOC_orig.µmol.,method = "spearman")

#May
formula_may<-subset(x=formula_norm2, subset=formula_norm2$Sampling_time=="18-May")
rownames(formula_may)<-(formula_may$Sample.Name)
formula_may<-formula_may[,-1]
formula_may<-formula_may[,-(4970:4979)]
formula_may<-formula_may[,colSums(formula_may)>0]


for_PCA<-(formula_may)
PCA_formula<-prcomp(for_PCA, scale=T)
plot(PCA_formula$x[,1:2])

PCA_formula$x <- cbind(rownames(PCA_formula$x), data.frame(PCA_formula$x, row.names=NULL))
names(PCA_formula$x)[names(PCA_formula$x) == "rownames(PCA_formula$x)"] <- "Sample.Name"
PCA_all<-merge.data.frame(x= PCA_formula$x, y=info, by.x = "Sample.Name", by.y = "Sample.Name")

summary_dsh<-data.frame(PCA_all[,1:3])
summary_dsh<-cbind(summary_dsh, PCA_all)
summary_dsh[,2:4]<-NULL

theme_nogrid <- function (base_size = 12, base_family = "") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace% 
    theme(
      panel.grid = element_blank()
    )   
}

pcpc<-ggplot(summary_dsh, aes(x=PC1, y=PC2, group=Site, col=Climate))+ ggtitle("Formula PCA May Climate")
pcpc<-pcpc+geom_text(data=summary_dsh,aes(x=PC1, y=PC2, label=Site), size=5, position=position_jitter(width=0, height=3))+theme_nogrid()+xlab("PC1")+ylab("PC2")+scale_colour_manual("Site",values=c("deepskyblue2", "darkorange1"))
pcpc+theme(axis.text=element_text(size=12),axis.title=element_text(size=10,face="bold"), title=element_text(size=16), legend.text=element_text(size=14,face="bold"))


pcpc<-ggplot(summary_dsh, aes(x=PC1, y=PC2, group=Site, col=Alteration))+ ggtitle("Formula PCA May Alteration")
pcpc<-pcpc+geom_text(data=summary_dsh,aes(x=PC1, y=PC2, label=Site), size=5, position=position_jitter(width=0, height=3))+theme_nogrid()+xlab("PC1")+ylab("PC2")+scale_colour_manual("Site",values=c("deepskyblue2", "darkorange1"))
pcpc+theme(axis.text=element_text(size=12),axis.title=element_text(size=10,face="bold"), title=element_text(size=16), legend.text=element_text(size=14,face="bold"))

#August
formula_aug<-subset(x=formula_norm2, subset=formula_norm2$Sampling_time=="18-Aug")
rownames(formula_aug)<-(formula_aug$Sample.Name)
formula_aug<-formula_aug[,-1]
formula_aug<-formula_aug[,-(4970:4979)]
formula_aug<-formula_aug[,colSums(formula_aug)>0]


for_PCA<-(formula_aug)
PCA_formula<-prcomp(for_PCA, scale=T)
plot(PCA_formula$x[,1:2])

PCA_formula$x <- cbind(rownames(PCA_formula$x), data.frame(PCA_formula$x, row.names=NULL))
names(PCA_formula$x)[names(PCA_formula$x) == "rownames(PCA_formula$x)"] <- "Sample.Name"
PCA_all<-merge.data.frame(x= PCA_formula$x, y=info, by.x = "Sample.Name", by.y = "Sample.Name")

summary_dsh<-data.frame(PCA_all[,1:3])
summary_dsh<-cbind(summary_dsh, PCA_all)
summary_dsh[,2:4]<-NULL

theme_nogrid <- function (base_size = 12, base_family = "") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace% 
    theme(
      panel.grid = element_blank()
    )   
}

pcpc<-ggplot(summary_dsh, aes(x=PC1, y=PC2, group=Site, col=Climate))+ ggtitle("Formula PCA August Climate")
pcpc<-pcpc+geom_text(data=summary_dsh,aes(x=PC1, y=PC2, label=Site), size=5, position=position_jitter(width=0, height=3))+theme_nogrid()+xlab("PC1")+ylab("PC2")+scale_colour_manual("Site",values=c("deepskyblue2", "darkorange1"))
pcpc+theme(axis.text=element_text(size=12),axis.title=element_text(size=10,face="bold"), title=element_text(size=16), legend.text=element_text(size=14,face="bold"))


pcpc<-ggplot(summary_dsh, aes(x=PC1, y=PC2, group=Site, col=Alteration))+ ggtitle("Formula PCA August Alteration")
pcpc<-pcpc+geom_text(data=summary_dsh,aes(x=PC1, y=PC2, label=Site), size=5, position=position_jitter(width=0, height=3))+theme_nogrid()+xlab("PC1")+ylab("PC2")+scale_colour_manual("Site",values=c("deepskyblue2", "darkorange1"))
pcpc+theme(axis.text=element_text(size=12),axis.title=element_text(size=10,face="bold"), title=element_text(size=16), legend.text=element_text(size=14,face="bold"))




