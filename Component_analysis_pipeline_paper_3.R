library(vegan)
library(geosphere)
library(gridExtra)
library(data.table)
library(ggplot2)
library(lubridate)
library(rcompanion)

####Getting the data####
sample_loadings <- as.data.table(read.csv("denormalised_sample_loadings.csv"))
sample_names <- as.data.table(read.csv("sample_name_site_match.csv"))
site_info <- as.data.table(read.csv("site_summary.csv"))

oct1<- as.data.table(cbind(read.table("abs_data/oct.txt", sep=",", header=TRUE), campaign="oct"))
feb1<- as.data.table(cbind(read.table("abs_data/feb.txt", sep=",", header=TRUE), campaign="feb"))
apr1<- as.data.table(cbind(read.table("abs_data/apr.txt", sep=",", header=TRUE), campaign="apr"))
may1<- as.data.table(cbind(read.table("abs_data/may.txt", sep=",", header=TRUE), campaign="may"))
aug1<- as.data.table(cbind(read.table("abs_data/aug.txt", sep=",", header=TRUE), campaign="aug"))
dec1<- as.data.table(cbind(read.table("abs_data/dec.txt", sep=",", header=TRUE), campaign="dec"))

oct1[,"site"] <- substr(oct1$site, 1, nchar(oct1$site)-4)
apr1[,"site"] <- substr(apr1$site, 1, nchar(apr1$site)-4)

abs_data = rbind(oct1,feb1, apr1, may1, aug1,dec1)
abs_data = subset(abs_data, subset = !grepl("^Blank", abs_data$site))

abs_data[abs_data[, HIX == "Inf"]]$HIX = NA
abs_data[abs_data[, HIX > 100]]$HIX = NA
abs_data[abs_data[, E2.to.E3 < 0]]$E2.to.E3 = NA
abs_data[abs_data[, E2.to.E3 > 30]]$E2.to.E3 = NA
abs_data[abs_data[, E4.to.E6 < 0]]$E4.to.E6 = NA


data_means = aggregate(x = abs_data[, c("FIX", "FIX2", "HIX", "HIX2", "beta.alpha","DecAbsCoeff254", "DecAbsCoeff420", "DecAbsCoeff430", "DecAbsCoeff436", "E2.to.E3", "E4.to.E6", 
                                    "slope_classic", "slope_lm", "slope_short_Loiselle", "slope_short_Helms", "SR_Loiselle", "SR_Helms")], 
                     by = abs_data[,c("site", "campaign")], FUN = mean, na.rm = T, trim = 0.25)

meta_sum_optical = as.data.table(cbind(read.csv("meta_file_2.csv"),data_means), key = c("site","campaign"))

# Converting the sample names without the index numbers - wouldn't need it if data were named with 000
sample_loadings[, c("campaign", "run", "sample2", "num") := tstrsplit(sample, "_")]
sample_loadings[,c("sample_code") := paste0(sample_loadings$campaign,"_", sample_loadings$run, "_", sample_loadings$sample2)]
setdiff(sample_names$sample_code, sample_loadings$sample_code) #Check why these are not matching. In my case, they are the samples removed during parafac because of weird looking eems.So it makes sense to remove them here too since there might be a contamination affecting all

##### Binding the data frames ####
meta_component_data = merge(sample_loadings, sample_names, by.x=c("sample_code", "campaign"), by.y=c("sample_code", "campaign"))
meta_component_data = meta_component_data[meta_component_data$site!="Carrion",]
meta_sum_optical = meta_sum_optical[meta_sum_optical$site!="Carrion",]
sample_names = sample_names[sample_names$site!="Carrion",]

##### Prepping for the DOC mean and variance plots ####
data_sum = aggregate(meta_component_data[,c("Comp.1","Comp.2","Comp.3","Comp.4","Comp.5","Comp.6","Comp.7","Comp.8")], by = meta_component_data[, c("campaign", "site")], FUN = mean, na.rm = T)
data_sum = merge(site_info, data_sum, by.x = c("site"), by.y = c("site"))
data_sum = merge(data_sum, meta_sum_optical[, -c(14, 15)], by = c("site","campaign"))
data_sum[, "SUVA254" := DecAbsCoeff254/NPOC]
data_sum[data_sum[, SUVA254 > 6]]$SUVA254 = NA #This is definitely a weird point. Sama sample as below. Seems like an outlier. Changes the PLSR table of now
data_sum[data_sum[, SR_Loiselle > 1.7]]$SR_Loiselle = NA

ggplot(data_sum)+
  geom_point(aes(E2.to.E3, site))

# Here I replace the under the detection limit values
below_detection_limit = function(data) {
  data = gsub(x = data, pattern = c("<0,1", "<0,01", "<0,02", "<0,06"), replacement = c("0.05","0.005", "0.01", "0.03") ) #Or replacement is c("0.1","0.01", "0.02", "0.06") for abslute detection limits
  return(as.numeric(data))
}

data_sum[, c("HMWS_C", "humic_like_substance_C", "LMWS_C", "HMWS_N", "HS", "SUVA_HS", "SUVA_ges")] = lapply(data_sum[,c("HMWS_C", "humic_like_substance_C", "LMWS_C", "HMWS_N", "HS", "SUVA_HS", "SUVA_ges")], FUN = below_detection_limit)

data_sum2 = data_sum[,c("campaign", "site", "groups.x","Class","alteration", "Comp.1", "Comp.2",  "Comp.3", "Comp.4", "Comp.5", "Comp.6", "Comp.7", "Comp.8", "BDOC", "NPOC", "SUVA254",
                         "FIX", "FIX2", "HIX2", "beta.alpha", "slope_classic", "slope_lm", "slope_short_Helms", "slope_short_Loiselle",
                         "SR_Helms", "SR_Loiselle", "E2.to.E3", "E4.to.E6", "DecAbsCoeff254")]
# DOC mean and variance plots and tests
theme_pca <- theme_bw() +
  theme(
    text = element_text(size=11),
    axis.title = element_text(size=11),
    axis.text = element_text(size=11, color="black"), 
    legend.position = "none", 
    strip.placement ="outside", 
    strip.background = element_blank(), 
    strip.text = element_blank(),
    axis.line =  element_line(color="black"),
    panel.border = element_blank(),
    axis.line.y.right = element_line(color="white"),
    panel.grid = element_blank())

CV = function(x){
  coeff = sd(x, na.rm = T)/mean(x, na.rm = T)
  print(coeff)
}

DOC_sum = merge(data_sum2[, .(mean_BDOC = mean(BDOC, na.rm = TRUE), mean_NPOC = mean(NPOC, na.rm = TRUE), 
                                    var_BDOC = CV(BDOC), var_NPOC = CV(NPOC)), by = .(site)], 
                                   site_info[site != "Carrion",c("site", "Class", "alteration", "groups")])

group_means = DOC_sum[, .(means = mean(mean_NPOC, na.rm = TRUE) , sd_NPOC = sd(mean_NPOC)), by = .(groups)]
group_CVs = DOC_sum[, .(means = mean(var_NPOC, na.rm = TRUE) , sd_NPOC = sd(var_NPOC)), by = .(groups)]



DOC_sum_DOC = melt(data_sum, id.vars = c("site", "campaign", "groups.x", "Class"), measure.vars = c("NPOC"))
DOC_sum_DOC$campaign = factor(x = DOC_sum_DOC$campaign, levels = c("oct", "dec", "feb", "apr", "may", "aug"))


DOC_sum_melt = melt(DOC_sum, id.vars = c("site", "Class", "groups", "alteration"), measure.vars = c("mean_NPOC", "var_NPOC"))

dir.create("plots") 

data_sum$groups.x=factor(data_sum$groups.x, levels=c("TempNat", "TempAlt", "MedNat", "MedAlt"))

mean_NPOC = ggplot(data_sum, aes(x=groups.x, y=(NPOC))) +
  theme(axis.line.x =  element_line(color="black"), 
        axis.line.y =  element_line(),
        panel.grid = element_blank(),
        panel.background = element_blank())+
  geom_boxplot(mapping=aes(fill=groups.x, group=reorder(site, NPOC, median)), varwidth = T, width=1, lwd=0.2, outlier.size = 0.5)+
  scale_shape_manual(values=c(23,22,25,24))+
  scale_fill_manual(values=c("#B4DCED", '#6996D1','#F5CB7D','#F09E41'))+
  scale_x_discrete(limits = c("TempNat", "TempAlt", "MedNat","MedAlt"),
                   labels = c("nA", "aA", "nM", "aM"))+
  theme(axis.text = element_text(color="black", size=11), axis.title.x=element_blank(),legend.position = "none")+
  ylab("DOC mg C/L")


pdf('plots/mean_NPOC.pdf', width = 4, height = 4)
plot(mean_NPOC)
dev.off()

monthly_means =  data_sum[, .(monthly_mean = mean(NPOC, na.rm = TRUE)), by = .(campaign, groups.x)]

mediterranean_DOC = ggplot(data_sum[Class == 'Mediterranean'], aes(x = campaign, y = NPOC)) +
  geom_line(data = monthly_means[groups.x == 'MedAlt' | groups.x == 'MedNat'], 
            aes(x=campaign, y = monthly_mean, group = groups.x, color = groups.x), lwd = 1.5) +
  geom_line(aes(group = site, color = groups.x), lwd = 0.5) +
  geom_point(aes(group = site, shape = groups.x, fill = groups.x), size = 2, color = 'black') +
  scale_shape_manual(values=c(25, 24)) +
  scale_color_manual(values=c('#F5CB7D','#F09E41')) +
  scale_fill_manual(values=c('#F5CB7D','#F09E41')) +
  scale_x_discrete(limits = c("feb", "may", "oct","apr", "aug", "dec"),
                   labels = c("Feb", "May", "Oct", "Apr", "Aug", "Dec")) +
  theme_pca+
  theme(axis.title.x = element_blank())+
  ylab("DOC mg C/L")

pdf('plots/mediterranean_DOC.pdf', width = 2.8, height = 2.8)
plot(mediterranean_DOC)
dev.off()


temperate_DOC = ggplot(data_sum[Class == 'Temperate'], aes(x = campaign, y = NPOC)) +
  geom_line(data = monthly_means[groups.x == 'TempAlt' | groups.x == 'TempNat'], 
            aes(x = campaign, y = monthly_mean, group = groups.x, color = groups.x), lwd = 1.5) +
  geom_line(aes(group = site, color = groups.x), lwd = 0.5) +
  geom_point(aes(group = site, shape = groups.x, fill = groups.x), size = 2, color = 'black') +
  scale_shape_manual(values=c(23, 22)) +
  scale_color_manual(values=c("#B4DCED", '#6996D1')) +
  scale_fill_manual(values=c("#B4DCED", '#6996D1')) +
  scale_x_discrete(limits = c("feb", "may", "oct","apr", "aug", "dec"),
                   labels = c("Feb", "May", "Oct", "Apr", "Aug", "Dec")) +
  theme_pca+
  theme(axis.title.x = element_blank())+
  ylab("DOC mg C/L")

pdf('plots/temperate_DOC.pdf', width = 2.8, height = 2.8)
plot(temperate_DOC)
dev.off()

data_sum$C_tot=data_sum[,Comp.1+Comp.2+Comp.3+Comp.4+Comp.5+Comp.6+Comp.7+Comp.8]

#### Pipeline 2:1-way flow regime strategy####
#Pipeline 2, aspect i - annaul mean analyis

shapiro.test((DOC_sum$mean_NPOC))
hist((DOC_sum$mean_NPOC))
boxplot(log((DOC_sum$mean_NPOC)))

bartlett.test(log(DOC_sum$mean_NPOC)~DOC_sum$groups)
oneway.test(log(DOC_sum$mean_NPOC)~DOC_sum$groups, var.equal = T)

var.test(log(mean_NPOC)~alteration, data = DOC_sum[Class == "Mediterranean"])
t.test(log(mean_NPOC)~alteration, data = DOC_sum[Class == "Mediterranean"], var.equal = T)

var.test(log(mean_NPOC)~alteration, data = DOC_sum[Class == "Temperate"])
t.test(log(mean_NPOC)~alteration, data = DOC_sum[Class == "Temperate"], var.equal = T)

var.test(log(mean_NPOC)~Class, data = DOC_sum[alteration == "Natural"])
t.test(log(mean_NPOC)~Class, data = DOC_sum[alteration == "Natural"], var.equal = T)

p.adjust(c( 0.2135,  0.02854, 0.08353), method = "bonferroni", n = 3)
m = pairwise.t.test((DOC_sum$mean_NPOC), DOC_sum$groups, p.adjust.method = "bonferroni", pool.sd = T)
multcompView::multcompLetters(fullPTable(m$p.value))


#Pipeline 2, aspect ii - VC analysis 
shapiro.test(log(DOC_sum$var_NPOC))
hist(log(DOC_sum$var_NPOC))

bartlett.test(log(DOC_sum$var_NPOC)~DOC_sum$groups)
oneway.test(log(DOC_sum$var_NPOC) ~ DOC_sum$groups, var.equal = F)

var.test(var_NPOC~alteration, data = DOC_sum[Class == "Mediterranean"])
t.test((var_NPOC)~alteration, data = DOC_sum[Class == "Mediterranean"], var.equal = T)

var.test(var_NPOC~alteration, data = DOC_sum[Class == "Temperate"])
t.test(var_NPOC~alteration, data = DOC_sum[Class == "Temperate"], var.equal = T)

var.test(var_NPOC~Class, data = DOC_sum[alteration == "Natural"])
t.test(var_NPOC~Class, data = DOC_sum[alteration == "Natural"], var.equal = F)

p.adjust(c( 0.01134,  0.5621, 0.5056), method = "bonferroni", n = 3)
m = pairwise.t.test(DOC_sum$var_NPOC, DOC_sum$groups, p.adjust.method = "bonferroni", pool.sd = F)
multcompView::multcompLetters(fullPTable(m$p.value),  threshold = 0.05)

#multcompLetters(m$p.value)

#Pipeline 2. aspect iii - variance tests
var.test(mean_NPOC~alteration, data = DOC_sum[Class == "Mediterranean"])
var.test(mean_NPOC~alteration, data = DOC_sum[Class == "Temperate"])
var.test(mean_NPOC~Class, data = DOC_sum[alteration =="Natural"])

var.test(var_NPOC~alteration, data = DOC_sum[Class == "Mediterranean"])
var.test(var_NPOC~alteration, data = DOC_sum[Class == "Temperate"])
var.test(var_NPOC~Class, data = DOC_sum[alteration == "Natural"])
#### End of DOC mean and variance 
#### DOM quality plots####
data_summary = data_sum[, .(FI2 = (FIX), HI2= (HIX2), 
                          beta.alpha2 = (beta.alpha ), SUVA254_2 = (SUVA254),
                          SR = (SR_Loiselle ), E2toE3 = (E2.to.E3 ),
                          C_HMWS = (HMWS_C/CDOC ), C_LMWS = (LMWS_C/CDOC ),
                          C_HS = (humic_like_substance_C/CDOC ),
                          C_terrestrial = ((Comp.1 + Comp.2 + Comp.3 + Comp.4 + Comp.5) / C_tot),
                          C_protein = ((Comp.6 + Comp.7 + Comp.8) / C_tot)),
                          by = .(site, alteration, Class,  groups.x, campaign)]

#data_summary[data_summary[, C_HMWS > 0.3]]$C_HMWS = NA
#data_summary[data_summary[, C_HS > 1.5]]$C_HS = NA

data_sum$C_tot=data_sum[,Comp.1+Comp.2+Comp.3+Comp.4+Comp.5+Comp.6+Comp.7+Comp.8]

#### Statistical tests for the quality parameters ####
DOM_averages = data_sum[,.(FI2 = mean(FIX, na.rm = T), HI2 = mean(HIX2, na.rm = T), 
                         beta.alpha2= mean(beta.alpha, na.rm = T), SUVA254_2 = mean(SUVA254, na.rm = T),
                         SR = mean(SR_Loiselle, na.rm = T), E2toE3 = mean(E2.to.E3, na.rm = T), 
                         C_humic = mean((Comp.1 + Comp.2 + Comp.3 + Comp.4 + Comp.5) / C_tot, na.rm = T),
                         C_protein = mean((Comp.6 + Comp.7 + Comp.8) / C_tot, na.rm = T), 
                         DOC = mean(NPOC, na.rm = T), 
                         percent_HMWS_C = mean(HMWS_C / CDOC, na.rm = T),
                         percent_Humic_C = mean(humic_like_substance_C / CDOC, na.rm = T),
                         percent_LMWS_C = mean(LMWS_C / CDOC, na.rm = T, trim = 5)),
                      by = .(site, groups.x, Class, alteration)]

shapiro.test((DOM_averages$beta.alpha2))
hist((DOM_averages$beta.alpha2))

bartlett.test((percent_Humic_C) ~groups.x, data = DOM_averages)
oneway.test((percent_Humic_C) ~groups.x, var.equal = T, data = DOM_averages)

var.test((percent_Humic_C) ~alteration, data = DOM_averages[Class == "Mediterranean"])
t.test((percent_Humic_C) ~alteration, data = DOM_averages[Class == "Mediterranean"], var.equal = T)

var.test((percent_Humic_C) ~alteration, data = DOM_averages[Class == "Temperate"])
t.test((percent_Humic_C)~alteration, data = DOM_averages[Class == "Temperate"], var.equal = T)

var.test((percent_Humic_C) ~Class, data = DOM_averages[alteration == "Natural"])
t.test((percent_Humic_C) ~Class, data = DOM_averages[alteration == "Natural"], var.equal = T)

p.adjust(c(0.02067, 0.9364, 0.1284), method = "bonferroni", n = 3)

m = pairwise.t.test(log(DOM_averages$percent_HMWS_C), DOM_averages$groups.x, p.adjust.method = "bonferroni", pool.sd = F, paired = F)
multcompView::multcompLetters(fullPTable(m$p.value))


#CV of PC2 aspect (ii)
DOM_CV_averages = data_sum[,.(FI2_cv = CV(FIX), HI2_cv = CV(HIX2), 
                            beta.alpha2_cv = CV(beta.alpha), SUVA254_2_cv = CV(SUVA254),
                            SR_cv = CV(SR_Loiselle), E2toE3_cv = CV(E2.to.E3), 
                            C_humic_cv = CV((Comp.1 + Comp.2 + Comp.3 + Comp.4 + Comp.5) / C_tot),
                            C_protein_cv = CV((Comp.6 + Comp.7 + Comp.8) / C_tot), 
                            DOC_cv = CV(NPOC),
                            percent_HMWS_C_cv = CV(HMWS_C / CDOC),
                            percent_Humic_C_cv = CV(humic_like_substance_C / CDOC),
                            percent_LMWS_C_cv = CV(LMWS_C / CDOC)),
                         by = .(site, groups.x, Class, alteration)]

signif((aggregate(DOM_averages$percent_Humic_C, by = DOM_averages[,"groups.x"], FUN = mean)$x), 2)
signif((aggregate(DOM_averages$percent_Humic_C, by = DOM_averages[,"groups.x"], FUN = sd))$x , 1)

signif((aggregate(DOM_CV_averages$percent_Humic_C, by=DOM_CV_averages[,"groups.x"], FUN = mean))$x, 1)
signif(aggregate(DOM_CV_averages$percent_Humic_C, by=DOM_CV_averages[,"groups.x"], FUN = sd)$x , 1)


shapiro.test(log(DOM_CV_averages$HI2_cv))
hist(log(DOM_CV_averages$HI2_cv))
boxplot(log(DOM_CV_averages$HI2_cv))
qqnorm(log(DOM_CV_averages$HI2_cv),main = "Normal Q-Q Plot");qqline(log(DOM_CV_averages$HI2_cv)) 


bartlett.test(log(percent_Humic_C_cv)~groups.x, data = DOM_CV_averages)
oneway.test(log(percent_Humic_C_cv)~groups.x, var.equal = F, data = DOM_CV_averages)

var.test(log(percent_Humic_C_cv)~alteration, data = DOM_CV_averages[Class == "Mediterranean"])
t.test(log(percent_Humic_C_cv)~alteration, data = DOM_CV_averages[Class == "Mediterranean"], var.equal = T)

var.test(log(percent_Humic_C_cv)~alteration, data = DOM_CV_averages[Class == "Temperate"])
t.test(log(percent_Humic_C_cv)~alteration, data = DOM_CV_averages[Class == "Temperate"], var.equal = T)

var.test(log(percent_Humic_C_cv)~Class, data = DOM_CV_averages[alteration == "Natural"])
t.test(log(percent_Humic_C_cv)~Class, data = DOM_CV_averages[alteration == "Natural"], var.equal = T)

p.adjust(c(0.05581, 0.003154, 0.08455), method = "bonferroni")

m = pairwise.t.test(log(DOM_CV_averages$C_protein_cv), DOM_CV_averages$groups.x, p.adjust.method = "bonferroni", 
                    pool.sd = F)
multcompView::multcompLetters(fullPTable(m$p.value))

#group variation aspect (iii)
### End of statistical analysisi on individual  indices ###

#This is an anova per river group and parameter, 
#I calculated the percentage of error that was 
#coming from within groups (so for the season, in my case). 
#For that I just divided the sum sq coming from season by the total one (coming from site and season) 

##### PCA data prep ####
data_sum$campaign <- factor(x = data_sum$campaign, levels = c("oct", "dec", "feb", "apr", "may", "aug"), labels = c("Oct", "Dec", "Feb", "Apr", "May", "Aug"))

pca_data <- data_sum[, .(alteration, Class,groups.x,
             HIX2, FIX, beta.alpha, SR_Loiselle,
             E2.to.E3, SUVA254,
             C_humic = (Comp.1 + Comp.2 + Comp.3 + Comp.4 + Comp.5) / C_tot, 
             C_protein = (Comp.6 + Comp.7 + Comp.8) / C_tot,
             C_LMWS = (LMWS_C / CDOC),
             C_HMWS = (HMWS_C / CDOC),
             C_HS = (humic_like_substance_C / CDOC)),
         by = .(site, campaign)]

pca_data[is.na(pca_data)] <- 0

wine.pca <- prcomp(pca_data[, -(1:5)], scale. = TRUE) 
summary(wine.pca)
group_types <- factor(c("TempAlt", "TempNat", "MedAlt", "MedNat"))

group_list <- vector("list")
for(i in seq_along(group_types)){
  group_list[[i]] =  site_info[groups == group_types[[i]], c("site", "Class", "alteration", "groups", "Alteration_type", "Catchment")]
}

#### PCA with automated prcomp calcualted rotations ####
PCAloadings <- data.frame(Variables = rownames(wine.pca$rotation), wine.pca$rotation)
ggplot(data_sum, aes(x = wine.pca$x[,1], y = wine.pca$x[,2])) +
  geom_point(aes(color = data_sum$groups.y, fill = data_sum$groups.x, shape = data_sum$groups.x),  size = 4) +
  scale_shape_manual(values = c(15,16,17,18)) +
  scale_fill_manual(values = c("#942D0A", "#E65525", "#043005","#4F9608", "red", "blue")) +
  geom_segment(data = PCAloadings, aes(x = 0, y = 0, xend = (PC1*10),
                                       yend = (PC2*10)), arrow = arrow(length = unit(1/2, "picas")),color = "black") +
  annotate("text", x = (PCAloadings$PC1*10.4), y = (PCAloadings$PC2*10.4), label = PCAloadings$Variables) +
  theme_classic() +
  guides(fill = "legend") +
  #theme(legend.position = c(-1,0))+
  labs(color = "Sites", x = "PC 1 (32%)", y = "PC 2 (20%)", title = "PCA of PARAFAC Components and Other Optical Parameters")

####  polygon plots #### 
xlim <- c(-7.5, 7.5)
ylim <- c(-7.5, 7.5)
par(mfrow = c(2,2), mai = c(0.5,0.5,0.3,0.3))

plot_TempNat <- plot(wine.pca$x[, c(1:2)], type = "n", ylim = ylim, xlim = xlim, cex.main = 1.5, cex.axis = 1.25) #, main = "Natural Atlantic"
ordi_TempNat <- ordihull(ord = wine.pca$x[, c(1:2)], groups = data_sum$site ,display = "sites",draw = "polygon", label = F, show.groups = group_list[[2]]$site,
                         alpha = 0.7, col = c("#B4DCED"))
abline(v = 0, h = 0, lty = 2)

#text(x=-5.5, y=5.5,labels="B", cex=1.7)

plot_TempAlt <- plot(wine.pca$x[, c(1:2)], type="n",  ylim = ylim, xlim = xlim, cex.main = 1.5, cex.axis = 1.25) #main = "Altered Atlantic",
ordi_TempAlt <- ordihull(ord = wine.pca$x[,c(1:2)],groups = data_sum$site ,display = "sites",draw="polygon", label = F, show.groups = group_list[[1]]$site,
                         alpha = 0.7, col = c('#6996D1'))
abline(v = 0, h = 0, lty = 2)

#text(x=-5.5, y=5.5, labels="A", cex=1.7)

plot_MedNat <- plot(wine.pca$x[, c(1:2)], type="n", ylim = ylim, xlim = xlim, cex.main = 1.5, cex.axis = 1.25) # main = "Natural Mediterranean",
ordi_MedNat <- ordihull(ord = wine.pca$x[, c(1:2)],groups = data_sum$site ,display = "sites",draw = "polygon", label = F, show.groups = group_list[[4]]$site,
                        alpha = 0.7, col = c('#F5CB7D'))
abline(v = 0, h = 0, lty = 2)

#text(x=-5.4, y=5.5, labels="D", cex=1.7)

plot_MedAlt <- plot(wine.pca$x[, c(1:2)], type = "n", ylim = ylim, xlim = xlim, cex.main = 1.5, cex.axis =1.25) #, main = "Altered Mediterranean"
ordi_MedAlt <- ordihull(ord = wine.pca$x[,c(1:2)],groups = data_sum$site ,display = "sites",draw = "polygon", label = F, show.groups = group_list[[3]]$site,
                        alpha = 0.7, col = c("#F09E41"))
abline(v = 0, h = 0, lty = 2)

#text(x=-5.5, y=5.5, labels="C", cex=1.7)

par(mfrow = c(1,1))

grDevices::pdf('PCA_ploygons2.png', width = 19, height = 20)
#PCA_polygons
dev.off()

summary(ordi_MedNat)
summary(ordi_MedAlt)

#### Pipeline 2 for the statistical tests ####
## Pipeline 2 aspect (i) - Annual mean of centroid location

assign(paste0("centroid_",  group_list[[1]][1, groups]), 
       do.call(rbind.data.frame, lapply(ordi_TempAlt, centroid)
))

assign(paste0("centroid_",  group_list[[2]][1,groups]), 
       do.call(rbind.data.frame, lapply(ordi_TempNat, centroid)
       ))

assign(paste0("centroid_",  group_list[[3]][1,groups]), 
       do.call(rbind.data.frame, lapply(ordi_MedAlt, centroid)
       ))

assign(paste0("centroid_",  group_list[[4]][1,groups]), 
       do.call(rbind.data.frame, lapply(ordi_MedNat, centroid)
       ))

centroids = data.table(rbind(centroid_MedAlt, centroid_MedNat, centroid_TempAlt, centroid_TempNat), keep.rownames="site")
centroids = data.table(merge(centroids, site_info, by="site")) # Can also do this whole thing with summary(ordi_MedNat)

#Pipeline 2 aspect (ii) - temporal distance to centroid location
# Alternative A: using the polygon area as a temporal dispersion proxy. 
polygon_data = data.table(t(cbind(summary(ordi_MedNat), summary(ordi_MedAlt), summary(ordi_TempAlt), summary(ordi_TempNat))), keep.rownames = "site")
polygon_data = merge(polygon_data, site_info, by = "site")

PCA_results = cbind(wine.pca$x[, c(1:2)], pca_data)

# var of PC axes
PC_info = merge(PCA_results[, .(mean_PC1 = mean(PC1, na.rm = TRUE), mean_PC2 = mean(PC2, na.rm = TRUE), 
                                var_PC1 = var(PC1), var_PC2 = var(PC2)), by = .(site)], 
                site_info[, c("site", "Class", "alteration", "groups")])

#Statistical results of the PC1 aspect (i)
shapiro.test((PC_info$mean_PC1))
hist((PC_info$mean_PC1))
boxplot((PC_info$mean_PC1))
qqnorm((PC_info$mean_PC1),main = "Normal Q-Q Plot");qqline((PC_info$mean_PC1)) 

bartlett.test(mean_PC2~groups, data = PC_info)
oneway.test(mean_PC2~groups, var.equal = F, data = PC_info)

var.test(mean_PC2~alteration, data = PC_info[Class == "Mediterranean"])
t.test(mean_PC2~alteration, data = PC_info[Class == "Mediterranean"], var.equal = T)

var.test(mean_PC2~alteration, data = PC_info[Class == "Temperate"])
t.test(mean_PC2~alteration, data = PC_info[Class == "Temperate"], var.equal = T)

var.test(mean_PC2~Class, data = PC_info[alteration == "Natural"])
t.test(mean_PC2~Class, data = PC_info[alteration == "Natural"], var.equal = F)

p.adjust(c( 0.03377, 0.0004782,  0.8223), method = "bonferroni", n = 3)


#CV of PC1 aspect (ii)
shapiro.test(log(PC_info$var_PC1))
hist(log(PC_info$var_PC1))

bartlett.test(log(var_PC1)~groups, data = PC_info)
oneway.test(log(var_PC1)~groups, var.equal = F, data = PC_info)

var.test((var_PC2)~alteration, data = PC_info[Class == "Mediterranean"])
t.test((var_PC2)~alteration, data = PC_info[Class == "Mediterranean"], var.equal = F)

var.test((var_PC2)~alteration, data = PC_info[Class == "Temperate"])
t.test((var_PC2)~alteration, data = PC_info[Class == "Temperate"], var.equal = F)

var.test((var_PC2)~Class, data = PC_info[alteration == "Natural"])
t.test((var_PC2)~Class, data = PC_info[alteration == "Natural"], var.equal = F)

p.adjust(c(0.4319, 0.00383, 0.7015), method="bonferroni", n = 3)

m = pairwise.t.test(log(PC_info$var_PC1), PC_info$groups, p.adjust.method = "bonferroni", pool.sd = F, paired = F)
multcompView::multcompLetters(fullPTable(m$p.value), , threshold = 0.07)


#Statistical results of the PC2 aspect (i)
shapiro.test((PC_info$mean_PC2))
hist((PC_info$mean_PC2))

bartlett.test(mean_PC2~groups, data = PC_info)
oneway.test(mean_PC2~groups, var.equal = T, data = PC_info)

var.test(mean_PC2~alteration, data = PC_info[Class == "Mediterranean"])
t.test((mean_PC2)~groups, data = PC_info[Class == "Mediterranean"], var.equal = F)

var.test(mean_PC2~alteration, data = PC_info[Class == "Temperate"])
t.test(mean_PC2~alteration, data = PC_info[Class == "Temperate"], var.equal = F)

var.test(mean_PC2~Class, data = PC_info[alteration == "Natural"])
t.test(mean_PC2~Class, data = PC_info[alteration == "Natural"], var.equal = T)

p.adjust(c(0.07537, 0.0008622, 0.8163), method = "bonferroni", n = 3)

m = pairwise.t.test(PC_info$mean_PC2, PC_info$groups, p.adjust.method = "bonferroni", pool.sd = F, paired = F)
multcompView::multcompLetters(fullPTable(m$p.value))


#CV of PC2 aspect (ii)
shapiro.test(log(PC_info$var_PC2))
hist(log(PC_info$var_PC2))

bartlett.test(log(var_PC2)~groups, data = PC_info)
oneway.test(log(var_PC2)~groups, var.equal = T, data=PC_info)

var.test( log(var_PC2)~alteration, data = PC_info[Class == "Mediterranean"])
t.test( log(var_PC2)~alteration, data = PC_info[Class == "Mediterranean"], var.equal = F)

var.test( log(var_PC2)~alteration, data = PC_info[Class == "Temperate"])
t.test(log(var_PC2)~alteration, data = PC_info[Class == "Temperate"], var.equal = T)

var.test( log(var_PC2)~Class, data = PC_info[alteration == "Natural"])
t.test( log(var_PC2)~Class, data = PC_info[alteration == "Natural"], var.equal = T)

p.adjust(c(0.7324, 0.00036,   0.1091), method = "bonferroni", n = 3)

m = pairwise.t.test(log(PC_info$var_PC2), PC_info$groups, p.adjust.method = "bonferroni", pool.sd = F, paired = F)
multcompView::multcompLetters(fullPTable(m$p.value))


# Alternative B: using the  average of individual sampling distances to the centroid of each polygon as a proxy for temporal dispersion
#ordi_all=c(ordi_MedAlt, ordi_MedNat, ordi_TempAlt, ordi_TempNat)

site_info = site_info[site_info$site != "Carrion",]

# Calculating the centroid for all the PC axis

betas_PC_all = betadisper(vegdist((wine.pca$x[,1:11]), method = "euclidean"),group = pca_data$site, type="centroid")
disp_PC_all = tapply(betas_PC_all[["distances"]], betas_PC_all[["group"]], mean)

#### Centroid with all PC axes ####
Multi_centroid = data.table(cbind(disp_PC_all, betas_PC_all[["centroids"]], match(rownames(betas_PC_all[["centroids"]]), rownames(disp_PC_all))), keep.rownames = "site")
Multi_centroid = Multi_centroid[site_info, on = .(site)]

adonis2(vegdist(Multi_centroid[,3:12], method = "euclidian")~Multi_centroid$groups, permutations = 100000)

adonis2(vegdist(Multi_centroid[Class == "Mediterranean", c(3:12)], method = "euclidian")~alteration, data = Multi_centroid[Class == "Mediterranean"], permutations = 100000)

adonis2(vegdist((Multi_centroid[Class == "Temperate", c(3:12)]), method = "euclidian")~alteration, data = Multi_centroid[Class == "Temperate"], permutations = 100000)

adonis2(vegdist((Multi_centroid[alteration == "Natural",c(3:12)]), method = "euclidian")~Class, data = Multi_centroid[alteration == "Natural"], permutations = 10000)

p.adjust(c(0.07102 ,   0.2437, 0.1165), method = "bonferroni", n = 3)


# Temporal dispersion with all PC axes 
shapiro.test(log(Multi_centroid$disp_PC_all))
hist(log(Multi_centroid[Class == "Temperate"]$disp_PC_all))
qqnorm(log(Multi_centroid$disp_PC_all), main="Normal Q-Q Plot of male");qqline(log(Multi_centroid$disp_PC_all))

bartlett.test(log(disp_PC_all)~groups, data = Multi_centroid)
oneway.test(log(disp_PC_all)~groups, var.equal = F, data = Multi_centroid)

var.test(log(disp_PC_all)~alteration, data = Multi_centroid[Class == "Mediterranean"])
t.test(log(disp_PC_all)~alteration, data = Multi_centroid[Class == "Mediterranean"], var.equal = F)

var.test(log(disp_PC_all)~alteration, data = Multi_centroid[Class == "Temperate"])
t.test(log(disp_PC_all)~alteration, data = Multi_centroid[Class == "Temperate"], var.equal = F)

var.test(log(disp_PC_all)~Class, data = Multi_centroid[alteration == "Natural"])
t.test(log(disp_PC_all)~Class, data = Multi_centroid[alteration == "Natural"], var.equal = F)

p.adjust(c(0.706, 0.005629, 0.1206), method = "bonferroni", n = 3)
m = pairwise.t.test(log(Multi_centroid$disp_PC_all), Multi_centroid$groups, p.adjust.method = "bonferroni", pool.sd = F)
multcompView::multcompLetters(fullPTable(m$p.value))

# Aspect (iii): mvd
summary(aov(betadisper(vegdist((Multi_centroid[, c(3:12)]), method = "euclidian"), Multi_centroid$groups)$distances~Multi_centroid$groups))

# Manuel calculations of multidimensional PCA centroid
pca_results = cbind(pca_data, wine.pca$x)
manuel_centroid = aggregate(x = pca_results[, c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "PC11")], by = pca_results[,c("site")], FUN = mean)
manuel_centroid = data.table(cbind(manuel_centroid, site_info, by = "site"))

#### Centroid with all PC axes end ####

disp_all = as.data.table(data.frame(disp_PC_all), keep.rownames = "site")
disp_all = merge(disp_all, unique(pca_data[, c("site", "groups.x", "Class", "alteration")]), by = "site")

disp_melted = melt(disp_all, id.vars = c("site", "groups.x", "Class", "alteration"), measure.vars = c("disp_PC_all"))

formula <- lm(formula = disp_PC2~groups.x, data = disp_all[disp_all$Class == "Temperate", ])

pl = ggplot(data = disp_melted[variable == "disp_PC_all"], aes(x = groups.x, y = value)) + 
  theme(axis.line.x =  element_line(color="black"),axis.line.y =  element_line(), panel.grid = element_blank(), panel.background = element_blank())+
  geom_boxplot(aes(x = groups.x, y = value, fill = groups.x), outlier.size = 0.5, lwd = 0.2)+
 # geom_point(aes(fill=groups.x, shape=groups.x),fill="black", size=3)+
  #facet_wrap(~variable,  labeller=as_labeller(c(disp_PC1="PC1", disp_PC2="PC2", disp_PC_all="Dispersion on the 2 PC space")))+
  labs(y = "Dispersion", x = NULL)+
  scale_fill_manual(values = c("#B4DCED", '#6996D1','#F5CB7D','#F09E41'), labels = c("aM", "nM", "aA", "nA")) +
  theme_pca() + theme(legend.title = element_blank(), legend.position = "bottom") + scale_x_discrete(labels = c("MedAlt" = "aM","MedNat" = "nM",
                                                                                                        "TempAlt" = "aA", "TempNat" = "nA"))+
  scale_shape_manual(values = c(23,22,25,24))+
  theme(panel.background = element_rect(fill = NA), strip.background.x = element_rect(fill="white"), strip.text = element_text(size = 11), axis.title = element_text(size = 11))+
  theme(text = element_text(size = 11),axis.text = element_text(size = 11, color = "black"), axis.text.x = element_text(hjust = 1),
        axis.title=element_blank(), legend.position = "none", strip.placement ="outside", strip.background = element_blank(), strip.text = element_text(size = 11, color = "black"))


pdf('plots/PCAll.pdf', width = 2.4, height = 2.5)
plot(pl)
dev.off()


pca_results$groups.x = factor(pca_results$groups.x, levels = c("TempNat", "TempAlt", "MedNat", "MedAlt"))

individual_pc1 = ggplot(data = pca_results, aes(x=groups.x, y=PC1)) + 
  theme(axis.line.x =  element_line(color="black"),axis.line.y =  element_line(), panel.grid = element_blank(), panel.background = element_blank())+
  geom_boxplot(aes(x=groups.x, y=PC1,  fill=groups.x, group=reorder(site, PC1, median)),  lwd = 0.2, outlier.size = 0.5)+
  labs(y="Axis locations", x=NULL)+scale_y_continuous(limits=c(-4,7.5))+
  scale_fill_manual(values =c("#B4DCED", '#6996D1','#F5CB7D','#F09E41'))+
  theme_pca()+theme(legend.title=element_blank(), legend.position = "bottom")+scale_x_discrete(labels=c("MedAlt" = "aM","MedNat" = "nM", "TempAlt" = "aA", "TempNat" = "nA"))+
  theme(panel.background = element_rect(fill = NA), axis.title = element_text(size=11))+
  theme(text = element_text(size=11),axis.text=element_text(size=11, color="black"),
        axis.title=element_text(size=11), legend.position = "none")+ylab("PC1")

individual_pc2 = ggplot(data = pca_results, aes(x=groups.x, y=PC2)) + 
  theme(axis.line.x =  element_line(color="black"),axis.line.y =  element_line(), panel.grid = element_blank(), panel.background = element_blank())+
  geom_boxplot(aes(x=groups.x, y=PC2, fill=groups.x, group=reorder(site, PC2, median)),lwd=0.2, outlier.size=0.5)+
  labs(y="Axis locations", x=NULL)+
  scale_fill_manual(values =c("#B4DCED", '#6996D1','#F5CB7D','#F09E41'))+
  theme_pca()+theme(legend.title=element_blank())+scale_x_discrete(labels=c("MedAlt" = "aM","MedNat" = "nM", "TempAlt" = "aA", "TempNat" = "nA"))+
  scale_shape_manual(values=c(23,22,25,24))+scale_y_continuous(limits=c(-5,7.5))+
  theme(panel.background = element_rect(fill = NA), axis.title = element_text(size=11))+
  theme(text = element_text(size=11),axis.text=element_text(size=11, color="black"),
        axis.title=element_text(size=11), legend.position = "none")+ylab("PC2")

pdf('plots/PC1.pdf', width = 4, height = 4)
plot(individual_pc1)
dev.off()
pdf('plots/PC2.pdf', width = 4, height = 4)
plot(individual_pc2)
dev.off()

#### data reflection onto the PCA space ####
PCA_scores <- data.table(pca_data, wine.pca$x)
PCA_rot <- data.table(t(cor(PCA_scores[,c("PC1", "PC2")], pca_data[,-(1:5)], method = "pearson")), keep.rownames = "Variables")

PCAloadings <- data.table(Variables = rownames(wine.pca$rotation), wine.pca$rotation)

pca.scores<-wine.pca$x
eigenvec12<-cbind(wine.pca$rotation[,1],wine.pca$rotation[,2])
PCAloadings<-data.frame(cor(pca_data[,-(1:5)],pca.scores))
PCAloadings$Variables=rownames(PCAloadings)

### PCA plots #### 
MedAlt <- PCA_scores[PCA_scores$groups == "MedAlt", ][chull(PCA_scores[PCA_scores$groups == "MedAlt", c("PC1", "PC2")]), ]  # hull values for grp A
MedNat <- PCA_scores[PCA_scores$groups == "MedNat", ][chull(PCA_scores[PCA_scores$groups == "MedNat", c("PC1", "PC2")]), ]  # hull values for grp A
TempAlt <- PCA_scores[PCA_scores$groups == "TempAlt", ][chull(PCA_scores[PCA_scores$groups == "TempAlt", c("PC1", "PC2")]), ]  # hull values for grp A
TempNat <- PCA_scores[PCA_scores$groups == "TempNat", ][chull(PCA_scores[PCA_scores$groups == "TempNat", c("PC1", "PC2")]), ]  # hull values for grp A

hull.data <- rbind(MedAlt, MedNat,TempAlt,TempNat)  #combine grp.a and grp.b
hull.data
# title="PCA of PARAFAC Components and Other Optical Parameters"
plot_all <- ggplot(data_sum, aes(x=wine.pca$x[,1], y=wine.pca$x[,2]))+
  geom_point(aes(color=data_sum$groups.x, fill=data_sum$groups.x, shape=data_sum$groups.x), size=2.5, colour="black")+
  scale_shape_manual(values=c(23,22,25,24))+
  scale_fill_manual(values=c("#942D0A","#E65525", "#043005","#4F9608"))+
  geom_segment(data = PCA_rot, aes(x = 0, y = 0, xend = (PC1*8), yend = (PC2*8)), arrow = arrow(length = unit(1/2, "picas")),color = "black") +
  annotate("text", x = (PCA_rot$PC1*8), y = (PCA_rot$PC2*8),label = PCA_rot$Variables, size=4, color="black")+
  labs(color="Sites", x="PC1 (44.2.8%)", y="PC2 (15.6%)", tag = "a. Loadings and sites")+
  theme_pca()+ 
  theme(plot.tag.position=c(0.24,0.97))+
  scale_x_continuous(limits=c(-8,8), n.breaks=10)+
  scale_y_continuous(limits=c(-8,8), n.breaks=10)+
  #guides(fill="legend")+
  #theme(legend.position = c(-1,0))+
  geom_vline(xintercept = 0, lty=2) + geom_hline(yintercept = 0, lty=2)

hull.data$groups=factor(hull.data$groups, levels = c("TempNat", "TempAlt", "MedNat", "MedAlt"))

plot_points <- ggplot(pca_results, aes(x = wine.pca$x[, 1], y = wine.pca$x[, 2]))+
  geom_polygon(data = hull.data, mapping=aes(x = PC1, y = PC2,color = groups, group = groups), fill = NA,lwd = 1, alpha = 0.7)+
  geom_point(aes(fill = pca_results$groups.x, shape=pca_results$groups.x), size=1, alpha=0.7,  color="black", stroke=0.25)+
  scale_shape_manual(values=c(23,22,25,24), labels=c("MedAlt" = "aM","MedNat" = "nM", "TempAlt" = "aA", "TempNat" = "nA"))+
  scale_fill_manual(values=c( "#B4DCED",'#6996D1','#F5CB7D','#F09E41'), labels=c("MedAlt" = "aM","MedNat" = "nM", "TempAlt" = "aA", "TempNat" = "nA"))+
  scale_color_manual(values =c( "#B4DCED",'#6996D1','#F5CB7D','#F09E41'), labels=c("MedAlt" = "aM","MedNat" = "nM", "TempAlt" = "aA", "TempNat" = "nA"))+
  labs(color="Sites", x="PC1 (44%)", y="PC2 (16%)")+ #, tag = "A"
  theme_pca()+
  theme(plot.tag.position=c(0.15,0.96))+
  scale_x_continuous(limits=c(-7.5,7.5), n.breaks=10)+
  scale_y_continuous(limits=c(-7.5,7.5), n.breaks=10)+
  geom_vline(xintercept = 0, lty=2) + geom_hline(yintercept = 0, lty=2)+
  theme(legend.title = element_blank(),legend.position = "none", text=element_text(size=11))+
  theme(axis.text =  element_text(size=11, color="black"),axis.title=element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size=5, alpha=1)))

pdf('plots/PCA_points.pdf', width = 2.7, height = 2.5)
plot(plot_points)
dev.off()

plot_optical <- ggplot(data_sum, aes(x=wine.pca$x[,1], y=wine.pca$x[,2]))+
  geom_segment(data = PCA_rot, aes(x = 0, y = 0, xend = (PC1), yend = (PC2)), arrow = arrow(length = unit(1/2, "picas")),color = "black") +
  annotate("text", x = (PCA_rot$PC1), y = (PCA_rot$PC2),label = PCA_rot$Variables, size=3.5, color="black")+
  labs(color="Sites", x="PC1 (47.1%)", y="PC2 (13.6%)")+
  theme_pca()+ 
  theme(plot.tag.position=c(0.20,0.95))+
  scale_x_continuous(limits=c(-1,1), n.breaks=10)+
  scale_y_continuous(limits=c(-1,1), n.breaks=10)+
  guides(fill="legend")+
  theme(legend.position = c(-1,0))+
  geom_vline(xintercept = 0, lty=2) + geom_hline(yintercept = 0, lty=2)+
  theme(axis.text=element_text(size=25),axis.title=element_text(size=25))+
  theme(text = element_text(size=11))

 PCA_optical = ggplot(data_sum, aes(x = wine.pca$x[,1], y = wine.pca$x[,2]))+
  geom_segment(data = PCA_rot, aes(x = 0, y = 0, xend = (PC1), yend = (PC2)), arrow = arrow(length = unit(1/2, "picas")),color = "black", lwd=0.7) +
 #geom_segment(data=LC_OCD_rot, aes(x = 0, y = 0, xend = (PC1), yend = (PC2)), arrow = arrow(length = unit(1/2, "picas")),color = "blue", lwd=0.7) + 
  #annotate("text", x = as.vector(labels_locator2[1:14,1]), y = as.vector(labels_locator[1:14,2]),label =var_names[1:14], size=5, color="black")+
  # annotate("text", as.vector(labels_locator2[15:17,1]), y = as.vector(labels_locator[15:17,2]), label =var_names[15:17], size=5, color="blue")+
  labs(color="Sites", x="PC1 (44%)", y="PC2 (16%)")+
  theme_pca()+ 
  theme(plot.tag.position=c(0.19,0.95))+
  scale_x_continuous(limits=c(-1,1), n.breaks=8)+
  scale_y_continuous(limits=c(-1,1), n.breaks=8)+
  guides(fill="legend")+
  theme(legend.position = c(-1,0) )+
  geom_vline(xintercept = 0, lty=2) + geom_hline(yintercept = 0, lty=2)+
  theme(text = element_text(size=20),axis.text=element_text(size=20, color = "black"))
 
 pdf('plots/PCA_optical.pdf', width = 6.7, height = 6)
 plot(PCA_optical)
 dev.off()
#### end ####

#### Hydrological Indices Analysis ####

source("indice_analysis.R")


#### Loading plots ####
parafac_em <- read.csv("Em.csv")
parafac_ex<- read.csv("Ex.csv")
ggplot()+
  geom_line(aes(x=parafac_em$Em,y=(1/3)*(parafac_em$Comp.8)))+
  geom_line(aes(x=parafac_ex$Ex,y=3*(parafac_ex$Comp.8)), linetype = "dashed")+
  theme_classic()+
  labs(title="Comp8")

# Hydrological plots
source("hydrological_analysis.R")
