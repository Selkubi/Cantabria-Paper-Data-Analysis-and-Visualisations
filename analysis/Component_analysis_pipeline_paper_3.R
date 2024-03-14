# library(vegan)
 library(geosphere)
# library(gridExtra)
# library(data.table)
# library(ggplot2)
# library(lubridate)

#### The functions to be used are here, but first you have to load the project package
roxygen2::roxygenize()

#### 1. Data Cleaning ####
source("analysis/Data_Cleaning.R")

##### 2. DOC Mean and CV plots #####
source("analysis/DOC_DOM_plots.R")

mean_NPOC 

if(!dir.exists("output/plots")) {dir.create("output/plots")}
pdf('output/plots/mean_NPOC.pdf', width = 4, height = 4)
plot(mean_NPOC)
dev.off()

mediterranean_DOC 

pdf('output/plots/mediterranean_DOC.pdf', width = 2.8, height = 2.8)
plot(mediterranean_DOC)
dev.off()

temperate_DOC

pdf('output/plots/temperate_DOC.pdf', width = 2.8, height = 2.8)
plot(temperate_DOC)
dev.off()

##### 3. Statistical tests for DOC concentration mean values #####
shapiro.test(log(DOC_sum$mean_NPOC))
hist(log(DOC_sum$mean_NPOC))
boxplot(log((DOC_sum$mean_NPOC)))

pairwise_results_mean <- lapply(DOC_sum[, c(3,5)], pairwise.t.test, DOC_sum$groups, p.adjust.methods = "bonferroni", pool.sd = T)

oneway_test_results(DOC_sum, "mean_NPOC", "groups", log_normalise = T)
statistics_pipeline_wrapper(DOC_sum, "mean_NPOC", "groups", log_normalise = T)
multcompView::multcompLetters(rcompanion::fullPTable(pairwise_results_mean[["mean_NPOC"]]$p.value))$Letters

####  4. Statistical tests for DOC concentration CV values ####
shapiro.test(log(DOC_sum$var_NPOC))
hist(log(DOC_sum$var_NPOC))
boxplot(log((DOC_sum$var_NPOC)))

oneway_test_results(DOC_sum, "var_NPOC", "groups", log_normalise = T)
statistics_pipeline_wrapper(DOC_sum, "var_NPOC", "groups", log_normalise = T)
multcompView::multcompLetters(rcompanion::fullPTable(pairwise_results_mean[["var_NPOC"]]$p.value))$Letters

#### End of DOC mean and CV ####

#### DOM composition ####
# data_summary = data_sum[, .(FI2 = (FIX), HI2= (HIX2), 
#                           beta.alpha2 = (beta.alpha ), SUVA254_2 = (SUVA254),
#                           SR = (SR_Loiselle ), E2toE3 = (E2.to.E3 ),
#                           C_HMWS = (HMWS_C/CDOC), C_LMWS = (LMWS_C/CDOC),
#                           C_HS = (humic_like_substance_C/CDOC ),
#                           C_terrestrial = ((Comp.1 + Comp.2 + Comp.3 + Comp.4 + Comp.5) / C_tot),
#                           C_protein = ((Comp.6 + Comp.7 + Comp.8) / C_tot)),
#                           by = .(site, alteration, Class,  groups.x, campaign)]

#data_summary[data_summary[, C_HMWS > 0.3]]$C_HMWS = NA
#data_summary[data_summary[, C_HS > 1.5]]$C_HS = NA

data_sum$C_tot=data_sum[,Comp.1+Comp.2+Comp.3+Comp.4+Comp.5+Comp.6+Comp.7+Comp.8]

##### 1. DOM composition Mean and sd Calculation ####
DOM_averages <- data_sum[,.(FI2 = mean(FIX, na.rm = T), HI2 = mean(HIX2, na.rm = T), 
                         beta.alpha2= mean(beta.alpha, na.rm = T), SUVA254_2 = mean(SUVA254, na.rm = T),
                         SR = mean(SR_Loiselle, na.rm = T), E2toE3 = mean(E2.to.E3, na.rm = T), 
                         C_humic = mean((Comp.1 + Comp.2 + Comp.3 + Comp.4 + Comp.5) / C_tot, na.rm = T),
                         C_protein = mean((Comp.6 + Comp.7 + Comp.8) / C_tot, na.rm = T), 
                         DOC = mean(NPOC, na.rm = T), 
                         percent_HMWS_C = mean(HMWS_C / CDOC, na.rm = T),
                         percent_Humic_C = mean(humic_like_substance_C / CDOC, na.rm = T),
                         percent_LMWS_C = mean(LMWS_C / CDOC, na.rm = T, trim = 5)),
                      by = .(site, groups.x, Class, alteration)]

mean_DOM_values <- data.table(get_mean_and_sd(DOM_averages, "percent_Humic_C", "groups.x")[,1])
for(i in names(DOM_averages[,-c(1:4)])) {
  mean_DOM_values <- merge(mean_DOM_values, get_mean_and_sd(DOM_averages,  names(DOM_averages[, ..i]), "groups.x"), by = "name")
}
# all the mean+-sd values for DOM composition variables
t(mean_DOM_values)

##### 2. Statistical tests for DOM composition mean values ####
# check the shapiro test results of all the values
shapiro_test_results <- sapply(DOM_averages[, -c(1:4)], function(x) shapiro.test(x)$p.value)

# plot the histograms to visually check the normality
par(mfrow = c(3, 4))
for (col in names(DOM_averages[, -c(1:4)])) {
  hist(DOM_averages[[col]], main = col, xlab = "Value")
}
par(mfrow = c(1, 1))

# Create a log_normalise vector of logical values looking at the above histograms and shapiro results
log_normalise <- rep(F, 12)
log_normalise[c(9, 10, 11)] <- T

# calculate the oneway tset results of all the values (ir the flow regime tests have p < 0.05)
pairwise_results_mean <- lapply(DOM_averages[,-c(1:4)], pairwise.t.test, DOM_averages$groups.x, p.adjust.methods = "bonferroni", pool.sd = T)

mean_results <- vector("list", length = length(names(DOM_averages[, -c(1:4)])))
names(mean_results) <-  names(DOM_averages[, -c(1:4)])
for(i in seq_along(DOM_averages[,-c(1:4)])) {
  names(mean_results[i]) <- names(DOM_averages[, -c(1:4)])[i]
  mean_results[[i]]<- statistics_pipeline_wrapper(DOM_averages,  names(DOM_averages[,-c(1:4)])[i], "groups.x", log_normalise = log_normalise[i])
  mean_results[[i]][["t_test_Letters"]] <- multcompView::multcompLetters(rcompanion::fullPTable(pairwise_results_mean[[i]]$p.value))$Letters
} 


# all the statistical test results for the mean values 
mean_results

##### 3. Printing F test and t-test results for DOM composition mean values #####
F_test_list <- list()
t_test_list <- list()
t_test_letters <- list()
for (variable_name in names(mean_results)) {
  F_test_list[[variable_name]] <- mean_results[[variable_name]]$F_test
  t_test_list[[variable_name]] <- mean_results[[variable_name]]$t_test
  t_test_letters[[variable_name]] <- mean_results[[variable_name]]$t_test_Letters
}
F_test_list
t_test_list
t_test_letters

#### End of mean +-sd results and the related statistical tests ####

##### 4. DOM composition CV and sd calculations #####
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

CV_DOM_values <- data.table(get_mean_and_sd(DOM_CV_averages, "FI2_cv", "groups.x")[,1])
for(i in names(DOM_CV_averages[,-c(1:4)])) {
  CV_DOM_values <- merge(CV_DOM_values, get_mean_and_sd(DOM_CV_averages,  names(DOM_CV_averages[, ..i]), "groups.x"), by = "name")
}
# all the mean+-sd values for DOM composition variables
t(CV_DOM_values)

##### 5. Statistical tests for DOM composition CV values ####
# check the shapiro test results of all the values
shapiro_test_results <- sapply(DOM_CV_averages[, -c(1:4)], function(x) shapiro.test(x)$p.value)

# plot the histograms to visually check the normality
par(mfrow = c(3, 4))
for (col in names(DOM_CV_averages[, -c(1:4)])) {
  hist(DOM_CV_averages[[col]], main = col, xlab = "Value")
}
par(mfrow = c(1, 1))

#  Create a log_normalise vector of logical values looking at the above histograms and shapiro results
log_normalise <- rep(F, 12)
log_normalise[c(6, 7, 8, 11)] <- T

# calculate the oneway test results of all the values (ir the flow regime tests have p < 0.05)
pairwise_results_CV <- lapply(DOM_CV_averages[,-c(1:4)], pairwise.t.test, DOM_CV_averages$groups.x, p.adjust.methods = "bonferroni", pool.sd = T)

CV_results <- vector("list", length = length(names(DOM_CV_averages[, -c(1:4)])))
names(CV_results) <-  names(DOM_CV_averages[, -c(1:4)])
for(i in seq_along(DOM_CV_averages[,-c(1:4)])) {
  names(CV_results[i]) <- names(DOM_CV_averages[, -c(1:4)])[i]
  CV_results[[i]]<- statistics_pipeline_wrapper(DOM_CV_averages,  names(DOM_CV_averages[,-c(1:4)])[i], "groups.x", log_normalise = log_normalise[i])
  CV_results[[i]][["t_test_Letters"]] <- multcompView::multcompLetters(rcompanion::fullPTable(pairwise_results_CV[[i]]$p.value))$Letters
} 

CV_results

##### 6. Printing F test and t-test results for DOM composition CV values #####
F_test_list <- list()
t_test_list <- list()
t_test_letters <- list()
for (variable_name in names(CV_results)) {
    F_test_list[[variable_name]] <- CV_results[[variable_name]]$F_test
    t_test_list[[variable_name]] <- CV_results[[variable_name]]$t_test
    t_test_letters[[variable_name]] <- CV_results[[variable_name]]$t_test_Letters
  }
F_test_list
t_test_list
t_test_letters
#### End of mean and CV calculations ####

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
source("analysis/multivariate_pca_plots.R")

# Centroid data (longitude and latitude) is extracted from here and named with centroid_<name>, and stored as a data.table within centroid
for(i in seq_along(group_list)){
assign(paste0("centroid_",  group_list[[i]][1, groups]), 
       do.call(rbind.data.frame, lapply(get(paste0("ordi_", (group_list[[i]][1, groups]))), centroid)
))
}

centroids = data.table(rbind(centroid_MedAlt, centroid_MedNat, centroid_TempAlt, centroid_TempNat), keep.rownames="site")
centroids = data.table(merge(centroids, site_info, by="site")) # Can also do this whole thing with summary(ordi_MedNat)

# The plots also produce ordihull objects that has the area which can be accessed with eg: summary(ordi_MedNat) 
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

bartlett.test(mean_PC1~groups, data = PC_info)
oneway.test(mean_PC1~groups, var.equal = T, data = PC_info)

var.test(mean_PC1~alteration, data = PC_info[Class == "Mediterranean"])
t.test(mean_PC1~alteration, data = PC_info[Class == "Mediterranean"], var.equal = T)

var.test(mean_PC1~alteration, data = PC_info[Class == "Temperate"])
t.test(mean_PC1~alteration, data = PC_info[Class == "Temperate"], var.equal = T)

var.test(mean_PC1~Class, data = PC_info[alteration == "Natural"])
t.test(mean_PC1~Class, data = PC_info[alteration == "Natural"], var.equal = T)

p.adjust(c( 0.03377, 0.0004782,  0.8223), method = "bonferroni", n = 3)


#CV of PC1 aspect (ii)
shapiro.test((PC_info$var_PC1))
hist((PC_info$var_PC1))

bartlett.test((var_PC1)~groups, data = PC_info)
oneway.test((var_PC1)~groups, var.equal = T, data = PC_info)

var.test((var_PC1)~alteration, data = PC_info[Class == "Mediterranean"])
t.test((var_PC1)~alteration, data = PC_info[Class == "Mediterranean"], var.equal = T)

var.test((var_PC1)~alteration, data = PC_info[Class == "Temperate"])
t.test((var_PC1)~alteration, data = PC_info[Class == "Temperate"], var.equal = T)

var.test((var_PC1)~Class, data = PC_info[alteration == "Natural"])
t.test((var_PC1)~Class, data = PC_info[alteration == "Natural"], var.equal = T)

p.adjust(c(0.1603, 0.01589,  0.0496), method="bonferroni", n = 3)

m = pairwise.t.test(log(PC_info$var_PC1), PC_info$groups, p.adjust.method = "bonferroni", pool.sd = F, paired = F)

multcompView::multcompLetters(fullPTable(m$p.value), , threshold = 0.07)

#Statistical results of the PC2 aspect (i)
shapiro.test((PC_info$mean_PC2))
hist((PC_info$mean_PC2))

bartlett.test(mean_PC2~groups, data = PC_info)
oneway.test(mean_PC2~groups, var.equal = F, data = PC_info)

var.test(mean_PC2~alteration, data = PC_info[Class == "Mediterranean"])
t.test((mean_PC2)~groups, data = PC_info[Class == "Mediterranean"], var.equal = T)

var.test(mean_PC2~alteration, data = PC_info[Class == "Temperate"])
t.test(mean_PC2~alteration, data = PC_info[Class == "Temperate"], var.equal = T)

var.test(mean_PC2~Class, data = PC_info[alteration == "Natural"])
t.test(mean_PC2~Class, data = PC_info[alteration == "Natural"], var.equal = T)

p.adjust(c(0.8163,  0.0004782, 0.03377), method = "bonferroni", n = 3)

m = pairwise.t.test(PC_info$mean_PC2, PC_info$groups, p.adjust.method = "bonferroni", pool.sd = F, paired = F)
multcompView::multcompLetters(fullPTable(m$p.value))


#CV of PC2 aspect (ii)
shapiro.test(log(PC_info$var_PC2))
hist(log(PC_info$var_PC2))

bartlett.test(log(var_PC2)~groups, data = PC_info)
oneway.test(log(var_PC2)~groups, var.equal = T, data=PC_info)

var.test(log(var_PC2)~alteration, data = PC_info[Class == "Mediterranean"])
t.test( log(var_PC2)~alteration, data = PC_info[Class == "Mediterranean"], var.equal = T)

var.test( log(var_PC2)~alteration, data = PC_info[Class == "Temperate"])
t.test(log(var_PC2)~alteration, data = PC_info[Class == "Temperate"], var.equal = F)

var.test( log(var_PC2)~Class, data = PC_info[alteration == "Natural"])
t.test( log(var_PC2)~Class, data = PC_info[alteration == "Natural"], var.equal = F)

p.adjust(c(0.6973, 0.0003826,    0.1542), method = "bonferroni", n = 3)

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

# These correspond to the "Mean DOM composotion" for the Global Tests on Dispersion table. 
gb_fr = adonis2(vegdist(Multi_centroid[,3:12], method = "euclidian")~Multi_centroid$groups, permutations = 100000)

gb_me = adonis2(vegdist(Multi_centroid[Class == "Mediterranean", c(3:12)], method = "euclidian")~alteration, data = Multi_centroid[Class == "Mediterranean"], permutations = 100000)

gb_te = adonis2(vegdist((Multi_centroid[Class == "Temperate", c(3:12)]), method = "euclidian")~alteration, data = Multi_centroid[Class == "Temperate"], permutations = 100000)

gb_nat = adonis2(vegdist((Multi_centroid[alteration == "Natural",c(3:12)]), method = "euclidian")~Class, data = Multi_centroid[alteration == "Natural"], permutations = 10000)

mean_DOM_composition_results = rbind(get_adonis_results(gb_fr), do_bonferroni_to_adonis(gb_nat, gb_te, gb_me))

colnames = c("Test", "flow_regime", "nA-nM", "nA-aA", "aM-aM")
mean_DOM = c("F" = gb_fr[["F"]][1])

#not putting the p adjustd values becaues these are F test-like p values
p.adjust(c(0.07102,  0.2437, 0.1165), method = "bonferroni", n = 3)

# Temporal dispersion with all PC axes 
shapiro.test(log(Multi_centroid$disp_PC_all))
hist(log(Multi_centroid[Class == "Temperate"]$disp_PC_all))
qqnorm(log(Multi_centroid$disp_PC_all), main="Normal Q-Q Plot of male");qqline(log(Multi_centroid$disp_PC_all))

bartlett.test(log(disp_PC_all)~groups, data = Multi_centroid)
oneway.test(log(disp_PC_all)~groups, var.equal = T, data = Multi_centroid)

var.test(log(disp_PC_all)~alteration, data = Multi_centroid[Class == "Mediterranean"])
t.test(log(disp_PC_all)~alteration, data = Multi_centroid[Class == "Mediterranean"], var.equal = T)

var.test(log(disp_PC_all)~alteration, data = Multi_centroid[Class == "Temperate"])
t.test(log(disp_PC_all)~alteration, data = Multi_centroid[Class == "Temperate"], var.equal = T)

var.test(log(disp_PC_all)~Class, data = Multi_centroid[alteration == "Natural"])
t.test(log(disp_PC_all)~Class, data = Multi_centroid[alteration == "Natural"], var.equal = T)

p.adjust(c(0.1079,  0.00273,  0.6857), method = "bonferroni", n = 3)
m = pairwise.t.test(log(Multi_centroid$disp_PC_all), Multi_centroid$groups, p.adjust.method = "bonferroni", pool.sd = F)
multcompView::multcompLetters(fullPTable(m$p.value))

# Aspect (iii): mvd
summary(aov(betadisper(vegdist((Multi_centroid[, c(3:12)]), method = "euclidian"), Multi_centroid$groups)$distances~Multi_centroid$groups))

# Aspect extra
river_centroid_dispersions = betadisper(vegdist(Multi_centroid[,3:12], method = "euclidian"), group = Multi_centroid$groups, type="centroid")
river_centroid_dispersions_sites = data.table(cbind(river_centroid_dispersions=river_centroid_dispersions[["distances"]], 
                                                                             site = Multi_centroid[["site"]],
                                                                             groups = Multi_centroid[["groups"]],
                                                                             alteration = Multi_centroid[["alteration"]],
                                                                             Class = Multi_centroid[["Class"]]))

summary(aov(river_centroid_dispersions_sites[Class == "Temperate"]$river_centroid_dispersions~river_centroid_dispersions_sites[Class == "Temperate"]$groups))
summary(aov(river_centroid_dispersions_sites[Class == "Mediterranean"]$river_centroid_dispersions~river_centroid_dispersions_sites[Class == "Mediterranean"]$groups))
summary(aov(river_centroid_dispersions_sites[alteration == "Natural"]$river_centroid_dispersions~river_centroid_dispersions_sites[alteration == "Natural"]$groups))

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
  theme_pca + theme(legend.title = element_blank(), legend.position = "bottom") + scale_x_discrete(labels = c("MedAlt" = "aM","MedNat" = "nM",

                                                                                                        "TempAlt" = "aA", "TempNat" = "nA"))+
  scale_shape_manual(values = c(23,22,25,24))+
  theme(panel.background = element_rect(fill = NA), strip.background.x = element_rect(fill="white"), strip.text = element_text(size = 11), axis.title = element_text(size = 11))+
  theme(text = element_text(size = 11),axis.text = element_text(size = 11, color = "black"), axis.text.x = element_text(hjust = 1),
        axis.title=element_blank(), legend.position = "none", strip.placement ="outside", strip.background = element_blank(), strip.text = element_text(size = 11, color = "black"))


pdf('output/plots/PCAll.pdf', width = 2.4, height = 2.5)
plot(pl)
dev.off()

pca_results$groups.x = factor(pca_results$groups.x, levels = c("TempNat", "TempAlt", "MedNat", "MedAlt"))

individual_pc1 = ggplot(data = pca_results, aes(x=groups.x, y=PC1)) + 
  theme(axis.line.x =  element_line(color="black"),axis.line.y =  element_line(), panel.grid = element_blank(), panel.background = element_blank())+
  geom_boxplot(aes(x=groups.x, y=PC1,  fill=groups.x, group=reorder(site, PC1, median)),  lwd = 0.2, outlier.size = 0.5)+
  labs(y="Axis locations", x=NULL)+scale_y_continuous(limits=c(-4,7.5))+
  scale_fill_manual(values =c("#B4DCED", '#6996D1','#F5CB7D','#F09E41'))+
  theme_pca +theme(legend.title=element_blank(), legend.position = "bottom")+scale_x_discrete(labels=c("MedAlt" = "aM","MedNat" = "nM", "TempAlt" = "aA", "TempNat" = "nA"))+
  theme(panel.background = element_rect(fill = NA), axis.title = element_text(size=11))+
  theme(text = element_text(size=11),axis.text=element_text(size=11, color="black"),
        axis.title=element_text(size=11), legend.position = "none")+ylab("PC1")

individual_pc2 = ggplot(data = pca_results, aes(x=groups.x, y=PC2)) + 
  theme(axis.line.x =  element_line(color="black"),axis.line.y =  element_line(), panel.grid = element_blank(), panel.background = element_blank())+
  geom_boxplot(aes(x=groups.x, y=PC2, fill=groups.x, group=reorder(site, PC2, median)),lwd=0.2, outlier.size=0.5)+
  labs(y="Axis locations", x=NULL)+
  scale_fill_manual(values =c("#B4DCED", '#6996D1','#F5CB7D','#F09E41'))+
  theme_pca +theme(legend.title=element_blank())+scale_x_discrete(labels=c("MedAlt" = "aM","MedNat" = "nM", "TempAlt" = "aA", "TempNat" = "nA"))+
  scale_shape_manual(values=c(23,22,25,24))+scale_y_continuous(limits=c(-5,7.5))+
  theme(panel.background = element_rect(fill = NA), axis.title = element_text(size=11))+
  theme(text = element_text(size=11),axis.text=element_text(size=11, color="black"),
        axis.title=element_text(size=11), legend.position = "none")+ylab("PC2")

pdf('output/plots/PC1.pdf', width = 4, height = 4)
plot(individual_pc1)
dev.off()
pdf('output/plots/PC2.pdf', width = 4, height = 4)
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

pdf('output/plots/PCA_points.pdf', width = 2.7, height = 2.5)
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
 
 pdf('output/plots/PCA_optical.pdf', width = 6.7, height = 6)
 plot(PCA_optical)
 dev.off()
#### end ####

#### Hydrological Indices Analysis ####

source("R/indice_analysis.R")


#### Loading plots ####
parafac_em <- read.csv("data/PARAFAC_data/Em.csv")
parafac_ex<- read.csv("data/PARAFAC_data/Ex.csv")
ggplot()+
  geom_line(aes(x=parafac_em$Em,y=(1/3)*(parafac_em$Comp.8)))+
  geom_line(aes(x=parafac_ex$Ex,y=3*(parafac_ex$Comp.8)), linetype = "dashed")+
  theme_classic()+
  labs(title="Comp8")

# Hydrological plots
source("R/hydrological_analysis.R")

