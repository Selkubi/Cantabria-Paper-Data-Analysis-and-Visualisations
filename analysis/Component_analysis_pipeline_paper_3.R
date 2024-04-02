# library(vegan)
 library(geosphere)
# library(gridExtra)
# library(data.table)
# library(ggplot2)
# library(lubridate)

#### The functions to be used are here, but first you have to load the project package
roxygen2::roxygenize()

##### 1. Data Cleaning #####
source("analysis/01_data_cleaning.R")

##### 2. DOC Mean and CV plots #####
source("analysis/02_DOC_DOM_plots.R")

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

#####  4. Statistical tests for DOC concentration CV values #####
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

##### 7. PCA data prep #####
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

##### 8. PCA polygon plots #####
PCAloadings <- data.frame(Variables = rownames(wine.pca$rotation), wine.pca$rotation)
source("analysis/03_PCA_polygon_plots.R")
dev.off()

##### 9. PCA PC1 and PC2 centroid calculations #####
# Centroid data (longitude and latitude) is extracted from here and named with centroid_<name>, and stored as a data.table within centroid
for(i in seq_along(group_list)){
assign(paste0("centroid_",  group_list[[i]][1, groups]), 
       do.call(rbind.data.frame, lapply(get(paste0("ordi_", (group_list[[i]][1, groups]))), centroid)
))
}

centroids <- data.table(rbind(centroid_MedAlt, centroid_MedNat, centroid_TempAlt, centroid_TempNat), keep.rownames = "site")
centroids <- data.table(merge(centroids, site_info, by = "site")) # Can also do this whole thing with summary(ordi_MedNat)

# The plots also produce ordihull objects that has the area which can be accessed with eg: summary(ordi_MedNat) 
polygon_data <- data.table(t(cbind(summary(ordi_MedNat), summary(ordi_MedAlt), summary(ordi_TempAlt), summary(ordi_TempNat))), keep.rownames = "site")
polygon_data <- merge(polygon_data, site_info, by = "site")

PCA_results <- cbind(wine.pca$x[, c(1:2)], pca_data)

# complete dataset of PCA information
PC_info <- merge(PCA_results[, .(mean_PC1 = mean(PC1, na.rm = TRUE), mean_PC2 = mean(PC2, na.rm = TRUE), 
                                var_PC1 = var(PC1), var_PC2 = var(PC2)), by = .(site)], 
                site_info[, c("site", "Class", "alteration", "groups")])

##### 10. Statistical tests on PCA polygons #####
# Mean of PC1
shapiro.test((PC_info$mean_PC1))
hist((PC_info$mean_PC1))
boxplot((PC_info$mean_PC1))

# Calculate the pairwise t-test for all 4 variables
pairwise_results_mean <- lapply(PC_info[, c(2:5)], pairwise.t.test, DOC_sum$groups, p.adjust.methods = "bonferroni", pool.sd = T)

oneway_test_results(PC_info, "mean_PC1", "groups", log_normalise = F)
statistics_pipeline_wrapper(PC_info, "mean_PC1", "groups", log_normalise = F)
#multcompView::multcompLetters(rcompanion::fullPTable(pairwise_results_mean[["mean_PC1"]]$p.value))$Letters

# CV of PC1
shapiro.test((PC_info$var_PC1))
hist((PC_info$var_PC1))

oneway_test_results(PC_info, "var_PC1", "groups", log_normalise = F)
statistics_pipeline_wrapper(PC_info, "var_PC1", "groups", log_normalise = F)

# Mean of PC2
shapiro.test((PC_info$mean_PC2))
hist((PC_info$mean_PC2))

oneway_test_results(PC_info, "mean_PC2", "groups", log_normalise = F)
statistics_pipeline_wrapper(PC_info, "mean_PC2", "groups", log_normalise = F)

# CV of PC2
shapiro.test(log(PC_info$var_PC2))
hist(log(PC_info$var_PC2))

# The var_PC2 value needs to be log transformed, so we re-do the pairwise test for the letters
pairwise_results_var_PC2 <- lapply(log(PC_info[, c(5)]), pairwise.t.test, DOC_sum$groups, p.adjust.methods = "bonferroni", pool.sd = T)

oneway_test_results(PC_info, "var_PC2", "groups", log_normalise = T)
statistics_pipeline_wrapper(PC_info, "var_PC2", "groups", log_normalise = T)
# Do for the exact p value which is < 0.05 var.test(log(var_PC2)~groups, data = PC_info[Class == "Temperate"])

##### 11. Global test results of the PCA including all of the PC axis #####
betas_PC_all <- betadisper(vegdist((wine.pca$x[,1:10]), method = "euclidean"),group = pca_data$site, type="centroid")
disp_PC_all <- tapply(betas_PC_all[["distances"]], betas_PC_all[["group"]], mean)
Multi_centroid <- data.table(cbind(disp_PC_all, betas_PC_all[["centroids"]], match(rownames(betas_PC_all[["centroids"]]), rownames(disp_PC_all))), keep.rownames = "site")
Multi_centroid <- Multi_centroid[site_info, on = .(site)]

# Mean centroid of all PC axis
gb_fr <- adonis2(vegdist(Multi_centroid[,3:12], method = "euclidian")~Multi_centroid$groups, permutations = 100000)
gb_me <- adonis2(vegdist(Multi_centroid[Class == "Mediterranean", c(3:12)], method = "euclidian")~alteration, data = Multi_centroid[Class == "Mediterranean"], permutations = 100000)
gb_te <- adonis2(vegdist((Multi_centroid[Class == "Temperate", c(3:12)]), method = "euclidian")~alteration, data = Multi_centroid[Class == "Temperate"], permutations = 100000)
gb_nat <- adonis2(vegdist((Multi_centroid[alteration == "Natural",c(3:12)]), method = "euclidian")~Class, data = Multi_centroid[alteration == "Natural"], permutations = 10000)

mean_DOM_composition_results <- rbind(get_adonis_results(gb_fr), do_bonferroni_to_adonis(gb_nat, gb_te, gb_me))

# Temporal dispersion (CV equivalent) of all PC axis
shapiro.test(log(Multi_centroid$disp_PC_all))
hist(log(Multi_centroid[Class == "Temperate"]$disp_PC_all))
qqnorm(log(Multi_centroid$disp_PC_all), main="Normal Q-Q Plot of male");qqline((Multi_centroid$disp_PC_all))

oneway_test_results(Multi_centroid, "disp_PC_all", "groups", log_normalise = T)
statistics_pipeline_wrapper(Multi_centroid, "disp_PC_all", "groups", log_normalise = T)

# Among river variation of mean DOM composition (dispersion of river centroids)
summary(aov(betadisper(vegdist((Multi_centroid[, c(3:12)]), method = "euclidian"), Multi_centroid$groups)$distances~Multi_centroid$groups))
##### 12. PCA Mean/Dispersion boxplots ####
source("analysis/04_PCA_axis_plots.R")
pl
individual_pc1 
individual_pc2 

pdf('output/plots/PCAll.pdf', width = 2.4, height = 2.5)
plot(pl)
dev.off()

pdf('output/plots/PC1.pdf', width = 4, height = 4)
plot(individual_pc1)
dev.off()

pdf('output/plots/PC2.pdf', width = 4, height = 4)
plot(individual_pc2)
dev.off()


##### 13. PCA ordihull calculations ####
PCA_scores <- data.table(pca_data, wine.pca$x)
PCA_rot <- data.table(t(cor(PCA_scores[,c("PC1", "PC2")], pca_data[,-(1:5)], method = "pearson")), keep.rownames = "Variables")

PCAloadings <- data.table(Variables = rownames(wine.pca$rotation), wine.pca$rotation)

pca.scores<-wine.pca$x
eigenvec12<-cbind(wine.pca$rotation[,1],wine.pca$rotation[,2])
PCAloadings<-data.frame(cor(pca_data[,-(1:5)],pca.scores))
PCAloadings$Variables=rownames(PCAloadings)

MedAlt <- PCA_scores[PCA_scores$groups == "MedAlt", ][chull(PCA_scores[PCA_scores$groups == "MedAlt", c("PC1", "PC2")]), ]  # hull values for grp A
MedNat <- PCA_scores[PCA_scores$groups == "MedNat", ][chull(PCA_scores[PCA_scores$groups == "MedNat", c("PC1", "PC2")]), ]  # hull values for grp A
TempAlt <- PCA_scores[PCA_scores$groups == "TempAlt", ][chull(PCA_scores[PCA_scores$groups == "TempAlt", c("PC1", "PC2")]), ]  # hull values for grp A
TempNat <- PCA_scores[PCA_scores$groups == "TempNat", ][chull(PCA_scores[PCA_scores$groups == "TempNat", c("PC1", "PC2")]), ]  # hull values for grp A

hull.data <- rbind(MedAlt, MedNat,TempAlt,TempNat)  #combine grp.a and grp.b
hull.data

# PCA with all the data
source("analysis/05_PCA_plots.R")
plot_all

hull.data$groups <- factor(hull.data$groups, levels = c("TempNat", "TempAlt", "MedNat", "MedAlt"))

plot_points 
pdf('output/plots/PCA_points.pdf', width = 2.7, height = 2.5)
plot(plot_points)
dev.off()

plot_optical # we took out just the loading arrows and put the variable names manually using Inkscape 
PCA_optical 
 pdf('output/plots/PCA_optical.pdf', width = 6.7, height = 6)
 plot(PCA_optical)
 dev.off()
#### End of PCA related calculations and plots ####

#### 14. Hydrological Indices Analysis ####
source("analysis/Indice_analysis.R")

#### 14. Loading plots ####
source("analysis/06_loading_plots.R")

# Hydrological plots
source("R/hydrological_analysis.R")
 