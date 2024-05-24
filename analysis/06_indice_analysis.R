#### PLSR trial ####
library(pls)
library(plsVarSel)
library(mdatools)
library(data.table)
library(ggplot2)

#### Hydrological Indices Analysis ####
magnitude_ind <- read.csv("data/hydrological_data/Hydra_Hydro_Ind_subs_group1.csv")
magnitude_duration_extreme_ind <- read.csv("data/hydrological_data/Hydra_Hydro_Ind_subs_group2.csv")
timing_extreme_ind <- read.csv("data/hydrological_data/Hydra_Hydro_Ind_subs_group3.csv")
freq_duration_pulses_ind <- read.csv("data/hydrological_data/Hydra_Hydro_Ind_subs_group4.csv")
rate_freq_ind <- read.csv("data/hydrological_data/Hydra_Hydro_Ind_subs_group5.csv")

meta_indices <- as.data.table(Reduce(function(x,y) merge(x=x, y=y, by="Stream"), x = list(magnitude_ind, magnitude_duration_extreme_ind, 
                                                                                     timing_extreme_ind,freq_duration_pulses_ind,rate_freq_ind)))
site_info <- as.data.table(read.csv("data/metafiles/site_summary.csv"))

# We remove the Carrion sample same as in other analysis due to low C concentration
site_info <- site_info[site_info$site!="Carrion",]
meta_indices <- meta_indices[Stream!="Carrion"]  

#Flow indice PCA
pca_data2 <- merge(meta_indices, site_info, by.x = "Stream", by.y = "site")
wine.pca2 <- prcomp(pca_data2[, -c(1,87:95)], scale. = TRUE) 
summary(wine.pca2)
screeplot(wine.pca2)

### PCA plots #### 
PCA_scores <- data.table(wine.pca2$x[,1:2],  pca_data2[,c("Stream", "Class", "alteration", "groups")])
MedAlt <- PCA_scores[PCA_scores$groups == "MedAlt", ][chull(PCA_scores[PCA_scores$groups == "MedAlt", c("PC1", "PC2")]), ]  # hull values for grp A
MedNat <- PCA_scores[PCA_scores$groups == "MedNat", ][chull(PCA_scores[PCA_scores$groups == "MedNat", c("PC1", "PC2")]), ]  # hull values for grp A
TempAlt <- PCA_scores[PCA_scores$groups == "TempAlt", ][chull(PCA_scores[PCA_scores$groups == "TempAlt", c("PC1", "PC2")]), ]  # hull values for grp A
TempNat <- PCA_scores[PCA_scores$groups == "TempNat", ][chull(PCA_scores[PCA_scores$groups == "TempNat", c("PC1", "PC2")]), ]  # hull values for grp A

hull.data <- rbind(MedAlt, MedNat,TempAlt,TempNat)
hull.data$groups <- factor(hull.data$groups, levels = c("TempNat", "TempAlt", "MedNat", "MedAlt"))
pca_data2$groups <- factor(pca_data2$groups, levels = c("TempNat", "TempAlt", "MedNat", "MedAlt"))

PCAloadings <- data.frame(Variables = rownames(wine.pca2$rotation), wine.pca2$rotation)
PCAloadings$Variables <- c("l2","lcv","lca","lkur", "M1","M2","M3", "M4" ,"M5" , "M6", "M7","M8" ,"M9" ,"M10" , "M11" ,"M12", "sdM1", "sdM2", "sdM3" , "sdM4", "sdM5", "sdM6","sdM7" , "sdM8" , "sdM9", "sdM10",      
                        "sdM11","sdM12","X5", "X25","X75", "X95","1LF", "1HF","3LF","3HF", "7LF","7HF", "30LF","30HF", "90LF","90HF", "sd1LF","sd1HF", "sd3LF", "sd3HF", "sd7LF","sd7HF","sd30LF","sd30HF","sd90LF","sd90HF" ,    
                        "ZFD","BFI", "sdZFD","sdBFI","JMin","JMAx", "sdJMIn","sdJMax","Pred","FRE1","FRE3","FRE7","sdFRE1", "sdFRE3","sdFRE7","nPLow","dPLow","nPHigh","dPHigh", "sdnPLow", "sddPLow", "sdnHigh","sddPHigh","nPos","Pos", "nNeg",       
                        "Neg","sdnPos","sdPos","sdnNeg","sdNeg","Rev","sdReversals")

pca.scores<-wine.pca2$x
eigenvec12<-cbind(wine.pca2$rotation[,1],wine.pca2$rotation[,2])
PCAloadings<-data.frame(cor(pca_data2[,-c(1,87:95)],pca.scores))
PCAloadings$Variables <- rownames(PCAloadings)

PCA_scores2 <- data.table(pca_data2, wine.pca2$x)
PCA_rot2 <- data.table(t(cor(PCA_scores2[,c("PC1", "PC2")], pca_data2[,-c(1,87:95)], method = "pearson")), keep.rownames = "Variables")

Indice_PCA1 <- ggplot(pca_data2, aes(x = wine.pca2$x[, 1], y=wine.pca2$x[, 2])) +
  geom_segment(data = PCA_rot2[!startsWith(PCAloadings$Variables, "sd"), ], aes(x = 0, y = 0, xend = (PC1*9.9), yend = (PC2*9.9)), 
               arrow = arrow(length = unit(1/2, "picas")),color = "black", alpha = 0.2) +
  annotate("text", x = (PCA_rot2[!startsWith(PCAloadings$Variables, "sd"), ]$PC1*10), y = (PCA_rot2[!startsWith(PCAloadings$Variables, "sd"), ]$PC2*10),
          label = PCA_rot2[!startsWith(PCAloadings$Variables, "sd"), ]$Variables, size=4, color="black")+
  theme_pca() +
  scale_fill_manual(values = c("#B4DCED", '#6996D1','#F5CB7D','#F09E41'), labels = c("nA", "aA", "nM", "aM")) +
  theme(legend.title = element_blank(), legend.position = "none") +
  labs(color = "Sites", x = "PC 1", y = "PC 2", title = NULL) +
  geom_vline(xintercept = 0, lty = 2) + geom_hline(yintercept = 0, lty = 2) + 
  xlim(-10, 10) + ylim(-10, 10)

pdf('output/plots/Indice_PCA1.pdf', width = 3, height = 3)
plot(Indice_PCA1)
dev.off()

Indice_PCA2 <- Indice_PCA1 <- ggplot(pca_data2, aes(x = wine.pca2$x[, 1], y=wine.pca2$x[, 2])) +
  geom_segment(data = PCA_rot2[startsWith(PCAloadings$Variables, "sd"), ], aes(x = 0, y = 0, xend = (PC1*9.9), yend = (PC2*9.9)), 
               arrow = arrow(length = unit(1/2, "picas")),color = "black", alpha = 0.2) +
  annotate("text", x = (PCA_rot2[startsWith(PCAloadings$Variables, "sd"), ]$PC1*10), y = (PCA_rot2[startsWith(PCAloadings$Variables, "sd"), ]$PC2*10),
           label = PCA_rot2[startsWith(PCAloadings$Variables, "sd"), ]$Variables, size=4, color="black")+
  theme_pca() +
  scale_fill_manual(values = c("#B4DCED", '#6996D1','#F5CB7D','#F09E41'), labels = c("nA", "aA", "nM", "aM")) +
  theme(legend.title = element_blank(), legend.position = "none") +
  labs(color = "Sites", x = "PC 1", y = "PC 2", title = NULL) +
  theme(text = element_text(size = 7), axis.title = element_text(size = 11, color = "black"), axis.text = element_text(size = 11, color = "black")) +
  geom_vline(xintercept = 0, lty = 2) + geom_hline(yintercept = 0, lty = 2) + 
  xlim(-10, 10) + ylim(-10, 10)

pdf('output/plots/Indice_PCA2.pdf', width = 3, height = 3)
plot(Indice_PCA2)
dev.off()

River_scores <- ggplot(pca_data2, aes(x = wine.pca2$x[, 1], y = wine.pca2$x[, 2])) +
  geom_polygon(data = hull.data, mapping = aes(x = PC1,y = PC2,fill = groups, group = groups), alpha = 0.8) +
  geom_point(aes(group = pca_data2$groups, fill = pca_data2$groups, shape = pca_data2$groups), size = 2.5) +
  theme_pca() +
  scale_fill_manual(values = c("#B4DCED", '#6996D1','#F5CB7D','#F09E41'), labels = c("nA", "aA", "nM", "aM")) +
  scale_shape_manual(values = c(23,22,25,24)) +
  theme(legend.title = element_blank(), legend.position = "none") +
  labs(color = "Sites", x = "PC 1", y = "PC 2", title = NULL) +
  theme(text = element_text(size = 7), axis.title = element_text(size = 11, color = "black"), axis.text = element_text(size = 11, color = "black")) +
  geom_vline(xintercept = 0, lty = 2) + geom_hline(yintercept = 0, lty = 2) + xlim(-10, 10) + ylim(-10, 10)

pdf('output/plots/River_scores.pdf', width = 3, height = 3)
plot(River_scores)
dev.off()

plot(wine.pca2$x[, c(1:2)], type = "n", ylim = c(-10, 10), xlim = c(-10,10), cex.main = 1.5, cex.axis = 1.25)
p <- ordihull(ord = wine.pca2$x[, c(1:2)], groups = pca_data2$groups ,display = "sites",draw = "polygon", label = F,
         alpha = 0.7, col = c("#B4DCED"))
summary(p)

#### The indice bpxplots selected from the previous PCA analyiss ####
data <- pca_data2[, c("Stream", "groups", "alteration","Class",  "nPLow", "dPLow", "nPHigh", "dPHigh")]
ggplot() +
  geom_point(data = data, aes(x = Stream, y = nPLow, col = groups), size = 8)

melted_data <- melt(data, id.vars = c("Stream", "groups", "alteration", "Class"), measure.vars = c("nPLow", "dPLow", "nPHigh", "dPHigh"))
melted_data$groups <- factor(melted_data$groups, levels = c("TempNat", "TempAlt", "MedNat", "MedAlt"))

hydro_indices <- ggplot(melted_data)+
  geom_boxplot(aes(x = groups, y = value, fill = groups), width = 0.5)+
  facet_wrap(~variable, scales = "free", strip.position = "left", labeller = as_labeller(c(nPLow="Low Flow Event Count", dPLow = "Low Flow Event Duration", 
                                                                                     nPHigh = "High Flow Event Count", dPHigh = "High Flow Event Duration"))) +
  scale_fill_manual(values = c("#B4DCED", '#6996D1','#F5CB7D','#F09E41'))+
  theme(axis.line.x = element_line(color = "black"),axis.line.y =  element_line(color = "black"), panel.grid = element_blank(), panel.background = element_blank()) +
  theme(legend.position = "none") + ylab("7-Day Maximum Flows (7HF)") +
  scale_x_discrete(limits = c( "TempNat","TempAlt", "MedNat", "MedAlt"),
                   labels = c("nA", "aA", "nM", "aM")) +
  theme(axis.text = element_text(size = 11, color = "black"), axis.title = element_blank(), legend.position = "none", strip.placement = "outside", 
        strip.background = element_blank(), strip.text = element_text(size = 11, color = "black"), panel.spacing = unit(1, "lines"))
  
pdf('output/plots/Indice_boxplots.pdf', width = 5, height = 5)
plot(hydro_indices)
dev.off()

data2 <- pca_data2[, c("Stream", "groups", "alteration","Class", "FRE1", "FRE3", "FRE7")]
melted_data2 = melt(data2, id.vars = c("Stream", "groups", "alteration", "Class"), measure.vars = c("FRE1", "FRE3", "FRE7"))
melted_data2$groups = factor(melted_data2$groups, levels = c("TempNat", "TempAlt", "MedNat", "MedAlt"))

fre_indice <- ggplot(melted_data2)+
  geom_boxplot(aes(x = groups, y = value, fill = groups), width = 0.5)+
  facet_wrap(~variable, scales = "free", strip.position = "left", labeller = as_labeller(c(FRE1 = "Frequency of 1-Day Events", 
                                                                                           FRE3 = "Frequency of 3-Day Events", 
                                                                                           FRE7 = "Frequency of 7-Day Events"))) +
  scale_fill_manual(values = c("#B4DCED", '#6996D1','#F5CB7D','#F09E41')) +
  theme(axis.line.x = element_line(color = "black"),axis.line.y =  element_line(color = "black"), panel.grid = element_blank(), panel.background = element_blank()) +
  theme(legend.position = "none") +
  scale_x_discrete(limits = c( "TempNat","TempAlt", "MedNat", "MedAlt"),
                   labels = c("nA", "aA", "nM", "aM")) +
  theme(axis.text = element_text(size = 11, color = "black"), axis.title = element_blank(), legend.position = "none", strip.placement ="outside", 
        strip.background = element_blank(), strip.text = element_text(size = 11, color = "black"), panel.spacing = unit(1, "lines"))

pdf('output/plots/fre_indice_boxplots.pdf', width = 6, height = 1.5)
plot(fre_indice)
dev.off()
#### End of indice boxplots

#### Statistical tests on the indices ####
shapiro.test((data$nPLow))
hist(data$nPLow)
boxplot(data$nPLow)

indice_results <- vector("list", length = length(names(data[, -c(1:4)])))
names(indice_results) <-  names(data[, -c(1:4)])
for(i in seq_along(data[, -c(1:4)])) {
  names(indice_results[i]) <- names(data[, -c(1:4)])[i]
  indice_results[[i]] <- statistics_pipeline_wrapper(data, names(data[, -c(1:4)])[i], "groups", log_normalise = T)
} 

bartlett.test((nPLow)~groups, data = data)
oneway.test((nPLow)~groups, var.equal = F, data = data)

var.test((nPLow)~alteration, data = data[Class == "Mediterranean"])
t.test((nPLow)~alteration, data=data[Class == "Mediterranean"], var.equal = T)

var.test((nPLow)~alteration, data=data[Class == "Temperate"])
t.test((nPLow)~alteration, data=data[Class == "Temperate"], var.equal = T)

var.test((nPLow)~Class, data = data[alteration == "Natural"])
t.test((nPLow)~Class, data = data[alteration == "Natural"], var.equal = F)

p.adjust(c( 0.001737,  0.4524,   0.603), method = "bonferroni", n = 3)

m <- pairwise.t.test((data$nPLow), data$groups, p.adjust.method = "bonferroni", pool.sd = F)
multcompView::multcompLetters(m$p.value)
### End of Statistical tests


#### PLSR model ####
# get VIP scores and reg coeffs

data <- merge(data.table(pca_data, wine.pca$x[, 1:2]), meta_indices, by.y = "Stream", by.x = "site")
data <- merge(unique(data[, c(1, 3:5, 19:103)]), PC_info, by = c("site", "Class", "alteration"))

shapiro.test(log(data$var_PC2))
hist(log(data$var_PC2))
boxplot(data$var_PC2)

model_varPC1 <- plsr(var_PC1~., data = data[, -c("site", "Class", "alteration", "groups.x", "mean_PC1", "mean_PC2", "var_PC2", "groups")], 
                  scale = T, validation = "CV")
summary(model_varPC1)

model_varPC2 <- plsr(var_PC2~., data = data[, -c("site", "Class", "alteration", "groups.x", "mean_PC1", "mean_PC2", "var_PC1", "groups")], scale = T, validation = "CV")
summary(model_varPC2)

data2 <- merge(Multi_centroid[,-c(3:13,17:22)], meta_indices, by.x = "site", by.y = "Stream")

model_dispPC_all <- plsr(disp_PC_all~., data = data2[, -c("site", "Class", "alteration", "groups")], scale = T, validation = "CV")
summary(model_dispPC_all)

VIP_PC1 <- plsVarSel::VIP(pls.object = model_varPC1, opt.comp = 1, p = dim(model_varPC1$coef)[1])
RC_PC1 <- model_varPC1[["coefficients"]][, , 1]
normalized_RC_PC1 <- RC_PC1 / (max(abs(RC_PC1)))

VIP_PC2 <- plsVarSel::VIP(pls.object = model_varPC2, opt.comp = 1, p = dim(model_varPC2$coef)[1])
RC_PC2 <- model_varPC2[["coefficients"]][, , 1]
normalized_RC_PC2 <- RC_PC2 / (max(abs(RC_PC2)))

VIP_PC_all <- plsVarSel::VIP(pls.object = model_dispPC_all, opt.comp = 1, p=dim(model_dispPC_all$coef)[1])
RC_PC_all <- model_dispPC_all[["coefficients"]][, , 1]
normalized_RC_PC_all <- RC_PC_all / (max(abs(RC_PC_all)))

rc_values = data.table(cbind(VIP_PC1, RC_PC1,VIP_PC2, RC_PC2, VIP_PC_all, RC_PC_all), keep.rownames = "index")
write.table(rc_values, file = "Regression_Coefficient_results.txt", sep=";", row.names = TRUE)


#### PLSR model plot 2 ####
melted_VIP <- rbind(melt(data.table(cbind("VIP" = VIP_PC1, normalized_RC_PC1), keep.rownames = "index"), id.vars = c("index", "VIP")), 
                 melt(data.table(cbind("VIP" = VIP_PC2, normalized_RC_PC2), keep.rownames = "index"), id.vars = c("index", "VIP")),
                 melt(data.table(cbind("VIP" = VIP_PC_all, normalized_RC_PC_all), keep.rownames = "index"), id.vars = c("index", "VIP")))  

used_indexes <- rc_values[VIP_PC1 > 1 | VIP_PC2 > 1 | VIP_PC_all > 1][, "index"]
used_indexes <- rbind(data.table(used_indexes[, 1], variable = "normalized_RC_PC1"), data.table(used_indexes[, 1], variable = "normalized_RC_PC2"), data.table(used_indexes[,1], variable = "normalized_RC_PC_all"))
plot_indexes <- merge(used_indexes,  melted_VIP[VIP > 1 ], all.x = TRUE, key = c("index", "variable"))


for (i in 1:length(plot_indexes$index)){
  
  if (length(grep(plot_indexes$index[[i]],  colnames(magnitude_ind)))>0){
    plot_indexes$groups[[i]]=c("Magnitude of Annual and Monthly Flows")
  } 
  else if (length(grep(plot_indexes$index[[i]],  colnames(magnitude_duration_extreme_ind)))>0){
    plot_indexes$groups[[i]]=c("Magnitude and Duration of Annual Extremes")}
  
  else if (length(grep(plot_indexes$index[[i]],  colnames(timing_extreme_ind)))>0){
    plot_indexes$groups[[i]]=c("Timing of Extreme Flow Events")}
  
  else if (length(grep(plot_indexes$index[[i]],  colnames(freq_duration_pulses_ind)))>0){
    plot_indexes$groups[[i]]=c("Frequency and Duration of High and Low Flow Pulses ")
  }
  else if (length(grep(plot_indexes$index[[i]],  colnames(rate_freq_ind)))>0){
    plot_indexes$groups[[i]]=c("Rate and frequency of flow changes ")}
  
  else  {plot_indexes$groups[[i]]="NA"}
  
}

plot_indexes$index <- factor(plot_indexes$index, levels = c("M1","M2","M3", "M4" ,"M5" , "M6", "M7","M8" ,"M9" ,"M10" , "M11" ,"M12", "sdM1", "sdM2", "sdM3" , "sdM4", "sdM5", "sdM6","sdM7" , "sdM8" , "sdM9", "sdM10",      
                                                        "sdM11","sdM12","l2","lcv","lca","lkur","X5", "X25","X75", "X95","X1LF", "X1HF","X3LF","X3HF", "X7LF","X7HF", "X30LF","X30HF", "X90LF","X90HF", "sd1LF","sd1HF", "sd3LF", "sd3HF", "sd7LF","sd7HF","sd30LF","sd30HF","sd90LF","sd90HF" ,    
                                                        "ZFD","BFI", "sdZFD","sdBFI","JMin","JMAx", "sdJMIn","sdJMax","Pred","FRE1","FRE3","FRE7","sdFRE1", "sdFRE3","sdFRE7","nPLow","dPLow","nPHigh","dPHigh", "sdnPLow", "sddPLow", "sdnHigh","sddPHigh","nPos","Pos", "nNeg",       
                                                        "Neg","sdnPos","sdPos","sdnNeg","sdNeg","Rev","sdReversals"), 
                           labels = c("January Flow","February Flow","March Flow", "April Flow" ,"May Flow" , "June Flow", "July Flow","August Flow" ,"September Flow" ,"October Flow" , "November Flow" ,"December Flow", 
                                      "January Flow SD ", "February Flow SD", "March Flow SD" , "April Flow SD", "May Flow SD", "June Flow SD","July Flow SD" , "August Flow SD" , "September Flow SD", "October Flow SD", "November Flow SD","December Flow SD",
                                      "Variance of the flow duration curve","Coefficient of variation of the flow duration curve",
                                      "Skewness of the flow duration curve","Kurtosis of the flow duration curve", 
                                      "High 5% Exceedence Flows", "High 25% Exceedence Flows","Low 75% Exceedence Flows", "Low 95% Exceedence Flows",
                                      "1-Day Minimum", "1-Day Maximum","3-Day Minimum","3-Day Maximum", "7-Day Minimum","7-Day Maximum", "30-Day Minimum","30-Day Maximum", "90-Day Minimum","90-Day Maximum", 
                                      "1-Day Minimum SD","1-Day Maximum SD", "3-Day Minimum SD", "3-Day Maximum SD", "7-Day Minimum SD","7-Day Maximum SD","30-Day Minimum SD","30-Day Maximum SD","90-Day Minimum SD","90-Day Maximum SD" ,    
                                      "Zero Flow Days","Base FLow Index", "Zero Flow Days SD","Base FLow Index SD","Date of Annual Minimum Flow","Date of Annual Maximum Flow", "Date of Annual Minimum Flow SD","Date of Annual Maximum Flow SD","Predictability",
                                      "1-Day Flood Frequency","3-Day Flood Frequency","7-Day Flood Frequency","1-Day Flood Frequency SD", "3-Day Flood Frequency SD","7-Day Flood Frequency SD",
                                      "Low Pulse Count","Low Pulse Duration","High Pulse Count","High Pulse Duration", "Low Pulse Count SD", "Low Pulse Duration SD", "High Pulse Count SD","High Pulse Duration SD",
                                      "No Increasing Flow Days","Rise Rate", "No Decreasing Flow Days",  "Fall Rate","No Increasing Flow Days SD","Rise Rate SD","No Decreasing Flow Days SD","Fall Rate SD","Reversal Count","Reversal Count SD"), )

plot_indexes[is.na(plot_indexes$value)]$value <- 0

PLSR_All1 <- ggplot(plot_indexes[!sapply(plot_indexes$index, FUN = grepl, pattern = "SD")]) +
  geom_col(aes(x = index, y = value, fill = variable), position = position_dodge(),
           color = "black", width = 0.7,  lwd = 0.2) +
  theme_pca() +
  theme(axis.text.x = element_text(angle = 65, hjust = 1),
        legend.position = "top", 
        panel.grid.major = element_line(color = "grey")) +
  scale_fill_manual(values = c("#eeca8e", "#9c8fb1", "#d1d1d1"),
                    limits = c("normalized_RC_PC_all", "normalized_RC_PC1", "normalized_RC_PC2"),
                    labels = c("PCA Dispersion", "PC1 Variance", "PC2 Variance")) +
  ylab("Normalized Regression Coefficient") + ylim(-1, 1) +
  guides(fill = guide_legend(title = "Model"))
  

PLSR_All2 <- ggplot(plot_indexes[sapply(plot_indexes$index, FUN = grepl, pattern = "SD")]) +
  geom_col(aes(x = index, y = value, fill = variable), position = position_dodge(),
           color = "black", width = 0.6,  lwd = 0.2) +
  theme_pca() +
  theme(axis.text.x = element_text(angle = 65, hjust = 1),
        legend.position = "top", 
        panel.grid.major = element_line(color = "grey")) +
  scale_fill_manual(values = c("#eeca8e", "#9c8fb1", "#d1d1d1"),
                    limits = c("normalized_RC_PC_all", "normalized_RC_PC1", "normalized_RC_PC2"),
                    labels = c("PCA Dispersion", "PC1 Variance", "PC2 Variance")) +
  ylab("Bormalized Regression Coefficient") + ylim(-1, 1) +
  guides(fill = guide_legend(title = "Model"))
# "#665191", "#dd5182", "#ffa600"

pdf('output/plots/PLSR_All1.pdf', width = 10, height = 6)
plot(PLSR_All1)
dev.off()
pdf('output/plots/PLSR_All2.pdf', width = 10, height = 6)
plot(PLSR_All2)
dev.off()

#### End of PLSR model plot 2
















