library(ggplot2)

##### Figure 6A -  All PCAs Dispersion #####

pl = ggplot(data = Multi_centroid, aes(x = groups, y = disp_PC_all)) + 
  geom_boxplot(aes(fill = groups), outlier.size = 0.5, lwd = 0.2) +
  labs(y = "Dispersion", x = NULL) +
  theme_pca() + 
  scale_x_discrete(limits = c("TempNat", "TempAlt", "MedNat","MedAlt"),
                   labels = c("nA", "aA", "nM", "aM")) +
  scale_fill_manual(values = c('#F09E41', '#F5CB7D', '#6996D1', "#B4DCED")) +
  scale_shape_manual(values = c(23,22,25,24)) +
  theme_boxplot()


##### Figure 6B -  PC1 Dispersion #####
PCA_results$alteration_type_grouping = factor(PCA_results$alteration_type_grouping, 
                              levels = c("TempNat", "TempAlt_irrigation", "TempAlt_hydropower", "MedNat", "MedAlt_irrigation"))

individual_pc1 = ggplot(data = PCA_results, aes(x = groups.x, y = PC1)) + 
  geom_boxplot(aes(x = groups.x, y = PC1,  fill = alteration_type_grouping, group = reorder(site, PC1, median)),  lwd = 0.2, outlier.size = 0.5) +
  scale_y_continuous(limits = c(-5, 7.5)) +
  scale_fill_manual(values =c("#B4DCED", '#6996D1', "#2B5FA2", '#F5CB7D','#F09E41'))+
  scale_x_discrete(labels = c("MedAlt" = "aM","MedNat" = "nM", "TempAlt" = "aA", "TempNat" = "nA")) +
  theme_pca() + theme(legend.position = "top")
  ylab("PC1") +
  theme_boxplot()

##### Figure 6C -  PC1 Dispersion #####
individual_pc2 = ggplot(data = PCA_results, aes(x = groups.x, y = PC2)) + 
  geom_boxplot(aes(x = groups.x, y = PC2, fill = groups.x, group = reorder(site, PC2, median)),lwd = 0.2, outlier.size = 0.5)+
  scale_fill_manual(values = c("#B4DCED", '#6996D1','#F5CB7D','#F09E41'))+
  scale_x_discrete(labels = c("MedAlt" = "aM","MedNat" = "nM", "TempAlt" = "aA", "TempNat" = "nA"))+
  scale_shape_manual(values = c(2.5,7.5))+
  theme_pca() +
  ylab("PC2") +
  theme_boxplot()


                                                                                                              
                                                                                                              