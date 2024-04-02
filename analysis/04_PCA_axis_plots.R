library(ggplot2)

##### Figure 6A -  All PCAs Dispersion #####

pl = ggplot(data = disp_melted[variable == "disp_PC_all"], aes(x = groups.x, y = value)) + 
  geom_boxplot(aes(x = groups.x, y = value, fill = groups.x), outlier.size = 0.5, lwd = 0.2)+
  # geom_point(aes(fill=groups.x, shape=groups.x),fill="black", size=3)+
  #facet_wrap(~variable,  labeller=as_labeller(c(disp_PC1="PC1", disp_PC2="PC2", disp_PC_all="Dispersion on the 2 PC space")))+
  labs(y = "Dispersion", x = NULL)+
  scale_fill_manual(values = c("#B4DCED", '#6996D1','#F5CB7D','#F09E41'), labels = c("aM", "nM", "aA", "nA")) +
  theme_pca + 
  scale_x_discrete(labels = c("MedAlt" = "aM","MedNat" = "nM", "TempAlt" = "aA", "TempNat" = "nA")) +
  scale_shape_manual(values = c(23,22,25,24))+
  theme(panel.background = element_rect(fill = NA),
        panel.grid = element_blank(),
        strip.background = element_blank(), 
        strip.text = element_text(size = 11, color = "black"), 
        strip.placement ="outside")+
  theme(text = element_text(size = 11),
        axis.line.x =  element_line(color = "black"),
        axis.line.y =  element_line(),
        axis.text = element_text(size = 11, color = "black"), 
        axis.text.x = element_text(hjust = 1),
        axis.title = element_blank(), 
        legend.position = "none") 


##### Figure 6B -  PC1 Dispersion #####
pca_results$groups.x = factor(pca_results$groups.x, 
                              levels = c("TempNat", "TempAlt", "MedNat", "MedAlt"))

individual_pc1 = ggplot(data = pca_results, aes(x=groups.x, y=PC1)) + 
  geom_boxplot(aes(x=groups.x, y=PC1,  fill=groups.x, group=reorder(site, PC1, median)),  lwd = 0.2, outlier.size = 0.5)+
  scale_y_continuous(limits=c(-4,7.5))+
  scale_fill_manual(values =c("#B4DCED", '#6996D1','#F5CB7D','#F09E41'))+
  scale_x_discrete(labels=c("MedAlt" = "aM","MedNat" = "nM", "TempAlt" = "aA", "TempNat" = "nA")) +
  theme_pca +
  ylab("PC1") +
  theme(panel.background = element_rect(fill = NA),
        panel.grid = element_blank(),
        strip.background = element_blank(), 
        strip.text = element_text(size = 11, color = "black"), 
        strip.placement ="outside")+
  theme(text = element_text(size = 11),
        axis.line.x =  element_line(color = "black"),
        axis.line.y =  element_line(),
        axis.text = element_text(size = 11, color = "black"), 
        axis.text.x = element_text(hjust = 1),
        axis.title = element_text(color = "black"), 
        legend.position = "none") 

##### Figure 6C -  PC1 Dispersion #####
individual_pc2 = ggplot(data = pca_results, aes(x = groups.x, y = PC2)) + 
  geom_boxplot(aes(x = groups.x, y = PC2, fill = groups.x, group = reorder(site, PC2, median)),lwd = 0.2, outlier.size=0.5)+
  scale_fill_manual(values = c("#B4DCED", '#6996D1','#F5CB7D','#F09E41'))+
  scale_x_discrete(labels=c("MedAlt" = "aM","MedNat" = "nM", "TempAlt" = "aA", "TempNat" = "nA"))+
  scale_shape_manual(values = c(23,22,25,24))+
  scale_y_continuous(limits = c(-5,7.5))+
  theme_pca +
  ylab("PC2") +
  theme(legend.title=element_blank())+
  theme(panel.background = element_rect(fill = NA),
        panel.grid = element_blank(),
        strip.background = element_blank(), 
        strip.text = element_text(size = 11, color = "black"), 
        strip.placement ="outside")+
  theme(text = element_text(size = 11),
        axis.line.x =  element_line(color = "black"),
        axis.line.y =  element_line(),
        axis.text = element_text(size = 11, color = "black"), 
        axis.text.x = element_text(hjust = 1),
        axis.title = element_text(color = "black"), 
        legend.position = "none") 

                                                                                                              
                                                                                                              