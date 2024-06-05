library(ggplot2)

plot_all <- ggplot(data_sum, aes(x = wine.pca$x[,1], y = wine.pca$x[,2])) +
  geom_point(aes(color = data_sum$groups.x, fill = data_sum$groups.x, shape = data_sum$groups.x), 
             size = 2.5, colour = "black") +
  scale_shape_manual(values = c(23,22,25,24)) +
  scale_fill_manual(values = c("#B4DCED", '#6996D1','#F5CB7D','#F09E41')) +
  geom_segment(data = PCA_rot, aes(x = 0, y = 0, xend = (PC1*8), yend = (PC2*8)), arrow = arrow(length = unit(1/2, "picas")),color = "black") +
  annotate("text", x = (PCA_rot$PC1*8), y = (PCA_rot$PC2*8),
           label = PCA_rot$Variables, size = 4, color = "black")+
  labs(color = "Sites", x = "PC1 (44.2.8%)", y = "PC2 (15.6%)", tag = "a. Loadings and sites")+
  theme_pca() + 
  theme(plot.tag.position = c(0.24,0.97)) +
  scale_x_continuous(limits = c(-8,8), n.breaks=10) +
  scale_y_continuous(limits = c(-8,8), n.breaks=10) +
  geom_vline(xintercept = 0, lty=2) + geom_hline(yintercept = 0, lty = 2)

plot_points <- ggplot(PCA_results, aes(x = wine.pca$x[, 1], y = wine.pca$x[, 2])) +
  geom_polygon(data = hull.data, 
               mapping = aes(x = PC1, y = PC2, color = groups.x, group = groups.x), 
               fill = NA,lwd = 1, alpha = 0.7) +
  geom_point(aes(fill = PCA_results$groups, shape = PCA_results$groups), 
             size = 1, alpha = 0.7,  color = "black", stroke = 0.25) +
  scale_shape_manual(values = c(23,22,25,24), 
                     labels = c("MedAlt" = "aM","MedNat" = "nM", "TempAlt" = "aA", "TempNat" = "nA")) +
  scale_fill_manual(values = c( "#B4DCED",'#6996D1','#F5CB7D','#F09E41'),
                    labels = c("MedAlt" = "aM","MedNat" = "nM", "TempAlt" = "aA", "TempNat" = "nA")) +
  scale_color_manual(values = c( "#B4DCED",'#6996D1','#F5CB7D','#F09E41'), 
                     labels = c("MedAlt" = "aM","MedNat" = "nM", "TempAlt" = "aA", "TempNat" = "nA")) +
  labs(color = "Sites", x = "PC1 (44%)", y = "PC2 (16%)")+ #, tag = "A"
  theme_pca() +
  theme(plot.tag.position = c(0.15, 0.96))+
  scale_x_continuous(limits = c(-7.5, 7.5), n.breaks = 10)+
  scale_y_continuous(limits = c(-7.5, 7.5), n.breaks = 10)+
  geom_vline(xintercept = 0, lty = 2) + geom_hline(yintercept = 0, lty=2)+
  theme(legend.title = element_blank(),
        legend.position = "none", 
        text=element_text(size = 11))+
  theme(axis.text =  element_text(size = 11, color = "black"),
        axis.title = element_text(size = 11)) +
  guides(shape = guide_legend(override.aes = list(size = 5, alpha = 1)))

plot_optical <- ggplot(data_sum, aes(x = wine.pca$x[,1], y = wine.pca$x[,2])) +
  geom_segment(data = PCA_rot, aes(x = 0, y = 0, xend = (PC1), yend = (PC2)), 
               arrow = arrow(length = unit(1/2, "picas")),color = "black") +
  annotate("text", x = (PCA_rot$PC1), y = (PCA_rot$PC2),
           label = PCA_rot$Variables, 
           size = 3.5, 
           color = "black")+
  labs(color = "Sites", x = "PC1 (47.1%)", y = "PC2 (13.6%)" )+
  theme_pca() + 
  theme(plot.tag.position = c(0.20,0.95)) +
  scale_x_continuous(limits = c(-1,1), n.breaks = 10) +
  scale_y_continuous(limits = c(-1,1), n.breaks = 10) +
  theme(legend.position = c(-1,0))+
  geom_vline(xintercept = 0, lty = 2) + 
  geom_hline(yintercept = 0, lty = 2) +
  theme(axis.text = element_text(size = 11), 
        axis.title = element_text(size = 11)) +
  theme(text = element_text(size = 11))

PCA_optical <- ggplot(data_sum, aes(x = wine.pca$x[,1], y = wine.pca$x[,2]))+
  geom_segment(data = PCA_rot, aes(x = 0, y = 0, xend = (PC1), yend = (PC2)), arrow = arrow(length = unit(1/2, "picas")),color = "black", lwd = 0.7) +
  #geom_segment(data=LC_OCD_rot, aes(x = 0, y = 0, xend = (PC1), yend = (PC2)), arrow = arrow(length = unit(1/2, "picas")),color = "blue", lwd=0.7) + 
  #annotate("text", x = as.vector(labels_locator2[1:14,1]), y = as.vector(labels_locator[1:14,2]),label =var_names[1:14], size=5, color="black")+
  # annotate("text", as.vector(labels_locator2[15:17,1]), y = as.vector(labels_locator[15:17,2]), label =var_names[15:17], size=5, color="blue")+
  labs(color = "Sites", x = "PC1 (44%)", y = "PC2 (16%)") +
  theme_pca() + 
  theme(plot.tag.position = c(0.19, 0.95)) +
  scale_x_continuous(limits = c(-1, 1), n.breaks = 8) +
  scale_y_continuous(limits = c(-1, 1), n.breaks = 8) +
  guides(fill = "legend") +
  theme(legend.position = c(-1, 0) ) +
  geom_vline(xintercept = 0, lty = 2) + geom_hline(yintercept = 0, lty = 2)+
  theme(text = element_text(size = 20),axis.text=element_text(size = 20, color = "black"))


