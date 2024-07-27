
##### Figure 5 -  Polygon plots #####

#png("PCA_polygons.pdf", width = 10, height = 12, units = "in", res = 300)  # Adjust resolution as needed

xlim <- c(-4, 7.5)
ylim <- c(-4, 6.5)
par(mfrow = c(2,2), mai = c(0.5, 0.5, 0.3, 0.3))

plot_TempNat <- plot(wine.pca$x[, c(1:2)], type = "n", ylim = ylim, xlim = xlim, cex.main = 1.5, cex.axis = 1.25) #, main = "Natural Atlantic"
ordi_TempNat <- ordihull(ord = wine.pca$x[, c(1:2)], groups = data_sum$site, display = "sites", draw = "polygon", label = F, show.groups = group_list[[2]]$site,
                         alpha = 0.7, col = c("#B4DCED"))
abline(v = 0, h = 0, lty = 2)

plot_TempAlt <- plot(wine.pca$x[, c(1:2)], type="n",  ylim = ylim, xlim = xlim, cex.main = 1.5, cex.axis = 1.25) #main = "Altered Atlantic",
ordi_TempAlt <- ordihull(ord = wine.pca$x[, c(1:2)], groups = data_sum$site, display = "sites", draw = "polygon", label = F, show.groups = group_list[[1]]$site,
                         alpha = 0.7, col = c('#6996D1'))
abline(v = 0, h = 0, lty = 2)

plot_MedNat <- plot(wine.pca$x[, c(1:2)], type="n", ylim = ylim, xlim = xlim, cex.main = 1.5, cex.axis = 1.25) # main = "Natural Mediterranean",
ordi_MedNat <- ordihull(ord = wine.pca$x[, c(1:2)], groups = data_sum$site, display = "sites", draw = "polygon", label = F, show.groups = group_list[[4]]$site,
                        alpha = 0.7, col = c('#F5CB7D'))
abline(v = 0, h = 0, lty = 2)

plot_MedAlt <- plot(wine.pca$x[, c(1:2)], type = "n", ylim = ylim, xlim = xlim, cex.main = 1.5, cex.axis =1.25) #, main = "Altered Mediterranean"
ordi_MedAlt <- ordihull(ord = wine.pca$x[, c(1:2)], groups = data_sum$site, display = "sites", draw = "polygon", label = F, show.groups = group_list[[3]]$site,
                        alpha = 0.7, col = c("#F09E41"))
abline(v = 0, h = 0, lty = 2)

par(mfrow = c(1,1))

#dev.off()

##### 