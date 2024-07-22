library(ggplot2)
library(data.table)
library(vegan)

##### Figure 4.A - mean DOC values of rivers #####
mean_NPOC = ggplot(data_sum, aes(x = groups.x, y = (NPOC))) +
  theme(axis.line.x =  element_line(color = "black"), 
        axis.line.y =  element_line(),
        panel.grid = element_blank(),
        panel.background = element_blank()) +
  geom_boxplot(mapping = aes(fill = alteration_type_grouping , group = reorder(site, NPOC, median)), varwidth = T, width = 1, lwd = 0.3, outlier.size = 0.5)+
  scale_fill_manual(values = c("#B4DCED", '#6996D1',"#2B5FA2",'#F5CB7D','#F09E41')) +
  scale_x_discrete(limits = c("TempNat", "TempAlt", "MedNat","MedAlt"),
                   labels = c("nA", "aA", "nM", "aM")) +
  theme(axis.text = element_text(color = "black", size = 11), axis.title.x = element_blank(), legend.position = "bottom") +
  ylab("DOC mg C/L")

##### Figure 4.B - mean Mediterranean DOC values #####

# first we create the montly_means data.table for the averaged values for both of the figures
monthly_means =  data_sum[, .(monthly_mean = mean(NPOC, na.rm = TRUE)), by = .(campaign, groups.x)]

temperate_DOC = ggplot2::ggplot(data_sum[Class == 'Temperate'], aes(x = campaign, y = NPOC)) +
  geom_line(data = monthly_means[groups.x == 'TempAlt' | groups.x == 'TempNat'], 
            aes(x = campaign, y = monthly_mean, group = groups.x, color = groups.x), lwd = 1.5) +
  scale_color_manual(values = c("#B4DCED", '#6996D1')) +
  ggnewscale::new_scale_color() +
  geom_line(aes(group = site, color = alteration_type_grouping), lwd = 0.5) +
  geom_point(aes(group = site, shape = alteration_type_grouping, fill = alteration_type_grouping), size = 2, color = 'black') +
  scale_shape_manual(values = c(23, 22, 21)) +
  scale_color_manual(values = c("#B4DCED", '#6996D1', "#2B5FA2")) +
  scale_fill_manual(values = c("#B4DCED", '#6996D1', "#2B5FA2")) +
  scale_x_discrete(limits = c("feb", "may", "oct","apr", "aug", "dec"),
                   labels = c("Feb", "May", "Oct", "Apr", "Aug", "Dec")) +
  theme_pca() +
  theme(axis.title.x = element_blank()) +
  ylab("DOC mg C/L")

##### Figure 4.C - mean Mediterranean DOC values #####
monthly_means =  data_sum[, .(monthly_mean = mean(NPOC, na.rm = TRUE)), by = .(campaign, groups.x)]

mediterranean_DOC = ggplot(data_sum[Class == 'Mediterranean'], aes(x = campaign, y = NPOC)) +
  geom_line(data = monthly_means[groups.x == 'MedAlt' | groups.x == 'MedNat'], 
            aes(x = campaign, y = monthly_mean, group = groups.x, color = groups.x), lwd = 1.5) +
  geom_line(aes(group = site, color = groups.x), lwd = 0.5) +
  geom_point(aes(group = site, shape = groups.x, fill = groups.x), size = 2, color = 'black') +
  scale_shape_manual(values = c(25, 24)) +
  scale_color_manual(values = c('#F5CB7D','#F09E41')) +
  scale_fill_manual(values = c('#F5CB7D','#F09E41')) +
  scale_x_discrete(limits = c("feb", "may", "oct","apr", "aug", "dec"),
                   labels = c("Feb", "May", "Oct", "Apr", "Aug", "Dec")) +
  theme_pca() + 
  theme(axis.title.x = element_blank()) +
  ylab("DOC mg C/L")

