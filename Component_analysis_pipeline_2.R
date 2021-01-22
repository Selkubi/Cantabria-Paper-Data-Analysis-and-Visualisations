library(tidyverse)
library(vegan)
library(geosphere)
library(lubridate)
library(tidyselect)
library(gridExtra)

#### Taking the data####
sample_loadings <- read_csv("denormalised_sample_loadings.csv")
sample_names <- read_csv("sample_name_site_match.csv")
site_info <- read_csv("site_summary.csv")

oct1<- read_delim("abs_data/oct.txt", delim=",") %>% cbind(campaign="oct") 
feb1<- read_delim("abs_data/feb.txt", delim=",") %>% cbind(campaign="feb")
apr1<- read_delim("abs_data/apr.txt", delim=",") %>% cbind(campaign="apr")
may1<- read_delim("abs_data/may.txt", delim=",") %>% cbind(campaign="may")
aug1<- read_delim("abs_data/aug.txt", delim=",") %>% cbind(campaign="aug")
dec1<- read_delim("abs_data/dec.txt", delim=",") %>% cbind(campaign="dec")

oct1$site <- as.character(lapply(oct1$site, str_sub, 1, -5))
apr1$site <- as.character(lapply(apr1$site, str_sub, 1, -5))

abs_data <- bind_rows(oct1,feb1, apr1, may1, aug1,dec1) %>%
  filter(!str_detect(site, "^Blank"))

abs_data <- abs_data %>%
  group_by(site,campaign) %>%
  summarise(DecAbsCoeff254=mean(DecAbsCoeff254, na.rm=T),
            FIX=mean(FIX, na.rm=T),
            FIX2=mean(FIX2, na.rm=T),
            HIX=mean(HIX, na.rm=T), #took (, trim=0.20) out for now cause porma december and tiron october results looks weird. 
            HIX2=mean(HIX2, na.rm=T),
            beta.alpha=mean(beta.alpha, na.rm=T),
            DecAbsCoeff420=mean(DecAbsCoeff420, na.rm=T),
            DecAbsCoeff430=mean(DecAbsCoeff430, na.rm=T),
            DecAbsCoeff436=mean(DecAbsCoeff436, na.rm=T),
            E2.to.E3=mean(E2.to.E3, na.rm=T), #I took out manually two data points of this (Tiron1C May and Sella1B 2B in August)
            E4.to.E6=mean(E4.to.E6, na.rm=T), 
            slope_classic=mean(slope_classic, na.rm=T),
            slope_lm=mean(slope_lm, na.rm=T),
            slope_short_Loiselle=mean(slope_short_Loiselle, na.rm=T),
            slope_short_Helms=mean(slope_short_Helms, na.rm=T),
            SR_Loiselle=mean(SR_Loiselle, na.rm=T),
            SR_Helms=mean(SR_Helms, na.rm=T))

LC_OCD <- read_csv("meta_file.csv")
meta_sum_optical <- left_join(LC_OCD,abs_data, by=c("site","campaign"))

# Converting the sample names without the index numbers - wouldn't need it if data were named with 000
sample_loadings$sample_code <- str_sub(sample_loadings$sample, 1, 13)
a <- str_split_fixed(sample_loadings$sample_code, "[_]", n=4)
sample_loadings$sample_code <- paste(a[,1], a[,2],a[,3], sep = "_")
sample_loadings$sample_code <- as.character(sample_loadings$sample_code)
sample_names$sample_code <- as.character(sample_names$sample_code)

##### Binding the data frames ####
meta_component_data <- left_join(sample_loadings, sample_names,by="sample_code")
meta_component_data <- meta_component_data %>%
  filter(site!="Carrion")

meta_sum_optical <- meta_sum_optical %>%
  filter(site!="Carrion")

sample_names <- sample_names %>%
  filter(site!="Carrion")

#mean averaged data of the component loadings
data_sum <- meta_component_data %>%
  group_by(campaign,site) %>%
  summarise(
    Comp1=mean(Comp.1, na.rm=T),
    Comp2=mean(Comp.2, na.rm=T),
    Comp3=mean(Comp.3, na.rm=T),
    Comp4=mean(Comp.4, na.rm=T),
    Comp5=mean(Comp.5, na.rm=T),
    Comp6=mean(Comp.6, na.rm=T),
    Comp7=mean(Comp.7, na.rm=T),
    Comp8=mean(Comp.8, na.rm=T)
  )
data_sum <- left_join(data_sum, site_info, by="site")
data_sum$campaign <- as.character(data_sum$campaign)
data_sum <- left_join(data_sum, meta_sum_optical, by=c("site","campaign")) %>%
  mutate(
    SUVA254=DecAbsCoeff254/NPOC
  )

# CHOOSE WISELY HERE ACCORDING TO THE ANALYSIS YOU WANT TO DO
# Here I replace the under the detection limit values to the absolute detection limit
data_sum[["HMWS_C"]] <- str_replace_all(data_sum[["HMWS_C"]], c("<0,1"="NA","<0,01"="NA", "<0,02"="NA", "<0,06"="NA")) %>% as.double()
data_sum[["humic_like_substance_C"]] <- str_replace_all(data_sum[["humic_like_substance_C"]],c("<0,1"="NA","<0,01"="NA", "<0,02"="NA", "<0,06"="NA")) %>% as.double()
data_sum[["LMWS_C"]] <- str_replace_all(data_sum[["LMWS_C"]], c("<0,1"="NA","<0,01"="NA", "<0,02"="NA", "<0,06"="NA")) %>% as.double()
data_sum[["HMWS_N"]] <- str_replace_all(data_sum[["HMWS_N"]], c("<0,1"="NA","<0,01"="NA", "<0,02"="NA", "<0,06"="NA"))  %>% as.double()
data_sum[["HS"]] <- str_replace_all(data_sum[["HS"]],c("<0,1"="NA","<0,01"="NA", "<0,02"="NA", "<0,06"="NA")) %>% as.double()
data_sum[["SUVA_HS"]] <- str_replace_all(data_sum[["SUVA_HS"]], c("<0,1"="NA","<0,01"="NA", "<0,02"="NA", "<0,06"="NA"))  %>% as.double()
data_sum[["SUVA_ges"]] <- str_replace_all(data_sum[["SUVA_ges"]], c("<0,1"="NA","<0,01"="NA", "<0,02"="NA", "<0,06"="NA"))  %>% as.double()

# Or this 

#data_sum[["HMWS_C"]] <- str_replace_all(data_sum[["HMWS_C"]], c("<0,1"="0.1","<0,01"="0.01", "<0,02"="0.02", "<0,06"="0.06")) %>% as.double()
#data_sum[["humic_like_substance_C"]] <- str_replace_all(data_sum[["humic_like_substance_C"]],c("<0,1"="0.1","<0,01"="0.01", "<0,02"="0.02", "<0,06"="0.06")) %>% as.double()
#data_sum[["LMWS_C"]] <- str_replace_all(data_sum[["LMWS_C"]], c("<0,1"="0.1","<0,01"="0.01", "<0,02"="0.02", "<0,06"="0.06")) %>% as.double()
#data_sum[["HMWS_N"]] <- str_replace_all(data_sum[["HMWS_N"]], c("<0,1"="0.1","<0,01"="0.01", "<0,02"="0.02", "<0,06"="0.06"))  %>% as.double()
#data_sum[["HS"]] <- str_replace_all(data_sum[["HS"]],c("<0,1"="0.1","<0,01"="0.01", "<0,02"="0.02", "<0,06"="0.06")) %>% as.double()
#data_sum[["SUVA_HS"]] <- str_replace_all(data_sum[["SUVA_HS"]], c("<0,1"="0.1","<0,01"="0.01", "<0,02"="0.02", "<0,06"="0.06"))  %>% as.double()
#data_sum[["SUVA_ges"]] <- str_replace_all(data_sum[["SUVA_ges"]], c("<0,1"="0.1","<0,01"="0.01", "<0,02"="0.02", "<0,06"="0.06"))  %>% as.double()

#### PCA analysis ####
data_sum2 <- data_sum %>%
  select(campaign, site, Comp1, Comp2, Comp3, Comp4, Comp5, Comp6, Comp7, Comp8, BDOC, NPOC, SUVA254, 
         FIX, FIX2, HIX2, beta.alpha, slope_classic, slope_lm, slope_short_Helms, slope_short_Loiselle,
         SR_Helms, SR_Loiselle, E2.to.E3, E4.to.E6, DecAbsCoeff254, groups.x)

BDOC_normalised_components<- data_sum2 %>%
  mutate(Comp1=Comp1/NPOC,
         Comp2=Comp2/NPOC,
         Comp3=Comp3/NPOC,
         Comp4=Comp4/NPOC,
         Comp5=Comp5/NPOC,
         Comp6=Comp6/NPOC,
         Comp7=Comp7/NPOC,
         Comp8=Comp8/NPOC
  ) %>%
  select(starts_with("Comp"),NPOC, SUVA254, 
         FIX, FIX2, HIX2, beta.alpha, slope_classic, slope_lm, slope_short_Helms, slope_short_Loiselle,
         SR_Helms, SR_Loiselle, E2.to.E3, E4.to.E6, DecAbsCoeff254, campaign, site, groups.x)

NPOC_data <- (data_sum2 %>% group_by(site) %>% summarise(NPOC_mean=mean(NPOC)))

BDOC_normalised_components[is.na(BDOC_normalised_components)] <- 0
data_sum2[is.na(data_sum2)] <- 0

DOC_variance_sum <- data_sum2 %>%
  group_by(site)%>%
  summarise(
    n=n(),
    variance_BDOC=var(BDOC, na.rm = T),
    variance_NPOC=var(NPOC, na.rm=T),
    groups=groups.x
  )


theme_pca <- function (base_size = 12, base_family = "") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace% 
    theme(
      panel.grid = element_blank()
    )   
}

NPOC_var_plot <- ggplot(DOC_variance_sum) +
  theme_pca()+
  geom_boxplot(aes(x=variance_NPOC, y=groups), fill=c(col=c('#942D0A','#E65525', '#043005', "#4F9608")))+
  coord_flip()+
  ylab("Groups")+xlab("DOC Concentration Variance")+
  theme(axis.text=element_text(size=11),axis.title=element_text(size=11))+
   scale_y_discrete(limits = c("MedAlt", "MedNat", "TempAlt","TempNat"),
                    labels = c("Altered \nMediterranean", "Natural \nMediterranean", "Altered \nTemperate", "Natural \nTemperate"))

#dev.new()
#png('Seasonal NPOC variance per hydrological class.png')
#NPOC_var_plot
#dev.off()

#pca_data <- data_sum2%>% #or put BDOC_normalised_components for the normalised comp analysis
#select(starts_with("Comp"),BDOC.x, SUVA254, 
#      FIX, HIX2, beta.alpha, DecAbsCoeff420, DecAbsCoeff430, DecAbsCoeff436, slope_classic, 
#     E2.to.E3, E4.to.E6)
pca_data <- BDOC_normalised_components%>% #or put data_sum2 
  select(starts_with("Comp"), site, campaign,
         SUVA254, FIX, HIX2, beta.alpha, slope_short_Helms,SR_Loiselle,
         E2.to.E3, E4.to.E6)

wine.pca <- prcomp(pca_data[,-(9:10)], scale. = TRUE) 
summary(wine.pca)
group_types <- factor(c("TempAlt", "TempNat", "MedAlt", "MedNat"))

group_list <- vector("list")
for(i in seq_along(group_types)){
  group_list[[i]] <- site_info %>%
    select(site, Class, alteration, groups, Alteration_type, Catchment)%>%
    filter(groups==group_types[[i]]) %>%
    unique()
}

#### PCA with automated prcomp calcualted rotations ####
PCAloadings <- data.frame(Variables = rownames(wine.pca$rotation), wine.pca$rotation)
ggplot(data_sum, aes(x=wine.pca$x[,1], y=wine.pca$x[,2]))+
  geom_point(aes(color=data_sum$groups.x, fill=data_sum$groups.x), shape=21, size=4, colour="black")+
  scale_fill_manual(values=c("#E65525", "#942D0A", "#043005","#4F9608"))+
  geom_segment(data = PCAloadings, aes(x = 0, y = 0, xend = (PC1*15),
                                       yend = (PC2*15)), arrow = arrow(length = unit(1/2, "picas")),color = "black") +
  annotate("text", x = (PCAloadings$PC1*15), y = (PCAloadings$PC2*15),
           label = PCAloadings$Variables)+
  theme_classic()+
  guides(fill="legend")+
  theme(legend.position = c(-1,0))+
  labs(color="Sites", x="PC 1 (32%)", y="PC 2 (20%)", title="PCA of PARAFAC Components and Other Optical Parameters")

# component loadings comparison boxplots

par(mfrow=c(2,2),mai=c(0.3,0.3,0.3,0.3))

plot_TempAlt <- plot(wine.pca$x[,c(1,2)], type="n", main = "Altered Temperate", ylim=c(-5,6), xlim=c(-5,6), cex.main=1.5, cex.axis=1.25)
ordi_TempAlt <- ordihull(ord=wine.pca$x[,c(1,2)],groups=data_sum$site ,display="sites",draw="polygon", label=F, show.groups = group_list[[1]]$site,
                         alpha=150, col=c('#043005'))

plot_TempNat <- plot(wine.pca$x[,c(1,2)], type="n", main = "Natural Temperate", ylim=c(-5,6), xlim=c(-5,6), cex.main=1.5, cex.axis=1.25)
ordi_TempNat <- ordihull(ord=wine.pca$x[,c(1,2)],groups=data_sum$site ,display="sites",draw="polygon", label=F, show.groups = group_list[[2]]$site,
                         alpha=150, col=c("#4F9608"))
plot_MedAlt <- plot(wine.pca$x[,c(1,2)], type="n", main = "Altered Mediterranean", ylim=c(-5,6), xlim=c(-5,6), cex.main=1.5, cex.axis=1.25)
<<<<<<< HEAD
ordi_MedAlt <- ordihull(ord=wine.pca$x[,c(1,2)],groups=data_sum$site ,display="sites",draw="polygon", label=T, show.groups = group_list[[3]]$site,
=======
ordi_MedAlt <- ordihull(ord=wine.pca$x[,c(1,2)],groups=data_sum$site ,display="sites",draw="polygon", label=F, show.groups = group_list[[3]]$site,
>>>>>>> d578d27ab27749c91ae2af9bf292b19f4b132310
                        alpha=150, col=c('#942D0A'))
plot_MedNat <- plot(wine.pca$x[,c(1,2)], type="n", main = "Natural Mediterranean", ylim=c(-5,6), xlim=c(-5,6), cex.main=1.5, cex.axis=1.25)
ordi_MedNat <- ordihull(ord=wine.pca$x[,c(1,2)],groups=data_sum$site ,display="sites",draw="polygon", label=F, show.groups = group_list[[4]]$site,
                        alpha=150, col=c('#E65525'))

par(mfrow=c(1,1))

#### Anova with betadisper(vegadist) ####
betas <- vegdist(wine.pca$x[,2],method="euclidean") %>% #if you want to do only PC1 and 2, indicate as such here
  betadisper(group = BDOC_normalised_components$site) 
disp <- betas[["distances"]] %>%
  tapply(BDOC_normalised_components$site, mean) %>%
  as_tibble(rownames="site")

par(mfrow=c(1,1),mai=c(0.6,0.5,0.3,0.3), mex=1.5)
boxplot(disp[["value"]]~BDOC_normalised_components$groups.x[match(disp[[1]],BDOC_normalised_components$site)], 
        xlab=NULL, ylab = "Dispersion of PCA axis", col=c('#942D0A','#E65525', '#043005', "#4F9608"), names=c("Altered \nMediterr.", "Natural\nMediterr.", "Altered \nTemperate", "Natural \nTemperate"))
text(x=4.4, y=3.1, labels="(C)")


#png('Seasonal quality boxplot of variance across sites')
#seasonal_variance_boxplot
#dev.off()
boxplot(betas, ylim=c(0,8))
unique_data<-site_info[,c("site", "groups")]
disp<-left_join(disp, as_tibble(unique_data), by=c("site"="site"))

disp_med<-disp[str_sub(disp$groups, 1,3) == "Med",]
disp_temp<-disp[str_sub(disp$groups, 1,3) =="Tem",]

anova(lm(disp_med[["value"]]~factor(disp_med$groups)))
anova(lm(disp_temp[["value"]]~factor(disp_temp$groups)))

# Average centroid distance boxplots 
#Distance to centroid of MEditerranean sites
distance_matrix <- tibble(distances=betas[["distances"]], site=betas[["group"]])
distance_matrix <- left_join(distance_matrix, unique_data, by="site")

#Distance to centroid of MEditerranean sites
distance_matrix %>% 
  filter(groups=="MedAlt" | groups=="MedNat") %>%
  ggplot() +
  geom_boxplot(aes(x=reorder(site, distances, mean), y=distances, fill=groups))+
  theme_pca()+ coord_flip()+
  scale_fill_manual(values= c("#942D0A", "#E65525"), name = "Class", labels = c("Altered \nMEditerranean", "Natural \nMediterranean"))+labs(color="Groups", x="Sites", y="Distance to centroid", tag =c("(B)\nR^2=0.937"))+ theme(plot.tag.position=c(0.70,0.96))+
  stat_summary(aes(x=reorder(site, distances, mean), y=distances, fill=groups),fun.y=mean, geom="point", shape=20, size=2, color="red", fill="red")
<<<<<<< HEAD
=======

#Distance to centroid of Temperate sites
>>>>>>> d578d27ab27749c91ae2af9bf292b19f4b132310

#Distance to centroid of Temperate sites
distance_matrix %>% 
  filter(groups=="TempAlt" | groups=="TempNat") %>%
  ggplot() +
  geom_boxplot(aes(x=reorder(site, distances, mean), y=distances, fill=groups))+
  theme_pca()+ coord_flip()+
  labs(color="Groups", x="Sites", y="Distance to centroid")+
  scale_fill_manual(values= c("#043005","#4F9608"), name = "Class", labels = c("Altered \nTemperate", "Natural \nTemperate"))+labs(color="Groups", x="Sites", y="Distance to centroid", tag =c("(A)\nR^2=0.001"))+theme(plot.tag.position=c(0.74,0.96))+
  stat_summary(aes(x=reorder(site, distances, mean), y=distances, fill=groups),fun.y=mean, geom="point", shape=20, size=2, color="red", fill="red")+
  theme(plot.tag.position=c(0.68,0.95))


#### LC-OCD data reflection onto the PCA space ####
PCA_scores <- cbind(pca_data, as_tibble(wine.pca$x)) %>% 
  select(site, campaign, starts_with("PC"))
LC_OCD_cor <- data_sum %>% select("site", "campaign","BDOC","CDOC","HMWS_C", 
                                  "humic_like_substance_C", "LMWS_C", "HMWS_N", "HS", "SUVA_HS", "SUVA_ges") %>%
  mutate("% HMWS C"=HMWS_C*100/BDOC,
         "% HS_C"=humic_like_substance_C*100/BDOC,
         "% LMWS_C"=LMWS_C*100/BDOC,
         "HMWS (C:N)"=HMWS_C/HMWS_N, 
         "humic (C:N)"=humic_like_substance_C/HS
  )

write_delim(x=(LC_OCD_cor %>% group_by(site) %>% summarise(
  HMWS_C=mean(HMWS_C, na.rm = T), 
  humic_like_substance_C=mean(humic_like_substance_C, na.rm=T), 
  LMWS_C= mean(LMWS_C, na.rm=T), 
  HMWS_N= mean(HMWS_N, na.rm=T), 
  HS=mean(HS, na.rm=T), 
  SUVA_HS= mean(SUVA_HS, na.rm=T), 
  SUVA_ges=mean(SUVA_ges, na.rm=T))), file="LC_OCD_Sum", delim=";"
)

LC_OCD_rot <- cor(PCA_scores[c("PC1", "PC2")], LC_OCD_cor[,-c(1:9)], method = "pearson", use="complete.obs") %>% 
  t() %>% as_tibble(rownames = "Variables")
PCA_rot <- cor(PCA_scores[c("PC1", "PC2")], pca_data[,-(9:10)], method = "pearson", use="complete.obs") %>% 
  t() %>% as_tibble(rownames = "Variables")

PCAloadings <- data.frame(Variables = rownames(wine.pca$rotation), wine.pca$rotation) %>% as_tibble
new_PCAloadings <- full_join(PCAloadings, LC_OCD_rot)


theme_pca <- function (base_size = 12, base_family = "") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace% 
    theme(
      panel.grid = element_blank()
    )   
}

### PCA plots #### 

# title="PCA of PARAFAC Components and Other Optical Parameters"
plot_all <- ggplot(data_sum, aes(x=wine.pca$x[,1], y=wine.pca$x[,2]))+
  geom_point(aes(color=data_sum$groups.x, fill=data_sum$groups.x, shape=data_sum$groups.x), size=2.5, colour="black")+
  scale_shape_manual(values=c(23,22,25,24))+
  scale_fill_manual(values=c("#E65525", "#942D0A", "#043005","#4F9608"))+
  geom_segment(data = LC_OCD_rot, aes(x = 0, y = 0, xend = (PC1*8), yend = (PC2*8)), arrow = arrow(length = unit(1/2, "picas")),color = "blue") +
  annotate("text", x = (LC_OCD_rot$PC1*8), y = (LC_OCD_rot$PC2*8), label = LC_OCD_rot$Variables, size=4, color="blue")+
  geom_segment(data = PCA_rot, aes(x = 0, y = 0, xend = (PC1*8), yend = (PC2*8)), arrow = arrow(length = unit(1/2, "picas")),color = "black") +
  annotate("text", x = (PCA_rot$PC1*8), y = (PCA_rot$PC2*8),label = PCA_rot$Variables, size=4, color="black")+
  labs(color="Sites", x="PC1 (34.8%)", y="PC2 (20.1%)", tag = "a. Loadings and sites")+
  theme_pca()+ 
  theme(plot.tag.position=c(0.24,0.97))+
  scale_x_continuous(limits=c(-8,8), n.breaks=10)+
  scale_y_continuous(limits=c(-8,8), n.breaks=10)+
  guides(fill="legend")+
  theme(legend.position = c(-1,0))+
  geom_vline(xintercept = 0, lty=2) + geom_hline(yintercept = 0, lty=2)


plot_points <- ggplot(data_sum, aes(x=wine.pca$x[,1], y=wine.pca$x[,2]))+
  geom_point(aes(color=data_sum$groups.x, fill=data_sum$groups.x, shape=data_sum$groups.x), size=2.5, colour="black")+
  scale_shape_manual(values=c(23,22,25,24))+
  scale_fill_manual(values=c("#E65525", "#942D0A", "#043005","#4F9608"))+
  labs(color="Sites", x="PC1 (34.8%)", y="PC2 (20.1%)", tag = "b. Sites")+
  theme_pca()+ 
  theme(plot.tag.position=c(0.15,0.97))+
  scale_x_continuous(limits=c(-8,8), n.breaks=10)+
  scale_y_continuous(limits=c(-8,8), n.breaks=10)+
  guides(fill="legend")+
  theme(legend.position = c(-1,0))+
  geom_vline(xintercept = 0, lty=2) + geom_hline(yintercept = 0, lty=2)

plot_optical <- ggplot(data_sum, aes(x=wine.pca$x[,1], y=wine.pca$x[,2]))+
  geom_segment(data = PCA_rot, aes(x = 0, y = 0, xend = (PC1), yend = (PC2)), arrow = arrow(length = unit(1/3, "picas")),color = "black") +
  annotate("text", x = (PCA_rot$PC1), y = (PCA_rot$PC2+0.05),label = PCA_rot$Variables, size=3.5, color="black")+
  labs(color="Sites", x="PC1 (34.8%)", y="PC2 (20.1%)", tag = "c. Optical parameters")+
  theme_pca()+ 
  theme(plot.tag.position=c(0.24,0.97))+
  scale_x_continuous(limits=c(-1,1), n.breaks=10)+
  scale_y_continuous(limits=c(-1,1), n.breaks=10)+
  guides(fill="legend")+
  theme(legend.position = c(-1,0))+
  geom_vline(xintercept = 0, lty=2) + geom_hline(yintercept = 0, lty=2)

plot_LCOCD <- ggplot(data_sum, aes(x=wine.pca$x[,1], y=wine.pca$x[,2]))+
  geom_segment(data = LC_OCD_rot, aes(x = 0, y = 0, xend = (PC1), yend = (PC2)), arrow = arrow(length = unit(1/2, "picas")),color = "blue") +
  annotate("text", x = (LC_OCD_rot$PC1), y = (LC_OCD_rot$PC2), label = LC_OCD_rot$Variables, size=4, color="blue")+
  labs(color="Sites", x="PC1 (34.8%)", y="PC2 (20.1%)", tag = "d. LC-OCD parameters")+
  theme_pca()+  
  theme(plot.tag.position=c(0.28,0.97))+ 
  scale_x_continuous(limits=c(-1,1), n.breaks=10)+
  scale_y_continuous(limits=c(-1,1), n.breaks=10)+
  guides(fill="legend")+
  theme(legend.position = c(-1,0))+
  geom_vline(xintercept = 0, lty=2) + geom_hline(yintercept = 0, lty=2)


#### end ####
# Treating the Lc-OCD data as new data and do a PCA on that so that later on you can do a procrustes anaylsis
# For this you can use the dataset without reduced information Choose the bottom option above in line 100

LC_OCD_cor[is.na(LC_OCD_cor)] <- 0
wine.pca3 <- prcomp(LC_OCD_cor[,-c(1:2)], scale=T)
pro <- procrustes(X = wine.pca, Y = wine.pca3, symmetric = T)
protest(X = wine.pca3, Y = wine.pca, scores = "sites", permutations = 999)

pr <- plot(pro, kind=1, type="text")
summary(pro)
biplot(wine.pca3)

#### Flow/Parameter Analysis ####
flow_data <- read_csv("Flow_HYDRA_2017_2018.csv")
flow_data$Date <- (mdy(flow_data$Date))

site_groups <- meta_sum_optical %>%
  select(site, groups) %>%
  unique()
group_types <- factor(c("TempAlt", "TempNat", "MedAlt", "MedNat"))
group_list <- vector("list")
for(i in seq_along(group_types)){
  group_list[[i]] <- site_groups %>%
    filter(groups==group_types[[i]]) %>%
    unique()
}

flow_data <- flow_data %>%
  filter(Date > dmy("01-10-2017"), Date < dmy("01-09-2018"))


# plot with normalized component PC's
data_sum4 <- cbind(wine.pca$x[,1:2], as.data.frame(pca_data[,c("site","campaign")]))
for(i in 1:nrow(data_sum4)){
  if (data_sum4[["campaign"]][[i]] == "apr") {
    data_sum4[["exact_date"]][[i]] <- as.character.Date(dmy("11/04/2018"))
  } else if (data_sum4[["campaign"]][[i]] == "oct") {
    data_sum4[["exact_date"]][[i]] <- as.character.Date(dmy("11/10/2017"))
  } else if (data_sum4[["campaign"]][[i]] == "dec"){
    data_sum4[["exact_date"]][[i]] <- as.character.Date(dmy("13/12/2017"))
  } else if (data_sum4[["campaign"]][[i]] == "feb"){
    data_sum4[["exact_date"]][[i]] <-  as.character.Date(dmy("14/02/2018"))
  } else if (data_sum4[["campaign"]][[i]] == "may"){
    data_sum4[["exact_date"]][[i]] <-  as.character.Date(dmy("13/06/2018"))
  } else if (data_sum4[["campaign"]][[i]] == "aug") {
    data_sum4[["exact_date"]][[i]] <- as.character.Date(dmy("14/08/2018"))
  } else {data_sum4[["exact_date"]][[i]] <- (0)}
}

data_sum4[["exact_date"]] <- as.Date(unlist(data_sum4[["exact_date"]]))
data_sum4 <- data_sum4 %>% as_tibble()
data_sum4 <- left_join(data_sum4, site_info, by=c("site"="site"))

data_sum5 <- data_sum4 %>%
  group_by(site) %>%
  mutate(meanPC1=mean(PC1),
         meanPC2=mean(PC2), 
         var_PC1=var(PC1),
         var_PC2=var(PC2)) %>%
  select(site, groups, meanPC1,meanPC2, starts_with("var")) %>% unique()

data_sum5 %>%
  ggplot() +
  geom_boxplot(aes(x=meanPC2, y=groups, fill=groups)) + coord_flip() +
  scale_fill_manual(values=c("#942D0A", "#E65525" , "#043005","#4F9608")) +
  ggtitle("PC2 mean")

#### Hydrological Indices Analysis ####
magnitude_ind <- read_csv("Hydra_Hydro_Ind_subs_group1.csv") %>% 
  mutate()
magnitude_duration_extreme_ind <- read_csv("Hydra_Hydro_Ind_subs_group2.csv")
timing_extreme_ind <- read_csv("Hydra_Hydro_Ind_subs_group3.csv")
freq_duration_pulses_ind <- read_csv("Hydra_Hydro_Ind_subs_group4.csv")
rate_freq_ind <- read_csv("Hydra_Hydro_Ind_subs_group5.csv")

meta_indices <-  left_join(magnitude_ind, magnitude_duration_extreme_ind, timing_extreme_ind, by="Stream") %>%
  left_join(timing_extreme_ind, by="Stream") %>%
  left_join(freq_duration_pulses_ind, by="Stream") %>%
  left_join(rate_freq_ind, by="Stream") %>%
  filter(Stream!="Carrion") #tkaing out carrion from the analysis

pca_data2 <- meta_indices %>%
  left_join(data_sum5[,1:2], by=c("Stream"="site"))

wine.pca2 <- prcomp(pca_data2[,-c(1,87)], scale. = TRUE) 
summary(wine.pca2)
screeplot(wine.pca2)
# The kaiser rule (with the eigenvalues at sdev) indicates the first 9 pc's are important
# The screeplot shows that the first 3,4 are the most important
# So first try modeling with only the 9 PC's

PCAloadings <- data.frame(Variables = rownames(wine.pca2$rotation), wine.pca2$rotation)
ggplot(pca_data2, aes(x=wine.pca2$x[,1], y=wine.pca2$x[,2]))+
  geom_point(aes(fill=data_sum5$groups), shape=21, size=2.5, colour="black")+
  scale_fill_manual(values=c("#E65525", "#942D0A", "#043005","#4F9608"))+
  geom_segment(data = PCAloadings, aes(x = 0, y = 0, xend = (PC3*2),
                                       yend = (PC4*2)), arrow = arrow(length = unit(1/10, "picas")),color = "black") +
  annotate("text", x = (PCAloadings$PC3*10), y = (PCAloadings$PC4*10),
           label = PCAloadings$Variables)+
  theme_classic()+
  guides(fill="legend")+
  theme(legend.position = c(-1,0))+
  labs(color="Sites", x="PC 3", y="PC 4", title="Principal Component Analysis of Flow Indices")


flow_indice_scores <- cbind("site"=pca_data2[["Stream"]], wine.pca2$x) %>% as_tibble

ordered_flow_scores <- flow_indice_scores %>%arrange(site) %>%
  left_join(data_sum5[c("var_PC1", "site")], by="site") %>%
  select(-site) %>% #we exclude PC19 bcs it is perfectly alined with something
  mutate_if(is.character, as.numeric)

#### Akaike ####
mlr <- lm(var_PC1~., data=ordered_flow_scores)
car::vif(mlr)
alias(mlr) # to check multicolinerarity, PC19 seems redundant
summary(mlr)
plot(ordered_flow_scores)  

simple.mod1<-lm(var_PC1~PC1,data=ordered_flow_scores)
simple.mod2 <- lm(var_PC1~PC4, data=ordered_flow_scores)
summary(simple.mod1)
summary(simple.mod2)

#### Step Akaike ####
step(simple.mod1,scope=~.+PC1+PC2+PC3+PC4,direction="both")
# same solutions! compare with selection by F!

# combined forward and backward selection based on AIC
step(simple.mod1,direction="both")
extractAIC(simple.mod1)

indice_factor_loadings <- wine.pca2[["rotation"]][,c("PC1", "PC2")]

#Grouped models for each hydrological class
simple.mod1<-lm(var_PC1~PC1,data=ordered_flow_scores)

ordered_flow_scores_Mediterranean <-  flow_indice_scores %>%arrange(site) %>%
  left_join(data_sum5[c("var_PC1", "site", "groups")], by="site") %>%
  filter(groups=="MedNat"|groups=="MedAlt") %>%
  select(-site, -groups) %>% #we exclude PC19 bcs it is perfectly alined with something
  mutate_if(is.character, as.numeric)

simple_mod_med<-lm(var_PC1~PC3,data=ordered_flow_scores_Mediterranean)
summary(simple_mod_med)
step(simple_mod_med, scope=~.+PC1+PC2+PC3+PC4, direction="both")
extractAIC(simple_mod_med)

ordered_flow_scores_Temperate <-  flow_indice_scores %>%arrange(site) %>%
  left_join(data_sum5[c("var_PC1", "site", "groups")], by="site") %>%
  filter(groups=="TempNat"|groups=="TempAlt") %>%
  select(-site, -groups) %>% #we exclude PC19 bcs it is perfectly alined with something
  mutate_if(is.character, as.numeric)

simple_mod_temp<-lm(var_PC1~PC3+PC1,data=ordered_flow_scores_Temperate)
summary(simple_mod_temp)
step(simple_mod_temp, scope=~.+PC1+PC2+PC3+PC4, direction="both")
extractAIC(simple_mod_temp)
plot(simple_mod_med)


#### PLSR trial ####
library(pls)
library(plsVarSel)
library(mdatools)

meta_indices %<>% arrange(Stream) %>%
  column_to_rownames("Stream")

data_sum5 %<>% arrange(site) %>%
  column_to_rownames("site")

# with mda tools
selected_variable <- data_sum5[5]
model <-  pls(meta_indices, selected_variable, ncomp=15,  cv=list("rand", 4, 4), ncomp.selcrit="min", scale = TRUE, info = "meta indice-")
show(model$ncomp.selected)
plotRMSE(model)
summary(model)
summary(model$res$cal)

model <-  selectCompNum(model, 7)
plot(model)
par(mfrow=c(2,2))
plotVariance(model, type="h", show.labels=T) # This helps pickng the right number of components in bulding a model. A big jump indicates overfitting
plotQDoF(model) # if this isn's in accord with the variance plots,then there is a problem.Probably overfitting. Check https://mdatools.com/docs/pca-model-complexity.html for more
plotT2DoF(model)
plotDistDoF(model)

par(mfrow = c(2, 2))
plotRegcoeffs(model)
plotRegcoeffs(model, ncomp = 2)
plot(model$coeffs, ncomp = 1, type = "b", show.labels = TRUE)
plot(model$coeffs, type = "b", show.labels = TRUE)
par(mfrow = c(1, 1))

plotVIPScores(model, show.labels=T)
plotVIPScores(model, ncomp = 1, type = "h", show.labels = TRUE)
vip <-  vipscores(model, ncomp = 1)
m3 <-  pls(meta_indices, selected_variable, 1, scale = T, cv = 1, exclcols = (vip < 0.5))

plotSelectivityRatio(model,  type = "h", show.labels = TRUE)
plotSelectivityRatio(model, ncomp = 2, type = "h", show.labels = TRUE)

c = categorize(model, newm$res$cal)
print(c) # No outliers so we keep all the sites

#This is chosen. Ncomp=1. Now we select the variables with VIP scores higher than 1 as the important ones to be plotted in the 
#Either select by looking at the VIP scores, or you can do a Regression coefficient analysis using jack.knifing technique (Mehmood et al. 2020) This is harsher in selecting but can be useful
model <-  selectCompNum(model, 1)
plotRegcoeffs(model, type = "h", show.ci = TRUE, show.labels = TRUE)
VIPscores <- plotVIPScores(model, ncomp =1, type = "h", show.labels = TRUE)
summary(model$coeffs, ncomp = 1)

exclcols <-  model$coeffs$p.values[, 1, 1] > 0.15
show(exclcols)

newm <-  pls(meta_indices, selected_variable, 1, scale = TRUE, cv = 1, exclcols = exclcols)
summary(newm)
plot(newm$coeffs, ncomp = 1, type = "b", show.labels = TRUE)
plotVIPScores(newm, ncomp = 1, type = "h", show.labels = TRUE)
show(getRegcoeffs(newm))

plotXResiduals(newm, ncomp = 1)
plotXResiduals(newm$res$cal)

#### End of RC perdictor selecting method ####

# WE move with manual selection of predictors by VIP scores. 

plotVIPScores(model, ncomp = 1, type = "h", show.labels = TRUE)
vip <-  vipscores(model, ncomp = 1)
model_selected <- pls(meta_indices, selected_variable, 1, scale = T, cv = 1, exclcols = (vip < 1))
plotVIPScores(model_selected, ncomp = 1, type = "h", show.labels = TRUE)
plot(model_selected$coeffs, ncomp = 1, type = "b", show.labels = TRUE)

#reg_coeffs_PC1_var  <- model_selected[["coeffs"]][["values"]] %>% as_tibble(rownames = NA) 
#reg_coeffs_PC1_var <- reg_coeffs_PC1_var# %>% filter(reg_coeffs_PC1_var [,1]!=0)
#vip_PC1_var <-  vipscores(model_selected, ncomp = 1)

# BE  CAREFUL HERE: CORRECT THIS FOR REPLICABILITY, OTHERWISE ONLY DO IT WHEN DATA_SUM5[[5]]
reg_coeffs_PC2_var <-  model_selected[["coeffs"]][["values"]] %>% as_tibble(rownames = NA) 
reg_coeffs_PC2_var <- reg_coeffs_PC2_var # %>% filter(reg_coeffs[,1]!=0)
vip_PC2_var <-  vipscores(model_selected, ncomp = 1)

reg_coeffs <- left_join(reg_coeffs_PC1_var%>% rownames_to_column("indice"), reg_coeffs_PC2_var%>% rownames_to_column("indice"), by="indice")
colnames(reg_coeffs) <- c("indice", "var_PC1", "var_PC2")
reg_coeffs_no_zero <- reg_coeffs %>%
  filter(var_PC1!="0" | var_PC2!="0") %>% t()


colnames(reg_coeffs_no_zero) <- reg_coeffs_no_zero[1,]
reg_coeffs2 <- reg_coeffs_no_zero[2:3,]
rownames(reg_coeffs2) <- c("var_PC1", "var_PC2")

reg_coeffs3 <-  reg_coeffs2 %>% as.array()
reg_coeffs4 <-  reg_coeffs2 %>% as_tibble(rownames=NA, )

for(i in 1:length(reg_coeffs4)) {
reg_coeffs4[i] <- as.numeric(unlist(reg_coeffs4[i]))
}
rownames(reg_coeffs4) <- c("var_PC1", "var_PC2")

#### Heatmap for indices and DOM variance ####
library(reshape2)
library(RColorBrewer)

head(meta_indices)
data_sum6 <- data_sum5

meta_indices2 <- meta_indices %>% arrange(rownames(meta_indices)) %>% 
  rownames_to_column() %>%
  left_join(rownames_to_column(data_sum5), by="rowname") %>%
  select(-var_PC1, -var_PC2, -meanPC1, -meanPC2) %>%
  column_to_rownames("rowname")

cormat <- round(cor(ordered_flow_scores[,-20], meta_indices2[,-86]),3)
cormat2 <- round(cor(data_sum6[,-1], meta_indices2[,-86]),3)
cormat_med_alt <- round(cor((data_sum6%>%filter(groups=="MedAlt")%>%select(-groups)), meta_indices2%>%filter(groups=="MedAlt")%>%select(-groups)),3)
cormat_temp_alt <- round(cor((data_sum6%>%filter(groups=="TempAlt")%>%select(-groups)), meta_indices2%>%filter(groups=="TempAlt")%>%select(-groups)),3)
cormat_med_nat <- round(cor((data_sum6%>%filter(groups=="MedNat")%>%select(-groups)), meta_indices2%>%filter(groups=="MedNat")%>%select(-groups)),3)
cormat_temp_nat <- round(cor((data_sum6%>%filter(groups=="TempNat")%>%select(-groups)), meta_indices2%>%filter(groups=="TempNat")%>%select(-groups)),3)


melted_indices <- reshape2::melt(cormat)
melted_indices2 <- reshape2::melt(cormat2)
melted_indices3 <- reshape2::melt(reg_coeffs3)
melted_indices3$value <- as.numeric(melted_indices3$value)


melted_indices_medalt <- reshape2::melt(cormat_med_alt)
melted_indices_tempalt <- reshape2::melt(cormat_temp_alt)
melted_indices_mednat <- reshape2::melt(cormat_med_nat)
melted_indices_tempnat <- reshape2::melt(cormat_temp_nat)

# Create a ggheatmap
ggplot(melted_indices2, aes(x=Var2, y= Var1  , fill = value))+ 
  geom_raster()+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-0.8,0.8), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()

# Print the heatmap
print(ggheatmap)
library(RColorBrewer)
library(pheatmap)

heatmaplabels <- c("I2","Icv", "Ica", "Ikur", "Mean Monthly Flow_January", "Mean Monthly Flow_February", "Mean Monthly Flow_March", "Mean Monthly Flow_April", "Mean Monthly Flow_May", "Mean Monthly Flow_June", "Mean Monthly Flow_July", "Mean Monthly Flow_August", "Mean Monthly Flow_September", "Mean Monthly Flow_October", "Mean Monthly Flow_November", "Mean Monthly Flow_December", 
                   "Standard Deviation January", "Standard Deviation February", "Standard Deviation March", "Standard Deviation April", "Standard Deviation May", "Standard Deviation June", "Standard Deviation July", "Standard Deviation August", "Standard Deviation September", "Standard Deviation October", "Standard Deviation November", "Standard Deviation December",
                   "Magnitude of flows exceeded 5% of the time", "Magnitude of flows exceeded 25% of the time", "Magnitude of flows exceeded 75% of the time", "Magnitude of flows exceeded 95% of the time", 
                   "Magnitude of Minimum annual flow of 1-day duration", "Magnitude of Maximum annual flow of 1-day duration", "Magnitude of Minimum annual flow of 3-day duration", "Magnitude of Maximum annual flow of 3-day duration","Magnitude of Minimum annual flow of 7-day duration", "Magnitude of Maximum annual flow of 7-day duration","Magnitude of Minimum annual flow of 30-day duration", "Magnitude of Maximum annual flow of 30-day duration","Magnitude of Minimum annual flow of 90-day duration", "Magnitude of Maximum annual flow of 90-day duration",
                   "sd1LF", "sd1HF","sd3LF", "sd3HF","sd7LF", "sd7HF","sd30LF", "sd30HF","sd90LF", "sd90HF", "Number of 0-flow days (ZFD)", "7-day min flow/mean annual daily flows (BFI)", "sdZFD", "sdBFI",
                   "Number of high flow events per year (upper threshold 1-time median flow overall years)","Number of high flow events per year (upper threshold 3-time median flow overall years)","Number of high flow events per year (upper threshold 7-time median flow overall years)","sdFRE1", "sdFRE3", "sdFRE7", 
                   "Number of low pulses per year","Duration of low pulses per year","Number of high pulses per year","Duration of high pulses per year","sdnPLow", "sddPLow", "sdnHigh", "sddPHigh")

coul <- rev(colorRampPalette(brewer.pal(7, "RdBu"))(9))

pheatmap::pheatmap((cormat), cluster_cols=F, cluster_rows=F, cellheight = 12, cellwidth = 12, display_numbers = F, number_format = "%.1f", fontsize_number=5,number_color = "black", gaps_col =c(4, 16, 28, 42, 52, 56, 62), border_color = "grey",
                   angle_col = 45, color = coul)
pheatmap::pheatmap((cormat2), cluster_cols=F, cluster_rows=F, cellheight = 12, cellwidth = 12, display_numbers = F, number_format = "%.1f", fontsize_number=5,number_color = "black", gaps_col =c(4, 16, 28, 42, 52, 56, 62), labels_col=heatmaplabels,border_color = "grey",
                   angle_col = 45, color = coul)

pheatmap::pheatmap((reg_coeffs4), cluster_cols=T, cluster_rows=F, cellheight = 15, cellwidth = 15, display_numbers = F, number_format = "%.3f", fontsize_number=5,number_color = "black",  border_color = "grey",
                   angle_col = 45, color = coul)
order <- c("nPos","BFI","sdJMIn", "sdReversals", "sdJMax",
            "1LF","3LF", "7LF", "30LF", "90LF",  "sd1LF", "sd90LF", "M6","M7", "M8",  "M9",  "sdM9", 
           "sdM6", "sdM7", "sdM8","sd3LF","sd7LF","sd30LF",  "sdBFI", "X75",
           "90HF","dPHigh","sdnNeg", "sdnPos","sddPHigh","Neg",
          "sdNeg",
          "Pos","sdPos", "FRE7","FRE3","nPHigh","sdM10", "sdM1","sdM11", "1HF","3HF" , "JMin","FRE1",
          "nNeg", "30HF",  "M1","M2","M3","M11", "M12", "sdM12",  "X5", "sd90HF", "l2", "lcv"
          )

reg_coeffs5 <- reg_coeffs4[,order]



ordered_heatmap <- pheatmap::pheatmap((reg_coeffs5), cluster_cols=F, cluster_rows=F, cellheight = 15, cellwidth =15, display_numbers = F, number_format = "%.3f", fontsize_number=5,number_color = "black",  border_color = "grey",
                  angle_col = 45, color = coul, main= "Explanatory Indices", fontsize = 12, labels_row=c("PC1 variation", "PC2 variation"))

#### Separate PC variation barplots ####
# Anova with betadisper(vegadist) 
betas <- vegdist(wine.pca$x[,1],method="euclidean") %>% #if you want to do only PC1 and 2, indicate as such here
  betadisper(group = BDOC_normalised_components$site) 
disp <- betas[["distances"]] %>%
  tapply(BDOC_normalised_components$site, mean) %>%
  as_tibble(rownames="site")
boxplot(disp[["value"]]~BDOC_normalised_components$groups.x[match(disp[[1]],BDOC_normalised_components$site)], 
        xlab=NULL, ylab = "Dispersion of PCA axis", col=c('#942D0A','#E65525', '#043005', "#4F9608"), main="Variance in DOM (quality) over the seasons")

unique_data<-site_info[,c("site", "groups")]
disp<-left_join(disp, as_tibble(unique_data), by=c("site"="site"))

disp_med<-disp[str_sub(disp$groups, 1,3) == "Med",]
disp_temp<-disp[str_sub(disp$groups, 1,3) =="Tem",]

anova(lm(disp_med[["value"]]~factor(disp_med$groups)))
anova(lm(disp_temp[["value"]]~factor(disp_temp$groups)))

distance_matrix <- tibble(distances=betas[["distances"]], site=betas[["group"]])
distance_matrix <- left_join(distance_matrix, unique_data, by="site")
distance_matrix %>% 
  filter(groups=="TempAlt" | groups=="TempNat") %>%
  ggplot() + 
  geom_boxplot(aes(x=reorder(site, distances, mean), y=distances, fill=groups))+
  theme_classic()+ coord_flip()+
  labs(color="Groups", x="Sites", y="Distance to centroid", title="Average Distance to Centroid")+
  #scale_fill_manual(values= c("#942D0A", "#E65525"))+
  scale_fill_manual(values= c("#043005","#4F9608"))+
  stat_summary(aes(x=reorder(site, distances, mean), y=distances, fill=groups),fun.y=mean, geom="point", shape=20, size=2, color="red", fill="red")+
  scale_x_discrete(expand = c(0,2))+
  geom_segment(aes(x = 0,  xend = 0,y=3, yend=6), arrow = arrow(length = unit(0.5, "cm")))
#######
