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
meta_sum_optical <- read_csv("meta_file.csv")
meta_sum_optical$campaign<- str_replace_all(meta_sum_optical$campaign, c("April"="apr", "February"="feb", "May"="may","August"="aug","October"="oct","December"="dec"))


# Converting the sample names without the index numbers - wouldn't need it if data were named with 000
sample_loadings$sample_code <- str_sub(sample_loadings$sample, 1, 13)
a <- str_split_fixed(sample_loadings$sample_code, "[_]", n=4)
sample_loadings$sample_code <- paste(a[,1], a[,2],a[,3], sep = "_")
sample_loadings$sample_code <- as.character(sample_loadings$sample_code)
sample_names$sample_code <- as.character(sample_names$sample_code)

##### Binding the data frames ####
meta_data <- left_join(sample_loadings, sample_names,by="sample_code")

#mean averged data of the component loadings
data_sum <- meta_data %>%
  group_by(campaign,site) %>%
  summarise(
    Comp1=mean(Comp.1, na.rm=T),
    Comp2=mean(Comp.2, na.rm=T),
    Comp3=mean(Comp.3, na.rm=T),
    Comp4=mean(Comp.4, na.rm=T),
    Comp5=mean(Comp.5, na.rm=T),
    Comp6=mean(Comp.6, na.rm=T),
    Comp7=mean(Comp.7, na.rm=T),
    Comp8=mean(Comp.8, na.rm=T), 
    BDOC=mean(BDOC, na.rm=T),
    NPOC=mean(NPOC, na.rm=T)
  )
data_sum <- left_join(data_sum, site_info, by="site")
data_sum$campaign <- as.character(data_sum$campaign)
data_sum <- left_join(data_sum, meta_sum_optical, by=c("site","campaign", "Catchment", "groups"))
data_sum$Alteration_type <- as_factor(data_sum$Alteration_type)

#### PCA analysis ####
data_sum2 <- data_sum %>%
  select(campaign, site, Comp1, Comp2, Comp3, Comp4, Comp5, Comp6, Comp7, Comp8, BDOC.x, SUVA254, 
         FIX, FIX2, HIX2, beta.alpha, DecAbsCoeff420, DecAbsCoeff430, DecAbsCoeff436, slope_classic, 
         E2.to.E3, E4.to.E6, DecAbsCoeff254)
BDOC_normalised_components<- data_sum2 %>%
  mutate(Comp1=Comp1/BDOC.x,
         Comp2=Comp2/BDOC.x,
         Comp3=Comp3/BDOC.x,
         Comp4=Comp4/BDOC.x,
         Comp5=Comp5/BDOC.x,
         Comp6=Comp6/BDOC.x,
         Comp7=Comp7/BDOC.x,
         Comp8=Comp8/BDOC.x
  ) %>%
  select(starts_with("Comp"),BDOC.x, SUVA254, 
         FIX, FIX2, HIX2, beta.alpha, DecAbsCoeff420, DecAbsCoeff430, DecAbsCoeff436, slope_classic, 
         E2.to.E3, E4.to.E6, DecAbsCoeff254, campaign, site)

BDOC_normalised_components[is.na(BDOC_normalised_components)] <- 0
data_sum2[is.na(data_sum2)] <- 0

#pca_data <- data_sum2%>% #or put BDOC_normalised_components for the normalised comp analysis
#select(starts_with("Comp"),BDOC.x, SUVA254, 
#      FIX, HIX2, beta.alpha, DecAbsCoeff420, DecAbsCoeff430, DecAbsCoeff436, slope_classic, 
#     E2.to.E3, E4.to.E6)
pca_data <- BDOC_normalised_components%>% #or put ddata_sum2 
  select(starts_with("Comp"), site, campaign ,SUVA254, FIX, HIX2, beta.alpha, DecAbsCoeff420, DecAbsCoeff430, DecAbsCoeff436, slope_classic, 
         E2.to.E3, E4.to.E6)

wine.pca <- prcomp(pca_data[,-(9:10)], scale. = TRUE) 

alteration_types <- factor(c("0", "1", "2", "3", "4", "5"))

alteration_list <- vector("list")
for(i in seq_along(alteration_types)){
  alteration_list[[i]] <- site_info %>%
    select(site, Class, alteration, groups, Alteration_type, Catchment)%>%
    filter(Alteration_type==alteration_types[[i]]) %>%
    unique()
}

#### Plots ####
PCAloadings <- data.frame(Variables = rownames(wine.pca$rotation), wine.pca$rotation)
ggplot(data_sum, aes(x=wine.pca$x[,1], y=wine.pca$x[,2]))+
  geom_point(aes(color=data_sum$groups, fill=data_sum$Alteration_type), shape=21, size=2.5, colour="black")+
  scale_fill_manual(values=c("#E65525", "#942D0A", "#043005","#4F9608", "blue", "blue"))+
  geom_segment(data = PCAloadings, aes(x = 0, y = 0, xend = (PC1*15),
                                       yend = (PC2*15)), arrow = arrow(length = unit(1/2, "picas")),color = "black") +
  annotate("text", x = (PCAloadings$PC1*15), y = (PCAloadings$PC2*15),
           label = PCAloadings$Variables)+
  theme_classic()+
  guides(fill="legend")+
  theme(legend.position = c(-1,0))+
  labs(color="Sites", x="PC 1", y="PC 2", title="Principal Component Analysis of PARAFAC Components")


# component loadings comparison boxplots

par(mfrow=c(2,3))

plot_TempNat <- plot(wine.pca$x[,1:2], type="n", main = "TempNat", ylim=c(-6,6), xlim=c(-6,6), cex.main=2, cex.axis=1.25)
ordi_TempNat <- ordihull(ord=wine.pca$x[,1:2],groups=data_sum$site ,display="sites",draw="polygon", label=F, show.groups = alteration_list[[1]]$site,
                         alpha=150, col=c("blue"))

plot_TempAlt_1 <- plot(wine.pca$x[,1:2], type="n", main = "TempAlt Alteration 1", ylim=c(-6,6), xlim=c(-6,6), cex.main=2, cex.axis=1.25)
ordi_TempAlt_1 <- ordihull(ord=wine.pca$x[,1:2],groups=data_sum$site ,display="sites",draw="polygon", label=F, show.groups = alteration_list[[2]]$site,
                         alpha=150, col=c("green"))

plot_TempAlt_2 <- plot(wine.pca$x[,1:2], type="n", main = "TempAlt Alteration 2", ylim=c(-6,6), xlim=c(-6,6), cex.main=2, cex.axis=1.25)
ordi_TempAlt_2 <- ordihull(ord=wine.pca$x[,1:2],groups=data_sum$site ,display="sites",draw="polygon", label=F, show.groups = alteration_list[[3]]$site,
                        alpha=150, col=c('green'))

plot_MedNat <- plot(wine.pca$x[,1:2], type="n", main = "MedNat", ylim=c(-6,6), xlim=c(-6,6), cex.main=2, cex.axis=1.25)
ordi_MedNat <- ordihull(ord=wine.pca$x[,1:2],groups=data_sum$site ,display="sites",draw="polygon", label=F, show.groups = alteration_list[[4]]$site,
                        alpha=150, col=c("blue"))

plot_MedAlt_4 <- plot(wine.pca$x[,1:2], type="n", main = "Med Alt & Esla Alteration 4", ylim=c(-6,6), xlim=c(-6,6), cex.main=2, cex.axis=1.25)
ordi_MedAlt_4 <- ordihull(ord=wine.pca$x[,1:2],groups=data_sum$site ,display="sites",draw="polygon", label=F, show.groups = alteration_list[[5]]$site,
                        alpha=150, col=c('orange'))

plot_MedAlt_5 <- plot(wine.pca$x[,1:2], type="n", main = "MedAlt Alteration 6", ylim=c(-6,6), xlim=c(-6,6), cex.main=2, cex.axis=1.25)
ordi_MedAlt_5 <- ordihull(ord=wine.pca$x[,1:2],groups=data_sum$site ,display="sites",draw="polygon", label=F, show.groups = alteration_list[[6]]$site,
                        alpha=150, col=c('orange'))

par(mfrow=c(1,1))

#### Saving the centroid values of each site per group. Can also do it with a nested loop within loop but looks like it's gonna take time ####
#TempAlt
centroid_values <- list()
for (i in 1:nrow(alteration_list[[1]])) {
  centroid_values [[i]] <- alteration_list[[1]][i,] %>%
    mutate(PC2=centroid(ordi_TempNat[[i]])[1],
           PC1=centroid(ordi_TempNat[[i]])[2])
}
centroid_TempNat <- do.call(rbind.data.frame, centroid_values)

#TempNat
centroid_values <- list()
for (i in 1:(nrow(alteration_list[[2]]))) {
  centroid_values [[i]] <- alteration_list[[2]][i,] %>%
    mutate(PC2=centroid(ordi_TempAlt_1[[i]])[1],
           PC1=centroid(ordi_TempAlt_1[[i]])[2])
}
centroid_TempAlt_1 <- do.call(rbind.data.frame, centroid_values)

#MedAlt
centroid_values <- list()
for (i in 1:nrow(alteration_list[[3]])) {
  centroid_values [[i]] <- alteration_list[[3]][i,] %>%
    mutate(PC2=centroid(ordi_TempAlt_2[[i]])[1],
           PC1=centroid(ordi_TempAlt_2[[i]])[2])
}
centroid_TempAlt_2 <- do.call(rbind.data.frame, centroid_values)

#MedNat
centroid_values <- list()
for (i in 1:(nrow(alteration_list[[4]]))) {
  centroid_values [[i]] <- alteration_list[[4]][i,] %>%
    mutate(PC2=centroid(ordi_MedNat[[i]])[1], # Here leaving the last point of the polygon doesn't change anything. I guess it's already taken into account
           PC1=centroid(ordi_MedNat[[i]])[2])
}

centroid_MedNat <- do.call(rbind.data.frame, centroid_values)

#MedNat
centroid_values <- list()
for (i in 1:(nrow(alteration_list[[5]]))) {
  centroid_values [[i]] <- alteration_list[[5]][i,] %>%
    mutate(PC2=centroid(ordi_MedAlt_4[[i]])[1], # Here leaving the last point of the polygon doesn't change anything. I guess it's already taken into account
           PC1=centroid(ordi_MedAlt_4[[i]])[2])
}

centroid_MedAlt_4 <- do.call(rbind.data.frame, centroid_values)

#MedNat
centroid_values <- list()
for (i in 1:(nrow(alteration_list[[6]]))) {
  centroid_values [[i]] <- alteration_list[[6]][i,] %>%
    mutate(PC2=centroid(ordi_MedAlt_5[[i]])[1], # Here leaving the last point of the polygon doesn't change anything. I guess it's already taken into account
           PC1=centroid(ordi_MedAlt_5[[i]])[2])
}

centroid_MedAlt_5 <- do.call(rbind.data.frame, centroid_values)


centoid_data <- rbind(centroid_TempNat,centroid_TempAlt_1,centroid_TempAlt_2, centroid_MedNat, centroid_MedAlt_4, centroid_MedAlt_5)


#### Distances of each sampling time to their respective centroid ####
# change the sum to mean hear if you want to get the mean of each point to the centroid instead of sum
ordi_data <- append(ordi_TempNat, values=c(ordi_TempAlt_1,ordi_TempAlt_2, ordi_MedNat, ordi_MedAlt_4, ordi_MedAlt_5)) #combines the ordi point lists as a 1 list to be used in the loop later
nsite <- ncol(site_info)+1
for(i in seq_along(ordi_data)){
  distances <-((ordi_data[[i]][,1]-as.numeric(centoid_data[centoid_data$site==names(ordi_data[i]),"PC1"]))^2 + (ordi_data[[i]][,2]-as.numeric(centoid_data[centoid_data$site==names(ordi_data[i]),"PC2"]))^2) %>% sqrt() 
  ordi_data[[i]] <- cbind(ordi_data[[i]],distances)
  x <- c(names(ordi_data[i]), mean(distances[1:length(distances)-1])) #hear the mean is the mean of centroid distances per site. you can also set it as sum etc. 
  site_info[site_info$site==x[1], nsite] <- x[2]
}
colnames(site_info) <- c("site" ,"Class" ,"alteration", "groups", "Alteration_type" ,"Catchment", "Upstream_Catchment_size", "Elevation_in_m" ,"mean_centroid_distance")

# Anova on the averaged site distances to centroid
data_anova <- site_info %>%
  filter(Class=="Temperate") %>%
  select(mean_centroid_distance, Alteration_type)

anova(lm(data_anova$mean_centroid_distance ~ factor(data_anova$Alteration_type)))

# Permanova on groups
adonis(as.numeric(data_anova$mean_centroid_distance)~ factor(data_anova$Alteration_type), data=data_anova, method="euclidian", perm=999)

# Average centroid distance boxplots 
ordi_data2 <- plyr::ldply(ordi_data) %>% unique() %>% as_tibble()
colnames(ordi_data2) <- c("site", "PC1","PC2","Distances")
ordi_data2 <- left_join(ordi_data2, site_info, by=c("site"="site")) %>% as_tibble()
ordi_data2$Alteration_type<-as_factor(ordi_data2$Alteration_type)

ordi_data2 %>%
  filter(Class=="Temperate") %>%
  ggplot() +
  geom_boxplot(aes(x=reorder(site, as.numeric(Alteration_type)), y=Distances, fill=Alteration_type))+
  theme_classic()+ coord_flip()+
  labs(color="Groups", x="Sites", y="Distance to centroid", title="Average Distance to Centroid")+
  #scale_fill_manual(values= c("blue", "orange", "orangered"))
  scale_fill_manual(values= c("blue","green", "darkgreen"))       

#### Flow/Parameter Analysis ####

#assign the real dates here
data_sum3 <- data_sum 
#for now, I'll use normalised componenents as qualitative measure
data_sum3 <- data_sum3 %>%
  mutate(Comp1=Comp1/BDOC.x,
         Comp2=Comp2/BDOC.x,
         Comp3=Comp3/BDOC.x,
         Comp4=Comp4/BDOC.x,
         Comp5=Comp5/BDOC.x,
         Comp6=Comp6/BDOC.x,
         Comp7=Comp7/BDOC.x,
         Comp8=Comp8/BDOC.x
  )

data_sum3["exact_date"] <- ''

for(i in 1:nrow(data_sum3)){
  if (data_sum3$campaign[[i]] == "apr") {
    data_sum3$exact_date[[i]] <- as.character(dmy("10/04/2018"))
  } else if (data_sum3$campaign[[i]] == "oct") {
    data_sum3$exact_date[[i]] <- as.character(dmy("28/09/2017"))
  } else if (data_sum3$campaign[[i]] == "dec"){
    data_sum3$exact_date[[i]] <- as.character(dmy("12/12/2017"))
  } else if (data_sum3$campaign[[i]] == "feb"){
    data_sum3$exact_date[[i]] <-  as.character(dmy("13/02/2018"))
  } else if (data_sum3$campaign[[i]] == "may"){
    data_sum3$exact_date[[i]] <-  as.character(dmy("12/06/2018"))
  } else if (data_sum3$campaign[[i]] == "aug") {
    data_sum3$exact_date[[i]] <- as.character(dmy("15/08/2018"))
  } else {data_sum3$exact_date[[i]] <- (0)}
}
data_sum3$exact_date <- as.Date(data_sum3$exact_date)

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

#for(j in seq_along(group_list)){
j <- 1
plot <- list()
scale_factor <- 100
for(i in seq_along(group_list[[j]][[1]])){
  plot[[i]] <- ggplot(flow_data) + 
    geom_area(aes_string(x="Date", y=as.name(vars_select(colnames(flow_data), starts_with(group_list[[j]][[1]][[i]])))))+
    geom_line(data=subset(data_sum3, site==group_list[[j]][[1]][[i]]), aes(x=exact_date, y=Comp1*scale_factor), color="blue")+ 
    geom_point(data=subset(data_sum3, site==group_list[[j]][[1]][[i]]), aes(x=exact_date, y=Comp1*scale_factor), color="blue")+ 
    geom_line(data=subset(data_sum3, site==group_list[[j]][[1]][[i]]), aes(x=exact_date, y=Comp2*scale_factor), color="blue")+
    geom_point(data=subset(data_sum3, site==group_list[[j]][[1]][[i]]), aes(x=exact_date, y=Comp2*scale_factor), color="blue")+
    geom_line(data=subset(data_sum3, site==group_list[[j]][[1]][[i]]), aes(x=exact_date, y=Comp3*scale_factor), color="blue")+
    geom_point(data=subset(data_sum3, site==group_list[[j]][[1]][[i]]), aes(x=exact_date, y=Comp3*scale_factor), color="blue")+
    geom_line(data=subset(data_sum3, site==group_list[[j]][[1]][[i]]), aes(x=exact_date, y=Comp4*scale_factor), color="blue")+
    geom_point(data=subset(data_sum3, site==group_list[[j]][[1]][[i]]), aes(x=exact_date, y=Comp4*scale_factor), color="blue")+
    geom_line(data=subset(data_sum3, site==group_list[[j]][[1]][[i]]), aes(x=exact_date, y=Comp5*scale_factor), color="blue")+
    geom_point(data=subset(data_sum3, site==group_list[[j]][[1]][[i]]), aes(x=exact_date, y=Comp5*scale_factor), color="blue")+
    geom_line(data=subset(data_sum3, site==group_list[[j]][[1]][[i]]), aes(x=exact_date, y=Comp6*scale_factor), color="blue")+
    geom_point(data=subset(data_sum3, site==group_list[[j]][[1]][[i]]), aes(x=exact_date, y=Comp6*scale_factor), color="red")+
    geom_line(data=subset(data_sum3, site==group_list[[j]][[1]][[i]]), aes(x=exact_date, y=Comp7*scale_factor), color="red")+
    geom_point(data=subset(data_sum3, site==group_list[[j]][[1]][[i]]), aes(x=exact_date, y=Comp7*scale_factor), color="red")+
    geom_line(data=subset(data_sum3, site==group_list[[j]][[1]][[i]]), aes(x=exact_date, y=Comp8*scale_factor), color="red")+
    geom_point(data=subset(data_sum3, site==group_list[[j]][[1]][[i]]), aes(x=exact_date, y=Comp8*scale_factor), color="red")+
    scale_y_continuous("Flow", sec.axis = sec_axis(~ (./scale_factor), name = "Normalised Components"))+
    scale_x_date("Day")+
    theme_linedraw()+
    ggtitle(paste(group_list[[j]][[1]][[i]]))
} 
print(do.call(grid.arrange, plot))
#} 

# plot with normalized component PC's
data_sum4 <- cbind(wine.pca$x[,1:2], as.data.frame(BDOC_normalised_components[,c("site","campaign")]))
for(i in 1:nrow(data_sum4)){
  if (data_sum4$campaign[[i]] == "apr") {
    data_sum4$exact_date[[i]] <- as.character(dmy("10/04/2018"))
  } else if (data_sum4$campaign[[i]] == "oct") {
    data_sum4$exact_date[[i]] <- as.character(dmy("28/09/2017"))
  } else if (data_sum4$campaign[[i]] == "dec"){
    data_sum4$exact_date[[i]] <- as.character(dmy("12/12/2017"))
  } else if (data_sum4$campaign[[i]] == "feb"){
    data_sum4$exact_date[[i]] <-  as.character(dmy("13/02/2018"))
  } else if (data_sum4$campaign[[i]] == "may"){
    data_sum4$exact_date[[i]] <-  as.character(dmy("12/06/2018"))
  } else if (data_sum4$campaign[[i]] == "aug") {
    data_sum4$exact_date[[i]] <- as.character(dmy("15/08/2018"))
  } else {data_sum4$exact_date[[i]] <- (0)}
}
data_sum4$exact_date <- as.Date(data_sum4$exact_date)

j <- 1
plot <- list()
scale_factor <- 10
for(i in seq_along(group_list[[j]][[1]])){
  plot[[i]] <- ggplot(flow_data) + 
    geom_area(aes_string(x="Date", y=as.name(vars_select(colnames(flow_data), starts_with(group_list[[j]][[1]][[i]])))))+
    geom_line(data=subset(data_sum4, site==group_list[[j]][[1]][[i]]), aes(x=exact_date, y=PC1*scale_factor), color="blue", size=1)+ 
    geom_point(data=subset(data_sum4, site==group_list[[j]][[1]][[i]]), aes(x=exact_date, y=PC1*scale_factor), color="blue")+ 
    geom_line(data=subset(data_sum4, site==group_list[[j]][[1]][[i]]), aes(x=exact_date, y=PC2*scale_factor), color="red", size=1)+
    geom_point(data=subset(data_sum4, site==group_list[[j]][[1]][[i]]), aes(x=exact_date, y=PC2*scale_factor), color="red")+
    scale_y_continuous("Flow",  sec.axis = sec_axis(~ (./scale_factor), name = "Normalised Components"))+
    scale_x_date("Day")+
    theme_linedraw()+
    ggtitle(paste(group_list[[j]][[1]][[i]]))
}
print(do.call(grid.arrange, plot))
print(do.call(grid.arrange, plot_PC))


#### Loading plots ####
parafac_em <- read.csv("Em.csv")
parafac_ex<- read.csv("Ex.csv")
ggplot()+
  geom_line(aes(x=parafac_em$Em,y=(1/3)*(parafac_em$Comp.1)))+
  geom_line(aes(x=parafac_ex$Ex,y=3*(parafac_ex$Comp.1)), linetype = "dashed")+
  theme_classic()+
  labs(title="Comp1")

