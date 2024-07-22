library(data.table)

#### Getting the data####
sample_loadings <- as.data.table(read.csv("data/PARAFAC_data/denormalised_sample_loadings.csv"))
sample_names <- as.data.table(read.csv("data/metafiles/sample_name_site_match.csv"))
site_info <- as.data.table(read.csv2("data/metafiles/site_summary.csv"))

oct1 <- as.data.table(cbind(read.table("data/abs_data/oct.txt", sep=",", header=TRUE), campaign="oct"))
feb1 <- as.data.table(cbind(read.table("data/abs_data/feb.txt", sep=",", header=TRUE), campaign="feb"))
apr1 <- as.data.table(cbind(read.table("data/abs_data/apr.txt", sep=",", header=TRUE), campaign="apr"))
may1 <- as.data.table(cbind(read.table("data/abs_data/may.txt", sep=",", header=TRUE), campaign="may"))
aug1 <- as.data.table(cbind(read.table("data/abs_data/aug.txt", sep=",", header=TRUE), campaign="aug"))
dec1 <- as.data.table(cbind(read.table("data/abs_data/dec.txt", sep=",", header=TRUE), campaign="dec"))

oct1[,"site"] <- substr(oct1$site, 1, nchar(oct1$site)-4)
apr1[,"site"] <- substr(apr1$site, 1, nchar(apr1$site)-4)

abs_data <- rbind(oct1,feb1, apr1, may1, aug1,dec1)
abs_data <- subset(abs_data, subset = !grepl("^Blank", abs_data$site))

abs_data[abs_data[, HIX == "Inf"]]$HIX <- NA
abs_data[abs_data[, HIX > 100]]$HIX <- NA
abs_data[abs_data[, E2.to.E3 < 0]]$E2.to.E3 <- NA
abs_data[abs_data[, E2.to.E3 > 30]]$E2.to.E3 <- NA
abs_data[abs_data[, E4.to.E6 < 0]]$E4.to.E6 <- NA

data_means <- aggregate(x = abs_data[, c("FIX", "FIX2", "HIX", "HIX2", "beta.alpha","DecAbsCoeff254", "DecAbsCoeff420", "DecAbsCoeff430", "DecAbsCoeff436", "E2.to.E3", "E4.to.E6", 
                                        "slope_classic", "slope_lm", "slope_short_Loiselle", "slope_short_Helms", "SR_Loiselle", "SR_Helms")], 
                       by = abs_data[,c("site", "campaign")], FUN = mean, na.rm = T, trim = 0.25)

meta_sum_optical <- as.data.table(cbind(read.csv("data/metafiles/meta_file_2.csv"), data_means), key = c("site","campaign"))

# Converting the sample names without the index numbers - wouldn't need it if data were named with 000
sample_loadings[, c("campaign", "run", "sample2", "num") := tstrsplit(sample, "_")]
sample_loadings[, c("sample_code") := paste0(sample_loadings$campaign,"_", sample_loadings$run, "_", sample_loadings$sample2)]
setdiff(sample_names$sample_code, sample_loadings$sample_code) 
#Check why these are not matching. 
#In my case, they are the samples removed during parafac because of weird looking eems.
#So it makes sense to remove them here too since there might be a contamination affecting all

##### Binding the data frames ####
meta_component_data <- merge(sample_loadings, sample_names, by.x = c("sample_code", "campaign"), by.y = c("sample_code", "campaign"))

# We excluded Carrion from here onwards due to very low concentrations on all seasons
meta_component_data <- meta_component_data[meta_component_data$site != "Carrion",]
meta_sum_optical <- meta_sum_optical[meta_sum_optical$site != "Carrion",]
sample_names <- sample_names[sample_names$site != "Carrion",]
site_info <- site_info[site_info$site != "Carrion",]

##### Prepping for the DOC mean and variance plots ####
data_sum <- aggregate(meta_component_data[,c("Comp.1","Comp.2","Comp.3","Comp.4","Comp.5","Comp.6","Comp.7","Comp.8")], 
                     by = meta_component_data[, c("campaign", "site")], FUN = mean, na.rm = T)

data_sum <- merge(site_info, data_sum, by.x = c("site"), by.y = c("site"))
data_sum <- merge(data_sum, meta_sum_optical[, -c(14, 15)], by = c("site","campaign"))
data_sum[, "SUVA254" := DecAbsCoeff254/NPOC]
data_sum[data_sum[, SUVA254 > 6]]$SUVA254 <- NA #This is definitely a weird point. Sama sample as below. Seems like an outlier. Changes the PLSR table of now
data_sum[data_sum[, SR_Loiselle > 1.7]]$SR_Loiselle <- NA

variables_to_convert <- c("HMWS_C", "humic_like_substance_C", "LMWS_C", "HMWS_N", "HS", "SUVA_HS", "SUVA_ges")
data_sum[, (variables_to_convert) := lapply(.SD, convert_below_detection_limit), .SDcols = variables_to_convert]

# Select the variables to be used for the rest of the analysis/paper
data_sum2 <- data_sum[,c("campaign", "site", "groups.x","Class","alteration", "alteration_type", "alteration_type_grouping", "storage_index","Comp.1", "Comp.2",  "Comp.3", "Comp.4", "Comp.5", "Comp.6", "Comp.7", "Comp.8", "BDOC", "NPOC", "SUVA254",
                        "FIX", "FIX2", "HIX2", "beta.alpha", "slope_classic", "slope_lm", "slope_short_Helms", "slope_short_Loiselle",
                        "SR_Helms", "SR_Loiselle", "E2.to.E3", "E4.to.E6", "DecAbsCoeff254")]

#### DOC Concentration #### 
##### 1. Data Mean and CV Calculation #####
DOC_sum <- merge(data_sum2[, .(mean_BDOC = mean(BDOC, na.rm = TRUE), mean_NPOC = mean(NPOC, na.rm = TRUE), 
                              var_BDOC = CV(BDOC), var_NPOC = CV(NPOC)), by = .(site)], 
                site_info[site != "Carrion", c("site", "Class", "alteration", "groups", "alteration_type", "alteration_type_grouping", "storage_index")])

group_means <- DOC_sum[, .(means = mean(mean_NPOC, na.rm = TRUE) , sd_NPOC = sd(mean_NPOC)), by = .(groups)]
group_CVs <- DOC_sum[, .(means = mean(var_NPOC, na.rm = TRUE) , sd_NPOC = sd(var_NPOC)), by = .(groups)]
data_sum$groups.x <- factor(data_sum$groups.x, levels = c("TempNat", "TempAlt", "MedNat", "MedAlt"))
data_sum$alteration_type_grouping <- factor(data_sum$alteration_type_grouping, levels = c("TempNat", "TempAlt_irrigation", "TempAlt_hydropower",
                                                                          "MedNat", "MedAlt_irrigation"))

data_sum$alteration_type  <- factor(data_sum$alteration_type, levels = c("irrigation", "hydropower", "natural"))
