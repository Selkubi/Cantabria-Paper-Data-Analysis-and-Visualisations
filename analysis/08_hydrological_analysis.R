library(data.table)
library(ggplot2)

site_info <- as.data.table(read.csv2("data/metafiles/site_summary.csv"))
flow_data <- data.table(read.csv("data/hydrological_data/Flow_HYDRA_2017_2018_cleanSK.csv", header = TRUE))
flow_data$date <- lubridate::mdy(flow_data$Date)
flow <- flow_data[flow_data$date < "2018-10-01" & flow_data$date >= "2017-10-01"]
#flow = flow_data use this if you want all years as part of the anaysis instead of the measurement years

flow[, c("year", "month", "day") := tstrsplit(date, "-")]
flow[, c("year", "month", "day") := lapply(flow[, c("year", "month", "day")], FUN = as.numeric)]

monthly_means <- data.table(aggregate(x = flow[, -c("month", "Date", "day", "date")],
                                      by = flow[, c("month")], FUN = mean, na.rm = TRUE, nan.rm = TRUE))
rivers <- c(colnames(monthly_means[, 2:22]))
monthly_max <- t(monthly_means[, lapply(.SD, max, na.rm = TRUE), .SDcols = rivers ])
annual_max_discharge <- data.table("variable" = row.names(monthly_max), "annual_max" = monthly_max[, 1])

means <- melt(monthly_means[, -22], id.vars = "month")
means <- means[annual_max_discharge, on = .(variable = variable)]
means <- data.table(merge(means, site_info, by.x = "variable", by.y = "site"))
means$month <- as.factor(means$month)
means <- means[means$variable != "Carrion", ]

means_summary <- means[, .(min = min(value, na.rm = TRUE),
                           max = max(value, na.rm = TRUE),
                           mean = mean(value, na.rm = TRUE),
                           normalized_mean = mean(value / annual_max, na.rm = TRUE)),
                      by = .(alteration_type_grouping, month, Class)]

means_summary$alteration_type_grouping <- factor(means_summary$alteration_type_grouping, levels = c("TempNat", "TempAlt_hydropower", "TempAlt_irrigation", "MedNat", "MedAlt_irrigation"))
means$alteration_type_grouping <- factor(means$alteration_type_grouping, levels = c("TempNat", "TempAlt_hydropower", "TempAlt_irrigation", "MedNat", "MedAlt_irrigation"))
means$groups <- factor(means$groups, levels = c("TempNat", "TempAlt", "MedNat", "MedAlt"))
labels_plot <- c('TempNat' = "nA", 
            'TempAlt_hydropower' = "aA - Hydropower", 
            'TempAlt_irrigation' = "aA - Irrigation", 
            'MedNat' = "nM",
            'MedAlt_irrigation' = "aM - Irrigation")


hydrographs <- ggplot(means) +
  facet_wrap(~ alteration_type_grouping, scales = "free",
             labeller = as_labeller(labels_plot)) +
  geom_vline(xintercept = c(2, 4, 6, 8, 10, 12), color = "grey", alpha = 0.35, size = 4.5) +
  geom_point(aes(x = month, y = value / annual_max, group = variable, shape = alteration_type_grouping,
                 fill = alteration_type_grouping,  
                 colour = alteration_type_grouping)) +
  geom_line(aes(x = month, y = value / annual_max, group = variable, color = alteration_type_grouping), size = 0.5, alpha = 0.7) +
  geom_line(data = means_summary, aes(x = month, y = normalized_mean, group = alteration_type_grouping, color = alteration_type_grouping), size = 1.5) +
  scale_color_manual(values = c("#B4DCED", "#6996D1", "#2B5FA2", "#F5CB7D", "#F09E41"), labels = labels_plot) +
  scale_fill_manual(values = c("#B4DCED", "#6996D1", "#2B5FA2", "#F5CB7D", "#F09E41"), labels = labels_plot) +
  scale_x_discrete(labels = c("1" = "J", "2" = "F", "3" = "M",
                              "4" = "A", "5" = "M", "6" = "J",
                              "7" = "J", "8" = "A", "9" = "S",
                              "10" = "O", "11" = "N", "12" = "D")) +
  scale_shape_manual(values = c(23, 22, 21, 25, 24), labels = labels_plot) +
  labs(x = "Months", y = paste("Normalised mean monthly flow")) +
  theme_pca() + theme(legend.position = c(0.85, 0.2), legend.title = element_blank())

pdf("output/plots/hydrographs.pdf", width = 6, height = 4)
plot(hydrographs)
dev.off()

