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

hydrographs <- ggplot(means) +
  facet_wrap(~ alteration_type_grouping, scales = "free",
             labeller = as_labeller(c(TempNat = "nA", TempAlt_hydropower = "aA - Hydropower", TempAlt_irrigation = "aA - Irrigation", MedNat = "nM", MedAlt_irrigation = "aM"))) +
  geom_line(aes(x = month, y = value / annual_max, group = variable, color = alteration_type_grouping), size = 0.5, alpha = 0.7) +
  #directlabels::geom_dl(aes(x = month, y = value/annual_max, group = variable, label = variable, color=groups), method = "first.qp") +
  geom_line(data = means_summary, aes(x = month, y = normalized_mean, group = alteration_type_grouping, color = alteration_type_grouping), size = 1.5) +
  geom_point(aes(x = month, y = value / annual_max, group = variable, color = alteration_type_grouping), size = 1.2, alpha = 0.7) +
  scale_color_manual(values = c("#B4DCED", "#6996D1", "#2B5FA2", "#F5CB7D", "#F09E41")) +
  scale_x_discrete(labels = c("1" = "J", "2" = "F", "3" = "M",
                              "4" = "A", "5" = "M", "6" = "J",
                              "7" = "J", "8" = "A", "9" = "S",
                              "10" = "O", "11" = "N", "12" = "D")) +
  labs(x = "Months", y = paste("Mean Montly flow (m3/s)")) +
  theme_pca()

pdf("output/plots/hydrographs.pdf", width = 5, height = 4)
plot(hydrographs)
dev.off()
