library(data.table)
data=as.data.table(read_excel("HYDRA_env_fluvialnetwork_input.xlsx",  sheet = "Climate"))
t.test(data[Climate=="Mediterranean", LC_PREC], data[Climate=="Temperate", LC_PREC])
t.test(data[Climate=="Mediterranean", LC_TEMP], data[Climate=="Temperate", LC_TEMP])

data=as.data.table(read_excel("HYDRA_env_fluvialnetwork_input.xlsx",  sheet = "Landuse"))
t.test(data[Climate=="Mediterranean", LC_WAE], data[Climate=="Temperate", LC_WAE])
t.test(data[Climate=="Mediterranean", MN_DEN], data[Climate=="Temperate", MN_DEN])
t.test(data[Climate=="Mediterranean", BF_UHD], data[Climate=="Temperate", BF_UHD])
t.test(data[Climate=="Mediterranean", MN_WAE], data[Climate=="Temperate", MN_WAE])


data=as.data.table(read_excel("HYDRA_env_fluvialnetwork_input.xlsx",  sheet = "Geo"))
t.test(data[Climate=="Mediterranean", MN_watr], data[Climate=="Temperate", MN_watr])
t.test(data[Climate=="Mediterranean", BF_PERM], data[Climate=="Temperate", BF_PERM])
