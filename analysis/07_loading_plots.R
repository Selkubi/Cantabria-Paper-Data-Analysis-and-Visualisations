library(ggplot2)
parafac_em <- read.csv("data/PARAFAC_data/Em.csv")
parafac_ex<- read.csv("data/PARAFAC_data/Ex.csv")
ggplot()+
  geom_line(aes(x=parafac_em$Em,y=(1/3)*(parafac_em$Comp.8)))+
  geom_line(aes(x=parafac_ex$Ex,y=3*(parafac_ex$Comp.8)), linetype = "dashed")+
  theme_classic()+
  labs(title="Comp8")