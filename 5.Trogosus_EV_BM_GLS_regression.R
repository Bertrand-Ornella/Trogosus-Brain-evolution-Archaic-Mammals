
##### script for linear regression and graph - Trogosus paper

library(RRPP) #pairwise comparisons

######################## ggplot - linear regression - Lambda ##########################

#open data for analyses # change rowname of X to Species
Trogosus.data<-read.csv("Trogosus_dataset_final_Eocene.csv", header=T)
names(Trogosus.data)[names(Trogosus.data) == "X"] <- "Species"
Trogosus.data<-subset(Trogosus.data, Trogosus.data[,6] =="middle Eocene")

#Transform data to log10
Trogosus.data$Brain_volume_cm3<-log10(Trogosus.data$Brain_volume_cm3)
names(Trogosus.data)[names(Trogosus.data) == "Brain_volume_cm3"] <- "Brain"
Trogosus.data$Body_mass_g<-log10(Trogosus.data$Body_mass_g)
names(Trogosus.data)[names(Trogosus.data) == "Body_mass_g"] <- "Body"

EocC<-Trogosus.data[ which(Trogosus.data$Group=='Eocene crown clades'), ]
EoclineC_Br_B <-gls(Brain ~ Body, data=EocC)
pgls.fit.EocC <- predict(EoclineC_Br_B) #predict values for brain size
predframe.EocC <- with(EocC, data.frame(Species, Group, Body, Brain = pgls.fit.EocC))

sink("EV_BM_ANOVA_MiddleEocene_Reg_Crown.txt")
summary(EoclineC_Br_B)
sink()

EocS<-Trogosus.data[ which(Trogosus.data$Group=='Eocene archaic taxa'), ]
EoclineS_Br_B <-gls(Brain ~ Body, data=EocS)
pgls.fit.EocS <- predict(EoclineS_Br_B) #predict values for brain size
predframe.EocS <- with(EocS, data.frame(Species, Group, Body, Brain = pgls.fit.EocS))

sink("EV_BM_ANOVA_MiddleEocene_Reg_Archaic.txt")
summary(EoclineS_Br_B)
sink()

#GLS to obtain residuals
Br_Bo<-Trogosus.data[ which(Trogosus.data$Group=='Eocene crown clades'| Trogosus.data$Group=='Eocene archaic taxa'), ]
GLS_Br_Bo<-gls(Brain ~ Body, data=Br_Bo)
pgls.fit.Br_Bo <- predict(GLS_Br_Bo) #predict values for brain size
predframe.Br_Bo <- with(Br_Bo, data.frame(Species, Group, Body, Brain = pgls.fit.Br_Bo))

sink("EV_BM_ANOVA_MiddleEocene_Reg_ALL.txt")
summary(GLS_Br_Bo) #Final model
sink()

#######
#Make graph GLS regressions
ggplot(Trogosus.data, aes(Body, Brain, color = Group)) +
  geom_point(data = dplyr::filter(Trogosus.data, Group == "Eocene archaic taxa"),
             size = 2, aes(color = "lightcoral")) +
  geom_point(data = dplyr::filter(Trogosus.data, Group == "Eocene crown clades"),
             size = 2, aes(color = "lightblue")) +
  theme_minimal() + 
  scale_color_manual(name = "", values = c("lightblue","lightcoral"),labels = c("Eocene crown clades","Eocene archaic taxa")) +
  geom_line(data = dplyr::filter(predframe.Br_Bo, Group == "Eocene archaic taxa"), color = "black",
            linetype = "dashed") +
  geom_line(data = dplyr::filter(predframe.Br_Bo, Group == "Eocene crown clades"), color = "black",
            linetype = "dashed") +
  geom_line(data = dplyr::filter(predframe.EocS, Group == "Eocene archaic taxa"), color = "lightcoral",
            linetype = 1.5) +
  geom_line(data = dplyr::filter(predframe.EocC, Group == "Eocene crown clades"), color = "lightblue",
            linetype = 1.5) +
  #theme(legend.position = "top") +
  #geom_text(data = dplyr::filter(Trogosus.data, Group == "Eocene archaic taxa"), color = "lightcoral",
  #aes(label = Abbreviations), hjust = -0.3, vjust = 1.1)+ 
  #geom_text(data = dplyr::filter(Trogosus.data, Group == "Eocene crown clades"), color = "lightblue",
  #aes(label = Abbreviations), hjust = -0.3, vjust = 1.1) +           
  labs(x = "log10(Body mass)", y = "log10(Endocranial volume)") +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12,face = "bold")) 

##### ANOVA

sink("EV_BM_ANOVA_MiddleEocene_regressions.txt")
fitOLS<-lm.rrpp(Brain~Body + Group, data=Trogosus.data)
summary(fitOLS) # significant 0.001 ***
anova(fitOLS)
coef(fitOLS, test=TRUE)
PT1 <- pairwise(fitOLS, groups = Trogosus.data$Group)
summary(PT1, confidence = 0.95)
sink()

## END!

