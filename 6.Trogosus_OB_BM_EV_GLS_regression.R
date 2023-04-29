
##### script for linear regression and graph - Trogosus paper

library(RRPP) #pairwise comparaisons

######################## 1. ggplot - linear regression - OB-BM ##########################

#open data for analyses # change rowname of X to Species
Trogosus.data<-read.csv("OB_BM_EV_Trogosus_Eocene.csv", header=T)
names(Trogosus.data)[names(Trogosus.data) == "X"] <- "Species"
Trogosus.data<-subset(Trogosus.data, Trogosus.data[,6] =="middle Eocene")

#Transform data to log10
Trogosus.data$Olfactory_bulbs_mm3<-log10(Trogosus.data$Olfactory_bulbs_mm3)
names(Trogosus.data)[names(Trogosus.data) == "Olfactory_bulbs_mm3"] <- "OB"
Trogosus.data$Body_mass_mg<-log10(Trogosus.data$Body_mass_mg)
names(Trogosus.data)[names(Trogosus.data) == "Body_mass_mg"] <- "Body"

EocC<-Trogosus.data[ which(Trogosus.data$Group=='Eocene crown clades'), ]
EoclineC_Br_B <-gls(OB ~ Body, data=EocC)
pgls.fit.EocC <- predict(EoclineC_Br_B) #predict values for OB size
predframe.EocC <- with(EocC, data.frame(Species, Group, Body, OB = pgls.fit.EocC))

sink("OB_BM_ANOVA_MiddleEocene_Reg_Crown.txt")
summary(EoclineC_Br_B)
sink()

EocS<-Trogosus.data[ which(Trogosus.data$Group=='Eocene archaic taxa'), ]
EoclineS_Br_B <-gls(OB ~ Body, data=EocS)
pgls.fit.EocS <- predict(EoclineS_Br_B) #predict values for OB size
predframe.EocS <- with(EocS, data.frame(Species, Group, Body, OB = pgls.fit.EocS))

sink("OB_BM_ANOVA_MiddleEocene_Reg_Archaic.txt")
summary(EoclineS_Br_B)
sink()

#GLS to obtain residuals
OB_BM<-Trogosus.data[ which(Trogosus.data$Group=='Eocene crown clades'| Trogosus.data$Group=='Eocene archaic taxa'), ]
GLS_OB_BM<-gls(OB ~ Body, data=OB_BM)
pgls.fit.OB_BM <- predict(GLS_OB_BM) #predict values for OB size
predframe.OB_BM <- with(OB_BM, data.frame(Species, Group, Body, OB = pgls.fit.OB_BM))

sink("OB_BM_ANOVA_MiddleEocene_Reg_ALL.txt")
summary(GLS_OB_BM) #Final model
sink()

#Make graph with PGLS corrected regressions - OB-Body
ggplot(Trogosus.data, aes(Body, OB, color = Group)) +
  geom_point(data = dplyr::filter(Trogosus.data, Group == "Eocene archaic taxa"),
             size = 2, aes(color = "lightcoral")) +
  geom_point(data = dplyr::filter(Trogosus.data, Group == "Eocene crown clades"),
             size = 2, aes(color = "lightblue")) +
  theme_minimal() + 
  scale_color_manual(name = "", values = c("lightblue","lightcoral"),labels = c("Eocene crown clades","Eocene archaic taxa")) +
  geom_line(data = dplyr::filter(predframe.OB_BM, Group == "Eocene archaic taxa"), color = "black",
            linetype = "dashed") +
  geom_line(data = dplyr::filter(predframe.OB_BM, Group == "Eocene crown clades"), color = "black",
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
  labs(x = "log10(Body mass)", y = "log10(Olfactory bulb volume)") +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12,face = "bold")) 

##### ANOVA
sink("OB_BM_ANOVA_MiddleEocene_regressions.txt")
fitOLS<-lm.rrpp(OB~Body + Group, data=Trogosus.data)
summary(fitOLS) # significant 0.001 ***
anova(fitOLS)
coef(fitOLS, test=TRUE)
PT1 <- pairwise(fitOLS, groups = Trogosus.data$Group)
summary(PT1, confidence = 0.95)
sink()

######################## 2. ggplot - linear regression - OB-EV ##########################

#open data for analyses # change rowname of X to Species
Trogosus.data<-read.csv("OB_BM_EV_Trogosus_Eocene.csv", header=T)
names(Trogosus.data)[names(Trogosus.data) == "X"] <- "Species"
Trogosus.data<-subset(Trogosus.data, Trogosus.data[,6] =="middle Eocene")

#Transform data to log10
Trogosus.data$Olfactory_bulbs_mm3<-log10(Trogosus.data$Olfactory_bulbs_mm3)
names(Trogosus.data)[names(Trogosus.data) == "Olfactory_bulbs_mm3"] <- "OB"
Trogosus.data$Brain_volume_mm3<-log10(Trogosus.data$Brain_volume_mm3)
names(Trogosus.data)[names(Trogosus.data) == "Brain_volume_mm3"] <- "Brain"

EocC<-Trogosus.data[ which(Trogosus.data$Group=='Eocene crown clades'), ]
EoclineC_Br_B <-gls(OB ~ Brain, data=EocC)
pgls.fit.EocC <- predict(EoclineC_Br_B) #predict values for OB size
predframe.EocC <- with(EocC, data.frame(Species, Group, Brain, OB = pgls.fit.EocC))

sink("OB_EV_ANOVA_MiddleEocene_Reg_Crown.txt")
summary(EoclineC_Br_B)
sink()

EocS<-Trogosus.data[ which(Trogosus.data$Group=='Eocene archaic taxa'), ]
EoclineS_Br_B <-gls(OB ~ Brain, data=EocS)
pgls.fit.EocS <- predict(EoclineS_Br_B) #predict values for OB size
predframe.EocS <- with(EocS, data.frame(Species, Group, Brain, OB = pgls.fit.EocS))

sink("OB_EV_ANOVA_MiddleEocene_Reg_Archaic.txt")
summary(EoclineS_Br_B)
sink()

#GLS to obtain residuals
OB_EV<-Trogosus.data[ which(Trogosus.data$Group=='Eocene crown clades'| Trogosus.data$Group=='Eocene archaic taxa'), ]
GLS_OB_EV<-gls(OB ~ Brain, data=OB_EV)
pgls.fit.OB_EV <- predict(GLS_OB_EV) #predict values for OB size
predframe.OB_EV <- with(OB_EV, data.frame(Species, Group, Brain, OB = pgls.fit.OB_EV))

sink("OB_EV_ANOVA_MiddleEocene_Reg_ALL.txt")
summary(GLS_OB_EV) #Final model
sink()

#Make graph with PGLS corrected regressions - OB-Brain
ggplot(Trogosus.data, aes(Brain, OB, color = Group)) +
  geom_point(data = dplyr::filter(Trogosus.data, Group == "Eocene archaic taxa"),
             size = 2, aes(color = "lightcoral")) +
  geom_point(data = dplyr::filter(Trogosus.data, Group == "Eocene crown clades"),
             size = 2, aes(color = "lightblue")) +
  theme_minimal() + 
  scale_color_manual(name = "", values = c("lightblue","lightcoral"),labels = c("Eocene crown clades","Eocene archaic taxa")) +
  geom_line(data = dplyr::filter(predframe.OB_EV, Group == "Eocene archaic taxa"), color = "black",
            linetype = "dashed") +
  geom_line(data = dplyr::filter(predframe.OB_EV, Group == "Eocene crown clades"), color = "black",
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
  labs(x = "log10(Endocranial volume)", y = "log10(Olfactory bulbs volume)") +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12,face = "bold")) 

##### ANOVA

sink("OB_EV_ANOVA_MiddleEocene_regressions.txt")
fitOLS<-lm.rrpp(OB~Brain + Group, data=Trogosus.data)
summary(fitOLS) # significant 0.001 ***
anova(fitOLS)
coef(fitOLS, test=TRUE)
PT1 <- pairwise(fitOLS, groups = Trogosus.data$Group)
summary(PT1, confidence = 0.95)
sink()
  
## END!

