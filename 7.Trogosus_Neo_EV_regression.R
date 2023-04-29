
##### script for linear regression and graph - Trogosus paper

library(RRPP) #pairwise comparaisons

######################## ggplot - linear regression - Neocortex ##########################

#open data for analyses # change rowname of X to Species
Trogosus.data<-read.csv("Neo_EV_Trogosus_Middle_Eocene.csv", header=T)
names(Trogosus.data)[names(Trogosus.data) == "X"] <- "Species"

#Transform data to log10
Trogosus.data$Neocortex_surface_mm2<-log10(Trogosus.data$Neocortex_surface_mm2)
names(Trogosus.data)[names(Trogosus.data) == "Neocortex_surface_mm2"] <- "Neo"
Trogosus.data$Brain_surface_mm2<-log10(Trogosus.data$Brain_surface_mm2)
names(Trogosus.data)[names(Trogosus.data) == "Brain_surface_mm2"] <- "BrainS"

EocC<-Trogosus.data[ which(Trogosus.data$Group=='Eocene crown clades'), ]
EoclineC_Br_B <-gls(Neo ~ BrainS, data=EocC)
pgls.fit.EocC <- predict(EoclineC_Br_B) #predict values for brain size
predframe.EocC <- with(EocC, data.frame(Species, Group, BrainS, Neo = pgls.fit.EocC))

sink("Neo_EV_ANOVA_MiddleEocene_Reg_Crown.txt")
summary(EoclineC_Br_B)
sink()

EocS<-Trogosus.data[ which(Trogosus.data$Group=='Eocene archaic taxa'), ]
EoclineS_Br_B <-gls(Neo ~ BrainS, data=EocS)
pgls.fit.EocS <- predict(EoclineS_Br_B) #predict values for brain size
predframe.EocS <- with(EocS, data.frame(Species, Group, BrainS, Neo = pgls.fit.EocS))

sink("Neo_EV_ANOVA_MiddleEocene_Reg_Archaic.txt")
summary(EoclineS_Br_B)
sink()

#GLS to obtain residuals
Neo_Br<-Trogosus.data[ which(Trogosus.data$Group=='Eocene crown clades'| Trogosus.data$Group=='Eocene archaic taxa'), ]
GLS_Neo_Br<-gls(Neo ~ BrainS, data=Neo_Br)
pgls.fit.Neo_Br <- predict(GLS_Neo_Br) #predict values for brain size
predframe.Neo_Br <- with(Neo_Br, data.frame(Species, Group, BrainS, Neo = pgls.fit.Neo_Br))

sink("Neo_EV_ANOVA_MiddleEocene_Reg_ALL.txt")
summary(GLS_Neo_Br) #Final model
sink()

#Make graph with PGLS corrected regressions 
ggplot(Trogosus.data, aes(BrainS, Neo, color = Group)) +
  geom_point(data = dplyr::filter(Trogosus.data, Group == "Eocene archaic taxa"),
             size = 2, aes(color = "lightcoral")) +
  geom_point(data = dplyr::filter(Trogosus.data, Group == "Eocene crown clades"),
             size = 2, aes(color = "lightblue")) +
  theme_minimal() + 
  scale_color_manual(name = "", values = c("lightblue","lightcoral"),labels = c("Eocene crown clades","Eocene archaic taxa")) +
  geom_line(data = dplyr::filter(predframe.Neo_Br, Group == "Eocene archaic taxa"), color = "black",
            linetype = "dashed") +
  geom_line(data = dplyr::filter(predframe.Neo_Br, Group == "Eocene crown clades"), color = "black",
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
  labs(x = "log10(Brain surface)", y = "log10(Neocortex surface)") +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12,face = "bold")) 

##### ANOVA

sink("Neo_EV_ANOVA_MiddleEocene_regressions.txt")
fitOLS<-lm.rrpp(Neo~BrainS + Group, data=Trogosus.data)
summary(fitOLS) # significant 0.001 ***
anova(fitOLS)
coef(fitOLS, test=TRUE)
PT1 <- pairwise(fitOLS, groups = Trogosus.data$Group)
summary(PT1, confidence = 0.95)
sink()
  
## END!

