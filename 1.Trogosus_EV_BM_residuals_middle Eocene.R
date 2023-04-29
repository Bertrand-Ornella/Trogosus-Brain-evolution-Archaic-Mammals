
#script for the residuals Brain volume vs. Body mass - Trogosus paper

### Used functions
library(nlme) # GLS analysis
library(tibble) #make dataframe
library(ggpubr) #ggboxplot
library(car) #Levene's test
library(coin) #one_way test
library(rcompanion) #pairwisePermutationTest
library(onewaytests) #Welsh test

setwd("~/Documents/Trogosus")

#################  Make regression + residuals for Middle Eocene taxa  ############

#Open dataset
Trogosus.data<-read.csv("Trogosus_dataset_final_Eocene.csv", header=T, row.names = 1)
Trogosus.data<-subset(Trogosus.data, Trogosus.data[,5] =="middle Eocene")

#Transform data to log10
Trogosus.data$Brain_volume_cm3<-log10(Trogosus.data$Brain_volume_cm3)
names(Trogosus.data)[names(Trogosus.data) == "Brain_volume_cm3"] <- "Brain"
Trogosus.data$Body_mass_g<-log10(Trogosus.data$Body_mass_g)
names(Trogosus.data)[names(Trogosus.data) == "Body_mass_g"] <- "Body"

#GLS to obtain EQ and residuals
GLS_Br_Bo<-gls(Brain ~ Body, data=Trogosus.data)
summary(GLS_Br_Bo) #Final model

#values obtained from model above
a <- 0.5898303 #slope
b <--0.9575396 #intercept

Exp_b<-exp(b)
Exp_b

#Open dataset without log10
Trogosus.data<-read.csv("Trogosus_dataset_final_Eocene.csv", header=T, row.names = 1)
Trogosus.data<-subset(Trogosus.data, Trogosus.data[,5] =="middle Eocene")

BM <-Trogosus.data$Body_mass_g
Ec <- Exp_b*(BM)^a #Expected brain size
Ei <-Trogosus.data$Brain_volume_cm3
PEQ <-Ei/Ec

#Save PEQ as column in dataset
PEQ_df<-tibble(PEQ)
colnames(PEQ_df)<-c("PEQ")
Trogosus.data$PEQ=PEQ

#### Now do the predicted and residuals
Trogosus.data$predicted <- predict(GLS_Br_Bo)
Trogosus.data$residuals <- residuals(GLS_Br_Bo)

#export dataframe
write.csv(Trogosus.data,'EV_BM_Trogosus_Middle_Eocene.csv')

############### Individual Archaic groups from the Middle Eocene ############

#Open data with PEQ & residuals
Trogosus.data<-read.csv("EV_BM_Trogosus_Middle_Eocene.csv", header=T, row.names = 1)
Trogosus.data<-subset(Trogosus.data, Trogosus.data[,3] =="Eocene archaic taxa")

##ggplot - boxplot - residuals
ggboxplot(Trogosus.data,x="Order2", y="residuals", fill="Order2", 
          palette=c("lavender","lavender","lavender","lavender","lavender",
                    "lavender","lavender","lavender"),
          order=c("Trogosus hillsii","Hyopsodus","Microsyops annectens",
                  "Dinocerata","Mesonyx obtusidens","Creodonta","Carnivoramorpha",
                  "Metacheiromys marshi"),
          xlab =FALSE,legend.title = "")+
  geom_point() + 
  geom_text(aes(label = Abbreviations), hjust = 0, nudge_x = 0.05)+
  font("x.text", size = 8)+
  rotate_x_text(90)+
  theme(legend.position = "right")+
  labs(x='Clades', y='Residuals Endocranial volume vs. Body mass')

#Only select categories that have more than 2 taxa to run test
Trogosus.data<-subset(Trogosus.data, Order2 =="Trogosus hillsii"|Order2 =="Hyopsodus"|
        Order2 =="Dinocerata"|Order2 =="Creodonta"|Order2 =="Carnivoramorpha")

sink("EV_BM_ANOVA_MiddleEocene_Archaic.txt")
shapiro.test(Trogosus.data$residuals) #p_value < 0.05 No normally distributed - levene
#leveneTest(residuals ~ Order2, data = Trogosus.data) #p-value < 0.05 No equality of variances - welsh
bartlett.test(residuals ~ Order2, data = Trogosus.data) 
oneway_test(residuals ~ Order2, data = Trogosus.data)
Trogosus.data$Order2 = factor(Trogosus.data$Order2, 
                              levels = c("Trogosus hillsii","Hyopsodus","Dinocerata","Creodonta","Carnivoramorpha")) 
PT_test<-pairwisePermutationTest(residuals ~ Order2, data = Trogosus.data,method="fdr")
PT_test
#Welsh.test<-welch.test(residuals~Order2,data=Trogosus.data) #p-value<0.05 - yes statistical diff found
##paircomp(Welsh.test) #p-value<0.05 - yes statistical diff
sink()

############# Individual Crown groups from the Middle Eocene ################

#Open data with PEQ & residuals
Trogosus.data<-read.csv("EV_BM_Trogosus_Middle_Eocene.csv", header=T, row.names = 1)
Trogosus.data<-subset(Trogosus.data, Trogosus.data[,8] !="Eocene archaic taxa")

##ggplot - boxplot - residuals
ggboxplot(Trogosus.data,x="Order2", y="residuals", fill="Order2", 
          palette=c("lavender","lavender","lavender","lavender","lavender",
                    "lavender","lavender"),
          order=c("Trogosus hillsii","Perissodactyla","Artiodactyla","Rodentia",
                  "Moeritherium lyonsi","Primates","Archaeoceti"),
          xlab =FALSE,legend.title = "")+
  geom_point() + 
  geom_text(aes(label = Abbreviations), hjust = 0, nudge_x = 0.05)+
  font("x.text", size = 8)+
  rotate_x_text(90)+
  theme(legend.position = "right")+
  labs(x='Clades', y='Residuals Brain volume vs. Body mass')

#Only select categories that have more than 2 taxa to run test
Trogosus.data<-subset(Trogosus.data, Order2 =="Trogosus hillsii"|Order2 =="Perissodactyla"|Order2 =="Artiodactyla"|Order2 =="Rodentia"|
                        Order2 =="Primates"| Order2 =="Archaeoceti")

sink("EV_BM_ANOVA_MiddleEocene_Crown.txt")
shapiro.test(Trogosus.data$residuals) #p_value < 0.05 No normally distributed - levene
leveneTest(residuals ~ Order2, data = Trogosus.data) #p-value < 0.05 No equality of variances - welsh
#bartlett.test(residuals ~ Order2, data = Trogosus.data) 
#oneway_test(residuals ~ Order2, data = Trogosus.data)
#Trogosus.data$Order2 = factor(Trogosus.data$Order2, levels = c("Trogosus hillsii","Perissodactyla","Artiodactyla","Rodentia","Primates","Archaeoceti")) 
#PT_test<-pairwisePermutationTest(residuals ~ Order2, data = Trogosus.data,method="fdr")
#PT_test
Welsh.test<-welch.test(residuals~Order2,data=Trogosus.data) #p-value<0.05 - yes statistical diff found
paircomp(Welsh.test) #p-value<0.05 - yes statistical diff
sink()

## END!


