
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

#Open dataset and delete NA values
data1<-read.csv("Trogosus_dataset_final_Eocene.csv", header=T, row.names = 1)
NEO1<-data1[c(1:9,18,19)]
Trogosus.data<-na.omit(NEO1)

#Save dataset
write.csv(Trogosus.data,'Neo_EV_Trogosus_Eocene.csv')

#################  Make regression + residuals for Middle Eocene taxa ############

#Open dataset
Trogosus.data<-read.csv("Neo_EV_Trogosus_Eocene.csv", header=T, row.names = 1)
Trogosus.data<-subset(Trogosus.data, Trogosus.data[,5] =="middle Eocene")

#Transform data to log10
Trogosus.data$Neocortex_surface_mm2<-log10(Trogosus.data$Neocortex_surface_mm2)
names(Trogosus.data)[names(Trogosus.data) == "Neocortex_surface_mm2"] <- "Neo"
Trogosus.data$Brain_surface_mm2<-log10(Trogosus.data$Brain_surface_mm2)
names(Trogosus.data)[names(Trogosus.data) == "Brain_surface_mm2"] <- "BrainS"

#GLS to obtain Neo-EQ if need and residuals
GLS_Neo_Br<-gls(Neo ~ BrainS, data=Trogosus.data)
summary(GLS_Neo_Br) #Final model

#values obtained from model above
a <- 0.7756580 #slope
b <- 0.1095542 #intercept

Exp_b<-exp(b)
Exp_b

#Open dataset without log10 -- dataset with just Neo and BrainS
Trogosus.data<-read.csv("Neo_EV_Trogosus_Eocene.csv", header=T, row.names = 1)
Trogosus.data<-subset(Trogosus.data, Trogosus.data[,5] =="middle Eocene")

BM <-Trogosus.data$Brain_surface_mm2
Ec <- Exp_b*(BM)^a #Expected brain size
Ei <-Trogosus.data$Neocortex_surface_mm2
PEQ <-Ei/Ec

#Save PEQ as column in dataset
PEQ_df<-tibble(PEQ)
colnames(PEQ_df)<-c("PEQ")
Trogosus.data$PEQ=PEQ

#### Now do the predicted and residuals
Trogosus.data$predicted <- predict(GLS_Neo_Br)
Trogosus.data$residuals <- residuals(GLS_Neo_Br)

#export dataframe
write.csv(Trogosus.data,'Neo_EV_Trogosus_Middle_Eocene.csv')

############### Individual Archaic groups from the Middle Eocene ############

#Open data with residuals
Trogosus.data<-read.csv("Neo_EV_Trogosus_Middle_Eocene.csv", header=T, row.names = 1)
Trogosus.data<-subset(Trogosus.data, Trogosus.data[,3] =="Eocene archaic taxa")

##ggplot - boxplot - residuals
ggboxplot(Trogosus.data,x="Order2", y="residuals", fill="Order2", 
          palette=c("lavender","lavender","lavender","lavender","lavender",
                    "lavender"),
          order=c("Trogosus hillsii","Hyopsodus","Microsyops annectens",
                  "Creodonta","Carnivoramorpha","Metacheiromys marshi"),
          xlab =FALSE,legend.title = "")+
  geom_point() + 
  geom_text(aes(label = Abbreviations), hjust = 0, nudge_x = 0.05)+
  font("x.text", size = 8)+
  rotate_x_text(90)+
  theme(legend.position = "right")+
  labs(x='Clades', y='Residuals Neocortical surface vs Endocranial surface')

############# 3.Individual Crown groups from the Middle Eocene ################

#Open data with residuals
Trogosus.data<-read.csv("Neo_EV_Trogosus_Middle_Eocene.csv", header=T, row.names = 1)
Trogosus.data<-subset(Trogosus.data, Trogosus.data[,8] !="Eocene archaic taxa")

##ggplot - boxplot - residuals
ggboxplot(Trogosus.data,x="Order2", y="residuals", fill="Order2", 
          palette=c("lavender","lavender","lavender","lavender","lavender"),
          order=c("Trogosus hillsii","Perissodactyla","Artiodactyla","Rodentia",
                  "Primates"),
          xlab =FALSE,legend.title = "")+
  geom_point() + 
  geom_text(aes(label = Abbreviations), hjust = 0, nudge_x = 0.05)+
  font("x.text", size = 8)+
  rotate_x_text(90)+
  theme(legend.position = "right")+
  labs(x='Clades', y='Residuals Neocortical surface vs Endocranial surface')

#Only select categories that have more than 2 taxa to run test
Trogosus.data<-subset(Trogosus.data, Order2 =="Rodentia"| Order2 =="Primates")

sink("Neo_EV_ANOVA_MiddleEocene_Crown.txt")
shapiro.test(Trogosus.data$residuals) #p_value < 0.05 No normally distributed - levene
#leveneTest(residuals ~ Order2, data = Trogosus.data) #p-value < 0.05 No equality of variances - welsh
bartlett.test(residuals ~ Order2, data = Trogosus.data) 
oneway_test(residuals ~ Order2, data = Trogosus.data)
Trogosus.data$Order2 = factor(Trogosus.data$Order2, levels = c("Rodentia","Primates")) 
PT_test<-pairwisePermutationTest(residuals ~ Order2, data = Trogosus.data,method="fdr")
PT_test
#Welsh.test<-welch.test(residuals~Order2,data=Trogosus.data) #p-value<0.05 - yes statistical diff found
#paircomp(Welsh.test) #p-value<0.05 - yes statistical diff
sink()

#END! :)


