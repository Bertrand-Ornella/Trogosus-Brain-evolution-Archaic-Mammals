
#script for the residuals OB volume vs. Endocranial volume - Trogosus paper

### Used functions
library(nlme) # GLS analysis
library(tibble) #make dataframe
library(ggpubr) #ggboxplot
library(car) #Levene's test
library(coin) #one_way test
library(rcompanion) #pairwisePermutationTest
library(onewaytests) #Welsh test

setwd("~/Documents/Trogosus")

############  Make regression + residuals for Middle Eocene taxa DIET  ############

#Open dataset
Trogosus.data<-read.csv("Diet_OB_BM_EV_Trogosus_Eocene.csv", header=T, row.names = 1)
Trogosus.data<-subset(Trogosus.data, Trogosus.data[,5] =="middle Eocene")

#Transform data to log10
Trogosus.data$Olfactory_bulbs_mm3<-log10(Trogosus.data$Olfactory_bulbs_mm3)
names(Trogosus.data)[names(Trogosus.data) == "Olfactory_bulbs_mm3"] <- "OB"
Trogosus.data$Brain_volume_mm3<-log10(Trogosus.data$Brain_volume_mm3)
names(Trogosus.data)[names(Trogosus.data) == "Brain_volume_mm3"] <- "Brain"

#GLS to obtain OB - EQ if need and residuals
GLS_Ob_Br<-gls(OB ~ Brain, data=Trogosus.data)
summary(GLS_Ob_Br) #Final model

#values obtained from model above
a <- 0.9870827 #slope
b <--1.1358720 #intercept

Exp_b<-exp(b)
Exp_b

#Open dataset without log10 -- dataset with just OB and Brain
Trogosus.data<-read.csv("Diet_OB_BM_EV_Trogosus_Eocene.csv", header=T, row.names = 1)
Trogosus.data<-subset(Trogosus.data, Trogosus.data[,5] =="middle Eocene")

BR <-Trogosus.data$Brain_volume_mm3
Ec <- Exp_b*(BR)^a #Expected OB size
Ei <-Trogosus.data$Olfactory_bulbs_mm3
PEQ <-Ei/Ec

#Save PEQ as column in dataset
PEQ_df<-tibble(PEQ)
colnames(PEQ_df)<-c("PEQ")
Trogosus.data$PEQ=PEQ

#### Now do the predicted and residuals
Trogosus.data$predicted <- predict(GLS_Ob_Br)
Trogosus.data$residuals <- residuals(GLS_Ob_Br)

#export dataframe
write.csv(Trogosus.data,'Diet_OB_EV_Trogosus_Middle_Eocene.csv')

############# Diet archaic vs crown groups from the Middle Eocene #################

#Open data residuals
Trogosus.data<-read.csv("Diet_OB_EV_Trogosus_Middle_Eocene.csv", header=T, row.names = 1)

####
ggplot(Trogosus.data, aes(x = Diet, y = residuals, color = Group)) + 
  geom_boxplot() +
  geom_point() +
  scale_color_manual(values=c("lightblue", "lightcoral")) +
  #scale_fill_manual(values=c("lightblue", "lightcoral")) +
  #geom_text(aes(label = Abbreviations), hjust = 0, nudge_x = 0.05)+
  facet_wrap(~Group)

###### ARCHAIC -- Open data residuals
Trogosus.data<-read.csv("Diet_OB_EV_Trogosus_Middle_Eocene.csv", header=T, row.names = 1)
Trogosus.data<-subset(Trogosus.data, Trogosus.data[,3] =="Eocene archaic taxa")

#Test normality and homogeneity of the variance
#Not normal but variances are equal - pairwisePermutation test ANOVA
sink("Diet_OB_EV_ANOVA_MiddleEocene_Archaic.txt")
shapiro.test(Trogosus.data$residuals) #p_value < 0.05 No normally distributed - levene
#leveneTest(residuals ~ Diet, data = Trogosus.data)
bartlett.test(residuals ~ Diet, data = Trogosus.data) #p-value < 0.05 No equality of variances - welsh
oneway_test(residuals ~ Diet, data = Trogosus.data)
Trogosus.data$Diet = factor(Trogosus.data$Diet, levels = c("Herbivorous","Omnivorous-Carnivorous")) 
PT_test<-pairwisePermutationTest(residuals ~ Diet, data = Trogosus.data,method="fdr")
PT_test
#Welsh.test<-welch.test(residuals~Diet,data=Trogosus.data) #p-value<0.05 - yes statistical diff found
#paircomp(Welsh.test) #p-value<0.05 - yes statistical diff
sink()

####### HERBIVOROUS -- Open data residuals
Trogosus.data<-read.csv("Diet_OB_EV_Trogosus_Middle_Eocene.csv", header=T, row.names = 1)
Trogosus.data<-subset(Trogosus.data, Trogosus.data[,10] !="Omnivorous-Carnivorous")

#Test normality and homogeneity of the variance
#Not normal but variances are equal - pairwisePermutation test ANOVA
sink("Diet_OB_EV_ANOVA_MiddleEocene_Herbivorous.txt")
shapiro.test(Trogosus.data$residuals) #p_value < 0.05 No normally distributed - levene
bartlett.test(residuals ~ Group, data = Trogosus.data) #p-value < 0.05 No equality of variances - welsh
oneway_test(residuals ~ Group, data = Trogosus.data)
Trogosus.data$Group = factor(Trogosus.data$Group, levels = c("Eocene archaic taxa","Eocene crown clades")) 
PT_test<-pairwisePermutationTest(residuals ~ Group, data = Trogosus.data,method="fdr")
PT_test
#Welsh.test<-welch.test(residuals~Group2,data=Trogosus.data) #p-value<0.05 - yes statistical diff found
#paircomp(Welsh.test) #p-value<0.05 - yes statistical diff
sink()

###### Crown herbivores vs. Archaic carnivorans -- Open data residuals
Trogosus.data<-read.csv("Diet_OB_EV_Trogosus_Middle_Eocene.csv", header=T, row.names = 1)
Trogosus.data<-subset(Trogosus.data, Order =="Perissodactyla" | Order =="Artiodactyla" 
                      |Order == "Carnivoramorpha" | Order == "Creodonta")

#Test normality and homogeneity of the variance
#Not normal but variances are equal - pairwisePermutation test ANOVA
sink("Diet_OB_EV_ANOVA_MiddleEocene_CrownHerb_ArchaicCarni.txt")
shapiro.test(Trogosus.data$residuals) #p_value < 0.05 No normally distributed - levene
bartlett.test(residuals ~ Diet, data = Trogosus.data) #p-value < 0.05 No equality of variances - welsh
oneway_test(residuals ~ Diet, data = Trogosus.data)
Trogosus.data$Diet = factor(Trogosus.data$Diet, levels = c("Herbivorous","Omnivorous-Carnivorous")) 
PT_test<-pairwisePermutationTest(residuals ~ Diet, data = Trogosus.data,method="fdr")
PT_test
#Welsh.test<-welch.test(residuals~Group2,data=Trogosus.data) #p-value<0.05 - yes statistical diff found
#paircomp(Welsh.test) #p-value<0.05 - yes statistical diff
sink()

# END!



