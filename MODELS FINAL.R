setwd("~/Desktop/Research/Thesis/ThesisAnalysis")

library(e1071)
library(foreign)
library(MASS)
library(Hmisc)
library(reshape2)
library(pracma)
library(ggplot2) ##H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016
library(MuMIn) #Kamil Barton (2019) #. MuMIn: Multi-Model Inference. R package version 1.43.6.
library(mice) #Stef van Buuren, Karin Groothuis-Oudshoorn (2011). mice: Multivariate Imputation by Chained Equations in.R. Journal of Statistical Software, 45(3), 1-67. URL https://www.jstatsoft.org/v45/i03/.
library(ggpubr) #Alboukadel Kassambara (2018). ggpubr: 'ggplot2' Based Publication Ready Plots. 
library(rstan) #Stan Development Team (2018). RStan: the R interface to Stan.
library(rstanarm) #Goodrich B, Gabry J, Ali I & Brilleman S. (2018). rstanarm: Bayesian applied regression modeling via
library(bayesplot) #Jonah Gabry and Tristan Mahr (2018). bayesplot: Plotting for Bayesian Models. 
library(tidyverse) #Hadley Wickham (2017). tidyverse: Easily Install and Load the 'Tidyverse'. 
library(broom) #David Robinson and Alex Hayes (2019). broom: Convert Statistical Analysis Objects into Tidy Tibbles.
library(coda) #Martyn Plummer, Nicky Best, Kate Cowles and Karen Vines (2006). CODA: Convergence Diagnosis and Output.Analysis for MCMC, R News, vol 6, 7-11
library(data.table) #Matt Dowle and Arun Srinivasan (2019). data.table: Extension of `data.frame`
library(lme4) #Douglas Bates, Martin Maechler, Ben Bolker, Steve Walker (2015) #. Fitting Linear Mixed-Effects Models.Using lme4. Journal of Statistical Software, 67(1), 1-48. doi:10.18637/jss.v067.i01.
library(brms) # Paul-Christian Bürkner (2017). brms: An R Package for Bayesian Multilevel Models Using Stan. Journal of StatisticalSoftware, 80(1), 1-28. doi:10.18637/jss.v080.i01
library(DHARMa)#Florian Hartig (2019). DHARMa: Residual Diagnostics for Hierarchical (Multi-Level / Mixed) Regression Models. R package version 0.2.4. https://CRAN.R-project.org/package=DHARMa
library(plyr)#Hadley Wickham (2011). The Split-Apply-Combine Strategy for Data Analysis. Journal of Statistical Software, 40(1), 1-29. URL http://www.jstatsoft.org/v40/i01/.
library(networkD3) # J.J. Allaire, Christopher Gandrud, Kenton Russell and CJ Yetman (2017). networkD3: D3 JavaScript Network Graphs from R. R package version 0.4. https://CRAN.R-project.org/package=networkD3
library(lme4) 
library(car)
library(ordinal)
library(lmtest)
library(forcats)


##Upload raw survey data
rawsurvey <- read.csv("RawSurvey_noNames.csv")


##Upload seafood perceptions
seafood <- read.csv("CleanedSurvey.csv", na.strings = "")


##Convert seafood responses to integers
seafood[seafood == "Very high"] <- "5"
seafood[seafood == "High"] <- "4"
seafood[seafood == "Medium"] <- "3"
seafood[seafood == "Low"] <- "2"
seafood[seafood == "Very low"] <- "1"
seafood[seafood == "None at all"] <- "0"




seafood <- seafood %>% mutate(agegroup = case_when(Age >= 60   ~ '3',
                                                   Age >= 40  & Age <= 59 ~ '2',
                                                   Age >= 18  & Age <= 39 ~ '1'))

seafood$MarineOccupation = ifelse(seafood$Fisher == 'Yes'|seafood$Gleaner == 'Yes', 'Fisher/Gleaner', 'Other')

seafood <- seafood %>% 
  mutate_at(c(18:301), as.numeric)

cols.num <- c("Age","Fisher_start", "Fisher_before_start", "Fisher_before_end", "Gleaner_start", "Gleaner_before_start", "Gleaner_before_end", "year_move")
seafood[cols.num] <- sapply(seafood[cols.num],as.numeric)


seafood[seafood == "9-Ifaty"] <- "Ifaty"
seafood[seafood == "10-Madiorano"] <- "Madiorano"
seafood[seafood == "1-Ambalaboy"] <- "Ambalaboy"

##Time Spent Fishing
seafood$start <- with(seafood,pmin(Fisher_start,Gleaner_start,Fisher_before_start, Gleaner_before_start, na.rm=TRUE))
seafood$end <- with(seafood,pmax(Fisher_before_end,Gleaner_before_end, na.rm=TRUE))
seafood$endreal <- ifelse(seafood$MarineOccupation == 'Fisher/Gleaner', '2022', NA)
seafood$endreal[is.na(seafood$endreal)] <- seafood$end[is.na(seafood$endreal)] 
seafood$zero <- "0"

seafood$YearsFishing <- as.numeric(seafood$endreal) - as.numeric(seafood$start)
seafood$YearsFishing[is.na(seafood$YearsFishing)] <- seafood$zero[is.na(seafood$YearsFishing)] 
seafood$YearsFishing <- sapply(seafood$YearsFishing,as.numeric)

##Time spent in Village 
seafood$VillageYears <- 2022 - as.numeric(seafood$year_move)
seafood$VillageYears[is.na(seafood$VillageYears)] <- seafood$Age[is.na(seafood$VillageYears)] 



write.csv(seafood, "seafoodnumbers.csv", row.names=FALSE)




seafood <- read.csv("seafoodnumbers.csv")



#Functions
#to draw histogram on pairs plot
panel_hist = function(x, ...)
{
  usr = par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h = hist(x, plot = FALSE)
  breaks = h$breaks; nB = length(breaks)
  y = h$counts; y = y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "grey", ...)
}
#Standardize function for continuous variables
standardize = function(x){(x-mean(x, na.rm=T))/(2*sd(x, na.rm=T))} 

#variance inflation factors  function for mixed models
vif.mer = function (fit) {
  ## adapted from rms::vif
  v <-vcov(fit)
  nam = names(fixef(fit))
  ## exclude intercepts
  ns = sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
  if (ns > 0) {
    v = v[-(1:ns), -(1:ns), drop = FALSE]
    nam = nam[-(1:ns)]
  }
  d = diag(v)^0.5
  v = diag(solve(v/(d %o% d)))
  names(v) = nam
  v
}

#correlation function
panel.cor = function(x, y, digits=2, prefix="", cex.cor, ...)
{
  usr = par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r = abs(cor(x, y, method = "pearson",use = "complete.obs"))
  txt = format(c(r, 0.123456789), digits=digits)[1]
  txt = paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex.cor = 0.9/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r*2)
}
#function to extract the probabilities from the predicted distribution
quantInv = function(distr, value) ecdf(distr)(value)

#Explanatory variables: relevel categorical variables and standardize continuous variables


seafood$Gender=relevel(factor(seafood$Gender),ref="Male")
seafood$Village=relevel(factor(seafood$Village),ref="Ifaty")
seafood$MarineOccupation=relevel(factor(seafood$MarineOccupation),ref="Fisher/Gleaner")
seafood$Age_standard=standardize(nthroot(seafood$Age, 2))
skewness(seafood$Age, na.rm=TRUE)
skewness(seafood$Age_standard, na.rm=TRUE)
seafood$YearsFishing_standard=standardize(nthroot(seafood$YearsFishing, 2))
skewness(seafood$YearsFishing, na.rm=TRUE)
skewness(seafood$YearsFishing_standard, na.rm=TRUE)
seafood$VillageYears_standard=standardize(seafood$VillageYears)
skewness(seafood$VillageYears, na.rm=TRUE)
skewness(seafood$VillageYears_standard, na.rm=TRUE)


##Reduced skew (normalized and standardized) age and years of fishing
ggplot(seafood, aes(x=Age_standard, y=..density..)) + 
  geom_histogram(bins = 7, color = "black", fill = "gray") +
  geom_density(aes(y = ..density..))

ggplot(seafood, aes(x=YearsFishing_standard, y=..density..)) + 
  geom_histogram(bins = 7, color = "black", fill = "gray") +
  geom_density(aes(y = ..density..))

ggplot(seafood, aes(x=VillageYears_standard, y=..density..)) + 
  geom_histogram(bins = 7, color = "black", fill = "gray") +
  geom_density(aes(y = ..density..))


##Check covariates
covariates <- c("Gender", "Village", "MarineOccupation", "Age_standard", "YearsFishing_standard", "VillageYears_standard")

pMiss <- function(x){sum(is.na(x))/length(x)*100}

covariatescolumns <- seafood[, covariates]

apply(covariatescolumns,2,pMiss)
mice::md.pattern(covariatescolumns)
VIM::aggr(covariatescolumns, col=c('navyblue','red'), numbers=TRUE, sortVars=TRUE, labels=names(data), cex.axis=.7, gap=3, ylab=c("Histogram of missing data","Pattern"))


#Pairs correlations among covariates

pairs(~Gender+
        Village+
        ##MarineOccupation+
        Age_standard+
        VillageYears_standard+
        YearsFishing_standard,
      data=seafood,lower.panel=panel.cor )


#variance inflation factor

model <- brm(Mullidae_cons_change ~ 
               Gender+
               VillageYears_standard+
               YearsFishing_standard+
               Age_standard +(1|Village),
             data=seafood[!is.na(seafood$Mullidae_cons_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))


VIF.table=as.data.frame(vif(model))


##Ordinal Logistic Regression 


##Mullidae
seafood$Mullidae_cons_change <- seafood$Mullidae_cons_now - seafood$Mullidae_cons_10
seafood$Mullidae_ea_change <- seafood$Mullidae_ea_now - seafood$Mullidae_ea_10
seafood$Mullidae_harv_change <- seafood$Mullidae_harv_now - seafood$Mullidae_harv_10

seafood$Mullidae_cons_change <- ordered(as.factor(seafood$Mullidae_cons_change), 
                                        levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

seafood$Mullidae_ea_change <- ordered(as.factor(seafood$Mullidae_ea_change), 
                                      levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

seafood$Mullidae_harv_change <- ordered(as.factor(seafood$Mullidae_harv_change), 
                                        levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

##Mullidae Cons
nullMullidaecons<-brm(Mullidae_cons_change ~ 1 +(1|Village),
                      data=seafood[!is.na(seafood$Mullidae_cons_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testMullidaecons <-brm(Mullidae_cons_change ~ 
                         Gender+
                         VillageYears_standard+
                         YearsFishing_standard+
                         Age_standard +(1|Village),
                       data=seafood[!is.na(seafood$Mullidae_cons_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

mullidaecons <- pp_check(testMullidaecons,ndraws = 50)+ggtitle("Model fit Mullidae Consumption")

loo(nullMullidaecons, testMullidaecons)

mullidaeconseffects<-as.data.frame(fixef(testMullidaecons, probs = c(0.1,0.9)))
mullidaeconseffects$community<-row.names(mullidaeconseffects)
colnames(mullidaeconseffects)<-c("estimate","se","conf.low","conf.high","community")
mullidaeconseff <- ggplot(mullidaeconseffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Mullidae Cons")




##Mullidae ea
nullMullidaeea<-brm(Mullidae_ea_change ~ 1 +(1|Village),
                    data=seafood[!is.na(seafood$Mullidae_ea_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testMullidaeea <-brm(Mullidae_ea_change ~ 
                       Gender+
                       VillageYears_standard+
                       YearsFishing_standard+
                       Age_standard +(1|Village),
                     data=seafood[!is.na(seafood$Mullidae_ea_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

mullidaeea <- pp_check(testMullidaeea,ndraws = 50)+ggtitle("Model fit Mullidae Abundance")

loo(nullMullidaeea, testMullidaeea)

mullidaeeaeffects<-as.data.frame(fixef(testMullidaeea, probs = c(0.1,0.9)))
mullidaeeaeffects$community<-row.names(mullidaeeaeffects)
colnames(mullidaeeaeffects)<-c("estimate","se","conf.low","conf.high","community")
mullidaeeaeff <- ggplot(mullidaeeaeffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Mullidae ea")


##Mullidae harv
nullMullidaeharv<-brm(Mullidae_harv_change ~ 1 +(1|Village),
                      data=seafood[!is.na(seafood$Mullidae_harv_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testMullidaeharv <-brm(Mullidae_harv_change ~ 
                         Gender+
                         VillageYears_standard+
                         YearsFishing_standard+
                         Age_standard +(1|Village),
                       data=seafood[!is.na(seafood$Mullidae_harv_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

mullidaeharv <- pp_check(testMullidaeharv,ndraws = 50)+ggtitle("Model fit Mullidae Harvest")

loo(nullMullidaeharv, testMullidaeharv)

mullidaeharveffects<-as.data.frame(fixef(testMullidaeharv, probs = c(0.1,0.9)))
mullidaeharveffects$community<-row.names(mullidaeharveffects)
colnames(mullidaeharveffects)<-c("estimate","se","conf.low","conf.high","community")
mullidaeharveff <- ggplot(mullidaeharveffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Mullidae harv")


##Export graphs
jpeg("PPcheck_Mullidae", width =3500, height = 2000, res=300)
ggarrange(mullidaecons, mullidaeea, mullidaeharv)
dev.off()

jpeg("FixedEffects_Mullidae", width =3500, height = 2000, res=300)
ggarrange(mullidaeconseff, mullidaeeaeff, mullidaeharveff)
dev.off()





##Pomacentridae
seafood$Pomacentridae_cons_change <- seafood$Pomacentridae_cons_now - seafood$Pomacentridae_cons_10
seafood$Pomacentridae_ea_change <- seafood$Pomacentridae_ea_now - seafood$Pomacentridae_ea_10
seafood$Pomacentridae_harv_change <- seafood$Pomacentridae_harv_now - seafood$Pomacentridae_harv_10

seafood$Pomacentridae_cons_change <- ordered(as.factor(seafood$Pomacentridae_cons_change), 
                                             levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

seafood$Pomacentridae_ea_change <- ordered(as.factor(seafood$Pomacentridae_ea_change), 
                                           levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

seafood$Pomacentridae_harv_change <- ordered(as.factor(seafood$Pomacentridae_harv_change), 
                                             levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))


##Pomacentridae Cons
nullPomacentridaecons<-brm(Pomacentridae_cons_change ~ 1 +(1|Village),
                           data=seafood[!is.na(seafood$Pomacentridae_cons_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testPomacentridaecons <-brm(Pomacentridae_cons_change ~ 
                              Gender+
                              VillageYears_standard+
                              YearsFishing_standard+
                              Age_standard +(1|Village),
                            data=seafood[!is.na(seafood$Pomacentridae_cons_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Pomacentridaecons <- pp_check(testPomacentridaecons,ndraws = 50)+ggtitle("Model fit Pomacentridae Consumption")

loo(nullPomacentridaecons, testPomacentridaecons)

Pomacentridaeconseffects<-as.data.frame(fixef(testPomacentridaecons, probs = c(0.1,0.9)))
Pomacentridaeconseffects$community<-row.names(Pomacentridaeconseffects)
colnames(Pomacentridaeconseffects)<-c("estimate","se","conf.low","conf.high","community")
Pomacentridaeconseff <- ggplot(Pomacentridaeconseffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Pomacentridae Cons")




##Pomacentridae ea
nullPomacentridaeea<-brm(Pomacentridae_ea_change ~ 1 +(1|Village),
                         data=seafood[!is.na(seafood$Pomacentridae_ea_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testPomacentridaeea <-brm(Pomacentridae_ea_change ~ 
                            Gender+
                            VillageYears_standard+
                            YearsFishing_standard+
                            Age_standard +(1|Village),
                          data=seafood[!is.na(seafood$Pomacentridae_ea_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Pomacentridaeea <- pp_check(testPomacentridaeea,ndraws = 50)+ggtitle("Model fit Pomacentridae Abundance")

loo(nullPomacentridaeea, testPomacentridaeea)

Pomacentridaeeaeffects<-as.data.frame(fixef(testPomacentridaeea, probs = c(0.1,0.9)))
Pomacentridaeeaeffects$community<-row.names(Pomacentridaeeaeffects)
colnames(Pomacentridaeeaeffects)<-c("estimate","se","conf.low","conf.high","community")
Pomacentridaeeaeff <- ggplot(Pomacentridaeeaeffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Pomacentridae ea")


##Pomacentridae harv
nullPomacentridaeharv<-brm(Pomacentridae_harv_change ~ 1 +(1|Village),
                           data=seafood[!is.na(seafood$Pomacentridae_harv_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testPomacentridaeharv <-brm(Pomacentridae_harv_change ~ 
                              Gender+
                              VillageYears_standard+
                              YearsFishing_standard+
                              Age_standard +(1|Village),
                            data=seafood[!is.na(seafood$Pomacentridae_harv_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Pomacentridaeharv <- pp_check(testPomacentridaeharv,ndraws = 50)+ggtitle("Model fit Pomacentridae Harvest")

loo(nullPomacentridaeharv, testPomacentridaeharv)

Pomacentridaeharveffects<-as.data.frame(fixef(testPomacentridaeharv, probs = c(0.1,0.9)))
Pomacentridaeharveffects$community<-row.names(Pomacentridaeharveffects)
colnames(Pomacentridaeharveffects)<-c("estimate","se","conf.low","conf.high","community")
Pomacentridaeharveff <- ggplot(Pomacentridaeharveffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Pomacentridae harv")


##Export graphs
jpeg("PPcheck_Pomacentridae", width =3500, height = 2000, res=300)
ggarrange(Pomacentridaecons, Pomacentridaeea, Pomacentridaeharv)
dev.off()

jpeg("FixedEffects_Pomacentridae", width =3500, height = 2000, res=300)
ggarrange(Pomacentridaeconseff, Pomacentridaeeaeff, Pomacentridaeharveff)
dev.off()






##Acanthuridae
seafood$Acanthuridae_cons_change <- seafood$Acanthuridae_cons_now - seafood$Acanthuridae_cons_10
seafood$Acanthuridae_ea_change <- seafood$Acanthuridae_ea_now - seafood$Acanthuridae_ea_10
seafood$Acanthuridae_harv_change <- seafood$Acanthuridae_harv_now - seafood$Acanthuridae_harv_10

seafood$Acanthuridae_cons_change <- ordered(as.factor(seafood$Acanthuridae_cons_change), 
                                            levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

seafood$Acanthuridae_ea_change <- ordered(as.factor(seafood$Acanthuridae_ea_change), 
                                          levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

seafood$Acanthuridae_harv_change <- ordered(as.factor(seafood$Acanthuridae_harv_change), 
                                            levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

##Acanthuridae Cons
nullAcanthuridaecons<-brm(Acanthuridae_cons_change ~ 1 +(1|Village),
                          data=seafood[!is.na(seafood$Acanthuridae_cons_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testAcanthuridaecons <-brm(Acanthuridae_cons_change ~ 
                             Gender+
                             VillageYears_standard+
                             YearsFishing_standard+
                             Age_standard +(1|Village),
                           data=seafood[!is.na(seafood$Acanthuridae_cons_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Acanthuridaecons <- pp_check(testAcanthuridaecons,ndraws = 50)+ggtitle("Model fit Acanthuridae Consumption")

loo(nullAcanthuridaecons, testAcanthuridaecons)

Acanthuridaeconseffects<-as.data.frame(fixef(testAcanthuridaecons, probs = c(0.1,0.9)))
Acanthuridaeconseffects$community<-row.names(Acanthuridaeconseffects)
colnames(Acanthuridaeconseffects)<-c("estimate","se","conf.low","conf.high","community")
Acanthuridaeconseff <- ggplot(Acanthuridaeconseffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Acanthuridae Cons")

##Acanthuridae ea
nullAcanthuridaeea<-brm(Acanthuridae_ea_change ~ 1 +(1|Village),
                        data=seafood[!is.na(seafood$Acanthuridae_ea_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testAcanthuridaeea <-brm(Acanthuridae_ea_change ~ 
                           Gender+
                           VillageYears_standard+
                           YearsFishing_standard+
                           Age_standard +(1|Village),
                         data=seafood[!is.na(seafood$Acanthuridae_ea_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Acanthuridaeea <- pp_check(testAcanthuridaeea,ndraws = 50)+ggtitle("Model fit Acanthuridae Abundance")

loo(nullAcanthuridaeea, testAcanthuridaeea)

Acanthuridaeeaeffects<-as.data.frame(fixef(testAcanthuridaeea, probs = c(0.1,0.9)))
Acanthuridaeeaeffects$community<-row.names(Acanthuridaeeaeffects)
colnames(Acanthuridaeeaeffects)<-c("estimate","se","conf.low","conf.high","community")
Acanthuridaeeaeff <- ggplot(Acanthuridaeeaeffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Acanthuridae ea")

##Acanthuridae harv
nullAcanthuridaeharv<-brm(Acanthuridae_harv_change ~ 1 +(1|Village),
                          data=seafood[!is.na(seafood$Acanthuridae_harv_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testAcanthuridaeharv <-brm(Acanthuridae_harv_change ~ 
                             Gender+
                             VillageYears_standard+
                             YearsFishing_standard+
                             Age_standard +(1|Village),
                           data=seafood[!is.na(seafood$Acanthuridae_harv_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Acanthuridaeharv <- pp_check(testAcanthuridaeharv,ndraws = 50)+ggtitle("Model fit Acanthuridae Harvest")

loo(nullAcanthuridaeharv, testAcanthuridaeharv)

Acanthuridaeharveffects<-as.data.frame(fixef(testAcanthuridaeharv, probs = c(0.1,0.9)))
Acanthuridaeharveffects$community<-row.names(Acanthuridaeharveffects)
colnames(Acanthuridaeharveffects)<-c("estimate","se","conf.low","conf.high","community")
Acanthuridaeharveff <- ggplot(Acanthuridaeharveffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Acanthuridae harv")

##Export graphs
jpeg("PPcheck_Acanthuridae", width =3500, height = 2000, res=300)
ggarrange(Acanthuridaecons, Acanthuridaeea, Acanthuridaeharv)
dev.off()

jpeg("FixedEffects_Acanthuridae", width =3500, height = 2000, res=300)
ggarrange(Acanthuridaeconseff, Acanthuridaeeaeff, Acanthuridaeharveff)
dev.off()







##Scaridae
seafood$Scaridae_cons_change <- seafood$Scaridae_cons_now - seafood$Scaridae_cons_10
seafood$Scaridae_ea_change <- seafood$Scaridae_ea_now - seafood$Scaridae_ea_10
seafood$Scaridae_harv_change <- seafood$Scaridae_harv_now - seafood$Scaridae_harv_10

seafood$Scaridae_cons_change <- ordered(as.factor(seafood$Scaridae_cons_change), 
                                        levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

seafood$Scaridae_ea_change <- ordered(as.factor(seafood$Scaridae_ea_change), 
                                      levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

seafood$Scaridae_harv_change <- ordered(as.factor(seafood$Scaridae_harv_change), 
                                        levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))


##Scaridae Cons
nullScaridaecons<-brm(Scaridae_cons_change ~ 1 +(1|Village),
                      data=seafood[!is.na(seafood$Scaridae_cons_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testScaridaecons <-brm(Scaridae_cons_change ~ 
                         Gender+
                         VillageYears_standard+
                         YearsFishing_standard+
                         Age_standard +(1|Village),
                       data=seafood[!is.na(seafood$Scaridae_cons_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Scaridaecons <- pp_check(testScaridaecons,ndraws = 50)+ggtitle("Model fit Scaridae Consumption")

loo(nullScaridaecons, testScaridaecons)

Scaridaeconseffects<-as.data.frame(fixef(testScaridaecons, probs = c(0.1,0.9)))
Scaridaeconseffects$community<-row.names(Scaridaeconseffects)
colnames(Scaridaeconseffects)<-c("estimate","se","conf.low","conf.high","community")
Scaridaeconseff <- ggplot(Scaridaeconseffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Scaridae Cons")

##Scaridae ea
nullScaridaeea<-brm(Scaridae_ea_change ~ 1 +(1|Village),
                    data=seafood[!is.na(seafood$Scaridae_ea_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testScaridaeea <-brm(Scaridae_ea_change ~ 
                       Gender+
                       VillageYears_standard+
                       YearsFishing_standard+
                       Age_standard +(1|Village),
                     data=seafood[!is.na(seafood$Scaridae_ea_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Scaridaeea <- pp_check(testScaridaeea,ndraws = 50)+ggtitle("Model fit Scaridae Abundance")

loo(nullScaridaeea, testScaridaeea)

Scaridaeeaeffects<-as.data.frame(fixef(testScaridaeea, probs = c(0.1,0.9)))
Scaridaeeaeffects$community<-row.names(Scaridaeeaeffects)
colnames(Scaridaeeaeffects)<-c("estimate","se","conf.low","conf.high","community")
Scaridaeeaeff <- ggplot(Scaridaeeaeffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Scaridae ea")


##Scaridae harv
nullScaridaeharv<-brm(Scaridae_harv_change ~ 1 +(1|Village),
                      data=seafood[!is.na(seafood$Scaridae_harv_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testScaridaeharv <-brm(Scaridae_harv_change ~ 
                         Gender+
                         VillageYears_standard+
                         YearsFishing_standard+
                         Age_standard +(1|Village),
                       data=seafood[!is.na(seafood$Scaridae_harv_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Scaridaeharv <- pp_check(testScaridaeharv,ndraws = 50)+ggtitle("Model fit Scaridae Harvest")

loo(nullScaridaeharv, testScaridaeharv)

Scaridaeharveffects<-as.data.frame(fixef(testScaridaeharv, probs = c(0.1,0.9)))
Scaridaeharveffects$community<-row.names(Scaridaeharveffects)
colnames(Scaridaeharveffects)<-c("estimate","se","conf.low","conf.high","community")
Scaridaeharveff <- ggplot(Scaridaeharveffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Scaridae harv")


##Export graphs
jpeg("PPcheck_Scaridae", width =3500, height = 2000, res=300)
ggarrange(Scaridaecons, Scaridaeea, Scaridaeharv)
dev.off()

jpeg("FixedEffects_Scaridae", width =3500, height = 2000, res=300)
ggarrange(Scaridaeconseff, Scaridaeeaeff, Scaridaeharveff)
dev.off()








##Scorpaenidae
seafood$Scorpaenidae_cons_change <- seafood$Scorpaenidae_cons_now - seafood$Scorpaenidae_cons_10
seafood$Scorpaenidae_ea_change <- seafood$Scorpaenidae_ea_now - seafood$Scorpaenidae_ea_10
seafood$Scorpaenidae_harv_change <- seafood$Scorpaenidae_harv_now - seafood$Scorpaenidae_harv_10

seafood$Scorpaenidae_cons_change <- ordered(as.factor(seafood$Scorpaenidae_cons_change), 
                                            levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

seafood$Scorpaenidae_ea_change <- ordered(as.factor(seafood$Scorpaenidae_ea_change), 
                                          levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

seafood$Scorpaenidae_harv_change <- ordered(as.factor(seafood$Scorpaenidae_harv_change), 
                                            levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

##Scorpaenidae Cons
nullScorpaenidaecons<-brm(Scorpaenidae_cons_change ~ 1 +(1|Village),
                          data=seafood[!is.na(seafood$Scorpaenidae_cons_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testScorpaenidaecons <-brm(Scorpaenidae_cons_change ~ 
                             Gender+
                             VillageYears_standard+
                             YearsFishing_standard+
                             Age_standard +(1|Village),
                           data=seafood[!is.na(seafood$Scorpaenidae_cons_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Scorpaenidaecons <- pp_check(testScorpaenidaecons,ndraws = 50)+ggtitle("Model fit Scorpaenidae Consumption")

loo(nullScorpaenidaecons, testScorpaenidaecons)

Scorpaenidaeconseffects<-as.data.frame(fixef(testScorpaenidaecons, probs = c(0.1,0.9)))
Scorpaenidaeconseffects$community<-row.names(Scorpaenidaeconseffects)
colnames(Scorpaenidaeconseffects)<-c("estimate","se","conf.low","conf.high","community")
Scorpaenidaeconseff <- ggplot(Scorpaenidaeconseffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Scorpaenidae Cons")

##Scorpaenidae ea
nullScorpaenidaeea<-brm(Scorpaenidae_ea_change ~ 1 +(1|Village),
                        data=seafood[!is.na(seafood$Scorpaenidae_ea_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testScorpaenidaeea <-brm(Scorpaenidae_ea_change ~ 
                           Gender+
                           VillageYears_standard+
                           YearsFishing_standard+
                           Age_standard +(1|Village),
                         data=seafood[!is.na(seafood$Scorpaenidae_ea_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Scorpaenidaeea <- pp_check(testScorpaenidaeea,ndraws = 50)+ggtitle("Model fit Scorpaenidae Abundance")

loo(nullScorpaenidaeea, testScorpaenidaeea)

Scorpaenidaeeaeffects<-as.data.frame(fixef(testScorpaenidaeea, probs = c(0.1,0.9)))
Scorpaenidaeeaeffects$community<-row.names(Scorpaenidaeeaeffects)
colnames(Scorpaenidaeeaeffects)<-c("estimate","se","conf.low","conf.high","community")
Scorpaenidaeeaeff <- ggplot(Scorpaenidaeeaeffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Scorpaenidae ea")

##Scorpaenidae harv
nullScorpaenidaeharv<-brm(Scorpaenidae_harv_change ~ 1 +(1|Village),
                          data=seafood[!is.na(seafood$Scorpaenidae_harv_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testScorpaenidaeharv <-brm(Scorpaenidae_harv_change ~ 
                             Gender+
                             VillageYears_standard+
                             YearsFishing_standard+
                             Age_standard +(1|Village),
                           data=seafood[!is.na(seafood$Scorpaenidae_harv_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Scorpaenidaeharv <- pp_check(testScorpaenidaeharv,ndraws = 50)+ggtitle("Model fit Scorpaenidae Harvest")

loo(nullScorpaenidaeharv, testScorpaenidaeharv)

Scorpaenidaeharveffects<-as.data.frame(fixef(testScorpaenidaeharv, probs = c(0.1,0.9)))
Scorpaenidaeharveffects$community<-row.names(Scorpaenidaeharveffects)
colnames(Scorpaenidaeharveffects)<-c("estimate","se","conf.low","conf.high","community")
Scorpaenidaeharveff <- ggplot(Scorpaenidaeharveffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Scorpaenidae harv")

##Export graphs
jpeg("PPcheck_Scorpaenidae", width =3500, height = 2000, res=300)
ggarrange(Scorpaenidaecons, Scorpaenidaeea, Scorpaenidaeharv)
dev.off()

jpeg("FixedEffects_Scorpaenidae", width =3500, height = 2000, res=300)
ggarrange(Scorpaenidaeconseff, Scorpaenidaeeaeff, Scorpaenidaeharveff)
dev.off()






##Siganidae
seafood$Siganidae_cons_change <- seafood$Siganidae_cons_now - seafood$Siganidae_cons_10
seafood$Siganidae_ea_change <- seafood$Siganidae_ea_now - seafood$Siganidae_ea_10
seafood$Siganidae_harv_change <- seafood$Siganidae_harv_now - seafood$Siganidae_harv_10

seafood$Siganidae_cons_change <- ordered(as.factor(seafood$Siganidae_cons_change), 
                                         levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

seafood$Siganidae_ea_change <- ordered(as.factor(seafood$Siganidae_ea_change), 
                                       levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

seafood$Siganidae_harv_change <- ordered(as.factor(seafood$Siganidae_harv_change), 
                                         levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

##Siganidae Cons
nullSiganidaecons<-brm(Siganidae_cons_change ~ 1 +(1|Village),
                       data=seafood[!is.na(seafood$Siganidae_cons_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testSiganidaecons <-brm(Siganidae_cons_change ~ 
                          Gender+
                          VillageYears_standard+
                          YearsFishing_standard+
                          Age_standard +(1|Village),
                        data=seafood[!is.na(seafood$Siganidae_cons_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Siganidaecons <- pp_check(testSiganidaecons,ndraws = 50)+ggtitle("Model fit Siganidae Consumption")

loo(nullSiganidaecons, testSiganidaecons)

Siganidaeconseffects<-as.data.frame(fixef(testSiganidaecons, probs = c(0.1,0.9)))
Siganidaeconseffects$community<-row.names(Siganidaeconseffects)
colnames(Siganidaeconseffects)<-c("estimate","se","conf.low","conf.high","community")
Siganidaeconseff <- ggplot(Siganidaeconseffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Siganidae Cons")

##Siganidae ea
nullSiganidaeea<-brm(Siganidae_ea_change ~ 1 +(1|Village),
                     data=seafood[!is.na(seafood$Siganidae_ea_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testSiganidaeea <-brm(Siganidae_ea_change ~ 
                        Gender+
                        VillageYears_standard+
                        YearsFishing_standard+
                        Age_standard +(1|Village),
                      data=seafood[!is.na(seafood$Siganidae_ea_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Siganidaeea <- pp_check(testSiganidaeea,ndraws = 50)+ggtitle("Model fit Siganidae Abundance")

loo(nullSiganidaeea, testSiganidaeea)

Siganidaeeaeffects<-as.data.frame(fixef(testSiganidaeea, probs = c(0.1,0.9)))
Siganidaeeaeffects$community<-row.names(Siganidaeeaeffects)
colnames(Siganidaeeaeffects)<-c("estimate","se","conf.low","conf.high","community")
Siganidaeeaeff <- ggplot(Siganidaeeaeffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Siganidae ea")

##Siganidae harv
nullSiganidaeharv<-brm(Siganidae_harv_change ~ 1 +(1|Village),
                       data=seafood[!is.na(seafood$Siganidae_harv_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testSiganidaeharv <-brm(Siganidae_harv_change ~ 
                          Gender+
                          VillageYears_standard+
                          YearsFishing_standard+
                          Age_standard +(1|Village),
                        data=seafood[!is.na(seafood$Siganidae_harv_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Siganidaeharv <- pp_check(testSiganidaeharv,ndraws = 50)+ggtitle("Model fit Siganidae Harvest")

loo(nullSiganidaeharv, testSiganidaeharv)

Siganidaeharveffects<-as.data.frame(fixef(testSiganidaeharv, probs = c(0.1,0.9)))
Siganidaeharveffects$community<-row.names(Siganidaeharveffects)
colnames(Siganidaeharveffects)<-c("estimate","se","conf.low","conf.high","community")
Siganidaeharveff <- ggplot(Siganidaeharveffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Siganidae harv")

##Export graphs
jpeg("PPcheck_Siganidae", width =3500, height = 2000, res=300)
ggarrange(Siganidaecons, Siganidaeea, Siganidaeharv)
dev.off()

jpeg("FixedEffects_Siganidae", width =3500, height = 2000, res=300)
ggarrange(Siganidaeconseff, Siganidaeeaeff, Siganidaeharveff)
dev.off()





##Sphyraenidae
seafood$Sphyraenidae_cons_change <- seafood$Sphyraenidae_cons_now - seafood$Sphyraenidae_cons_10
seafood$Sphyraenidae_ea_change <- seafood$Sphyraenidae_ea_now - seafood$Sphyraenidae_ea_10
seafood$Sphyraenidae_harv_change <- seafood$Sphyraenidae_harv_now - seafood$Sphyraenidae_harv_10

seafood$Sphyraenidae_cons_change <- ordered(as.factor(seafood$Sphyraenidae_cons_change), 
                                            levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

seafood$Sphyraenidae_ea_change <- ordered(as.factor(seafood$Sphyraenidae_ea_change), 
                                          levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

seafood$Sphyraenidae_harv_change <- ordered(as.factor(seafood$Sphyraenidae_harv_change), 
                                            levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

##Sphyraenidae Cons
nullSphyraenidaecons<-brm(Sphyraenidae_cons_change ~ 1 +(1|Village),
                          data=seafood[!is.na(seafood$Sphyraenidae_cons_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testSphyraenidaecons <-brm(Sphyraenidae_cons_change ~ 
                             Gender+
                             VillageYears_standard+
                             YearsFishing_standard+
                             Age_standard +(1|Village),
                           data=seafood[!is.na(seafood$Sphyraenidae_cons_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Sphyraenidaecons <- pp_check(testSphyraenidaecons,ndraws = 50)+ggtitle("Model fit Sphyraenidae Consumption")

loo(nullSphyraenidaecons, testSphyraenidaecons)

Sphyraenidaeconseffects<-as.data.frame(fixef(testSphyraenidaecons, probs = c(0.1,0.9)))
Sphyraenidaeconseffects$community<-row.names(Sphyraenidaeconseffects)
colnames(Sphyraenidaeconseffects)<-c("estimate","se","conf.low","conf.high","community")
Sphyraenidaeconseff <- ggplot(Sphyraenidaeconseffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Sphyraenidae Cons")

##Sphyraenidae ea
nullSphyraenidaeea<-brm(Sphyraenidae_ea_change ~ 1 +(1|Village),
                        data=seafood[!is.na(seafood$Sphyraenidae_ea_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testSphyraenidaeea <-brm(Sphyraenidae_ea_change ~ 
                           Gender+
                           VillageYears_standard+
                           YearsFishing_standard+
                           Age_standard +(1|Village),
                         data=seafood[!is.na(seafood$Sphyraenidae_ea_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Sphyraenidaeea <- pp_check(testSphyraenidaeea,ndraws = 50)+ggtitle("Model fit Sphyraenidae Abundance")

loo(nullSphyraenidaeea, testSphyraenidaeea)

Sphyraenidaeeaeffects<-as.data.frame(fixef(testSphyraenidaeea, probs = c(0.1,0.9)))
Sphyraenidaeeaeffects$community<-row.names(Sphyraenidaeeaeffects)
colnames(Sphyraenidaeeaeffects)<-c("estimate","se","conf.low","conf.high","community")
Sphyraenidaeeaeff <- ggplot(Sphyraenidaeeaeffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Sphyraenidae ea")

##Sphyraenidae harv
nullSphyraenidaeharv<-brm(Sphyraenidae_harv_change ~ 1 +(1|Village),
                          data=seafood[!is.na(seafood$Sphyraenidae_harv_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testSphyraenidaeharv <-brm(Sphyraenidae_harv_change ~ 
                             Gender+
                             VillageYears_standard+
                             YearsFishing_standard+
                             Age_standard +(1|Village),
                           data=seafood[!is.na(seafood$Sphyraenidae_harv_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Sphyraenidaeharv <- pp_check(testSphyraenidaeharv,ndraws = 50)+ggtitle("Model fit Sphyraenidae Harvest")

loo(nullSphyraenidaeharv, testSphyraenidaeharv)

Sphyraenidaeharveffects<-as.data.frame(fixef(testSphyraenidaeharv, probs = c(0.1,0.9)))
Sphyraenidaeharveffects$community<-row.names(Sphyraenidaeharveffects)
colnames(Sphyraenidaeharveffects)<-c("estimate","se","conf.low","conf.high","community")
Sphyraenidaeharveff <- ggplot(Sphyraenidaeharveffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Sphyraenidae harv")

##Export graphs
jpeg("PPcheck_Sphyraenidae", width =3500, height = 2000, res=300)
ggarrange(Sphyraenidaecons, Sphyraenidaeea, Sphyraenidaeharv)
dev.off()

jpeg("FixedEffects_Sphyraenidae", width =3500, height = 2000, res=300)
ggarrange(Sphyraenidaeconseff, Sphyraenidaeconseff, Sphyraenidaeconseff)
dev.off()





##Labridae
seafood$Labridae_cons_change <- seafood$Labridae_cons_now - seafood$Labridae_cons_10
seafood$Labridae_ea_change <- seafood$Labridae_ea_now - seafood$Labridae_ea_10
seafood$Labridae_harv_change <- seafood$Labridae_harv_now - seafood$Labridae_harv_10

seafood$Labridae_cons_change <- ordered(as.factor(seafood$Labridae_cons_change), 
                                        levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

seafood$Labridae_ea_change <- ordered(as.factor(seafood$Labridae_ea_change), 
                                      levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

seafood$Labridae_harv_change <- ordered(as.factor(seafood$Labridae_harv_change), 
                                        levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

##Labridae Cons
nullLabridaecons<-brm(Labridae_cons_change ~ 1 +(1|Village),
                      data=seafood[!is.na(seafood$Labridae_cons_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testLabridaecons <-brm(Labridae_cons_change ~ 
                         Gender+
                         VillageYears_standard+
                         YearsFishing_standard+
                         Age_standard +(1|Village),
                       data=seafood[!is.na(seafood$Labridae_cons_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Labridaecons <- pp_check(testLabridaecons,ndraws = 50)+ggtitle("Model fit Labridae Consumption")

loo(nullLabridaecons, testLabridaecons)

Labridaeconseffects<-as.data.frame(fixef(testLabridaecons, probs = c(0.1,0.9)))
Labridaeconseffects$community<-row.names(Labridaeconseffects)
colnames(Labridaeconseffects)<-c("estimate","se","conf.low","conf.high","community")
Labridaeconseff <- ggplot(Labridaeconseffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Labridae Cons")

##Labridae ea
nullLabridaeea<-brm(Labridae_ea_change ~ 1 +(1|Village),
                    data=seafood[!is.na(seafood$Labridae_ea_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testLabridaeea <-brm(Labridae_ea_change ~ 
                       Gender+
                       VillageYears_standard+
                       YearsFishing_standard+
                       Age_standard +(1|Village),
                     data=seafood[!is.na(seafood$Labridae_ea_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Labridaeea <- pp_check(testLabridaeea,ndraws = 50)+ggtitle("Model fit Labridae Abundance")

loo(nullLabridaeea, testLabridaeea)

Labridaeeaeffects<-as.data.frame(fixef(testLabridaeea, probs = c(0.1,0.9)))
Labridaeeaeffects$community<-row.names(Labridaeeaeffects)
colnames(Labridaeeaeffects)<-c("estimate","se","conf.low","conf.high","community")
Labridaeeaeff <- ggplot(Labridaeeaeffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Labridae ea")

##Labridae harv
nullLabridaeharv<-brm(Labridae_harv_change ~ 1 +(1|Village),
                      data=seafood[!is.na(seafood$Labridae_harv_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testLabridaeharv <-brm(Labridae_harv_change ~ 
                         Gender+
                         VillageYears_standard+
                         YearsFishing_standard+
                         Age_standard +(1|Village),
                       data=seafood[!is.na(seafood$Labridae_harv_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Labridaeharv <- pp_check(testLabridaeharv,ndraws = 50)+ggtitle("Model fit Labridae Harvest")

loo(nullLabridaeharv, testLabridaeharv)

Labridaeharveffects<-as.data.frame(fixef(testLabridaeharv, probs = c(0.1,0.9)))
Labridaeharveffects$community<-row.names(Labridaeharveffects)
colnames(Labridaeharveffects)<-c("estimate","se","conf.low","conf.high","community")
Labridaeharveff <- ggplot(Labridaeharveffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Labridae harv")

##Export graphs
jpeg("PPcheck_Labridae", width =3500, height = 2000, res=300)
ggarrange(Labridaecons, Labridaeea, Labridaeharv)
dev.off()

jpeg("FixedEffects_Labridae", width =3500, height = 2000, res=300)
ggarrange(Labridaeconseff, Labridaeeaeff, Labridaeharveff)
dev.off()








##Priacanthidae
seafood$Priacanthidae_cons_change <- seafood$Priacanthidae_cons_now - seafood$Priacanthidae_cons_10
seafood$Priacanthidae_ea_change <- seafood$Priacanthidae_ea_now - seafood$Priacanthidae_ea_10
seafood$Priacanthidae_harv_change <- seafood$Priacanthidae_harv_now - seafood$Priacanthidae_harv_10

seafood$Priacanthidae_cons_change <- ordered(as.factor(seafood$Priacanthidae_cons_change), 
                                             levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

seafood$Priacanthidae_ea_change <- ordered(as.factor(seafood$Priacanthidae_ea_change), 
                                           levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

seafood$Priacanthidae_harv_change <- ordered(as.factor(seafood$Priacanthidae_harv_change), 
                                             levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

##Priacanthidae Cons
nullPriacanthidaecons<-brm(Priacanthidae_cons_change ~ 1 +(1|Village),
                           data=seafood[!is.na(seafood$Priacanthidae_cons_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testPriacanthidaecons <-brm(Priacanthidae_cons_change ~ 
                              Gender+
                              VillageYears_standard+
                              YearsFishing_standard+
                              Age_standard +(1|Village),
                            data=seafood[!is.na(seafood$Priacanthidae_cons_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Priacanthidaecons <- pp_check(testPriacanthidaecons,ndraws = 50)+ggtitle("Model fit Priacanthidae Consumption")

loo(nullPriacanthidaecons, testPriacanthidaecons)

Priacanthidaeconseffects<-as.data.frame(fixef(testPriacanthidaecons, probs = c(0.1,0.9)))
Priacanthidaeconseffects$community<-row.names(Priacanthidaeconseffects)
colnames(Priacanthidaeconseffects)<-c("estimate","se","conf.low","conf.high","community")
Priacanthidaeconseff <- ggplot(Priacanthidaeconseffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Priacanthidae Cons")

##Priacanthidae ea
nullPriacanthidaeea<-brm(Priacanthidae_ea_change ~ 1 +(1|Village),
                         data=seafood[!is.na(seafood$Priacanthidae_ea_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testPriacanthidaeea <-brm(Priacanthidae_ea_change ~ 
                            Gender+
                            VillageYears_standard+
                            YearsFishing_standard+
                            Age_standard +(1|Village),
                          data=seafood[!is.na(seafood$Priacanthidae_ea_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Priacanthidaeea <- pp_check(testPriacanthidaeea,ndraws = 50)+ggtitle("Model fit Priacanthidae Abundance")

loo(nullPriacanthidaeea, testPriacanthidaeea)

Priacanthidaeeaeffects<-as.data.frame(fixef(testPriacanthidaeea, probs = c(0.1,0.9)))
Priacanthidaeeaeffects$community<-row.names(Priacanthidaeeaeffects)
colnames(Priacanthidaeeaeffects)<-c("estimate","se","conf.low","conf.high","community")
Priacanthidaeeaeff <- ggplot(Priacanthidaeeaeffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Priacanthidae ea")

##Priacanthidae harv
nullPriacanthidaeharv<-brm(Priacanthidae_harv_change ~ 1 +(1|Village),
                           data=seafood[!is.na(seafood$Priacanthidae_harv_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testPriacanthidaeharv <-brm(Priacanthidae_harv_change ~ 
                              Gender+
                              VillageYears_standard+
                              YearsFishing_standard+
                              Age_standard +(1|Village),
                            data=seafood[!is.na(seafood$Priacanthidae_harv_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Priacanthidaeharv <- pp_check(testPriacanthidaeharv,ndraws = 50)+ggtitle("Model fit Priacanthidae Harvest")

loo(nullPriacanthidaeharv, testPriacanthidaeharv)

Priacanthidaeharveffects<-as.data.frame(fixef(testPriacanthidaeharv, probs = c(0.1,0.9)))
Priacanthidaeharveffects$community<-row.names(Priacanthidaeharveffects)
colnames(Priacanthidaeharveffects)<-c("estimate","se","conf.low","conf.high","community")
Priacanthidaeharveff <- ggplot(Priacanthidaeharveffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Priacanthidae harv")

##Export graphs
jpeg("PPcheck_Priacanthidae", width =3500, height = 2000, res=300)
ggarrange(Priacanthidaecons, Priacanthidaeea, Priacanthidaeharv)
dev.off()

jpeg("FixedEffects_Priacanthidae", width =3500, height = 2000, res=300)
ggarrange(Priacanthidaeconseff, Priacanthidaeeaeff, Priacanthidaeharveff)
dev.off()






##Chaetodontidae
seafood$Chaetodontidae_cons_change <- seafood$Chaetodontidae_cons_now - seafood$Chaetodontidae_cons_10
seafood$Chaetodontidae_ea_change <- seafood$Chaetodontidae_ea_now - seafood$Chaetodontidae_ea_10
seafood$Chaetodontidae_harv_change <- seafood$Chaetodontidae_harv_now - seafood$Chaetodontidae_harv_10

seafood$Chaetodontidae_cons_change <- ordered(as.factor(seafood$Chaetodontidae_cons_change), 
                                              levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

seafood$Chaetodontidae_ea_change <- ordered(as.factor(seafood$Chaetodontidae_ea_change), 
                                            levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

seafood$Chaetodontidae_harv_change <- ordered(as.factor(seafood$Chaetodontidae_harv_change), 
                                              levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

##Chaetodontidae Cons
nullChaetodontidaecons<-brm(Chaetodontidae_cons_change ~ 1 +(1|Village),
                            data=seafood[!is.na(seafood$Chaetodontidae_cons_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testChaetodontidaecons <-brm(Chaetodontidae_cons_change ~ 
                               Gender+
                               VillageYears_standard+
                               YearsFishing_standard+
                               Age_standard +(1|Village),
                             data=seafood[!is.na(seafood$Chaetodontidae_cons_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Chaetodontidaecons <- pp_check(testChaetodontidaecons,ndraws = 50)+ggtitle("Model fit Chaetodontidae Consumption")

loo(nullChaetodontidaecons, testChaetodontidaecons)

Chaetodontidaeconseffects<-as.data.frame(fixef(testChaetodontidaecons, probs = c(0.1,0.9)))
Chaetodontidaeconseffects$community<-row.names(Chaetodontidaeconseffects)
colnames(Chaetodontidaeconseffects)<-c("estimate","se","conf.low","conf.high","community")
Chaetodontidaeconseff <- ggplot(Chaetodontidaeconseffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Chaetodontidae Cons")

##Chaetodontidae ea
nullChaetodontidaeea<-brm(Chaetodontidae_ea_change ~ 1 +(1|Village),
                          data=seafood[!is.na(seafood$Chaetodontidae_ea_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testChaetodontidaeea <-brm(Chaetodontidae_ea_change ~ 
                             Gender+
                             VillageYears_standard+
                             YearsFishing_standard+
                             Age_standard +(1|Village),
                           data=seafood[!is.na(seafood$Chaetodontidae_ea_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Chaetodontidaeea <- pp_check(testChaetodontidaeea,ndraws = 50)+ggtitle("Model fit Chaetodontidae Abundance")

loo(nullChaetodontidaeea, testChaetodontidaeea)

Chaetodontidaeeaeffects<-as.data.frame(fixef(testChaetodontidaeea, probs = c(0.1,0.9)))
Chaetodontidaeeaeffects$community<-row.names(Chaetodontidaeeaeffects)
colnames(Chaetodontidaeeaeffects)<-c("estimate","se","conf.low","conf.high","community")
Chaetodontidaeeaeff <- ggplot(Chaetodontidaeeaeffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Chaetodontidae ea")

##Chaetodontidae harv
nullChaetodontidaeharv<-brm(Chaetodontidae_harv_change ~ 1 +(1|Village),
                            data=seafood[!is.na(seafood$Chaetodontidae_harv_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testChaetodontidaeharv <-brm(Chaetodontidae_harv_change ~ 
                               Gender+
                               VillageYears_standard+
                               YearsFishing_standard+
                               Age_standard +(1|Village),
                             data=seafood[!is.na(seafood$Chaetodontidae_harv_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Chaetodontidaeharv <- pp_check(testChaetodontidaeharv,ndraws = 50)+ggtitle("Model fit Chaetodontidae Harvest")

loo(nullChaetodontidaeharv, testChaetodontidaeharv)

Chaetodontidaeharveffects<-as.data.frame(fixef(testChaetodontidaeharv, probs = c(0.1,0.9)))
Chaetodontidaeharveffects$community<-row.names(Chaetodontidaeharveffects)
colnames(Chaetodontidaeharveffects)<-c("estimate","se","conf.low","conf.high","community")
Chaetodontidaeharveff <- ggplot(Chaetodontidaeharveffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Chaetodontidae harv")

##Export graphs
jpeg("PPcheck_Chaetodontidae", width =3500, height = 2000, res=300)
ggarrange(Chaetodontidaecons, Chaetodontidaeea, Chaetodontidaeharv)
dev.off()

jpeg("FixedEffects_Chaetodontidae", width =3500, height = 2000, res=300)
ggarrange(Chaetodontidaeconseff, Chaetodontidaeeaeff, Chaetodontidaeharveff)
dev.off()






##Fistularidae
seafood$Fistularidae_cons_change <- seafood$Fistularidae_cons_now - seafood$Fistularidae_cons_10
seafood$Fistularidae_ea_change <- seafood$Fistularidae_ea_now - seafood$Fistularidae_ea_10
seafood$Fistularidae_harv_change <- seafood$Fistularidae_harv_now - seafood$Fistularidae_harv_10

seafood$Fistularidae_cons_change <- ordered(as.factor(seafood$Fistularidae_cons_change), 
                                            levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

seafood$Fistularidae_ea_change <- ordered(as.factor(seafood$Fistularidae_ea_change), 
                                          levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

seafood$Fistularidae_harv_change <- ordered(as.factor(seafood$Fistularidae_harv_change), 
                                            levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

##Fistularidae Cons
nullFistularidaecons<-brm(Fistularidae_cons_change ~ 1 +(1|Village),
                          data=seafood[!is.na(seafood$Fistularidae_cons_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testFistularidaecons <-brm(Fistularidae_cons_change ~ 
                             Gender+
                             VillageYears_standard+
                             YearsFishing_standard+
                             Age_standard +(1|Village),
                           data=seafood[!is.na(seafood$Fistularidae_cons_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Fistularidaecons <- pp_check(testFistularidaecons,ndraws = 50)+ggtitle("Model fit Fistularidae Consumption")

loo(nullFistularidaecons, testFistularidaecons)

Fistularidaeconseffects<-as.data.frame(fixef(testFistularidaecons, probs = c(0.1,0.9)))
Fistularidaeconseffects$community<-row.names(Fistularidaeconseffects)
colnames(Fistularidaeconseffects)<-c("estimate","se","conf.low","conf.high","community")
Fistularidaeconseff <- ggplot(Fistularidaeconseffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Fistularidae Cons")

##Fistularidae ea
nullFistularidaeea<-brm(Fistularidae_ea_change ~ 1 +(1|Village),
                        data=seafood[!is.na(seafood$Fistularidae_ea_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testFistularidaeea <-brm(Fistularidae_ea_change ~ 
                           Gender+
                           VillageYears_standard+
                           YearsFishing_standard+
                           Age_standard +(1|Village),
                         data=seafood[!is.na(seafood$Fistularidae_ea_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Fistularidaeea <- pp_check(testFistularidaeea,ndraws = 50)+ggtitle("Model fit Fistularidae Abundance")

loo(nullFistularidaeea, testFistularidaeea)

Fistularidaeeaeffects<-as.data.frame(fixef(testFistularidaeea, probs = c(0.1,0.9)))
Fistularidaeeaeffects$community<-row.names(Fistularidaeeaeffects)
colnames(Fistularidaeeaeffects)<-c("estimate","se","conf.low","conf.high","community")
Fistularidaeeaeff <- ggplot(Fistularidaeeaeffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Fistularidae ea")

##Fistularidae harv
nullFistularidaeharv<-brm(Fistularidae_harv_change ~ 1 +(1|Village),
                          data=seafood[!is.na(seafood$Fistularidae_harv_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testFistularidaeharv <-brm(Fistularidae_harv_change ~ 
                             Gender+
                             VillageYears_standard+
                             YearsFishing_standard+
                             Age_standard +(1|Village),
                           data=seafood[!is.na(seafood$Fistularidae_harv_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Fistularidaeharv <- pp_check(testFistularidaeharv,ndraws = 50)+ggtitle("Model fit Fistularidae Harvest")

loo(nullFistularidaeharv, testFistularidaeharv)

Fistularidaeharveffects<-as.data.frame(fixef(testFistularidaeharv, probs = c(0.1,0.9)))
Fistularidaeharveffects$community<-row.names(Fistularidaeharveffects)
colnames(Fistularidaeharveffects)<-c("estimate","se","conf.low","conf.high","community")
Fistularidaeharveff <- ggplot(Fistularidaeharveffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Fistularidae harv")

##Export graphs
jpeg("PPcheck_Fistularidae", width =3500, height = 2000, res=300)
ggarrange(Fistularidaecons, Fistularidaeea, Fistularidaeharv)
dev.off()

jpeg("FixedEffects_Fistularidae", width =3500, height = 2000, res=300)
ggarrange(Fistularidaeconseff, Fistularidaeeaeff, Fistularidaeharveff)
dev.off()






##Gerreidae
seafood$Gerreidae_cons_change <- seafood$Gerreidae_cons_now - seafood$Gerreidae_cons_10
seafood$Gerreidae_ea_change <- seafood$Gerreidae_ea_now - seafood$Gerreidae_ea_10
seafood$Gerreidae_harv_change <- seafood$Gerreidae_harv_now - seafood$Gerreidae_harv_10

seafood$Gerreidae_cons_change <- ordered(as.factor(seafood$Gerreidae_cons_change), 
                                         levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

seafood$Gerreidae_ea_change <- ordered(as.factor(seafood$Gerreidae_ea_change), 
                                       levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

seafood$Gerreidae_harv_change <- ordered(as.factor(seafood$Gerreidae_harv_change), 
                                         levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

##Gerreidae Cons
nullGerreidaecons<-brm(Gerreidae_cons_change ~ 1 +(1|Village),
                       data=seafood[!is.na(seafood$Gerreidae_cons_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testGerreidaecons <-brm(Gerreidae_cons_change ~ 
                          Gender+
                          VillageYears_standard+
                          YearsFishing_standard+
                          Age_standard +(1|Village),
                        data=seafood[!is.na(seafood$Gerreidae_cons_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Gerreidaecons <- pp_check(testGerreidaecons,ndraws = 50)+ggtitle("Model fit Gerreidae Consumption")

loo(nullGerreidaecons, testGerreidaecons)

Gerreidaeconseffects<-as.data.frame(fixef(testGerreidaecons, probs = c(0.1,0.9)))
Gerreidaeconseffects$community<-row.names(Gerreidaeconseffects)
colnames(Gerreidaeconseffects)<-c("estimate","se","conf.low","conf.high","community")
Gerreidaeconseff <- ggplot(Gerreidaeconseffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Gerreidae Cons")

##Gerreidae ea
nullGerreidaeea<-brm(Gerreidae_ea_change ~ 1 +(1|Village),
                     data=seafood[!is.na(seafood$Gerreidae_ea_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testGerreidaeea <-brm(Gerreidae_ea_change ~ 
                        Gender+
                        VillageYears_standard+
                        YearsFishing_standard+
                        Age_standard +(1|Village),
                      data=seafood[!is.na(seafood$Gerreidae_ea_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Gerreidaeea <- pp_check(testGerreidaeea,ndraws = 50)+ggtitle("Model fit Gerreidae Abundance")

loo(nullGerreidaeea, testGerreidaeea)

Gerreidaeeaeffects<-as.data.frame(fixef(testGerreidaeea, probs = c(0.1,0.9)))
Gerreidaeeaeffects$community<-row.names(Gerreidaeeaeffects)
colnames(Gerreidaeeaeffects)<-c("estimate","se","conf.low","conf.high","community")
Gerreidaeeaeff <- ggplot(Gerreidaeeaeffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Gerreidae ea")

##Gerreidae harv
nullGerreidaeharv<-brm(Gerreidae_harv_change ~ 1 +(1|Village),
                       data=seafood[!is.na(seafood$Gerreidae_harv_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testGerreidaeharv <-brm(Gerreidae_harv_change ~ 
                          Gender+
                          VillageYears_standard+
                          YearsFishing_standard+
                          Age_standard +(1|Village),
                        data=seafood[!is.na(seafood$Gerreidae_harv_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Gerreidaeharv <- pp_check(testGerreidaeharv,ndraws = 50)+ggtitle("Model fit Gerreidae Harvest")

loo(nullGerreidaeharv, testGerreidaeharv)

Gerreidaeharveffects<-as.data.frame(fixef(testGerreidaeharv, probs = c(0.1,0.9)))
Gerreidaeharveffects$community<-row.names(Gerreidaeharveffects)
colnames(Gerreidaeharveffects)<-c("estimate","se","conf.low","conf.high","community")
Gerreidaeharveff <- ggplot(Gerreidaeharveffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Gerreidae harv")

##Export graphs
jpeg("PPcheck_Gerreidae", width =3500, height = 2000, res=300)
ggarrange(Gerreidaecons, Gerreidaeea, Gerreidaeharv)
dev.off()

jpeg("FixedEffects_Gerreidae", width =3500, height = 2000, res=300)
ggarrange(Gerreidaeconseff, Gerreidaeeaeff, Gerreidaeharveff)
dev.off()






##Haemulidae
seafood$Haemulidae_cons_change <- seafood$Haemulidae_cons_now - seafood$Haemulidae_cons_10
seafood$Haemulidae_ea_change <- seafood$Haemulidae_ea_now - seafood$Haemulidae_ea_10
seafood$Haemulidae_harv_change <- seafood$Haemulidae_harv_now - seafood$Haemulidae_harv_10

seafood$Haemulidae_cons_change <- ordered(as.factor(seafood$Haemulidae_cons_change), 
                                          levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

seafood$Haemulidae_ea_change <- ordered(as.factor(seafood$Haemulidae_ea_change), 
                                        levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

seafood$Haemulidae_harv_change <- ordered(as.factor(seafood$Haemulidae_harv_change), 
                                          levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

##Haemulidae Cons
nullHaemulidaecons<-brm(Haemulidae_cons_change ~ 1 +(1|Village),
                        data=seafood[!is.na(seafood$Haemulidae_cons_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testHaemulidaecons <-brm(Haemulidae_cons_change ~ 
                           Gender+
                           VillageYears_standard+
                           YearsFishing_standard+
                           Age_standard +(1|Village),
                         data=seafood[!is.na(seafood$Haemulidae_cons_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Haemulidaecons <- pp_check(testHaemulidaecons,ndraws = 50)+ggtitle("Model fit Haemulidae Consumption")

loo(nullHaemulidaecons, testHaemulidaecons)

Haemulidaeconseffects<-as.data.frame(fixef(testHaemulidaecons, probs = c(0.1,0.9)))
Haemulidaeconseffects$community<-row.names(Haemulidaeconseffects)
colnames(Haemulidaeconseffects)<-c("estimate","se","conf.low","conf.high","community")
Haemulidaeconseff <- ggplot(Haemulidaeconseffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Haemulidae Cons")

##Haemulidae ea
nullHaemulidaeea<-brm(Haemulidae_ea_change ~ 1 +(1|Village),
                      data=seafood[!is.na(seafood$Haemulidae_ea_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testHaemulidaeea <-brm(Haemulidae_ea_change ~ 
                         Gender+
                         VillageYears_standard+
                         YearsFishing_standard+
                         Age_standard +(1|Village),
                       data=seafood[!is.na(seafood$Haemulidae_ea_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Haemulidaeea <- pp_check(testHaemulidaeea,ndraws = 50)+ggtitle("Model fit Haemulidae Abundance")

loo(nullHaemulidaeea, testHaemulidaeea)

Haemulidaeeaeffects<-as.data.frame(fixef(testHaemulidaeea, probs = c(0.1,0.9)))
Haemulidaeeaeffects$community<-row.names(Haemulidaeeaeffects)
colnames(Haemulidaeeaeffects)<-c("estimate","se","conf.low","conf.high","community")
Haemulidaeeaeff <- ggplot(Haemulidaeeaeffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Haemulidae ea")

##Haemulidae harv
nullHaemulidaeharv<-brm(Haemulidae_harv_change ~ 1 +(1|Village),
                        data=seafood[!is.na(seafood$Haemulidae_harv_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testHaemulidaeharv <-brm(Haemulidae_harv_change ~ 
                           Gender+
                           VillageYears_standard+
                           YearsFishing_standard+
                           Age_standard +(1|Village),
                         data=seafood[!is.na(seafood$Haemulidae_harv_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Haemulidaeharv <- pp_check(testHaemulidaeharv,ndraws = 50)+ggtitle("Model fit Haemulidae Harvest")

loo(nullHaemulidaeharv, testHaemulidaeharv)

Haemulidaeharveffects<-as.data.frame(fixef(testHaemulidaeharv, probs = c(0.1,0.9)))
Haemulidaeharveffects$community<-row.names(Haemulidaeharveffects)
colnames(Haemulidaeharveffects)<-c("estimate","se","conf.low","conf.high","community")
Haemulidaeharveff <- ggplot(Haemulidaeharveffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Haemulidae harv")

##Export graphs
jpeg("PPcheck_Haemulidae", width =3500, height = 2000, res=300)
ggarrange(Haemulidaecons, Haemulidaeea, Haemulidaeharv)
dev.off()

jpeg("FixedEffects_Haemulidae", width =3500, height = 2000, res=300)
ggarrange(Haemulidaeconseff, Haemulidaeeaeff, Haemulidaeharveff)
dev.off()





##Lethrinidae
seafood$Lethrinidae_cons_change <- seafood$Lethrinidae_cons_now - seafood$Lethrinidae_cons_10
seafood$Lethrinidae_ea_change <- seafood$Lethrinidae_ea_now - seafood$Lethrinidae_ea_10
seafood$Lethrinidae_harv_change <- seafood$Lethrinidae_harv_now - seafood$Lethrinidae_harv_10

seafood$Lethrinidae_cons_change <- ordered(as.factor(seafood$Lethrinidae_cons_change), 
                                           levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

seafood$Lethrinidae_ea_change <- ordered(as.factor(seafood$Lethrinidae_ea_change), 
                                         levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

seafood$Lethrinidae_harv_change <- ordered(as.factor(seafood$Lethrinidae_harv_change), 
                                           levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

##Lethrinidae Cons
nullLethrinidaecons<-brm(Lethrinidae_cons_change ~ 1 +(1|Village),
                         data=seafood[!is.na(seafood$Lethrinidae_cons_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testLethrinidaecons <-brm(Lethrinidae_cons_change ~ 
                            Gender+
                            VillageYears_standard+
                            YearsFishing_standard+
                            Age_standard +(1|Village),
                          data=seafood[!is.na(seafood$Lethrinidae_cons_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Lethrinidaecons <- pp_check(testLethrinidaecons,ndraws = 50)+ggtitle("Model fit Lethrinidae Consumption")

loo(nullLethrinidaecons, testLethrinidaecons)

Lethrinidaeconseffects<-as.data.frame(fixef(testLethrinidaecons, probs = c(0.1,0.9)))
Lethrinidaeconseffects$community<-row.names(Lethrinidaeconseffects)
colnames(Lethrinidaeconseffects)<-c("estimate","se","conf.low","conf.high","community")
Lethrinidaeconseff <- ggplot(Lethrinidaeconseffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Lethrinidae Cons")

##Lethrinidae ea
nullLethrinidaeea<-brm(Lethrinidae_ea_change ~ 1 +(1|Village),
                       data=seafood[!is.na(seafood$Lethrinidae_ea_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testLethrinidaeea <-brm(Lethrinidae_ea_change ~ 
                          Gender+
                          VillageYears_standard+
                          YearsFishing_standard+
                          Age_standard +(1|Village),
                        data=seafood[!is.na(seafood$Lethrinidae_ea_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Lethrinidaeea <- pp_check(testLethrinidaeea,ndraws = 50)+ggtitle("Model fit Lethrinidae Abundance")

loo(nullLethrinidaeea, testLethrinidaeea)

Lethrinidaeeaeffects<-as.data.frame(fixef(testLethrinidaeea, probs = c(0.1,0.9)))
Lethrinidaeeaeffects$community<-row.names(Lethrinidaeeaeffects)
colnames(Lethrinidaeeaeffects)<-c("estimate","se","conf.low","conf.high","community")
Lethrinidaeeaeff <- ggplot(Lethrinidaeeaeffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Lethrinidae ea")

##Lethrinidae harv
nullLethrinidaeharv<-brm(Lethrinidae_harv_change ~ 1 +(1|Village),
                         data=seafood[!is.na(seafood$Lethrinidae_harv_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testLethrinidaeharv <-brm(Lethrinidae_harv_change ~ 
                            Gender+
                            VillageYears_standard+
                            YearsFishing_standard+
                            Age_standard +(1|Village),
                          data=seafood[!is.na(seafood$Lethrinidae_harv_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Lethrinidaeharv <- pp_check(testLethrinidaeharv,ndraws = 50)+ggtitle("Model fit Lethrinidae Harvest")

loo(nullLethrinidaeharv, testLethrinidaeharv)

Lethrinidaeharveffects<-as.data.frame(fixef(testLethrinidaeharv, probs = c(0.1,0.9)))
Lethrinidaeharveffects$community<-row.names(Lethrinidaeharveffects)
colnames(Lethrinidaeharveffects)<-c("estimate","se","conf.low","conf.high","community")
Lethrinidaeharveff <- ggplot(Lethrinidaeharveffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Lethrinidae harv")

##Export graphs
jpeg("PPcheck_Lethrinidae", width =3500, height = 2000, res=300)
ggarrange(Lethrinidaecons, Lethrinidaeea, Lethrinidaeharv)
dev.off()

jpeg("FixedEffects_Lethrinidae", width =3500, height = 2000, res=300)
ggarrange(Lethrinidaeconseff, Lethrinidaeeaeff, Lethrinidaeharveff)
dev.off()





##Ostraciidae
seafood$Ostraciidae_cons_change <- seafood$Ostraciidae_cons_now - seafood$Ostraciidae_cons_10
seafood$Ostraciidae_ea_change <- seafood$Ostraciidae_ea_now - seafood$Ostraciidae_ea_10
seafood$Ostraciidae_harv_change <- seafood$Ostraciidae_harv_now - seafood$Ostraciidae_harv_10

seafood$Ostraciidae_cons_change <- ordered(as.factor(seafood$Ostraciidae_cons_change), 
                                           levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

seafood$Ostraciidae_ea_change <- ordered(as.factor(seafood$Ostraciidae_ea_change), 
                                         levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

seafood$Ostraciidae_harv_change <- ordered(as.factor(seafood$Ostraciidae_harv_change), 
                                           levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

##Ostraciidae Cons
nullOstraciidaecons<-brm(Ostraciidae_cons_change ~ 1 +(1|Village),
                         data=seafood[!is.na(seafood$Ostraciidae_cons_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testOstraciidaecons <-brm(Ostraciidae_cons_change ~ 
                            Gender+
                            VillageYears_standard+
                            YearsFishing_standard+
                            Age_standard +(1|Village),
                          data=seafood[!is.na(seafood$Ostraciidae_cons_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Ostraciidaecons <- pp_check(testOstraciidaecons,ndraws = 50)+ggtitle("Model fit Ostraciidae Consumption")

loo(nullOstraciidaecons, testOstraciidaecons)

Ostraciidaeconseffects<-as.data.frame(fixef(testOstraciidaecons, probs = c(0.1,0.9)))
Ostraciidaeconseffects$community<-row.names(Ostraciidaeconseffects)
colnames(Ostraciidaeconseffects)<-c("estimate","se","conf.low","conf.high","community")
Ostraciidaeconseff <- ggplot(Ostraciidaeconseffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Ostraciidae Cons")

##Ostraciidae ea
nullOstraciidaeea<-brm(Ostraciidae_ea_change ~ 1 +(1|Village),
                       data=seafood[!is.na(seafood$Ostraciidae_ea_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testOstraciidaeea <-brm(Ostraciidae_ea_change ~ 
                          Gender+
                          VillageYears_standard+
                          YearsFishing_standard+
                          Age_standard +(1|Village),
                        data=seafood[!is.na(seafood$Ostraciidae_ea_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Ostraciidaeea <- pp_check(testOstraciidaeea,ndraws = 50)+ggtitle("Model fit Ostraciidae Abundance")

loo(nullOstraciidaeea, testOstraciidaeea)

Ostraciidaeeaeffects<-as.data.frame(fixef(testOstraciidaeea, probs = c(0.1,0.9)))
Ostraciidaeeaeffects$community<-row.names(Ostraciidaeeaeffects)
colnames(Ostraciidaeeaeffects)<-c("estimate","se","conf.low","conf.high","community")
Ostraciidaeeaeff <- ggplot(Ostraciidaeeaeffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Ostraciidae ea")

##Ostraciidae harv
nullOstraciidaeharv<-brm(Ostraciidae_harv_change ~ 1 +(1|Village),
                         data=seafood[!is.na(seafood$Ostraciidae_harv_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testOstraciidaeharv <-brm(Ostraciidae_harv_change ~ 
                            Gender+
                            VillageYears_standard+
                            YearsFishing_standard+
                            Age_standard +(1|Village),
                          data=seafood[!is.na(seafood$Ostraciidae_harv_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Ostraciidaeharv <- pp_check(testOstraciidaeharv,ndraws = 50)+ggtitle("Model fit Ostraciidae Harvest")

loo(nullOstraciidaeharv, testOstraciidaeharv)

Ostraciidaeharveffects<-as.data.frame(fixef(testOstraciidaeharv, probs = c(0.1,0.9)))
Ostraciidaeharveffects$community<-row.names(Ostraciidaeharveffects)
colnames(Ostraciidaeharveffects)<-c("estimate","se","conf.low","conf.high","community")
Ostraciidaeharveff <- ggplot(Ostraciidaeharveffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Ostraciidae harv")

##Export graphs
jpeg("PPcheck_Ostraciidae", width =3500, height = 2000, res=300)
ggarrange(Ostraciidaecons, Ostraciidaeea, Ostraciidaeharv)
dev.off()

jpeg("FixedEffects_Ostraciidae", width =3500, height = 2000, res=300)
ggarrange(Ostraciidaeconseff, Ostraciidaeeaeff, Ostraciidaeharveff)
dev.off()










##Pomacanthidae
seafood$Pomacanthidae_cons_change <- seafood$Pomacanthidae_cons_now - seafood$Pomacanthidae_cons_10
seafood$Pomacanthidae_ea_change <- seafood$Pomacanthidae_ea_now - seafood$Pomacanthidae_ea_10
seafood$Pomacanthidae_harv_change <- seafood$Pomacanthidae_harv_now - seafood$Pomacanthidae_harv_10

seafood$Pomacanthidae_cons_change <- ordered(as.factor(seafood$Pomacanthidae_cons_change), 
                                             levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

seafood$Pomacanthidae_ea_change <- ordered(as.factor(seafood$Pomacanthidae_ea_change), 
                                           levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

seafood$Pomacanthidae_harv_change <- ordered(as.factor(seafood$Pomacanthidae_harv_change), 
                                             levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

##Pomacanthidae Cons
nullPomacanthidaecons<-brm(Pomacanthidae_cons_change ~ 1 +(1|Village),
                           data=seafood[!is.na(seafood$Pomacanthidae_cons_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testPomacanthidaecons <-brm(Pomacanthidae_cons_change ~ 
                              Gender+
                              VillageYears_standard+
                              YearsFishing_standard+
                              Age_standard +(1|Village),
                            data=seafood[!is.na(seafood$Pomacanthidae_cons_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Pomacanthidaecons <- pp_check(testPomacanthidaecons,ndraws = 50)+ggtitle("Model fit Pomacanthidae Consumption")

loo(nullPomacanthidaecons, testPomacanthidaecons)

Pomacanthidaeconseffects<-as.data.frame(fixef(testPomacanthidaecons, probs = c(0.1,0.9)))
Pomacanthidaeconseffects$community<-row.names(Pomacanthidaeconseffects)
colnames(Pomacanthidaeconseffects)<-c("estimate","se","conf.low","conf.high","community")
Pomacanthidaeconseff <- ggplot(Pomacanthidaeconseffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Pomacanthidae Cons")

##Pomacanthidae ea
nullPomacanthidaeea<-brm(Pomacanthidae_ea_change ~ 1 +(1|Village),
                         data=seafood[!is.na(seafood$Pomacanthidae_ea_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testPomacanthidaeea <-brm(Pomacanthidae_ea_change ~ 
                            Gender+
                            VillageYears_standard+
                            YearsFishing_standard+
                            Age_standard +(1|Village),
                          data=seafood[!is.na(seafood$Pomacanthidae_ea_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Pomacanthidaeea <- pp_check(testPomacanthidaeea,ndraws = 50)+ggtitle("Model fit Pomacanthidae Abundance")

loo(nullPomacanthidaeea, testPomacanthidaeea)

Pomacanthidaeeaeffects<-as.data.frame(fixef(testPomacanthidaeea, probs = c(0.1,0.9)))
Pomacanthidaeeaeffects$community<-row.names(Pomacanthidaeeaeffects)
colnames(Pomacanthidaeeaeffects)<-c("estimate","se","conf.low","conf.high","community")
Pomacanthidaeeaeff <- ggplot(Pomacanthidaeeaeffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Pomacanthidae ea")

##Pomacanthidae harv
nullPomacanthidaeharv<-brm(Pomacanthidae_harv_change ~ 1 +(1|Village),
                           data=seafood[!is.na(seafood$Pomacanthidae_harv_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testPomacanthidaeharv <-brm(Pomacanthidae_harv_change ~ 
                              Gender+
                              VillageYears_standard+
                              YearsFishing_standard+
                              Age_standard +(1|Village),
                            data=seafood[!is.na(seafood$Pomacanthidae_harv_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Pomacanthidaeharv <- pp_check(testPomacanthidaeharv,ndraws = 50)+ggtitle("Model fit Pomacanthidae Harvest")

loo(nullPomacanthidaeharv, testPomacanthidaeharv)

Pomacanthidaeharveffects<-as.data.frame(fixef(testPomacanthidaeharv, probs = c(0.1,0.9)))
Pomacanthidaeharveffects$community<-row.names(Pomacanthidaeharveffects)
colnames(Pomacanthidaeharveffects)<-c("estimate","se","conf.low","conf.high","community")
Pomacanthidaeharveff <- ggplot(Pomacanthidaeharveffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Pomacanthidae harv")

##Export graphs
jpeg("PPcheck_Pomacanthidae", width =3500, height = 2000, res=300)
ggarrange(Pomacanthidaecons, Pomacanthidaeea, Pomacanthidaeharv)
dev.off()

jpeg("FixedEffects_Pomacanthidae", width =3500, height = 2000, res=300)
ggarrange(Pomacanthidaeconseff, Pomacanthidaeeaeff, Pomacanthidaeharveff)
dev.off()








##Rajidae
seafood$Rajidae_cons_change <- seafood$Rajidae_cons_now - seafood$Rajidae_cons_10
seafood$Rajidae_ea_change <- seafood$Rajidae_ea_now - seafood$Rajidae_ea_10
seafood$Rajidae_harv_change <- seafood$Rajidae_harv_now - seafood$Rajidae_harv_10

seafood$Rajidae_cons_change <- ordered(as.factor(seafood$Rajidae_cons_change), 
                                       levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

seafood$Rajidae_ea_change <- ordered(as.factor(seafood$Rajidae_ea_change), 
                                     levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

seafood$Rajidae_harv_change <- ordered(as.factor(seafood$Rajidae_harv_change), 
                                       levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

##Rajidae Cons
nullRajidaecons<-brm(Rajidae_cons_change ~ 1 +(1|Village),
                     data=seafood[!is.na(seafood$Rajidae_cons_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testRajidaecons <-brm(Rajidae_cons_change ~ 
                        Gender+
                        VillageYears_standard+
                        YearsFishing_standard+
                        Age_standard +(1|Village),
                      data=seafood[!is.na(seafood$Rajidae_cons_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Rajidaecons <- pp_check(testRajidaecons,ndraws = 50)+ggtitle("Model fit Rajidae Consumption")

loo(nullRajidaecons, testRajidaecons)

Rajidaeconseffects<-as.data.frame(fixef(testRajidaecons, probs = c(0.1,0.9)))
Rajidaeconseffects$community<-row.names(Rajidaeconseffects)
colnames(Rajidaeconseffects)<-c("estimate","se","conf.low","conf.high","community")
Rajidaeconseff <- ggplot(Rajidaeconseffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Rajidae Cons")

##Rajidae ea
nullRajidaeea<-brm(Rajidae_ea_change ~ 1 +(1|Village),
                   data=seafood[!is.na(seafood$Rajidae_ea_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testRajidaeea <-brm(Rajidae_ea_change ~ 
                      Gender+
                      VillageYears_standard+
                      YearsFishing_standard+
                      Age_standard +(1|Village),
                    data=seafood[!is.na(seafood$Rajidae_ea_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Rajidaeea <- pp_check(testRajidaeea,ndraws = 50)+ggtitle("Model fit Rajidae Abundance")

loo(nullRajidaeea, testRajidaeea)

Rajidaeeaeffects<-as.data.frame(fixef(testRajidaeea, probs = c(0.1,0.9)))
Rajidaeeaeffects$community<-row.names(Rajidaeeaeffects)
colnames(Rajidaeeaeffects)<-c("estimate","se","conf.low","conf.high","community")
Rajidaeeaeff <- ggplot(Rajidaeeaeffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Rajidae ea")

##Rajidae harv
nullRajidaeharv<-brm(Rajidae_harv_change ~ 1 +(1|Village),
                     data=seafood[!is.na(seafood$Rajidae_harv_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testRajidaeharv <-brm(Rajidae_harv_change ~ 
                        Gender+
                        VillageYears_standard+
                        YearsFishing_standard+
                        Age_standard +(1|Village),
                      data=seafood[!is.na(seafood$Rajidae_harv_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Rajidaeharv <- pp_check(testRajidaeharv,ndraws = 50)+ggtitle("Model fit Rajidae Harvest")

loo(nullRajidaeharv, testRajidaeharv)

Rajidaeharveffects<-as.data.frame(fixef(testRajidaeharv, probs = c(0.1,0.9)))
Rajidaeharveffects$community<-row.names(Rajidaeharveffects)
colnames(Rajidaeharveffects)<-c("estimate","se","conf.low","conf.high","community")
Rajidaeharveff <- ggplot(Rajidaeharveffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Rajidae harv")

##Export graphs
jpeg("PPcheck_Rajidae", width =3500, height = 2000, res=300)
ggarrange(Rajidaecons, Rajidaeea, Rajidaeharv)
dev.off()

jpeg("FixedEffects_Rajidae", width =3500, height = 2000, res=300)
ggarrange(Rajidaeconseff, Rajidaeeaeff, Rajidaeharveff)
dev.off()







##Rhinobatidae
seafood$Rhinobatidae_cons_change <- seafood$Rhinobatidae_cons_now - seafood$Rhinobatidae_cons_10
seafood$Rhinobatidae_ea_change <- seafood$Rhinobatidae_ea_now - seafood$Rhinobatidae_ea_10
seafood$Rhinobatidae_harv_change <- seafood$Rhinobatidae_harv_now - seafood$Rhinobatidae_harv_10

seafood$Rhinobatidae_cons_change <- ordered(as.factor(seafood$Rhinobatidae_cons_change), 
                                            levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

seafood$Rhinobatidae_ea_change <- ordered(as.factor(seafood$Rhinobatidae_ea_change), 
                                          levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

seafood$Rhinobatidae_harv_change <- ordered(as.factor(seafood$Rhinobatidae_harv_change), 
                                            levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

##Rhinobatidae Cons
nullRhinobatidaecons<-brm(Rhinobatidae_cons_change ~ 1 +(1|Village),
                          data=seafood[!is.na(seafood$Rhinobatidae_cons_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testRhinobatidaecons <-brm(Rhinobatidae_cons_change ~ 
                             Gender+
                             VillageYears_standard+
                             YearsFishing_standard+
                             Age_standard +(1|Village),
                           data=seafood[!is.na(seafood$Rhinobatidae_cons_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Rhinobatidaecons <- pp_check(testRhinobatidaecons,ndraws = 50)+ggtitle("Model fit Rhinobatidae Consumption")

loo(nullRhinobatidaecons, testRhinobatidaecons)

Rhinobatidaeconseffects<-as.data.frame(fixef(testRhinobatidaecons, probs = c(0.1,0.9)))
Rhinobatidaeconseffects$community<-row.names(Rhinobatidaeconseffects)
colnames(Rhinobatidaeconseffects)<-c("estimate","se","conf.low","conf.high","community")
Rhinobatidaeconseff <- ggplot(Rhinobatidaeconseffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Rhinobatidae Cons")

##Rhinobatidae ea
nullRhinobatidaeea<-brm(Rhinobatidae_ea_change ~ 1 +(1|Village),
                        data=seafood[!is.na(seafood$Rhinobatidae_ea_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testRhinobatidaeea <-brm(Rhinobatidae_ea_change ~ 
                           Gender+
                           VillageYears_standard+
                           YearsFishing_standard+
                           Age_standard +(1|Village),
                         data=seafood[!is.na(seafood$Rhinobatidae_ea_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Rhinobatidaeea <- pp_check(testRhinobatidaeea,ndraws = 50)+ggtitle("Model fit Rhinobatidae Abundance")

loo(nullRhinobatidaeea, testRhinobatidaeea)

Rhinobatidaeeaeffects<-as.data.frame(fixef(testRhinobatidaeea, probs = c(0.1,0.9)))
Rhinobatidaeeaeffects$community<-row.names(Rhinobatidaeeaeffects)
colnames(Rhinobatidaeeaeffects)<-c("estimate","se","conf.low","conf.high","community")
Rhinobatidaeeaeff <- ggplot(Rhinobatidaeeaeffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Rhinobatidae ea")

##Rhinobatidae harv
nullRhinobatidaeharv<-brm(Rhinobatidae_harv_change ~ 1 +(1|Village),
                          data=seafood[!is.na(seafood$Rhinobatidae_harv_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testRhinobatidaeharv <-brm(Rhinobatidae_harv_change ~ 
                             Gender+
                             VillageYears_standard+
                             YearsFishing_standard+
                             Age_standard +(1|Village),
                           data=seafood[!is.na(seafood$Rhinobatidae_harv_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Rhinobatidaeharv <- pp_check(testRhinobatidaeharv,ndraws = 50)+ggtitle("Model fit Rhinobatidae Harvest")

loo(nullRhinobatidaeharv, testRhinobatidaeharv)

Rhinobatidaeharveffects<-as.data.frame(fixef(testRhinobatidaeharv, probs = c(0.1,0.9)))
Rhinobatidaeharveffects$community<-row.names(Rhinobatidaeharveffects)
colnames(Rhinobatidaeharveffects)<-c("estimate","se","conf.low","conf.high","community")
Rhinobatidaeharveff <- ggplot(Rhinobatidaeharveffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Rhinobatidae harv")

##Export graphs
jpeg("PPcheck_Rhinobatidae", width =3500, height = 2000, res=300)
ggarrange(Rhinobatidaecons, Rhinobatidaeea, Rhinobatidaeharv)
dev.off()

jpeg("FixedEffects_Rhinobatidae", width =3500, height = 2000, res=300)
ggarrange(Rhinobatidaeconseff, Rhinobatidaeeaeff, Rhinobatidaeharveff)
dev.off()






##Soleidae
seafood$Soleidae_cons_change <- seafood$Soleidae_cons_now - seafood$Soleidae_cons_10
seafood$Soleidae_ea_change <- seafood$Soleidae_ea_now - seafood$Soleidae_ea_10
seafood$Soleidae_harv_change <- seafood$Soleidae_harv_now - seafood$Soleidae_harv_10

seafood$Soleidae_cons_change <- ordered(as.factor(seafood$Soleidae_cons_change), 
                                        levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

seafood$Soleidae_ea_change <- ordered(as.factor(seafood$Soleidae_ea_change), 
                                      levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

seafood$Soleidae_harv_change <- ordered(as.factor(seafood$Soleidae_harv_change), 
                                        levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

##Soleidae Cons
nullSoleidaecons<-brm(Soleidae_cons_change ~ 1 +(1|Village),
                      data=seafood[!is.na(seafood$Soleidae_cons_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testSoleidaecons <-brm(Soleidae_cons_change ~ 
                         Gender+
                         VillageYears_standard+
                         YearsFishing_standard+
                         Age_standard +(1|Village),
                       data=seafood[!is.na(seafood$Soleidae_cons_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Soleidaecons <- pp_check(testSoleidaecons,ndraws = 50)+ggtitle("Model fit Soleidae Consumption")

loo(nullSoleidaecons, testSoleidaecons)

Soleidaeconseffects<-as.data.frame(fixef(testSoleidaecons, probs = c(0.1,0.9)))
Soleidaeconseffects$community<-row.names(Soleidaeconseffects)
colnames(Soleidaeconseffects)<-c("estimate","se","conf.low","conf.high","community")
Soleidaeconseff <- ggplot(Soleidaeconseffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Soleidae Cons")

##Soleidae ea
nullSoleidaeea<-brm(Soleidae_ea_change ~ 1 +(1|Village),
                    data=seafood[!is.na(seafood$Soleidae_ea_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testSoleidaeea <-brm(Soleidae_ea_change ~ 
                       Gender+
                       VillageYears_standard+
                       YearsFishing_standard+
                       Age_standard +(1|Village),
                     data=seafood[!is.na(seafood$Soleidae_ea_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Soleidaeea <- pp_check(testSoleidaeea,ndraws = 50)+ggtitle("Model fit Soleidae Abundance")

loo(nullSoleidaeea, testSoleidaeea)

Soleidaeeaeffects<-as.data.frame(fixef(testSoleidaeea, probs = c(0.1,0.9)))
Soleidaeeaeffects$community<-row.names(Soleidaeeaeffects)
colnames(Soleidaeeaeffects)<-c("estimate","se","conf.low","conf.high","community")
Soleidaeeaeff <- ggplot(Soleidaeeaeffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Soleidae ea")

##Soleidae harv
nullSoleidaeharv<-brm(Soleidae_harv_change ~ 1 +(1|Village),
                      data=seafood[!is.na(seafood$Soleidae_harv_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testSoleidaeharv <-brm(Soleidae_harv_change ~ 
                         Gender+
                         VillageYears_standard+
                         YearsFishing_standard+
                         Age_standard +(1|Village),
                       data=seafood[!is.na(seafood$Soleidae_harv_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Soleidaeharv <- pp_check(testSoleidaeharv,ndraws = 50)+ggtitle("Model fit Soleidae Harvest")

loo(nullSoleidaeharv, testSoleidaeharv)

Soleidaeharveffects<-as.data.frame(fixef(testSoleidaeharv, probs = c(0.1,0.9)))
Soleidaeharveffects$community<-row.names(Soleidaeharveffects)
colnames(Soleidaeharveffects)<-c("estimate","se","conf.low","conf.high","community")
Soleidaeharveff <- ggplot(Soleidaeharveffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Soleidae harv")

##Export graphs
jpeg("PPcheck_Soleidae", width =3500, height = 2000, res=300)
ggarrange(Soleidaecons, Soleidaeea, Soleidaeharv)
dev.off()

jpeg("FixedEffects_Soleidae", width =3500, height = 2000, res=300)
ggarrange(Soleidaeconseff, Soleidaeeaeff, Soleidaeharveff)
dev.off()







##Tetraodontidae
seafood$Tetraodontidae_cons_change <- seafood$Tetraodontidae_cons_now - seafood$Tetraodontidae_cons_10
seafood$Tetraodontidae_ea_change <- seafood$Tetraodontidae_ea_now - seafood$Tetraodontidae_ea_10
seafood$Tetraodontidae_harv_change <- seafood$Tetraodontidae_harv_now - seafood$Tetraodontidae_harv_10

seafood$Tetraodontidae_cons_change <- ordered(as.factor(seafood$Tetraodontidae_cons_change), 
                                              levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

seafood$Tetraodontidae_ea_change <- ordered(as.factor(seafood$Tetraodontidae_ea_change), 
                                            levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

seafood$Tetraodontidae_harv_change <- ordered(as.factor(seafood$Tetraodontidae_harv_change), 
                                              levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

##Tetraodontidae Cons
nullTetraodontidaecons<-brm(Tetraodontidae_cons_change ~ 1 +(1|Village),
                            data=seafood[!is.na(seafood$Tetraodontidae_cons_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testTetraodontidaecons <-brm(Tetraodontidae_cons_change ~ 
                               Gender+
                               VillageYears_standard+
                               YearsFishing_standard+
                               Age_standard +(1|Village),
                             data=seafood[!is.na(seafood$Tetraodontidae_cons_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Tetraodontidaecons <- pp_check(testTetraodontidaecons,ndraws = 50)+ggtitle("Model fit Tetraodontidae Consumption")

loo(nullTetraodontidaecons, testTetraodontidaecons)

Tetraodontidaeconseffects<-as.data.frame(fixef(testTetraodontidaecons, probs = c(0.1,0.9)))
Tetraodontidaeconseffects$community<-row.names(Tetraodontidaeconseffects)
colnames(Tetraodontidaeconseffects)<-c("estimate","se","conf.low","conf.high","community")
Tetraodontidaeconseff <- ggplot(Tetraodontidaeconseffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Tetraodontidae Cons")

##Tetraodontidae ea
nullTetraodontidaeea<-brm(Tetraodontidae_ea_change ~ 1 +(1|Village),
                          data=seafood[!is.na(seafood$Tetraodontidae_ea_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testTetraodontidaeea <-brm(Tetraodontidae_ea_change ~ 
                             Gender+
                             VillageYears_standard+
                             YearsFishing_standard+
                             Age_standard +(1|Village),
                           data=seafood[!is.na(seafood$Tetraodontidae_ea_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Tetraodontidaeea <- pp_check(testTetraodontidaeea,ndraws = 50)+ggtitle("Model fit Tetraodontidae Abundance")

loo(nullTetraodontidaeea, testTetraodontidaeea)

Tetraodontidaeeaeffects<-as.data.frame(fixef(testTetraodontidaeea, probs = c(0.1,0.9)))
Tetraodontidaeeaeffects$community<-row.names(Tetraodontidaeeaeffects)
colnames(Tetraodontidaeeaeffects)<-c("estimate","se","conf.low","conf.high","community")
Tetraodontidaeeaeff <- ggplot(Tetraodontidaeeaeffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Tetraodontidae ea")

##Tetraodontidae harv
nullTetraodontidaeharv<-brm(Tetraodontidae_harv_change ~ 1 +(1|Village),
                            data=seafood[!is.na(seafood$Tetraodontidae_harv_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testTetraodontidaeharv <-brm(Tetraodontidae_harv_change ~ 
                               Gender+
                               VillageYears_standard+
                               YearsFishing_standard+
                               Age_standard +(1|Village),
                             data=seafood[!is.na(seafood$Tetraodontidae_harv_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Tetraodontidaeharv <- pp_check(testTetraodontidaeharv,ndraws = 50)+ggtitle("Model fit Tetraodontidae Harvest")

loo(nullTetraodontidaeharv, testTetraodontidaeharv)

Tetraodontidaeharveffects<-as.data.frame(fixef(testTetraodontidaeharv, probs = c(0.1,0.9)))
Tetraodontidaeharveffects$community<-row.names(Tetraodontidaeharveffects)
colnames(Tetraodontidaeharveffects)<-c("estimate","se","conf.low","conf.high","community")
Tetraodontidaeharveff <- ggplot(Tetraodontidaeharveffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Tetraodontidae harv")

##Export graphs
jpeg("PPcheck_Tetraodontidae", width =3500, height = 2000, res=300)
ggarrange(Tetraodontidaecons, Tetraodontidaeea, Tetraodontidaeharv)
dev.off()

jpeg("FixedEffects_Tetraodontidae", width =3500, height = 2000, res=300)
ggarrange(Tetraodontidaeconseff, Tetraodontidaeeaeff, Tetraodontidaeharveff)
dev.off()









##Zanclidae
seafood$Zanclidae_cons_change <- seafood$Zanclidae_cons_now - seafood$Zanclidae_cons_10
seafood$Zanclidae_ea_change <- seafood$Zanclidae_ea_now - seafood$Zanclidae_ea_10
seafood$Zanclidae_harv_change <- seafood$Zanclidae_harv_now - seafood$Zanclidae_harv_10

seafood$Zanclidae_cons_change <- ordered(as.factor(seafood$Zanclidae_cons_change), 
                                         levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

seafood$Zanclidae_ea_change <- ordered(as.factor(seafood$Zanclidae_ea_change), 
                                       levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

seafood$Zanclidae_harv_change <- ordered(as.factor(seafood$Zanclidae_harv_change), 
                                         levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

##Zanclidae Cons
nullZanclidaecons<-brm(Zanclidae_cons_change ~ 1 +(1|Village),
                       data=seafood[!is.na(seafood$Zanclidae_cons_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testZanclidaecons <-brm(Zanclidae_cons_change ~ 
                          Gender+
                          VillageYears_standard+
                          YearsFishing_standard+
                          Age_standard +(1|Village),
                        data=seafood[!is.na(seafood$Zanclidae_cons_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Zanclidaecons <- pp_check(testZanclidaecons,ndraws = 50)+ggtitle("Model fit Zanclidae Consumption")

loo(nullZanclidaecons, testZanclidaecons)

Zanclidaeconseffects<-as.data.frame(fixef(testZanclidaecons, probs = c(0.1,0.9)))
Zanclidaeconseffects$community<-row.names(Zanclidaeconseffects)
colnames(Zanclidaeconseffects)<-c("estimate","se","conf.low","conf.high","community")
Zanclidaeconseff <- ggplot(Zanclidaeconseffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Zanclidae Cons")

##Zanclidae ea
nullZanclidaeea<-brm(Zanclidae_ea_change ~ 1 +(1|Village),
                     data=seafood[!is.na(seafood$Zanclidae_ea_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testZanclidaeea <-brm(Zanclidae_ea_change ~ 
                        Gender+
                        VillageYears_standard+
                        YearsFishing_standard+
                        Age_standard +(1|Village),
                      data=seafood[!is.na(seafood$Zanclidae_ea_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Zanclidaeea <- pp_check(testZanclidaeea,ndraws = 50)+ggtitle("Model fit Zanclidae Abundance")

loo(nullZanclidaeea, testZanclidaeea)

Zanclidaeeaeffects<-as.data.frame(fixef(testZanclidaeea, probs = c(0.1,0.9)))
Zanclidaeeaeffects$community<-row.names(Zanclidaeeaeffects)
colnames(Zanclidaeeaeffects)<-c("estimate","se","conf.low","conf.high","community")
Zanclidaeeaeff <- ggplot(Zanclidaeeaeffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Zanclidae ea")

##Zanclidae harv
nullZanclidaeharv<-brm(Zanclidae_harv_change ~ 1 +(1|Village),
                       data=seafood[!is.na(seafood$Zanclidae_harv_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testZanclidaeharv <-brm(Zanclidae_harv_change ~ 
                          Gender+
                          VillageYears_standard+
                          YearsFishing_standard+
                          Age_standard +(1|Village),
                        data=seafood[!is.na(seafood$Zanclidae_harv_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Zanclidaeharv <- pp_check(testZanclidaeharv,ndraws = 50)+ggtitle("Model fit Zanclidae Harvest")

loo(nullZanclidaeharv, testZanclidaeharv)

Zanclidaeharveffects<-as.data.frame(fixef(testZanclidaeharv, probs = c(0.1,0.9)))
Zanclidaeharveffects$community<-row.names(Zanclidaeharveffects)
colnames(Zanclidaeharveffects)<-c("estimate","se","conf.low","conf.high","community")
Zanclidaeharveff <- ggplot(Zanclidaeharveffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Zanclidae harv")

##Export graphs
jpeg("PPcheck_Zanclidae", width =3500, height = 2000, res=300)
ggarrange(Zanclidaecons, Zanclidaeea, Zanclidaeharv)
dev.off()

jpeg("FixedEffects_Zanclidae", width =3500, height = 2000, res=300)
ggarrange(Zanclidaeconseff, Zanclidaeeaeff, Zanclidaeharveff)
dev.off()







##Octopus
seafood$Octopus_cons_change <- seafood$Octopus_cons_now - seafood$Octopus_cons_10
seafood$Octopus_ea_change <- seafood$Octopus_ea_now - seafood$Octopus_ea_10
seafood$Octopus_harv_change <- seafood$Octopus_harv_now - seafood$Octopus_harv_10

seafood$Octopus_cons_change <- ordered(as.factor(seafood$Octopus_cons_change), 
                                       levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

seafood$Octopus_ea_change <- ordered(as.factor(seafood$Octopus_ea_change), 
                                     levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

seafood$Octopus_harv_change <- ordered(as.factor(seafood$Octopus_harv_change), 
                                       levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

##Octopus Cons
nullOctopuscons<-brm(Octopus_cons_change ~ 1 +(1|Village),
                     data=seafood[!is.na(seafood$Octopus_cons_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testOctopuscons <-brm(Octopus_cons_change ~ 
                        Gender+
                        VillageYears_standard+
                        YearsFishing_standard+
                        Age_standard +(1|Village),
                      data=seafood[!is.na(seafood$Octopus_cons_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Octopuscons <- pp_check(testOctopuscons,ndraws = 50)+ggtitle("Model fit Octopus Consumption")

loo(nullOctopuscons, testOctopuscons)

Octopusconseffects<-as.data.frame(fixef(testOctopuscons, probs = c(0.1,0.9)))
Octopusconseffects$community<-row.names(Octopusconseffects)
colnames(Octopusconseffects)<-c("estimate","se","conf.low","conf.high","community")
Octopusconseff <- ggplot(Octopusconseffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Octopus Cons")

##Octopus ea
nullOctopusea<-brm(Octopus_ea_change ~ 1 +(1|Village),
                   data=seafood[!is.na(seafood$Octopus_ea_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testOctopusea <-brm(Octopus_ea_change ~ 
                      Gender+
                      VillageYears_standard+
                      YearsFishing_standard+
                      Age_standard +(1|Village),
                    data=seafood[!is.na(seafood$Octopus_ea_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Octopusea <- pp_check(testOctopusea,ndraws = 50)+ggtitle("Model fit Octopus Abundance")

loo(nullOctopusea, testOctopusea)

Octopuseaeffects<-as.data.frame(fixef(testOctopusea, probs = c(0.1,0.9)))
Octopuseaeffects$community<-row.names(Octopuseaeffects)
colnames(Octopuseaeffects)<-c("estimate","se","conf.low","conf.high","community")
Octopuseaeff <- ggplot(Octopuseaeffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Octopus ea")

##Octopus harv
nullOctopusharv<-brm(Octopus_harv_change ~ 1 +(1|Village),
                     data=seafood[!is.na(seafood$Octopus_harv_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testOctopusharv <-brm(Octopus_harv_change ~ 
                        Gender+
                        VillageYears_standard+
                        YearsFishing_standard+
                        Age_standard +(1|Village),
                      data=seafood[!is.na(seafood$Octopus_harv_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Octopusharv <- pp_check(testOctopusharv,ndraws = 50)+ggtitle("Model fit Octopus Harvest")

loo(nullOctopusharv, testOctopusharv)

Octopusharveffects<-as.data.frame(fixef(testOctopusharv, probs = c(0.1,0.9)))
Octopusharveffects$community<-row.names(Octopusharveffects)
colnames(Octopusharveffects)<-c("estimate","se","conf.low","conf.high","community")
Octopusharveff <- ggplot(Octopusharveffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Octopus harv")

##Export graphs
jpeg("PPcheck_Octopus", width =3500, height = 2000, res=300)
ggarrange(Octopuscons, Octopusea, Octopusharv)
dev.off()

jpeg("FixedEffects_Octopus", width =3500, height = 2000, res=300)
ggarrange(Octopusconseff, Octopuseaeff, Octopusharveff)
dev.off()









##Loligo
seafood$Loligo_cons_change <- seafood$Loligo_cons_now - seafood$Loligo_cons_10
seafood$Loligo_ea_change <- seafood$Loligo_ea_now - seafood$Loligo_ea_10
seafood$Loligo_harv_change <- seafood$Loligo_harv_now - seafood$Loligo_harv_10

seafood$Loligo_cons_change <- ordered(as.factor(seafood$Loligo_cons_change), 
                                      levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

seafood$Loligo_ea_change <- ordered(as.factor(seafood$Loligo_ea_change), 
                                    levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

seafood$Loligo_harv_change <- ordered(as.factor(seafood$Loligo_harv_change), 
                                      levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

##Loligo Cons
nullLoligocons<-brm(Loligo_cons_change ~ 1 +(1|Village),
                    data=seafood[!is.na(seafood$Loligo_cons_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testLoligocons <-brm(Loligo_cons_change ~ 
                       Gender+
                       VillageYears_standard+
                       YearsFishing_standard+
                       Age_standard +(1|Village),
                     data=seafood[!is.na(seafood$Loligo_cons_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Loligocons <- pp_check(testLoligocons,ndraws = 50)+ggtitle("Model fit Loligo Consumption")

loo(nullLoligocons, testLoligocons)

Loligoconseffects<-as.data.frame(fixef(testLoligocons, probs = c(0.1,0.9)))
Loligoconseffects$community<-row.names(Loligoconseffects)
colnames(Loligoconseffects)<-c("estimate","se","conf.low","conf.high","community")
Loligoconseff <- ggplot(Loligoconseffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Loligo Cons")

##Loligo ea
nullLoligoea<-brm(Loligo_ea_change ~ 1 +(1|Village),
                  data=seafood[!is.na(seafood$Loligo_ea_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testLoligoea <-brm(Loligo_ea_change ~ 
                     Gender+
                     VillageYears_standard+
                     YearsFishing_standard+
                     Age_standard +(1|Village),
                   data=seafood[!is.na(seafood$Loligo_ea_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Loligoea <- pp_check(testLoligoea,ndraws = 50)+ggtitle("Model fit Loligo Abundance")

loo(nullLoligoea, testLoligoea)

Loligoeaeffects<-as.data.frame(fixef(testLoligoea, probs = c(0.1,0.9)))
Loligoeaeffects$community<-row.names(Loligoeaeffects)
colnames(Loligoeaeffects)<-c("estimate","se","conf.low","conf.high","community")
Loligoeaeff <- ggplot(Loligoeaeffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Loligo ea")

##Loligo harv
nullLoligoharv<-brm(Loligo_harv_change ~ 1 +(1|Village),
                    data=seafood[!is.na(seafood$Loligo_harv_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testLoligoharv <-brm(Loligo_harv_change ~ 
                       Gender+
                       VillageYears_standard+
                       YearsFishing_standard+
                       Age_standard +(1|Village),
                     data=seafood[!is.na(seafood$Loligo_harv_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Loligoharv <- pp_check(testLoligoharv,ndraws = 50)+ggtitle("Model fit Loligo Harvest")

loo(nullLoligoharv, testLoligoharv)

Loligoharveffects<-as.data.frame(fixef(testLoligoharv, probs = c(0.1,0.9)))
Loligoharveffects$community<-row.names(Loligoharveffects)
colnames(Loligoharveffects)<-c("estimate","se","conf.low","conf.high","community")
Loligoharveff <- ggplot(Loligoharveffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Loligo harv")

##Export graphs
jpeg("PPcheck_Loligo", width =3500, height = 2000, res=300)
ggarrange(Loligocons, Loligoea, Loligoharv)
dev.off()

jpeg("FixedEffects_Loligo", width =3500, height = 2000, res=300)
ggarrange(Loligoconseff, Loligoeaeff, Loligoharveff)
dev.off()






##Sepia
seafood$Sepia_cons_change <- seafood$Sepia_cons_now - seafood$Sepia_cons_10
seafood$Sepia_ea_change <- seafood$Sepia_ea_now - seafood$Sepia_ea_10
seafood$Sepia_harv_change <- seafood$Sepia_harv_now - seafood$Sepia_harv_10

seafood$Sepia_cons_change <- ordered(as.factor(seafood$Sepia_cons_change), 
                                     levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

seafood$Sepia_ea_change <- ordered(as.factor(seafood$Sepia_ea_change), 
                                   levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

seafood$Sepia_harv_change <- ordered(as.factor(seafood$Sepia_harv_change), 
                                     levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

##Sepia Cons
nullSepiacons<-brm(Sepia_cons_change ~ 1 +(1|Village),
                   data=seafood[!is.na(seafood$Sepia_cons_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testSepiacons <-brm(Sepia_cons_change ~ 
                      Gender+
                      VillageYears_standard+
                      YearsFishing_standard+
                      Age_standard +(1|Village),
                    data=seafood[!is.na(seafood$Sepia_cons_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Sepiacons <- pp_check(testSepiacons,ndraws = 50)+ggtitle("Model fit Sepia Consumption")

loo(nullSepiacons, testSepiacons)

Sepiaconseffects<-as.data.frame(fixef(testSepiacons, probs = c(0.1,0.9)))
Sepiaconseffects$community<-row.names(Sepiaconseffects)
colnames(Sepiaconseffects)<-c("estimate","se","conf.low","conf.high","community")
Sepiaconseff <- ggplot(Sepiaconseffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Sepia Cons")

##Sepia ea
nullSepiaea<-brm(Sepia_ea_change ~ 1 +(1|Village),
                 data=seafood[!is.na(seafood$Sepia_ea_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testSepiaea <-brm(Sepia_ea_change ~ 
                    Gender+
                    VillageYears_standard+
                    YearsFishing_standard+
                    Age_standard +(1|Village),
                  data=seafood[!is.na(seafood$Sepia_ea_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Sepiaea <- pp_check(testSepiaea,ndraws = 50)+ggtitle("Model fit Sepia Abundance")

loo(nullSepiaea, testSepiaea)

Sepiaeaeffects<-as.data.frame(fixef(testSepiaea, probs = c(0.1,0.9)))
Sepiaeaeffects$community<-row.names(Sepiaeaeffects)
colnames(Sepiaeaeffects)<-c("estimate","se","conf.low","conf.high","community")
Sepiaeaeff <- ggplot(Sepiaeaeffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Sepia ea")

##Sepia harv
nullSepiaharv<-brm(Sepia_harv_change ~ 1 +(1|Village),
                   data=seafood[!is.na(seafood$Sepia_harv_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testSepiaharv <-brm(Sepia_harv_change ~ 
                      Gender+
                      VillageYears_standard+
                      YearsFishing_standard+
                      Age_standard +(1|Village),
                    data=seafood[!is.na(seafood$Sepia_harv_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Sepiaharv <- pp_check(testSepiaharv,ndraws = 50)+ggtitle("Model fit Sepia Harvest")

loo(nullSepiaharv, testSepiaharv)

Sepiaharveffects<-as.data.frame(fixef(testSepiaharv, probs = c(0.1,0.9)))
Sepiaharveffects$community<-row.names(Sepiaharveffects)
colnames(Sepiaharveffects)<-c("estimate","se","conf.low","conf.high","community")
Sepiaharveff <- ggplot(Sepiaharveffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Sepia harv")

##Export graphs
jpeg("PPcheck_Sepia", width =3500, height = 2000, res=300)
ggarrange(Sepiacons, Sepiaea, Sepiaharv)
dev.off()

jpeg("FixedEffects_Sepia", width =3500, height = 2000, res=300)
ggarrange(Sepiaconseff, Sepiaeaeff, Sepiaharveff)
dev.off()






##Cypraea.Conus
seafood$Cypraea.Conus_cons_change <- seafood$Cypraea.Conus_cons_now - seafood$Cypraea.Conus_cons_10
seafood$Cypraea.Conus_ea_change <- seafood$Cypraea.Conus_ea_now - seafood$Cypraea.Conus_ea_10
seafood$Cypraea.Conus_harv_change <- seafood$Cypraea.Conus_harv_now - seafood$Cypraea.Conus_harv_10

seafood$Cypraea.Conus_cons_change <- ordered(as.factor(seafood$Cypraea.Conus_cons_change), 
                                             levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

seafood$Cypraea.Conus_ea_change <- ordered(as.factor(seafood$Cypraea.Conus_ea_change), 
                                           levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

seafood$Cypraea.Conus_harv_change <- ordered(as.factor(seafood$Cypraea.Conus_harv_change), 
                                             levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

##Cypraea.Conus Cons
nullCypraea.Conuscons<-brm(Cypraea.Conus_cons_change ~ 1 +(1|Village),
                           data=seafood[!is.na(seafood$Cypraea.Conus_cons_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testCypraea.Conuscons <-brm(Cypraea.Conus_cons_change ~ 
                              Gender+
                              VillageYears_standard+
                              YearsFishing_standard+
                              Age_standard +(1|Village),
                            data=seafood[!is.na(seafood$Cypraea.Conus_cons_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Cypraea.Conuscons <- pp_check(testCypraea.Conuscons,ndraws = 50)+ggtitle("Model fit Cypraea.Conus Consumption")

loo(nullCypraea.Conuscons, testCypraea.Conuscons)

Cypraea.Conusconseffects<-as.data.frame(fixef(testCypraea.Conuscons, probs = c(0.1,0.9)))
Cypraea.Conusconseffects$community<-row.names(Cypraea.Conusconseffects)
colnames(Cypraea.Conusconseffects)<-c("estimate","se","conf.low","conf.high","community")
Cypraea.Conusconseff <- ggplot(Cypraea.Conusconseffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Cypraea.Conus Cons")

##Cypraea.Conus ea
nullCypraea.Conusea<-brm(Cypraea.Conus_ea_change ~ 1 +(1|Village),
                         data=seafood[!is.na(seafood$Cypraea.Conus_ea_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testCypraea.Conusea <-brm(Cypraea.Conus_ea_change ~ 
                            Gender+
                            VillageYears_standard+
                            YearsFishing_standard+
                            Age_standard +(1|Village),
                          data=seafood[!is.na(seafood$Cypraea.Conus_ea_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Cypraea.Conusea <- pp_check(testCypraea.Conusea,ndraws = 50)+ggtitle("Model fit Cypraea.Conus Abundance")

loo(nullCypraea.Conusea, testCypraea.Conusea)

Cypraea.Conuseaeffects<-as.data.frame(fixef(testCypraea.Conusea, probs = c(0.1,0.9)))
Cypraea.Conuseaeffects$community<-row.names(Cypraea.Conuseaeffects)
colnames(Cypraea.Conuseaeffects)<-c("estimate","se","conf.low","conf.high","community")
Cypraea.Conuseaeff <- ggplot(Cypraea.Conuseaeffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Cypraea.Conus ea")

##Cypraea.Conus harv
nullCypraea.Conusharv<-brm(Cypraea.Conus_harv_change ~ 1 +(1|Village),
                           data=seafood[!is.na(seafood$Cypraea.Conus_harv_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testCypraea.Conusharv <-brm(Cypraea.Conus_harv_change ~ 
                              Gender+
                              VillageYears_standard+
                              YearsFishing_standard+
                              Age_standard +(1|Village),
                            data=seafood[!is.na(seafood$Cypraea.Conus_harv_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Cypraea.Conusharv <- pp_check(testCypraea.Conusharv,ndraws = 50)+ggtitle("Model fit Cypraea.Conus Harvest")

loo(nullCypraea.Conusharv, testCypraea.Conusharv)

Cypraea.Conusharveffects<-as.data.frame(fixef(testCypraea.Conusharv, probs = c(0.1,0.9)))
Cypraea.Conusharveffects$community<-row.names(Cypraea.Conusharveffects)
colnames(Cypraea.Conusharveffects)<-c("estimate","se","conf.low","conf.high","community")
Cypraea.Conusharveff <- ggplot(Cypraea.Conusharveffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Cypraea.Conus harv")

##Export graphs
jpeg("PPcheck_Cypraea.Conus", width =3500, height = 2000, res=300)
ggarrange(Cypraea.Conuscons, Cypraea.Conusea, Cypraea.Conusharv)
dev.off()

jpeg("FixedEffects_Cypraea.Conus", width =3500, height = 2000, res=300)
ggarrange(Cypraea.Conusconseff, Cypraea.Conuseaeff, Cypraea.Conusharveff)
dev.off()








##Pyrasus
seafood$Pyrasus_cons_change <- seafood$Pyrasus_cons_now - seafood$Pyrasus_cons_10
seafood$Pyrasus_ea_change <- seafood$Pyrasus_ea_now - seafood$Pyrasus_ea_10
seafood$Pyrasus_harv_change <- seafood$Pyrasus_harv_now - seafood$Pyrasus_harv_10

seafood$Pyrasus_cons_change <- ordered(as.factor(seafood$Pyrasus_cons_change), 
                                       levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

seafood$Pyrasus_ea_change <- ordered(as.factor(seafood$Pyrasus_ea_change), 
                                     levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

seafood$Pyrasus_harv_change <- ordered(as.factor(seafood$Pyrasus_harv_change), 
                                       levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

##Pyrasus Cons
nullPyrasuscons<-brm(Pyrasus_cons_change ~ 1 +(1|Village),
                     data=seafood[!is.na(seafood$Pyrasus_cons_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testPyrasuscons <-brm(Pyrasus_cons_change ~ 
                        Gender+
                        VillageYears_standard+
                        YearsFishing_standard+
                        Age_standard +(1|Village),
                      data=seafood[!is.na(seafood$Pyrasus_cons_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Pyrasuscons <- pp_check(testPyrasuscons,ndraws = 50)+ggtitle("Model fit Pyrasus Consumption")

loo(nullPyrasuscons, testPyrasuscons)

Pyrasusconseffects<-as.data.frame(fixef(testPyrasuscons, probs = c(0.1,0.9)))
Pyrasusconseffects$community<-row.names(Pyrasusconseffects)
colnames(Pyrasusconseffects)<-c("estimate","se","conf.low","conf.high","community")
Pyrasusconseff <- ggplot(Pyrasusconseffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Pyrasus Cons")

##Pyrasus ea
nullPyrasusea<-brm(Pyrasus_ea_change ~ 1 +(1|Village),
                   data=seafood[!is.na(seafood$Pyrasus_ea_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testPyrasusea <-brm(Pyrasus_ea_change ~ 
                      Gender+
                      VillageYears_standard+
                      YearsFishing_standard+
                      Age_standard +(1|Village),
                    data=seafood[!is.na(seafood$Pyrasus_ea_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Pyrasusea <- pp_check(testPyrasusea,ndraws = 50)+ggtitle("Model fit Pyrasus Abundance")

loo(nullPyrasusea, testPyrasusea)

Pyrasuseaeffects<-as.data.frame(fixef(testPyrasusea, probs = c(0.1,0.9)))
Pyrasuseaeffects$community<-row.names(Pyrasuseaeffects)
colnames(Pyrasuseaeffects)<-c("estimate","se","conf.low","conf.high","community")
Pyrasuseaeff <- ggplot(Pyrasuseaeffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Pyrasus ea")

##Pyrasus harv
nullPyrasusharv<-brm(Pyrasus_harv_change ~ 1 +(1|Village),
                     data=seafood[!is.na(seafood$Pyrasus_harv_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testPyrasusharv <-brm(Pyrasus_harv_change ~ 
                        Gender+
                        VillageYears_standard+
                        YearsFishing_standard+
                        Age_standard +(1|Village),
                      data=seafood[!is.na(seafood$Pyrasus_harv_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Pyrasusharv <- pp_check(testPyrasusharv,ndraws = 50)+ggtitle("Model fit Pyrasus Harvest")

loo(nullPyrasusharv, testPyrasusharv)

Pyrasusharveffects<-as.data.frame(fixef(testPyrasusharv, probs = c(0.1,0.9)))
Pyrasusharveffects$community<-row.names(Pyrasusharveffects)
colnames(Pyrasusharveffects)<-c("estimate","se","conf.low","conf.high","community")
Pyrasusharveff <- ggplot(Pyrasusharveffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Pyrasus harv")

##Export graphs
jpeg("PPcheck_Pyrasus", width =3500, height = 2000, res=300)
ggarrange(Pyrasuscons, Pyrasusea, Pyrasusharv)
dev.off()

jpeg("FixedEffects_Pyrasus", width =3500, height = 2000, res=300)
ggarrange(Pyrasusconseff, Pyrasuseaeff, Pyrasusharveff)
dev.off()




##Charonia
seafood$Charonia_cons_change <- seafood$Charonia_cons_now - seafood$Charonia_cons_10
seafood$Charonia_ea_change <- seafood$Charonia_ea_now - seafood$Charonia_ea_10
seafood$Charonia_harv_change <- seafood$Charonia_harv_now - seafood$Charonia_harv_10

seafood$Charonia_cons_change <- ordered(as.factor(seafood$Charonia_cons_change), 
                                        levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

seafood$Charonia_ea_change <- ordered(as.factor(seafood$Charonia_ea_change), 
                                      levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

seafood$Charonia_harv_change <- ordered(as.factor(seafood$Charonia_harv_change), 
                                        levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

##Charonia Cons
nullCharoniacons<-brm(Charonia_cons_change ~ 1 +(1|Village),
                      data=seafood[!is.na(seafood$Charonia_cons_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testCharoniacons <-brm(Charonia_cons_change ~ 
                         Gender+
                         VillageYears_standard+
                         YearsFishing_standard+
                         Age_standard +(1|Village),
                       data=seafood[!is.na(seafood$Charonia_cons_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Charoniacons <- pp_check(testCharoniacons,ndraws = 50)+ggtitle("Model fit Charonia Consumption")

loo(nullCharoniacons, testCharoniacons)

Charoniaconseffects<-as.data.frame(fixef(testCharoniacons, probs = c(0.1,0.9)))
Charoniaconseffects$community<-row.names(Charoniaconseffects)
colnames(Charoniaconseffects)<-c("estimate","se","conf.low","conf.high","community")
Charoniaconseff <- ggplot(Charoniaconseffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Charonia Cons")

##Charonia ea
nullCharoniaea<-brm(Charonia_ea_change ~ 1 +(1|Village),
                    data=seafood[!is.na(seafood$Charonia_ea_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testCharoniaea <-brm(Charonia_ea_change ~ 
                       Gender+
                       VillageYears_standard+
                       YearsFishing_standard+
                       Age_standard +(1|Village),
                     data=seafood[!is.na(seafood$Charonia_ea_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Charoniaea <- pp_check(testCharoniaea,ndraws = 50)+ggtitle("Model fit Charonia Abundance")

loo(nullCharoniaea, testCharoniaea)

Charoniaeaeffects<-as.data.frame(fixef(testCharoniaea, probs = c(0.1,0.9)))
Charoniaeaeffects$community<-row.names(Charoniaeaeffects)
colnames(Charoniaeaeffects)<-c("estimate","se","conf.low","conf.high","community")
Charoniaeaeff <- ggplot(Charoniaeaeffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Charonia ea")

##Charonia harv
nullCharoniaharv<-brm(Charonia_harv_change ~ 1 +(1|Village),
                      data=seafood[!is.na(seafood$Charonia_harv_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testCharoniaharv <-brm(Charonia_harv_change ~ 
                         Gender+
                         VillageYears_standard+
                         YearsFishing_standard+
                         Age_standard +(1|Village),
                       data=seafood[!is.na(seafood$Charonia_harv_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Charoniaharv <- pp_check(testCharoniaharv,ndraws = 50)+ggtitle("Model fit Charonia Harvest")

loo(nullCharoniaharv, testCharoniaharv)

Charoniaharveffects<-as.data.frame(fixef(testCharoniaharv, probs = c(0.1,0.9)))
Charoniaharveffects$community<-row.names(Charoniaharveffects)
colnames(Charoniaharveffects)<-c("estimate","se","conf.low","conf.high","community")
Charoniaharveff <- ggplot(Charoniaharveffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Charonia harv")

##Export graphs
jpeg("PPcheck_Charonia", width =3500, height = 2000, res=300)
ggarrange(Charoniacons, Charoniaea, Charoniaharv)
dev.off()

jpeg("FixedEffects_Charonia", width =3500, height = 2000, res=300)
ggarrange(Charoniaconseff, Charoniaeaeff, Charoniaharveff)
dev.off()





##Murex.Fasciolaria
seafood$Murex.Fasciolaria_cons_change <- seafood$Murex.Fasciolaria_cons_now - seafood$Murex.Fasciolaria_cons_10
seafood$Murex.Fasciolaria_ea_change <- seafood$Murex.Fasciolaria_ea_now - seafood$Murex.Fasciolaria_ea_10
seafood$Murex.Fasciolaria_harv_change <- seafood$Murex.Fasciolaria_harv_now - seafood$Murex.Fasciolaria_harv_10

seafood$Murex.Fasciolaria_cons_change <- ordered(as.factor(seafood$Murex.Fasciolaria_cons_change), 
                                                 levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

seafood$Murex.Fasciolaria_ea_change <- ordered(as.factor(seafood$Murex.Fasciolaria_ea_change), 
                                               levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

seafood$Murex.Fasciolaria_harv_change <- ordered(as.factor(seafood$Murex.Fasciolaria_harv_change), 
                                                 levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

##Murex.Fasciolaria Cons
nullMurex.Fasciolariacons<-brm(Murex.Fasciolaria_cons_change ~ 1 +(1|Village),
                               data=seafood[!is.na(seafood$Murex.Fasciolaria_cons_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testMurex.Fasciolariacons <-brm(Murex.Fasciolaria_cons_change ~ 
                                  Gender+
                                  VillageYears_standard+
                                  YearsFishing_standard+
                                  Age_standard +(1|Village),
                                data=seafood[!is.na(seafood$Murex.Fasciolaria_cons_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Murex.Fasciolariacons <- pp_check(testMurex.Fasciolariacons,ndraws = 50)+ggtitle("Model fit Murex.Fasciolaria Consumption")

loo(nullMurex.Fasciolariacons, testMurex.Fasciolariacons)

Murex.Fasciolariaconseffects<-as.data.frame(fixef(testMurex.Fasciolariacons, probs = c(0.1,0.9)))
Murex.Fasciolariaconseffects$community<-row.names(Murex.Fasciolariaconseffects)
colnames(Murex.Fasciolariaconseffects)<-c("estimate","se","conf.low","conf.high","community")
Murex.Fasciolariaconseff <- ggplot(Murex.Fasciolariaconseffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Murex.Fasciolaria Cons")

##Murex.Fasciolaria ea
nullMurex.Fasciolariaea<-brm(Murex.Fasciolaria_ea_change ~ 1 +(1|Village),
                             data=seafood[!is.na(seafood$Murex.Fasciolaria_ea_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testMurex.Fasciolariaea <-brm(Murex.Fasciolaria_ea_change ~ 
                                Gender+
                                VillageYears_standard+
                                YearsFishing_standard+
                                Age_standard +(1|Village),
                              data=seafood[!is.na(seafood$Murex.Fasciolaria_ea_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Murex.Fasciolariaea <- pp_check(testMurex.Fasciolariaea,ndraws = 50)+ggtitle("Model fit Murex.Fasciolaria Abundance")

loo(nullMurex.Fasciolariaea, testMurex.Fasciolariaea)

Murex.Fasciolariaeaeffects<-as.data.frame(fixef(testMurex.Fasciolariaea, probs = c(0.1,0.9)))
Murex.Fasciolariaeaeffects$community<-row.names(Murex.Fasciolariaeaeffects)
colnames(Murex.Fasciolariaeaeffects)<-c("estimate","se","conf.low","conf.high","community")
Murex.Fasciolariaeaeff <- ggplot(Murex.Fasciolariaeaeffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Murex.Fasciolaria ea")

##Murex.Fasciolaria harv
nullMurex.Fasciolariaharv<-brm(Murex.Fasciolaria_harv_change ~ 1 +(1|Village),
                               data=seafood[!is.na(seafood$Murex.Fasciolaria_harv_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testMurex.Fasciolariaharv <-brm(Murex.Fasciolaria_harv_change ~ 
                                  Gender+
                                  VillageYears_standard+
                                  YearsFishing_standard+
                                  Age_standard +(1|Village),
                                data=seafood[!is.na(seafood$Murex.Fasciolaria_harv_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Murex.Fasciolariaharv <- pp_check(testMurex.Fasciolariaharv,ndraws = 50)+ggtitle("Model fit Murex.Fasciolaria Harvest")

loo(nullMurex.Fasciolariaharv, testMurex.Fasciolariaharv)

Murex.Fasciolariaharveffects<-as.data.frame(fixef(testMurex.Fasciolariaharv, probs = c(0.1,0.9)))
Murex.Fasciolariaharveffects$community<-row.names(Murex.Fasciolariaharveffects)
colnames(Murex.Fasciolariaharveffects)<-c("estimate","se","conf.low","conf.high","community")
Murex.Fasciolariaharveff <- ggplot(Murex.Fasciolariaharveffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Murex.Fasciolaria harv")

##Export graphs
jpeg("PPcheck_Murex.Fasciolaria", width =3500, height = 2000, res=300)
ggarrange(Murex.Fasciolariacons, Murex.Fasciolariaea, Murex.Fasciolariaharv)
dev.off()

jpeg("FixedEffects_Murex.Fasciolaria", width =3500, height = 2000, res=300)
ggarrange(Murex.Fasciolariaconseff, Murex.Fasciolariaeaeff, Murex.Fasciolariaharveff)
dev.off()








##Lambis
seafood$Lambis_cons_change <- seafood$Lambis_cons_now - seafood$Lambis_cons_10
seafood$Lambis_ea_change <- seafood$Lambis_ea_now - seafood$Lambis_ea_10
seafood$Lambis_harv_change <- seafood$Lambis_harv_now - seafood$Lambis_harv_10

seafood$Lambis_cons_change <- ordered(as.factor(seafood$Lambis_cons_change), 
                                      levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

seafood$Lambis_ea_change <- ordered(as.factor(seafood$Lambis_ea_change), 
                                    levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

seafood$Lambis_harv_change <- ordered(as.factor(seafood$Lambis_harv_change), 
                                      levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

##Lambis Cons
nullLambiscons<-brm(Lambis_cons_change ~ 1 +(1|Village),
                    data=seafood[!is.na(seafood$Lambis_cons_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testLambiscons <-brm(Lambis_cons_change ~ 
                       Gender+
                       VillageYears_standard+
                       YearsFishing_standard+
                       Age_standard +(1|Village),
                     data=seafood[!is.na(seafood$Lambis_cons_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Lambiscons <- pp_check(testLambiscons,ndraws = 50)+ggtitle("Model fit Lambis Consumption")

loo(nullLambiscons, testLambiscons)

Lambisconseffects<-as.data.frame(fixef(testLambiscons, probs = c(0.1,0.9)))
Lambisconseffects$community<-row.names(Lambisconseffects)
colnames(Lambisconseffects)<-c("estimate","se","conf.low","conf.high","community")
Lambisconseff <- ggplot(Lambisconseffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Lambis Cons")

##Lambis ea
nullLambisea<-brm(Lambis_ea_change ~ 1 +(1|Village),
                  data=seafood[!is.na(seafood$Lambis_ea_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testLambisea <-brm(Lambis_ea_change ~ 
                     Gender+
                     VillageYears_standard+
                     YearsFishing_standard+
                     Age_standard +(1|Village),
                   data=seafood[!is.na(seafood$Lambis_ea_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Lambisea <- pp_check(testLambisea,ndraws = 50)+ggtitle("Model fit Lambis Abundance")

loo(nullLambisea, testLambisea)

Lambiseaeffects<-as.data.frame(fixef(testLambisea, probs = c(0.1,0.9)))
Lambiseaeffects$community<-row.names(Lambiseaeffects)
colnames(Lambiseaeffects)<-c("estimate","se","conf.low","conf.high","community")
Lambiseaeff <- ggplot(Lambiseaeffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Lambis ea")

##Lambis harv
nullLambisharv<-brm(Lambis_harv_change ~ 1 +(1|Village),
                    data=seafood[!is.na(seafood$Lambis_harv_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testLambisharv <-brm(Lambis_harv_change ~ 
                       Gender+
                       VillageYears_standard+
                       YearsFishing_standard+
                       Age_standard +(1|Village),
                     data=seafood[!is.na(seafood$Lambis_harv_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Lambisharv <- pp_check(testLambisharv,ndraws = 50)+ggtitle("Model fit Lambis Harvest")

loo(nullLambisharv, testLambisharv)

Lambisharveffects<-as.data.frame(fixef(testLambisharv, probs = c(0.1,0.9)))
Lambisharveffects$community<-row.names(Lambisharveffects)
colnames(Lambisharveffects)<-c("estimate","se","conf.low","conf.high","community")
Lambisharveff <- ggplot(Lambisharveffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Lambis harv")

##Export graphs
jpeg("PPcheck_Lambis", width =3500, height = 2000, res=300)
ggarrange(Lambiscons, Lambisea, Lambisharv)
dev.off()

jpeg("FixedEffects_Lambis", width =3500, height = 2000, res=300)
ggarrange(Lambisconseff, Lambiseaeff, Lambisharveff)
dev.off()





##Anadara
seafood$Anadara_cons_change <- seafood$Anadara_cons_now - seafood$Anadara_cons_10
seafood$Anadara_ea_change <- seafood$Anadara_ea_now - seafood$Anadara_ea_10
seafood$Anadara_harv_change <- seafood$Anadara_harv_now - seafood$Anadara_harv_10

seafood$Anadara_cons_change <- ordered(as.factor(seafood$Anadara_cons_change), 
                                       levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

seafood$Anadara_ea_change <- ordered(as.factor(seafood$Anadara_ea_change), 
                                     levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

seafood$Anadara_harv_change <- ordered(as.factor(seafood$Anadara_harv_change), 
                                       levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

##Anadara Cons
nullAnadaracons<-brm(Anadara_cons_change ~ 1 +(1|Village),
                     data=seafood[!is.na(seafood$Anadara_cons_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testAnadaracons <-brm(Anadara_cons_change ~ 
                        Gender+
                        VillageYears_standard+
                        YearsFishing_standard+
                        Age_standard +(1|Village),
                      data=seafood[!is.na(seafood$Anadara_cons_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Anadaracons <- pp_check(testAnadaracons,ndraws = 50)+ggtitle("Model fit Anadara Consumption")

loo(nullAnadaracons, testAnadaracons)

Anadaraconseffects<-as.data.frame(fixef(testAnadaracons, probs = c(0.1,0.9)))
Anadaraconseffects$community<-row.names(Anadaraconseffects)
colnames(Anadaraconseffects)<-c("estimate","se","conf.low","conf.high","community")
Anadaraconseff <- ggplot(Anadaraconseffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Anadara Cons")

##Anadara ea
nullAnadaraea<-brm(Anadara_ea_change ~ 1 +(1|Village),
                   data=seafood[!is.na(seafood$Anadara_ea_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testAnadaraea <-brm(Anadara_ea_change ~ 
                      Gender+
                      VillageYears_standard+
                      YearsFishing_standard+
                      Age_standard +(1|Village),
                    data=seafood[!is.na(seafood$Anadara_ea_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Anadaraea <- pp_check(testAnadaraea,ndraws = 50)+ggtitle("Model fit Anadara Abundance")

loo(nullAnadaraea, testAnadaraea)

Anadaraeaeffects<-as.data.frame(fixef(testAnadaraea, probs = c(0.1,0.9)))
Anadaraeaeffects$community<-row.names(Anadaraeaeffects)
colnames(Anadaraeaeffects)<-c("estimate","se","conf.low","conf.high","community")
Anadaraeaeff <- ggplot(Anadaraeaeffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Anadara ea")

##Anadara harv
nullAnadaraharv<-brm(Anadara_harv_change ~ 1 +(1|Village),
                     data=seafood[!is.na(seafood$Anadara_harv_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testAnadaraharv <-brm(Anadara_harv_change ~ 
                        Gender+
                        VillageYears_standard+
                        YearsFishing_standard+
                        Age_standard +(1|Village),
                      data=seafood[!is.na(seafood$Anadara_harv_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Anadaraharv <- pp_check(testAnadaraharv,ndraws = 50)+ggtitle("Model fit Anadara Harvest")

loo(nullAnadaraharv, testAnadaraharv)

Anadaraharveffects<-as.data.frame(fixef(testAnadaraharv, probs = c(0.1,0.9)))
Anadaraharveffects$community<-row.names(Anadaraharveffects)
colnames(Anadaraharveffects)<-c("estimate","se","conf.low","conf.high","community")
Anadaraharveff <- ggplot(Anadaraharveffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Anadara harv")

##Export graphs
jpeg("PPcheck_Anadara", width =3500, height = 2000, res=300)
ggarrange(Anadaracons, Anadaraea, Anadaraharv)
dev.off()

jpeg("FixedEffects_Anadara", width =3500, height = 2000, res=300)
ggarrange(Anadaraconseff, Anadaraeaeff, Anadaraharveff)
dev.off()






##Pinctada.Isognomon.Atrina.Pinna
seafood$Pinctada.Isognomon.Atrina.Pinna_cons_change <- seafood$Pinctada.Isognomon.Atrina.Pinna_cons_now - seafood$Pinctada.Isognomon.Atrina.Pinna_cons_10
seafood$Pinctada.Isognomon.Atrina.Pinna_ea_change <- seafood$Pinctada.Isognomon.Atrina.Pinna_ea_now - seafood$Pinctada.Isognomon.Atrina.Pinna_ea_10
seafood$Pinctada.Isognomon.Atrina.Pinna_harv_change <- seafood$Pinctada.Isognomon.Atrina.Pinna_harv_now - seafood$Pinctada.Isognomon.Atrina.Pinna_harv_10

seafood$Pinctada.Isognomon.Atrina.Pinna_cons_change <- ordered(as.factor(seafood$Pinctada.Isognomon.Atrina.Pinna_cons_change), 
                                                               levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

seafood$Pinctada.Isognomon.Atrina.Pinna_ea_change <- ordered(as.factor(seafood$Pinctada.Isognomon.Atrina.Pinna_ea_change), 
                                                             levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

seafood$Pinctada.Isognomon.Atrina.Pinna_harv_change <- ordered(as.factor(seafood$Pinctada.Isognomon.Atrina.Pinna_harv_change), 
                                                               levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

##Pinctada.Isognomon.Atrina.Pinna Cons
nullPinctada.Isognomon.Atrina.Pinnacons<-brm(Pinctada.Isognomon.Atrina.Pinna_cons_change ~ 1 +(1|Village),
                                             data=seafood[!is.na(seafood$Pinctada.Isognomon.Atrina.Pinna_cons_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testPinctada.Isognomon.Atrina.Pinnacons <-brm(Pinctada.Isognomon.Atrina.Pinna_cons_change ~ 
                                                Gender+
                                                VillageYears_standard+
                                                YearsFishing_standard+
                                                Age_standard +(1|Village),
                                              data=seafood[!is.na(seafood$Pinctada.Isognomon.Atrina.Pinna_cons_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Pinctada.Isognomon.Atrina.Pinnacons <- pp_check(testPinctada.Isognomon.Atrina.Pinnacons,ndraws = 50)+ggtitle("Model fit Pinctada.Isognomon.Atrina.Pinna Consumption")

loo(nullPinctada.Isognomon.Atrina.Pinnacons, testPinctada.Isognomon.Atrina.Pinnacons)

Pinctada.Isognomon.Atrina.Pinnaconseffects<-as.data.frame(fixef(testPinctada.Isognomon.Atrina.Pinnacons, probs = c(0.1,0.9)))
Pinctada.Isognomon.Atrina.Pinnaconseffects$community<-row.names(Pinctada.Isognomon.Atrina.Pinnaconseffects)
colnames(Pinctada.Isognomon.Atrina.Pinnaconseffects)<-c("estimate","se","conf.low","conf.high","community")
Pinctada.Isognomon.Atrina.Pinnaconseff <- ggplot(Pinctada.Isognomon.Atrina.Pinnaconseffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Pinctada.Isognomon.Atrina.Pinna Cons")

##Pinctada.Isognomon.Atrina.Pinna ea
nullPinctada.Isognomon.Atrina.Pinnaea<-brm(Pinctada.Isognomon.Atrina.Pinna_ea_change ~ 1 +(1|Village),
                                           data=seafood[!is.na(seafood$Pinctada.Isognomon.Atrina.Pinna_ea_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testPinctada.Isognomon.Atrina.Pinnaea <-brm(Pinctada.Isognomon.Atrina.Pinna_ea_change ~ 
                                              Gender+
                                              VillageYears_standard+
                                              YearsFishing_standard+
                                              Age_standard +(1|Village),
                                            data=seafood[!is.na(seafood$Pinctada.Isognomon.Atrina.Pinna_ea_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Pinctada.Isognomon.Atrina.Pinnaea <- pp_check(testPinctada.Isognomon.Atrina.Pinnaea,ndraws = 50)+ggtitle("Model fit Pinctada.Isognomon.Atrina.Pinna Abundance")

loo(nullPinctada.Isognomon.Atrina.Pinnaea, testPinctada.Isognomon.Atrina.Pinnaea)

Pinctada.Isognomon.Atrina.Pinnaeaeffects<-as.data.frame(fixef(testPinctada.Isognomon.Atrina.Pinnaea, probs = c(0.1,0.9)))
Pinctada.Isognomon.Atrina.Pinnaeaeffects$community<-row.names(Pinctada.Isognomon.Atrina.Pinnaeaeffects)
colnames(Pinctada.Isognomon.Atrina.Pinnaeaeffects)<-c("estimate","se","conf.low","conf.high","community")
Pinctada.Isognomon.Atrina.Pinnaeaeff <- ggplot(Pinctada.Isognomon.Atrina.Pinnaeaeffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Pinctada.Isognomon.Atrina.Pinna ea")

##Pinctada.Isognomon.Atrina.Pinna harv
nullPinctada.Isognomon.Atrina.Pinnaharv<-brm(Pinctada.Isognomon.Atrina.Pinna_harv_change ~ 1 +(1|Village),
                                             data=seafood[!is.na(seafood$Pinctada.Isognomon.Atrina.Pinna_harv_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testPinctada.Isognomon.Atrina.Pinnaharv <-brm(Pinctada.Isognomon.Atrina.Pinna_harv_change ~ 
                                                Gender+
                                                VillageYears_standard+
                                                YearsFishing_standard+
                                                Age_standard +(1|Village),
                                              data=seafood[!is.na(seafood$Pinctada.Isognomon.Atrina.Pinna_harv_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Pinctada.Isognomon.Atrina.Pinnaharv <- pp_check(testPinctada.Isognomon.Atrina.Pinnaharv,ndraws = 50)+ggtitle("Model fit Pinctada.Isognomon.Atrina.Pinna Harvest")

loo(nullPinctada.Isognomon.Atrina.Pinnaharv, testPinctada.Isognomon.Atrina.Pinnaharv)

Pinctada.Isognomon.Atrina.Pinnaharveffects<-as.data.frame(fixef(testPinctada.Isognomon.Atrina.Pinnaharv, probs = c(0.1,0.9)))
Pinctada.Isognomon.Atrina.Pinnaharveffects$community<-row.names(Pinctada.Isognomon.Atrina.Pinnaharveffects)
colnames(Pinctada.Isognomon.Atrina.Pinnaharveffects)<-c("estimate","se","conf.low","conf.high","community")
Pinctada.Isognomon.Atrina.Pinnaharveff <- ggplot(Pinctada.Isognomon.Atrina.Pinnaharveffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Pinctada.Isognomon.Atrina.Pinna harv")

##Export graphs
jpeg("PPcheck_Pinctada.Isognomon.Atrina.Pinna", width =3500, height = 2000, res=300)
ggarrange(Pinctada.Isognomon.Atrina.Pinnacons, Pinctada.Isognomon.Atrina.Pinnaea, Pinctada.Isognomon.Atrina.Pinnaharv)
dev.off()

jpeg("FixedEffects_Pinctada.Isognomon.Atrina.Pinna", width =3500, height = 2000, res=300)
ggarrange(Pinctada.Isognomon.Atrina.Pinnaconseff, Pinctada.Isognomon.Atrina.Pinnaeaeff, Pinctada.Isognomon.Atrina.Pinnaharveff)
dev.off()







##Tridacna
seafood$Tridacna_cons_change <- seafood$Tridacna_cons_now - seafood$Tridacna_cons_10
seafood$Tridacna_ea_change <- seafood$Tridacna_ea_now - seafood$Tridacna_ea_10
seafood$Tridacna_harv_change <- seafood$Tridacna_harv_now - seafood$Tridacna_harv_10

seafood$Tridacna_cons_change <- ordered(as.factor(seafood$Tridacna_cons_change), 
                                        levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

seafood$Tridacna_ea_change <- ordered(as.factor(seafood$Tridacna_ea_change), 
                                      levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

seafood$Tridacna_harv_change <- ordered(as.factor(seafood$Tridacna_harv_change), 
                                        levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

##Tridacna Cons
nullTridacnacons<-brm(Tridacna_cons_change ~ 1 +(1|Village),
                      data=seafood[!is.na(seafood$Tridacna_cons_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testTridacnacons <-brm(Tridacna_cons_change ~ 
                         Gender+
                         VillageYears_standard+
                         YearsFishing_standard+
                         Age_standard +(1|Village),
                       data=seafood[!is.na(seafood$Tridacna_cons_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Tridacnacons <- pp_check(testTridacnacons,ndraws = 50)+ggtitle("Model fit Tridacna Consumption")

loo(nullTridacnacons, testTridacnacons)

Tridacnaconseffects<-as.data.frame(fixef(testTridacnacons, probs = c(0.1,0.9)))
Tridacnaconseffects$community<-row.names(Tridacnaconseffects)
colnames(Tridacnaconseffects)<-c("estimate","se","conf.low","conf.high","community")
Tridacnaconseff <- ggplot(Tridacnaconseffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Tridacna Cons")

##Tridacna ea
nullTridacnaea<-brm(Tridacna_ea_change ~ 1 +(1|Village),
                    data=seafood[!is.na(seafood$Tridacna_ea_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testTridacnaea <-brm(Tridacna_ea_change ~ 
                       Gender+
                       VillageYears_standard+
                       YearsFishing_standard+
                       Age_standard +(1|Village),
                     data=seafood[!is.na(seafood$Tridacna_ea_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Tridacnaea <- pp_check(testTridacnaea,ndraws = 50)+ggtitle("Model fit Tridacna Abundance")

loo(nullTridacnaea, testTridacnaea)

Tridacnaeaeffects<-as.data.frame(fixef(testTridacnaea, probs = c(0.1,0.9)))
Tridacnaeaeffects$community<-row.names(Tridacnaeaeffects)
colnames(Tridacnaeaeffects)<-c("estimate","se","conf.low","conf.high","community")
Tridacnaeaeff <- ggplot(Tridacnaeaeffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Tridacna ea")

##Tridacna harv
nullTridacnaharv<-brm(Tridacna_harv_change ~ 1 +(1|Village),
                      data=seafood[!is.na(seafood$Tridacna_harv_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testTridacnaharv <-brm(Tridacna_harv_change ~ 
                         Gender+
                         VillageYears_standard+
                         YearsFishing_standard+
                         Age_standard +(1|Village),
                       data=seafood[!is.na(seafood$Tridacna_harv_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Tridacnaharv <- pp_check(testTridacnaharv,ndraws = 50)+ggtitle("Model fit Tridacna Harvest")

loo(nullTridacnaharv, testTridacnaharv)

Tridacnaharveffects<-as.data.frame(fixef(testTridacnaharv, probs = c(0.1,0.9)))
Tridacnaharveffects$community<-row.names(Tridacnaharveffects)
colnames(Tridacnaharveffects)<-c("estimate","se","conf.low","conf.high","community")
Tridacnaharveff <- ggplot(Tridacnaharveffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Tridacna harv")

##Export graphs
jpeg("PPcheck_Tridacna", width =3500, height = 2000, res=300)
ggarrange(Tridacnacons, Tridacnaea, Tridacnaharv)
dev.off()

jpeg("FixedEffects_Tridacna", width =3500, height = 2000, res=300)
ggarrange(Tridacnaconseff, Tridacnaeaeff, Tridacnaharveff)
dev.off()





##Tripneustes
seafood$Tripneustes_cons_change <- seafood$Tripneustes_cons_now - seafood$Tripneustes_cons_10
seafood$Tripneustes_ea_change <- seafood$Tripneustes_ea_now - seafood$Tripneustes_ea_10
seafood$Tripneustes_harv_change <- seafood$Tripneustes_harv_now - seafood$Tripneustes_harv_10

seafood$Tripneustes_cons_change <- ordered(as.factor(seafood$Tripneustes_cons_change), 
                                           levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

seafood$Tripneustes_ea_change <- ordered(as.factor(seafood$Tripneustes_ea_change), 
                                         levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

seafood$Tripneustes_harv_change <- ordered(as.factor(seafood$Tripneustes_harv_change), 
                                           levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

##Tripneustes Cons
nullTripneustescons<-brm(Tripneustes_cons_change ~ 1 +(1|Village),
                         data=seafood[!is.na(seafood$Tripneustes_cons_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testTripneustescons <-brm(Tripneustes_cons_change ~ 
                            Gender+
                            VillageYears_standard+
                            YearsFishing_standard+
                            Age_standard +(1|Village),
                          data=seafood[!is.na(seafood$Tripneustes_cons_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Tripneustescons <- pp_check(testTripneustescons,ndraws = 50)+ggtitle("Model fit Tripneustes Consumption")

loo(nullTripneustescons, testTripneustescons)

Tripneustesconseffects<-as.data.frame(fixef(testTripneustescons, probs = c(0.1,0.9)))
Tripneustesconseffects$community<-row.names(Tripneustesconseffects)
colnames(Tripneustesconseffects)<-c("estimate","se","conf.low","conf.high","community")
Tripneustesconseff <- ggplot(Tripneustesconseffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Tripneustes Cons")

##Tripneustes ea
nullTripneustesea<-brm(Tripneustes_ea_change ~ 1 +(1|Village),
                       data=seafood[!is.na(seafood$Tripneustes_ea_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testTripneustesea <-brm(Tripneustes_ea_change ~ 
                          Gender+
                          VillageYears_standard+
                          YearsFishing_standard+
                          Age_standard +(1|Village),
                        data=seafood[!is.na(seafood$Tripneustes_ea_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Tripneustesea <- pp_check(testTripneustesea,ndraws = 50)+ggtitle("Model fit Tripneustes Abundance")

loo(nullTripneustesea, testTripneustesea)

Tripneusteseaeffects<-as.data.frame(fixef(testTripneustesea, probs = c(0.1,0.9)))
Tripneusteseaeffects$community<-row.names(Tripneusteseaeffects)
colnames(Tripneusteseaeffects)<-c("estimate","se","conf.low","conf.high","community")
Tripneusteseaeff <- ggplot(Tripneusteseaeffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Tripneustes ea")

##Tripneustes harv
nullTripneustesharv<-brm(Tripneustes_harv_change ~ 1 +(1|Village),
                         data=seafood[!is.na(seafood$Tripneustes_harv_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testTripneustesharv <-brm(Tripneustes_harv_change ~ 
                            Gender+
                            VillageYears_standard+
                            YearsFishing_standard+
                            Age_standard +(1|Village),
                          data=seafood[!is.na(seafood$Tripneustes_harv_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Tripneustesharv <- pp_check(testTripneustesharv,ndraws = 50)+ggtitle("Model fit Tripneustes Harvest")

loo(nullTripneustesharv, testTripneustesharv)

Tripneustesharveffects<-as.data.frame(fixef(testTripneustesharv, probs = c(0.1,0.9)))
Tripneustesharveffects$community<-row.names(Tripneustesharveffects)
colnames(Tripneustesharveffects)<-c("estimate","se","conf.low","conf.high","community")
Tripneustesharveff <- ggplot(Tripneustesharveffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Tripneustes harv")

##Export graphs
jpeg("PPcheck_Tripneustes", width =3500, height = 2000, res=300)
ggarrange(Tripneustescons, Tripneustesea, Tripneustesharv)
dev.off()

jpeg("FixedEffects_Tripneustes", width =3500, height = 2000, res=300)
ggarrange(Tripneustesconseff, Tripneusteseaeff, Tripneustesharveff)
dev.off()








##Aristeidae
seafood$Aristeidae_cons_change <- seafood$Aristeidae_cons_now - seafood$Aristeidae_cons_10
seafood$Aristeidae_ea_change <- seafood$Aristeidae_ea_now - seafood$Aristeidae_ea_10
seafood$Aristeidae_harv_change <- seafood$Aristeidae_harv_now - seafood$Aristeidae_harv_10

seafood$Aristeidae_cons_change <- ordered(as.factor(seafood$Aristeidae_cons_change), 
                                          levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

seafood$Aristeidae_ea_change <- ordered(as.factor(seafood$Aristeidae_ea_change), 
                                        levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

seafood$Aristeidae_harv_change <- ordered(as.factor(seafood$Aristeidae_harv_change), 
                                          levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

##Aristeidae Cons
nullAristeidaecons<-brm(Aristeidae_cons_change ~ 1 +(1|Village),
                        data=seafood[!is.na(seafood$Aristeidae_cons_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testAristeidaecons <-brm(Aristeidae_cons_change ~ 
                           Gender+
                           VillageYears_standard+
                           YearsFishing_standard+
                           Age_standard +(1|Village),
                         data=seafood[!is.na(seafood$Aristeidae_cons_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Aristeidaecons <- pp_check(testAristeidaecons,ndraws = 50)+ggtitle("Model fit Aristeidae Consumption")

loo(nullAristeidaecons, testAristeidaecons)

Aristeidaeconseffects<-as.data.frame(fixef(testAristeidaecons, probs = c(0.1,0.9)))
Aristeidaeconseffects$community<-row.names(Aristeidaeconseffects)
colnames(Aristeidaeconseffects)<-c("estimate","se","conf.low","conf.high","community")
Aristeidaeconseff <- ggplot(Aristeidaeconseffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Aristeidae Cons")

##Aristeidae ea
nullAristeidaeea<-brm(Aristeidae_ea_change ~ 1 +(1|Village),
                      data=seafood[!is.na(seafood$Aristeidae_ea_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testAristeidaeea <-brm(Aristeidae_ea_change ~ 
                         Gender+
                         VillageYears_standard+
                         YearsFishing_standard+
                         Age_standard +(1|Village),
                       data=seafood[!is.na(seafood$Aristeidae_ea_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Aristeidaeea <- pp_check(testAristeidaeea,ndraws = 50)+ggtitle("Model fit Aristeidae Abundance")

loo(nullAristeidaeea, testAristeidaeea)

Aristeidaeeaeffects<-as.data.frame(fixef(testAristeidaeea, probs = c(0.1,0.9)))
Aristeidaeeaeffects$community<-row.names(Aristeidaeeaeffects)
colnames(Aristeidaeeaeffects)<-c("estimate","se","conf.low","conf.high","community")
Aristeidaeeaeff <- ggplot(Aristeidaeeaeffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Aristeidae ea")

##Aristeidae harv
nullAristeidaeharv<-brm(Aristeidae_harv_change ~ 1 +(1|Village),
                        data=seafood[!is.na(seafood$Aristeidae_harv_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testAristeidaeharv <-brm(Aristeidae_harv_change ~ 
                           Gender+
                           VillageYears_standard+
                           YearsFishing_standard+
                           Age_standard +(1|Village),
                         data=seafood[!is.na(seafood$Aristeidae_harv_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Aristeidaeharv <- pp_check(testAristeidaeharv,ndraws = 50)+ggtitle("Model fit Aristeidae Harvest")

loo(nullAristeidaeharv, testAristeidaeharv)

Aristeidaeharveffects<-as.data.frame(fixef(testAristeidaeharv, probs = c(0.1,0.9)))
Aristeidaeharveffects$community<-row.names(Aristeidaeharveffects)
colnames(Aristeidaeharveffects)<-c("estimate","se","conf.low","conf.high","community")
Aristeidaeharveff <- ggplot(Aristeidaeharveffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Aristeidae harv")

##Export graphs
jpeg("PPcheck_Aristeidae", width =3500, height = 2000, res=300)
ggarrange(Aristeidaecons, Aristeidaeea, Aristeidaeharv)
dev.off()

jpeg("FixedEffects_Aristeidae", width =3500, height = 2000, res=300)
ggarrange(Aristeidaeconseff, Aristeidaeeaeff, Aristeidaeharveff)
dev.off()









##Scylla
seafood$Scylla_cons_change <- seafood$Scylla_cons_now - seafood$Scylla_cons_10
seafood$Scylla_ea_change <- seafood$Scylla_ea_now - seafood$Scylla_ea_10
seafood$Scylla_harv_change <- seafood$Scylla_harv_now - seafood$Scylla_harv_10

seafood$Scylla_cons_change <- ordered(as.factor(seafood$Scylla_cons_change), 
                                      levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

seafood$Scylla_ea_change <- ordered(as.factor(seafood$Scylla_ea_change), 
                                    levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

seafood$Scylla_harv_change <- ordered(as.factor(seafood$Scylla_harv_change), 
                                      levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

##Scylla Cons
nullScyllacons<-brm(Scylla_cons_change ~ 1 +(1|Village),
                    data=seafood[!is.na(seafood$Scylla_cons_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testScyllacons <-brm(Scylla_cons_change ~ 
                       Gender+
                       VillageYears_standard+
                       YearsFishing_standard+
                       Age_standard +(1|Village),
                     data=seafood[!is.na(seafood$Scylla_cons_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Scyllacons <- pp_check(testScyllacons,ndraws = 50)+ggtitle("Model fit Scylla Consumption")

loo(nullScyllacons, testScyllacons)

Scyllaconseffects<-as.data.frame(fixef(testScyllacons, probs = c(0.1,0.9)))
Scyllaconseffects$community<-row.names(Scyllaconseffects)
colnames(Scyllaconseffects)<-c("estimate","se","conf.low","conf.high","community")
Scyllaconseff <- ggplot(Scyllaconseffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Scylla Cons")

##Scylla ea
nullScyllaea<-brm(Scylla_ea_change ~ 1 +(1|Village),
                  data=seafood[!is.na(seafood$Scylla_ea_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testScyllaea <-brm(Scylla_ea_change ~ 
                     Gender+
                     VillageYears_standard+
                     YearsFishing_standard+
                     Age_standard +(1|Village),
                   data=seafood[!is.na(seafood$Scylla_ea_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Scyllaea <- pp_check(testScyllaea,ndraws = 50)+ggtitle("Model fit Scylla Abundance")

loo(nullScyllaea, testScyllaea)

Scyllaeaeffects<-as.data.frame(fixef(testScyllaea, probs = c(0.1,0.9)))
Scyllaeaeffects$community<-row.names(Scyllaeaeffects)
colnames(Scyllaeaeffects)<-c("estimate","se","conf.low","conf.high","community")
Scyllaeaeff <- ggplot(Scyllaeaeffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Scylla ea")

##Scylla harv
nullScyllaharv<-brm(Scylla_harv_change ~ 1 +(1|Village),
                    data=seafood[!is.na(seafood$Scylla_harv_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testScyllaharv <-brm(Scylla_harv_change ~ 
                       Gender+
                       VillageYears_standard+
                       YearsFishing_standard+
                       Age_standard +(1|Village),
                     data=seafood[!is.na(seafood$Scylla_harv_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Scyllaharv <- pp_check(testScyllaharv,ndraws = 50)+ggtitle("Model fit Scylla Harvest")

loo(nullScyllaharv, testScyllaharv)

Scyllaharveffects<-as.data.frame(fixef(testScyllaharv, probs = c(0.1,0.9)))
Scyllaharveffects$community<-row.names(Scyllaharveffects)
colnames(Scyllaharveffects)<-c("estimate","se","conf.low","conf.high","community")
Scyllaharveff <- ggplot(Scyllaharveffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Scylla harv")

##Export graphs
jpeg("PPcheck_Scylla", width =3500, height = 2000, res=300)
ggarrange(Scyllacons, Scyllaea, Scyllaharv)
dev.off()

jpeg("FixedEffects_Scylla", width =3500, height = 2000, res=300)
ggarrange(Scyllaconseff, Scyllaeaeff, Scyllaharveff)
dev.off()






##Palunirus
seafood$Palunirus_cons_change <- seafood$Palunirus_cons_now - seafood$Palunirus_cons_10
seafood$Palunirus_ea_change <- seafood$Palunirus_ea_now - seafood$Palunirus_ea_10
seafood$Palunirus_harv_change <- seafood$Palunirus_harv_now - seafood$Palunirus_harv_10

seafood$Palunirus_cons_change <- ordered(as.factor(seafood$Palunirus_cons_change), 
                                         levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

seafood$Palunirus_ea_change <- ordered(as.factor(seafood$Palunirus_ea_change), 
                                       levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

seafood$Palunirus_harv_change <- ordered(as.factor(seafood$Palunirus_harv_change), 
                                         levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

##Palunirus Cons
nullPaluniruscons<-brm(Palunirus_cons_change ~ 1 +(1|Village),
                       data=seafood[!is.na(seafood$Palunirus_cons_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testPaluniruscons <-brm(Palunirus_cons_change ~ 
                          Gender+
                          VillageYears_standard+
                          YearsFishing_standard+
                          Age_standard +(1|Village),
                        data=seafood[!is.na(seafood$Palunirus_cons_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Paluniruscons <- pp_check(testPaluniruscons,ndraws = 50)+ggtitle("Model fit Palunirus Consumption")

loo(nullPaluniruscons, testPaluniruscons)

Palunirusconseffects<-as.data.frame(fixef(testPaluniruscons, probs = c(0.1,0.9)))
Palunirusconseffects$community<-row.names(Palunirusconseffects)
colnames(Palunirusconseffects)<-c("estimate","se","conf.low","conf.high","community")
Palunirusconseff <- ggplot(Palunirusconseffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Palunirus Cons")

##Palunirus ea
nullPalunirusea<-brm(Palunirus_ea_change ~ 1 +(1|Village),
                     data=seafood[!is.na(seafood$Palunirus_ea_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testPalunirusea <-brm(Palunirus_ea_change ~ 
                        Gender+
                        VillageYears_standard+
                        YearsFishing_standard+
                        Age_standard +(1|Village),
                      data=seafood[!is.na(seafood$Palunirus_ea_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Palunirusea <- pp_check(testPalunirusea,ndraws = 50)+ggtitle("Model fit Palunirus Abundance")

loo(nullPalunirusea, testPalunirusea)

Paluniruseaeffects<-as.data.frame(fixef(testPalunirusea, probs = c(0.1,0.9)))
Paluniruseaeffects$community<-row.names(Paluniruseaeffects)
colnames(Paluniruseaeffects)<-c("estimate","se","conf.low","conf.high","community")
Paluniruseaeff <- ggplot(Paluniruseaeffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Palunirus ea")

##Palunirus harv
nullPalunirusharv<-brm(Palunirus_harv_change ~ 1 +(1|Village),
                       data=seafood[!is.na(seafood$Palunirus_harv_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testPalunirusharv <-brm(Palunirus_harv_change ~ 
                          Gender+
                          VillageYears_standard+
                          YearsFishing_standard+
                          Age_standard +(1|Village),
                        data=seafood[!is.na(seafood$Palunirus_harv_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Palunirusharv <- pp_check(testPalunirusharv,ndraws = 50)+ggtitle("Model fit Palunirus Harvest")

loo(nullPalunirusharv, testPalunirusharv)

Palunirusharveffects<-as.data.frame(fixef(testPalunirusharv, probs = c(0.1,0.9)))
Palunirusharveffects$community<-row.names(Palunirusharveffects)
colnames(Palunirusharveffects)<-c("estimate","se","conf.low","conf.high","community")
Palunirusharveff <- ggplot(Palunirusharveffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Palunirus harv")

##Export graphs
jpeg("PPcheck_Palunirus", width =3500, height = 2000, res=300)
ggarrange(Paluniruscons, Palunirusea, Palunirusharv)
dev.off()

jpeg("FixedEffects_Palunirus", width =3500, height = 2000, res=300)
ggarrange(Palunirusconseff, Paluniruseaeff, Palunirusharveff)
dev.off()







##Carangidae
seafood$Carangidae_cons_change <- seafood$Carangidae_cons_now - seafood$Carangidae_cons_10
seafood$Carangidae_ea_change <- seafood$Carangidae_ea_now - seafood$Carangidae_ea_10
seafood$Carangidae_harv_change <- seafood$Carangidae_harv_now - seafood$Carangidae_harv_10

seafood$Carangidae_cons_change <- ordered(as.factor(seafood$Carangidae_cons_change), 
                                          levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

seafood$Carangidae_ea_change <- ordered(as.factor(seafood$Carangidae_ea_change), 
                                        levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

seafood$Carangidae_harv_change <- ordered(as.factor(seafood$Carangidae_harv_change), 
                                          levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

##Carangidae Cons
nullCarangidaecons<-brm(Carangidae_cons_change ~ 1 +(1|Village),
                        data=seafood[!is.na(seafood$Carangidae_cons_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testCarangidaecons <-brm(Carangidae_cons_change ~ 
                           Gender+
                           VillageYears_standard+
                           YearsFishing_standard+
                           Age_standard +(1|Village),
                         data=seafood[!is.na(seafood$Carangidae_cons_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Carangidaecons <- pp_check(testCarangidaecons,ndraws = 50)+ggtitle("Model fit Carangidae Consumption")

loo(nullCarangidaecons, testCarangidaecons)

Carangidaeconseffects<-as.data.frame(fixef(testCarangidaecons, probs = c(0.1,0.9)))
Carangidaeconseffects$community<-row.names(Carangidaeconseffects)
colnames(Carangidaeconseffects)<-c("estimate","se","conf.low","conf.high","community")
Carangidaeconseff <- ggplot(Carangidaeconseffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Carangidae Cons")

##Carangidae ea
nullCarangidaeea<-brm(Carangidae_ea_change ~ 1 +(1|Village),
                      data=seafood[!is.na(seafood$Carangidae_ea_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testCarangidaeea <-brm(Carangidae_ea_change ~ 
                         Gender+
                         VillageYears_standard+
                         YearsFishing_standard+
                         Age_standard +(1|Village),
                       data=seafood[!is.na(seafood$Carangidae_ea_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Carangidaeea <- pp_check(testCarangidaeea,ndraws = 50)+ggtitle("Model fit Carangidae Abundance")

loo(nullCarangidaeea, testCarangidaeea)

Carangidaeeaeffects<-as.data.frame(fixef(testCarangidaeea, probs = c(0.1,0.9)))
Carangidaeeaeffects$community<-row.names(Carangidaeeaeffects)
colnames(Carangidaeeaeffects)<-c("estimate","se","conf.low","conf.high","community")
Carangidaeeaeff <- ggplot(Carangidaeeaeffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Carangidae ea")

##Carangidae harv
nullCarangidaeharv<-brm(Carangidae_harv_change ~ 1 +(1|Village),
                        data=seafood[!is.na(seafood$Carangidae_harv_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testCarangidaeharv <-brm(Carangidae_harv_change ~ 
                           Gender+
                           VillageYears_standard+
                           YearsFishing_standard+
                           Age_standard +(1|Village),
                         data=seafood[!is.na(seafood$Carangidae_harv_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Carangidaeharv <- pp_check(testCarangidaeharv,ndraws = 50)+ggtitle("Model fit Carangidae Harvest")

loo(nullCarangidaeharv, testCarangidaeharv)

Carangidaeharveffects<-as.data.frame(fixef(testCarangidaeharv, probs = c(0.1,0.9)))
Carangidaeharveffects$community<-row.names(Carangidaeharveffects)
colnames(Carangidaeharveffects)<-c("estimate","se","conf.low","conf.high","community")
Carangidaeharveff <- ggplot(Carangidaeharveffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Carangidae harv")

##Export graphs
jpeg("PPcheck_Carangidae", width =3500, height = 2000, res=300)
ggarrange(Carangidaecons, Carangidaeea, Carangidaeharv)
dev.off()

jpeg("FixedEffects_Carangidae", width =3500, height = 2000, res=300)
ggarrange(Carangidaeconseff, Carangidaeeaeff, Carangidaeharveff)
dev.off()







##Balistidae
seafood$Balistidae_cons_change <- seafood$Balistidae_cons_now - seafood$Balistidae_cons_10
seafood$Balistidae_ea_change <- seafood$Balistidae_ea_now - seafood$Balistidae_ea_10
seafood$Balistidae_harv_change <- seafood$Balistidae_harv_now - seafood$Balistidae_harv_10

seafood$Balistidae_cons_change <- ordered(as.factor(seafood$Balistidae_cons_change), 
                                          levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

seafood$Balistidae_ea_change <- ordered(as.factor(seafood$Balistidae_ea_change), 
                                        levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

seafood$Balistidae_harv_change <- ordered(as.factor(seafood$Balistidae_harv_change), 
                                          levels = c("-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4","5", "NA"))

##Balistidae Cons
nullBalistidaecons<-brm(Balistidae_cons_change ~ 1 +(1|Village),
                        data=seafood[!is.na(seafood$Balistidae_cons_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testBalistidaecons <-brm(Balistidae_cons_change ~ 
                           Gender+
                           VillageYears_standard+
                           YearsFishing_standard+
                           Age_standard +(1|Village),
                         data=seafood[!is.na(seafood$Balistidae_cons_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Balistidaecons <- pp_check(testBalistidaecons,ndraws = 50)+ggtitle("Model fit Balistidae Consumption")

loo(nullBalistidaecons, testBalistidaecons)

Balistidaeconseffects<-as.data.frame(fixef(testBalistidaecons, probs = c(0.1,0.9)))
Balistidaeconseffects$community<-row.names(Balistidaeconseffects)
colnames(Balistidaeconseffects)<-c("estimate","se","conf.low","conf.high","community")
Balistidaeconseff <- ggplot(Balistidaeconseffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Balistidae Cons")

##Balistidae ea
nullBalistidaeea<-brm(Balistidae_ea_change ~ 1 +(1|Village),
                      data=seafood[!is.na(seafood$Balistidae_ea_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testBalistidaeea <-brm(Balistidae_ea_change ~ 
                         Gender+
                         VillageYears_standard+
                         YearsFishing_standard+
                         Age_standard +(1|Village),
                       data=seafood[!is.na(seafood$Balistidae_ea_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Balistidaeea <- pp_check(testBalistidaeea,ndraws = 50)+ggtitle("Model fit Balistidae Abundance")

loo(nullBalistidaeea, testBalistidaeea)

Balistidaeeaeffects<-as.data.frame(fixef(testBalistidaeea, probs = c(0.1,0.9)))
Balistidaeeaeffects$community<-row.names(Balistidaeeaeffects)
colnames(Balistidaeeaeffects)<-c("estimate","se","conf.low","conf.high","community")
Balistidaeeaeff <- ggplot(Balistidaeeaeffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Balistidae ea")

##Balistidae harv
nullBalistidaeharv<-brm(Balistidae_harv_change ~ 1 +(1|Village),
                        data=seafood[!is.na(seafood$Balistidae_harv_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

testBalistidaeharv <-brm(Balistidae_harv_change ~ 
                           Gender+
                           VillageYears_standard+
                           YearsFishing_standard+
                           Age_standard +(1|Village),
                         data=seafood[!is.na(seafood$Balistidae_harv_change),],family=cumulative("logit"), iter=3000, control = list(adapt_delta=0.99))

Balistidaeharv <- pp_check(testBalistidaeharv,ndraws = 50)+ggtitle("Model fit Balistidae Harvest")

loo(nullBalistidaeharv, testBalistidaeharv)

Balistidaeharveffects<-as.data.frame(fixef(testBalistidaeharv, probs = c(0.1,0.9)))
Balistidaeharveffects$community<-row.names(Balistidaeharveffects)
colnames(Balistidaeharveffects)<-c("estimate","se","conf.low","conf.high","community")
Balistidaeharveff <- ggplot(Balistidaeharveffects,aes(x=estimate,y=community,xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  ggtitle("Offset from intercept_Balistidae harv")

##Export graphs
jpeg("PPcheck_Balistidae", width =3500, height = 2000, res=300)
ggarrange(Balistidaecons, Balistidaeea, Balistidaeharv)
dev.off()

jpeg("FixedEffects_Balistidae", width =3500, height = 2000, res=300)
ggarrange(Balistidaeconseff, Balistidaeeaeff, Balistidaeharveff)
dev.off()





##Generating Effect Size Graphs for Consumption
mullidaeconseffects$Family <- "Mullidae"
Pomacentridaeconseffects$Family <- "Pomacentridae"
Acanthuridaeconseffects$Family <- "Acanthuridae"
Scaridaeconseffects$Family <- "Scaridae"
Scorpaenidaeconseffects$Family <- "Scorpaenidae"
Siganidaeconseffects$Family <- "Siganidae"
Sphyraenidaeconseffects$Family <- "Sphyraenidae"
Labridaeconseffects$Family <- "Labridae"
Priacanthidaeconseffects$Family <- "Priacanthidae"
Chaetodontidaeconseffects$Family <- "Chaetodontidae"
Fistularidaeconseffects$Family <- "Fistularidae"
Gerreidaeconseffects$Family <- "Gerreidae"
Haemulidaeconseffects$Family <- "Haemulidae"
Lethrinidaeconseffects$Family <- "Lethrinidae"
Ostraciidaeconseffects$Family <- "Ostraciidae"
Pomacanthidaeconseffects$Family <- "Pomacanthidae"
Rajidaeconseffects$Family <- "Rajidae"
Rhinobatidaeconseffects$Family <- "Rhinobatidae"
Soleidaeconseffects$Family <- "Soleidae"
Tetraodontidaeconseffects$Family <- "Tetraodontidae"
Carangidaeconseffects$Family <- "Carangidae"
Balistidaeconseffects$Family <- "Balistidae"
Carangidaeconseffects$Family <- "Carangidae"
Balistidaeconseffects$Family<- "Balistidae"
Zanclidaeconseffects$Family <- "Zanclidae"
Octopusconseffects$Family <- "Octopus"
Loligoconseffects$Family <- "Loligo"
Sepiaconseffects$Family <- "Sepia"
Cypraea.Conusconseffects$Family <- "Cypraea.Conus"
Pyrasusconseffects$Family <- "Pyrasus"
Charoniaconseffects$Family <- "Charonia"
Murex.Fasciolariaconseffects$Family <- "Murex.Fasciolaria"
Lambisconseffects$Family <- "Lambis"
Anadaraconseffects$Family <- "Anadara"
Pinctada.Isognomon.Atrina.Pinnaconseffects$Family <- "Pinctada.Isognomon.Atrina.Pinna"
Tridacnaconseffects$Family <- "Tridacna"
Tripneustesconseffects$Family <- "Tripneustes"
Aristeidaeconseffects$Family <- "Aristeidae"
Scyllaconseffects$Family <- "Scylla"
Palunirusconseffects$Family <- "Palunirus"



EffectsizeConsumption <- rbind(mullidaeconseffects,
                               Pomacentridaeconseffects,
                               Acanthuridaeconseffects,
                               Scaridaeconseffects,
                               Scorpaenidaeconseffects,
                               Siganidaeconseffects,
                               Sphyraenidaeconseffects,
                               Labridaeconseffects,
                               Priacanthidaeconseffects,
                               Chaetodontidaeconseffects,
                               Fistularidaeconseffects,
                               Gerreidaeconseffects,
                               Haemulidaeconseffects,
                               Lethrinidaeconseffects,
                               Ostraciidaeconseffects,
                               Pomacanthidaeconseffects,
                               Rajidaeconseffects,
                               Rhinobatidaeconseffects,
                               Soleidaeconseffects,
                               Tetraodontidaeconseffects,
                               Carangidaeconseffects,
                               Balistidaeconseffects,
                               Zanclidaeconseffects,
                               Octopusconseffects,
                               Loligoconseffects,
                               Sepiaconseffects,
                               Cypraea.Conusconseffects,
                               Pyrasusconseffects,
                               Charoniaconseffects,
                               Murex.Fasciolariaconseffects,
                               Lambisconseffects,
                               Anadaraconseffects,
                               Pinctada.Isognomon.Atrina.Pinnaconseffects,
                               Tridacnaconseffects,
                               Tripneustesconseffects,
                               Aristeidaeconseffects,
                               Scyllaconseffects,
                               Palunirusconseffects)

EffectsizeConsumption[EffectsizeConsumption == "Cypraea.Conus"] <- "Cypraeidae/Conus"
EffectsizeConsumption[EffectsizeConsumption == "Pinctada.Isognomon.Atrina.Pinna"] <- "Pteriidae/Isognomonidae"
EffectsizeConsumption[EffectsizeConsumption == "Murex.Fasciolaria"] <- "Muricidae"
EffectsizeConsumption[EffectsizeConsumption == "Octopus"] <- "Octopoda"


write.csv(EffectsizeConsumption, "ConsumptionEffectSizes.csv", row.names=FALSE)

EffectsizeConsumption <- read.csv("ConsumptionEffectSizes.csv")


EffectsizeConsumptionGender <- EffectsizeConsumption[EffectsizeConsumption$community == 'GenderFemale',]

jpeg("ConsumptionGenderEffectSize", width =2000, height = 3500, res=300)
ggplot(EffectsizeConsumptionGender,aes(x=estimate,y=fct_rev(fct_inorder(Family)),xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  geom_smooth() +
  labs(title = str_wrap("Effect of Gender on Perception of Change in Consumption of each Marine Resource Group", 45))+
  labs(x="Effect Size of Female", y="Marine Resource Group")+
  theme(plot.title = element_text(hjust = 0.5, margin=margin(b=20)))+
  theme(axis.title.x = element_text(margin=margin(t=20)), axis.title.y = element_text(margin=margin(r=20)))
dev.off()

EffectsizeConsumptionAge <- EffectsizeConsumption[EffectsizeConsumption$community == 'Age_standard',]

jpeg("ConsumptionAgeEffectSize", width =2000, height = 3500, res=300)
ggplot(EffectsizeConsumptionAge,aes(x=estimate,y=fct_rev(fct_inorder(Family)),xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  geom_smooth() +
  labs(title = str_wrap("Effect of Age on Perception of Change in Consumption of each Marine Resource Group", 45))+
  labs(x="Effect Size", y="Marine Resource Group")+
  theme(plot.title = element_text(hjust = 0.5, margin=margin(b=20)))+
  theme(axis.title.x = element_text(margin=margin(t=20)), axis.title.y = element_text(margin=margin(r=20)))
dev.off()

EffectsizeConsumptionLived <- EffectsizeConsumption[EffectsizeConsumption$community == 'VillageYears_standard',]

jpeg("ConsumptionYearsVillageEffectSize", width =2000, height = 3500, res=300)
ggplot(EffectsizeConsumptionLived,aes(x=estimate,y=fct_rev(fct_inorder(Family)),xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  geom_smooth() +
  labs(title = str_wrap("Effect of Years Spent in Current Village on Perception of Change in Consumption of each Marine Resource Group", 45))+
  labs(x="Effect Size", y="Marine Resource Group")+
  theme(plot.title = element_text(hjust = 0.5, margin=margin(b=20)))+
  theme(axis.title.x = element_text(margin=margin(t=20)), axis.title.y = element_text(margin=margin(r=20)))
dev.off()

EffectsizeConsumptionFishing <- EffectsizeConsumption[EffectsizeConsumption$community == 'YearsFishing_standard',]

jpeg("ConsumptionYearsFishingEffectSize", width =2000, height = 3500, res=300)
ggplot(EffectsizeConsumptionFishing,aes(x=estimate,y=fct_rev(fct_inorder(Family)),xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="darkred",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  geom_smooth() +
  labs(title = str_wrap("Effect of Years Spent Fishing/Gleaning on Perception of Change in Consumption of each Marine Resource Group", 45))+
  labs(x="Effect Size", y="Marine Resource Group")+
  theme(plot.title = element_text(hjust = 0.5, margin=margin(b=20)))+
  theme(axis.title.x = element_text(margin=margin(t=20)), axis.title.y = element_text(margin=margin(r=20)))
dev.off()










##Generating Effect Size Graphs for Ecological Abundance
mullidaeeaeffects$Family <- "Mullidae"
Pomacentridaeeaeffects$Family <- "Pomacentridae"
Acanthuridaeeaeffects$Family <- "Acanthuridae"
Scaridaeeaeffects$Family <- "Scaridae"
Scorpaenidaeeaeffects$Family <- "Scorpaenidae"
Siganidaeeaeffects$Family <- "Siganidae"
Sphyraenidaeeaeffects$Family <- "Sphyraenidae"
Labridaeeaeffects$Family <- "Labridae"
Priacanthidaeeaeffects$Family <- "Priacanthidae"
Chaetodontidaeeaeffects$Family <- "Chaetodontidae"
Fistularidaeeaeffects$Family <- "Fistularidae"
Gerreidaeeaeffects$Family <- "Gerreidae"
Haemulidaeeaeffects$Family <- "Haemulidae"
Lethrinidaeeaeffects$Family <- "Lethrinidae"
Ostraciidaeeaeffects$Family <- "Ostraciidae"
Pomacanthidaeeaeffects$Family <- "Pomacanthidae"
Rajidaeeaeffects$Family <- "Rajidae"
Rhinobatidaeeaeffects$Family <- "Rhinobatidae"
Soleidaeeaeffects$Family <- "Soleidae"
Tetraodontidaeeaeffects$Family <- "Tetraodontidae"
Carangidaeeaeffects$Family <- "Carangidae"
Balistidaeeaeffects$Family <- "Balistidae"
Carangidaeeaeffects$Family <- "Carangidae"
Balistidaeeaeffects$Family<- "Balistidae"
Zanclidaeeaeffects$Family <- "Zanclidae"
Octopuseaeffects$Family <- "Octopus"
Loligoeaeffects$Family <- "Loligo"
Sepiaeaeffects$Family <- "Sepia"
Cypraea.Conuseaeffects$Family <- "Cypraea.Conus"
Pyrasuseaeffects$Family <- "Pyrasus"
Charoniaeaeffects$Family <- "Charonia"
Murex.Fasciolariaeaeffects$Family <- "Murex.Fasciolaria"
Lambiseaeffects$Family <- "Lambis"
Anadaraeaeffects$Family <- "Anadara"
Pinctada.Isognomon.Atrina.Pinnaeaeffects$Family <- "Pinctada.Isognomon.Atrina.Pinna"
Tridacnaeaeffects$Family <- "Tridacna"
Tripneusteseaeffects$Family <- "Tripneustes"
Aristeidaeeaeffects$Family <- "Aristeidae"
Scyllaeaeffects$Family <- "Scylla"
Paluniruseaeffects$Family <- "Palunirus"



EffectsizeEcologicalAbundance <- rbind(mullidaeeaeffects,
                                       Pomacentridaeeaeffects,
                                       Acanthuridaeeaeffects,
                                       Scaridaeeaeffects,
                                       Scorpaenidaeeaeffects,
                                       Siganidaeeaeffects,
                                       Sphyraenidaeeaeffects,
                                       Labridaeeaeffects,
                                       Priacanthidaeeaeffects,
                                       Chaetodontidaeeaeffects,
                                       Fistularidaeeaeffects,
                                       Gerreidaeeaeffects,
                                       Haemulidaeeaeffects,
                                       Lethrinidaeeaeffects,
                                       Ostraciidaeeaeffects,
                                       Pomacanthidaeeaeffects,
                                       Rajidaeeaeffects,
                                       Rhinobatidaeeaeffects,
                                       Soleidaeeaeffects,
                                       Tetraodontidaeeaeffects,
                                       Carangidaeeaeffects,
                                       Balistidaeeaeffects,
                                       Zanclidaeeaeffects,
                                       Octopuseaeffects,
                                       Loligoeaeffects,
                                       Sepiaeaeffects,
                                       Cypraea.Conuseaeffects,
                                       Pyrasuseaeffects,
                                       Charoniaeaeffects,
                                       Murex.Fasciolariaeaeffects,
                                       Lambiseaeffects,
                                       Anadaraeaeffects,
                                       Pinctada.Isognomon.Atrina.Pinnaeaeffects,
                                       Tridacnaeaeffects,
                                       Tripneusteseaeffects,
                                       Aristeidaeeaeffects,
                                       Scyllaeaeffects,
                                       Paluniruseaeffects)

EffectsizeEcologicalAbundance[EffectsizeEcologicalAbundance == "Cypraea.Conus"] <- "Cypraeidae/Conus"
EffectsizeEcologicalAbundance[EffectsizeEcologicalAbundance == "Pinctada.Isognomon.Atrina.Pinna"] <- "Pteriidae/Isognomonidae"
EffectsizeEcologicalAbundance[EffectsizeEcologicalAbundance == "Murex.Fasciolaria"] <- "Muricidae"
EffectsizeEcologicalAbundance[EffectsizeEcologicalAbundance == "Octopus"] <- "Octopoda"

write.csv(EffectsizeEcologicalAbundance, "EcologicalAbundanceEffectSizes.csv", row.names=FALSE)

EffectsizeEcologicalAbundance <- read.csv("EcologicalAbundanceEffectSizes.csv")


EffectsizeEcologicalAbundanceGender <- EffectsizeEcologicalAbundance[EffectsizeEcologicalAbundance$community == 'GenderFemale',]

jpeg("EcologicalAbundanceGenderEffectSize", width =2000, height = 3500, res=300)
ggplot(EffectsizeEcologicalAbundanceGender,aes(x=estimate,y=fct_rev(fct_inorder(Family)),xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="forestgreen",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  geom_smooth() +
  labs(title = str_wrap("Effect of Gender on Perception of Change in Ecological Abundance of each Marine Resource Group", 45))+
  labs(x="Effect Size of Female", y="Marine Resource Group")+
  theme(plot.title = element_text(hjust = 0.5, margin=margin(b=20)))+
  theme(axis.title.x = element_text(margin=margin(t=20)), axis.title.y = element_text(margin=margin(r=20)))
dev.off()

EffectsizeEcologicalAbundanceAge <- EffectsizeEcologicalAbundance[EffectsizeEcologicalAbundance$community == 'Age_standard',]

jpeg("EcologicalAbundanceAgeEffectSize", width =2000, height = 3500, res=300)
ggplot(EffectsizeEcologicalAbundanceAge,aes(x=estimate,y=fct_rev(fct_inorder(Family)),xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="forestgreen",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  geom_smooth() +
  labs(title = str_wrap("Effect of Age on Perception of Change in Ecological Abundance of each Marine Resource Group", 45))+
  labs(x="Effect Size", y="Marine Resource Group")+
  theme(plot.title = element_text(hjust = 0.5, margin=margin(b=20)))+
  theme(axis.title.x = element_text(margin=margin(t=20)), axis.title.y = element_text(margin=margin(r=20)))
dev.off()

EffectsizeEcologicalAbundanceLived <- EffectsizeEcologicalAbundance[EffectsizeEcologicalAbundance$community == 'VillageYears_standard',]

jpeg("EcologicalAbundanceYearsVillageEffectSize", width =2000, height = 3500, res=300)
ggplot(EffectsizeEcologicalAbundanceLived,aes(x=estimate,y=fct_rev(fct_inorder(Family)),xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="forestgreen",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  geom_smooth() +
  labs(title = str_wrap("Effect of Years Spent in Current Village on Perception of Change in Ecological Abundance of each Marine Resource Group", 45))+
  labs(x="Effect Size", y="Marine Resource Group")+
  theme(plot.title = element_text(hjust = 0.5, margin=margin(b=20)))+
  theme(axis.title.x = element_text(margin=margin(t=20)), axis.title.y = element_text(margin=margin(r=20)))
dev.off()

EffectsizeEcologicalAbundanceFishing <- EffectsizeEcologicalAbundance[EffectsizeEcologicalAbundance$community == 'YearsFishing_standard',]

jpeg("EcologicalAbundanceYearsFishingEffectSize", width =2000, height = 3500, res=300)
ggplot(EffectsizeEcologicalAbundanceFishing,aes(x=estimate,y=fct_rev(fct_inorder(Family)),xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="forestgreen",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  geom_smooth() +
  labs(title = str_wrap("Effect of Years Spent Fishing/Gleaning on Perception of Change in Ecological Abundance of each Marine Resource Group", 45))+
  labs(x="Effect Size", y="Marine Resource Group")+
  theme(plot.title = element_text(hjust = 0.5, margin=margin(b=20)))+
  theme(axis.title.x = element_text(margin=margin(t=20)), axis.title.y = element_text(margin=margin(r=20)))
dev.off()






##Generating Effect Size Graphs for Harvest
mullidaeharveffects$Family <- "Mullidae"
Pomacentridaeharveffects$Family <- "Pomacentridae"
Acanthuridaeharveffects$Family <- "Acanthuridae"
Scaridaeharveffects$Family <- "Scaridae"
Scorpaenidaeharveffects$Family <- "Scorpaenidae"
Siganidaeharveffects$Family <- "Siganidae"
Sphyraenidaeharveffects$Family <- "Sphyraenidae"
Labridaeharveffects$Family <- "Labridae"
Priacanthidaeharveffects$Family <- "Priacanthidae"
Chaetodontidaeharveffects$Family <- "Chaetodontidae"
Fistularidaeharveffects$Family <- "Fistularidae"
Gerreidaeharveffects$Family <- "Gerreidae"
Haemulidaeharveffects$Family <- "Haemulidae"
Lethrinidaeharveffects$Family <- "Lethrinidae"
Ostraciidaeharveffects$Family <- "Ostraciidae"
Pomacanthidaeharveffects$Family <- "Pomacanthidae"
Rajidaeharveffects$Family <- "Rajidae"
Rhinobatidaeharveffects$Family <- "Rhinobatidae"
Soleidaeharveffects$Family <- "Soleidae"
Tetraodontidaeharveffects$Family <- "Tetraodontidae"
Carangidaeharveffects$Family <- "Carangidae"
Balistidaeharveffects$Family <- "Balistidae"
Carangidaeharveffects$Family <- "Carangidae"
Balistidaeharveffects$Family<- "Balistidae"
Zanclidaeharveffects$Family <- "Zanclidae"
Octopusharveffects$Family <- "Octopus"
Loligoharveffects$Family <- "Loligo"
Sepiaharveffects$Family <- "Sepia"
Cypraea.Conusharveffects$Family <- "Cypraea.Conus"
Pyrasusharveffects$Family <- "Pyrasus"
Charoniaharveffects$Family <- "Charonia"
Murex.Fasciolariaharveffects$Family <- "Murex.Fasciolaria"
Lambisharveffects$Family <- "Lambis"
Anadaraharveffects$Family <- "Anadara"
Pinctada.Isognomon.Atrina.Pinnaharveffects$Family <- "Pinctada.Isognomon.Atrina.Pinna"
Tridacnaharveffects$Family <- "Tridacna"
Tripneustesharveffects$Family <- "Tripneustes"
Aristeidaeharveffects$Family <- "Aristeidae"
Scyllaharveffects$Family <- "Scylla"
Palunirusharveffects$Family <- "Palunirus"



EffectsizeHarvest <- rbind(mullidaeharveffects,
                           Pomacentridaeharveffects,
                           Acanthuridaeharveffects,
                           Scaridaeharveffects,
                           Scorpaenidaeharveffects,
                           Siganidaeharveffects,
                           Sphyraenidaeharveffects,
                           Labridaeharveffects,
                           Priacanthidaeharveffects,
                           Chaetodontidaeharveffects,
                           Fistularidaeharveffects,
                           Gerreidaeharveffects,
                           Haemulidaeharveffects,
                           Lethrinidaeharveffects,
                           Ostraciidaeharveffects,
                           Pomacanthidaeharveffects,
                           Rajidaeharveffects,
                           Rhinobatidaeharveffects,
                           Soleidaeharveffects,
                           Tetraodontidaeharveffects,
                           Carangidaeharveffects,
                           Balistidaeharveffects,
                           Zanclidaeharveffects,
                           Octopusharveffects,
                           Loligoharveffects,
                           Sepiaharveffects,
                           Cypraea.Conusharveffects,
                           Pyrasusharveffects,
                           Charoniaharveffects,
                           Murex.Fasciolariaharveffects,
                           Lambisharveffects,
                           Anadaraharveffects,
                           Pinctada.Isognomon.Atrina.Pinnaharveffects,
                           Tridacnaharveffects,
                           Tripneustesharveffects,
                           Aristeidaeharveffects,
                           Scyllaharveffects,
                           Palunirusharveffects)

EffectsizeHarvest[EffectsizeHarvest == "Cypraea.Conus"] <- "Cypraeidae/Conus"
EffectsizeHarvest[EffectsizeHarvest == "Pinctada.Isognomon.Atrina.Pinna"] <- "Pteriidae/Isognomonidae"
EffectsizeHarvest[EffectsizeHarvest == "Murex.Fasciolaria"] <- "Muricidae"
EffectsizeHarvest[EffectsizeHarvest == "Octopus"] <- "Octopoda"

write.csv(EffectsizeHarvest, "HarvestEffectSizes.csv", row.names=FALSE)

EffectsizeHarvest <- read.csv("HarvestEffectSizes.csv")


EffectsizeHarvestGender <- EffectsizeHarvest[EffectsizeHarvest$community == 'GenderFemale',]

jpeg("HarvestGenderEffectSize", width =2000, height = 3500, res=300)
ggplot(EffectsizeHarvestGender,aes(x=estimate,y=fct_rev(fct_inorder(Family)),xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="blue3",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  geom_smooth() +
  labs(title = str_wrap("Effect of Gender on Perception of Change in Harvest of each Marine Resource Group", 45))+
  labs(x="Effect Size of Female", y="Marine Resource Group")+
  theme(plot.title = element_text(hjust = 0.5, margin=margin(b=20)))+
  theme(axis.title.x = element_text(margin=margin(t=20)), axis.title.y = element_text(margin=margin(r=20)))
dev.off()

EffectsizeHarvestAge <- EffectsizeHarvest[EffectsizeHarvest$community == 'Age_standard',]

jpeg("HarvestAgeEffectSize", width =2000, height = 3500, res=300)
ggplot(EffectsizeHarvestAge,aes(x=estimate,y=fct_rev(fct_inorder(Family)),xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="blue3",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  geom_smooth() +
  labs(title = str_wrap("Effect of Age on Perception of Change in Harvest of each Marine Resource Group", 45))+
  labs(x="Effect Size", y="Marine Resource Group")+
  theme(plot.title = element_text(hjust = 0.5, margin=margin(b=20)))+
  theme(axis.title.x = element_text(margin=margin(t=20)), axis.title.y = element_text(margin=margin(r=20)))
dev.off()

EffectsizeHarvestLived <- EffectsizeHarvest[EffectsizeHarvest$community == 'VillageYears_standard',]

jpeg("HarvestYearsVillageEffectSize", width =2000, height = 3500, res=300)
ggplot(EffectsizeHarvestLived,aes(x=estimate,y=fct_rev(fct_inorder(Family)),xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="blue3",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  geom_smooth() +
  labs(title = str_wrap("Effect of Years Spent in Current Village on Perception of Change in Harvest of each Marine Resource Group", 45))+
  labs(x="Effect Size", y="Marine Resource Group")+
  theme(plot.title = element_text(hjust = 0.5, margin=margin(b=20)))+
  theme(axis.title.x = element_text(margin=margin(t=20)), axis.title.y = element_text(margin=margin(r=20)))
dev.off()

EffectsizeHarvestFishing <- EffectsizeHarvest[EffectsizeHarvest$community == 'YearsFishing_standard',]

jpeg("HarvestYearsFishingEffectSize", width =2000, height = 3500, res=300)
ggplot(EffectsizeHarvestFishing,aes(x=estimate,y=fct_rev(fct_inorder(Family)),xmin=conf.low,xmax=conf.high))+geom_errorbar()+geom_point(fill="blue3",pch=21,size=3)+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  geom_smooth() +
  labs(title = str_wrap("Effect of Years Spent Fishing/Gleaning on Perception of Change in Harvest of each Marine Resource Group", 45))+
  labs(x="Effect Size", y="Marine Resource Group")+
  theme(plot.title = element_text(hjust = 0.5, margin=margin(b=20)))+
  theme(axis.title.x = element_text(margin=margin(t=20)), axis.title.y = element_text(margin=margin(r=20)))
dev.off()




##Generate Summary Averages Graphs



averages <- select(seafood, c("Mullidae_cons_change", 
                              "Pomacentridae_cons_change",
                              "Acanthuridae_cons_change", 
                              "Scaridae_cons_change", 
                              "Scorpaenidae_cons_change", 
                              "Siganidae_cons_change", 
                              "Sphyraenidae_cons_change", 
                              "Labridae_cons_change", 
                              "Priacanthidae_cons_change", 
                              "Chaetodontidae_cons_change", 
                              "Fistularidae_cons_change", 
                              "Gerreidae_cons_change", 
                              "Haemulidae_cons_change", 
                              "Lethrinidae_cons_change", 
                              "Ostraciidae_cons_change", 
                              "Pomacanthidae_cons_change", 
                              "Rajidae_cons_change", 
                              "Rhinobatidae_cons_change", 
                              "Soleidae_cons_change", 
                              "Tetraodontidae_cons_change", 
                              "Carangidae_cons_change", 
                              "Balistidae_cons_change", 
                              "Carangidae_cons_change",
                              "Zanclidae_cons_change",
                              "Octopus_cons_change", 
                              "Loligo_cons_change", 
                              "Sepia_cons_change", 
                              "Cypraea.Conus_cons_change", 
                              "Pyrasus_cons_change", 
                              "Charonia_cons_change", 
                              "Murex.Fasciolaria_cons_change", 
                              "Lambis_cons_change", 
                              "Anadara_cons_change", 
                              "Pinctada.Isognomon.Atrina.Pinna_cons_change", 
                              "Tridacna_cons_change", 
                              "Tripneustes_cons_change", 
                              "Aristeidae_cons_change", 
                              "Scylla_cons_change", 
                              "Palunirus_cons_change",
                              "Mullidae_ea_change", 
                              "Pomacentridae_ea_change", 
                              "Acanthuridae_ea_change", 
                              "Scaridae_ea_change", 
                              "Scorpaenidae_ea_change", 
                              "Siganidae_ea_change", 
                              "Sphyraenidae_ea_change", 
                              "Labridae_ea_change", 
                              "Priacanthidae_ea_change", 
                              "Chaetodontidae_ea_change", 
                              "Fistularidae_ea_change", 
                              "Gerreidae_ea_change", 
                              "Haemulidae_ea_change", 
                              "Lethrinidae_ea_change", 
                              "Ostraciidae_ea_change", 
                              "Pomacanthidae_ea_change", 
                              "Rajidae_ea_change", 
                              "Rhinobatidae_ea_change", 
                              "Soleidae_ea_change", 
                              "Tetraodontidae_ea_change", 
                              "Carangidae_ea_change", 
                              "Balistidae_ea_change",
                              "Zanclidae_ea_change", 
                              "Octopus_ea_change", 
                              "Loligo_ea_change", 
                              "Sepia_ea_change", 
                              "Cypraea.Conus_ea_change", 
                              "Pyrasus_ea_change", 
                              "Charonia_ea_change", 
                              "Murex.Fasciolaria_ea_change", 
                              "Lambis_ea_change", 
                              "Anadara_ea_change", 
                              "Pinctada.Isognomon.Atrina.Pinna_ea_change", 
                              "Tridacna_ea_change", 
                              "Tripneustes_ea_change", 
                              "Aristeidae_ea_change", 
                              "Scylla_ea_change", 
                              "Palunirus_ea_change",
                              "Mullidae_harv_change", 
                              "Pomacentridae_harv_change", 
                              "Acanthuridae_harv_change", 
                              "Scaridae_harv_change", 
                              "Scorpaenidae_harv_change", 
                              "Siganidae_harv_change", 
                              "Sphyraenidae_harv_change", 
                              "Labridae_harv_change", 
                              "Priacanthidae_harv_change", 
                              "Chaetodontidae_harv_change", 
                              "Fistularidae_harv_change", 
                              "Gerreidae_harv_change", 
                              "Haemulidae_harv_change", 
                              "Lethrinidae_harv_change", 
                              "Ostraciidae_harv_change", 
                              "Pomacanthidae_harv_change", 
                              "Rajidae_harv_change", 
                              "Rhinobatidae_harv_change", 
                              "Soleidae_harv_change", 
                              "Tetraodontidae_harv_change", 
                              "Carangidae_harv_change", 
                              "Balistidae_harv_change",
                              "Zanclidae_harv_change", 
                              "Octopus_harv_change", 
                              "Loligo_harv_change", 
                              "Sepia_harv_change", 
                              "Cypraea.Conus_harv_change", 
                              "Pyrasus_harv_change", 
                              "Charonia_harv_change", 
                              "Murex.Fasciolaria_harv_change", 
                              "Lambis_harv_change", 
                              "Anadara_harv_change", 
                              "Pinctada.Isognomon.Atrina.Pinna_harv_change", 
                              "Tridacna_harv_change", 
                              "Tripneustes_harv_change", 
                              "Aristeidae_harv_change", 
                              "Scylla_harv_change", 
                              "Palunirus_harv_change"))


averages = as.data.frame(sapply(averages, as.character))
averages = as.data.frame(sapply(averages, as.numeric))


tableofmeans <- as.data.frame(colMeans(averages, na.rm=TRUE))

consumptionaverages <- select(averages, c("Mullidae_cons_change", 
                                          "Pomacentridae_cons_change",
                                          "Acanthuridae_cons_change", 
                                          "Scaridae_cons_change", 
                                          "Scorpaenidae_cons_change", 
                                          "Siganidae_cons_change", 
                                          "Sphyraenidae_cons_change", 
                                          "Labridae_cons_change", 
                                          "Priacanthidae_cons_change", 
                                          "Chaetodontidae_cons_change", 
                                          "Fistularidae_cons_change", 
                                          "Gerreidae_cons_change", 
                                          "Haemulidae_cons_change", 
                                          "Lethrinidae_cons_change", 
                                          "Ostraciidae_cons_change", 
                                          "Pomacanthidae_cons_change", 
                                          "Rajidae_cons_change", 
                                          "Rhinobatidae_cons_change", 
                                          "Soleidae_cons_change", 
                                          "Tetraodontidae_cons_change", 
                                          "Carangidae_cons_change", 
                                          "Balistidae_cons_change", 
                                          "Zanclidae_cons_change", 
                                          "Octopus_cons_change", 
                                          "Loligo_cons_change", 
                                          "Sepia_cons_change", 
                                          "Cypraea.Conus_cons_change", 
                                          "Pyrasus_cons_change", 
                                          "Charonia_cons_change", 
                                          "Murex.Fasciolaria_cons_change", 
                                          "Lambis_cons_change", 
                                          "Anadara_cons_change", 
                                          "Pinctada.Isognomon.Atrina.Pinna_cons_change", 
                                          "Tridacna_cons_change", 
                                          "Tripneustes_cons_change", 
                                          "Aristeidae_cons_change", 
                                          "Scylla_cons_change", 
                                          "Palunirus_cons_change"))

tableofmeansconsumption <- as.data.frame(colMeans(consumptionaverages, na.rm=TRUE))



rownames(tableofmeansconsumption)[rownames(tableofmeansconsumption) == "Mullidae_cons_change"] <- "Mullidae"
rownames(tableofmeansconsumption)[rownames(tableofmeansconsumption) == "Pomacentridae_cons_change"] <- "Pomacentridae"
rownames(tableofmeansconsumption)[rownames(tableofmeansconsumption) == "Acanthuridae_cons_change"] <- "Acanthuridae"
rownames(tableofmeansconsumption)[rownames(tableofmeansconsumption) == "Scaridae_cons_change"] <- "Scaridae"
rownames(tableofmeansconsumption)[rownames(tableofmeansconsumption) == "Scorpaenidae_cons_change"] <- "Scorpaenidae"
rownames(tableofmeansconsumption)[rownames(tableofmeansconsumption) == "Siganidae_cons_change"] <- "Siganidae"
rownames(tableofmeansconsumption)[rownames(tableofmeansconsumption) == "Sphyraenidae_cons_change"] <- "Sphyraenidae"
rownames(tableofmeansconsumption)[rownames(tableofmeansconsumption) == "Labridae_cons_change"] <- "Labridae"
rownames(tableofmeansconsumption)[rownames(tableofmeansconsumption) == "Priacanthidae_cons_change"] <- "Priacanthidae"
rownames(tableofmeansconsumption)[rownames(tableofmeansconsumption) == "Chaetodontidae_cons_change"] <- "Chaetodontidae"
rownames(tableofmeansconsumption)[rownames(tableofmeansconsumption) == "Fistularidae_cons_change"] <- "Fistularidae"
rownames(tableofmeansconsumption)[rownames(tableofmeansconsumption) == "Gerreidae_cons_change"] <- "Gerreidae"
rownames(tableofmeansconsumption)[rownames(tableofmeansconsumption) == "Haemulidae_cons_change"] <- "Haemulidae"
rownames(tableofmeansconsumption)[rownames(tableofmeansconsumption) == "Lethrinidae_cons_change"] <- "Lethrinidae"
rownames(tableofmeansconsumption)[rownames(tableofmeansconsumption) == "Ostraciidae_cons_change"] <- "Ostraciidae"
rownames(tableofmeansconsumption)[rownames(tableofmeansconsumption) == "Pomacanthidae_cons_change"] <- "Pomacanthidae"
rownames(tableofmeansconsumption)[rownames(tableofmeansconsumption) == "Rajidae_cons_change"] <- "Rajidae"
rownames(tableofmeansconsumption)[rownames(tableofmeansconsumption) == "Rhinobatidae_cons_change"] <- "Rhinobatidae"
rownames(tableofmeansconsumption)[rownames(tableofmeansconsumption) == "Soleidae_cons_change"] <- "Soleidae"
rownames(tableofmeansconsumption)[rownames(tableofmeansconsumption) == "Tetraodontidae_cons_change"] <- "Tetraodontidae"
rownames(tableofmeansconsumption)[rownames(tableofmeansconsumption) == "Carangidae_cons_change"] <- "Carangidae"
rownames(tableofmeansconsumption)[rownames(tableofmeansconsumption) == "Balistidae_cons_change"] <- "Balistidae"
rownames(tableofmeansconsumption)[rownames(tableofmeansconsumption) == "Zanclidae_cons_change"] <- "Zanclidae"
rownames(tableofmeansconsumption)[rownames(tableofmeansconsumption) == "Octopus_cons_change"] <- "Octopoda"
rownames(tableofmeansconsumption)[rownames(tableofmeansconsumption) == "Loligo_cons_change"] <- "Loligo"
rownames(tableofmeansconsumption)[rownames(tableofmeansconsumption) == "Sepia_cons_change"] <- "Sepia"
rownames(tableofmeansconsumption)[rownames(tableofmeansconsumption) == "Cypraea.Conus_cons_change"] <- "Cypraeidae/Conus"
rownames(tableofmeansconsumption)[rownames(tableofmeansconsumption) == "Pyrasus_cons_change"] <- "Pyrasus"
rownames(tableofmeansconsumption)[rownames(tableofmeansconsumption) == "Charonia_cons_change"] <- "Charonia"
rownames(tableofmeansconsumption)[rownames(tableofmeansconsumption) == "Murex.Fasciolaria_cons_change"] <- "Muricidae"
rownames(tableofmeansconsumption)[rownames(tableofmeansconsumption) == "Lambis_cons_change"] <- "Lambis"
rownames(tableofmeansconsumption)[rownames(tableofmeansconsumption) == "Anadara_cons_change"] <- "Anadara"
rownames(tableofmeansconsumption)[rownames(tableofmeansconsumption) == "Pinctada.Isognomon.Atrina.Pinna_cons_change"] <- "Pteriidae/Isognomonidae"
rownames(tableofmeansconsumption)[rownames(tableofmeansconsumption) == "Tridacna_cons_change"] <- "Tridacna"
rownames(tableofmeansconsumption)[rownames(tableofmeansconsumption) == "Tripneustes_cons_change"] <- "Tripneustes"
rownames(tableofmeansconsumption)[rownames(tableofmeansconsumption) == "Aristeidae_cons_change"] <- "Aristeidae"
rownames(tableofmeansconsumption)[rownames(tableofmeansconsumption) == "Scylla_cons_change"] <- "Scylla"
rownames(tableofmeansconsumption)[rownames(tableofmeansconsumption) == "Palunirus_cons_change"] <- "Palunirus"

colnames(tableofmeansconsumption)[colnames(tableofmeansconsumption) == "colMeans(consumptionaverages, na.rm = TRUE)"] <- "Average"

tableofmeansconsumption$row_names <- row.names(tableofmeansconsumption) 

write.csv(tableofmeansconsumption, "AverageConsumption.csv", row.names=FALSE)



jpeg("AverageConsumption", width =4000, height = 2000, res=300)
ggplot(tableofmeansconsumption, aes(x=fct_inorder(row_names), y=Average), )+
  geom_bar(position = "dodge", stat = "identity", width = 0.80, fill="darkred") +
  ggtitle("Average Perceived Change in Consumption over 10 Years by Marine Resource Group") +
  labs(x="Marine Resource Group", y="Perceived Change") +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_text(aes(y = Average + 0.1* sign(Average), label = round(Average, digits = 2)), position = position_dodge(width=0.75), vjust=1, check_overlap = TRUE, size=2.75) + 
  theme(plot.title = element_text(hjust = 0.5, margin=margin(b=20)))+
  theme(axis.title.x = element_text(margin=margin(t=20)), axis.title.y = element_text(margin=margin(r=20)))
dev.off()











abundanceaverages <- select(averages, c("Mullidae_ea_change", 
                                        "Pomacentridae_ea_change",
                                        "Acanthuridae_ea_change", 
                                        "Scaridae_ea_change", 
                                        "Scorpaenidae_ea_change", 
                                        "Siganidae_ea_change", 
                                        "Sphyraenidae_ea_change", 
                                        "Labridae_ea_change", 
                                        "Priacanthidae_ea_change", 
                                        "Chaetodontidae_ea_change", 
                                        "Fistularidae_ea_change", 
                                        "Gerreidae_ea_change", 
                                        "Haemulidae_ea_change", 
                                        "Lethrinidae_ea_change", 
                                        "Ostraciidae_ea_change", 
                                        "Pomacanthidae_ea_change", 
                                        "Rajidae_ea_change", 
                                        "Rhinobatidae_ea_change", 
                                        "Soleidae_ea_change", 
                                        "Tetraodontidae_ea_change", 
                                        "Carangidae_ea_change", 
                                        "Balistidae_ea_change", 
                                        "Zanclidae_ea_change", 
                                        "Octopus_ea_change", 
                                        "Loligo_ea_change", 
                                        "Sepia_ea_change", 
                                        "Cypraea.Conus_ea_change", 
                                        "Pyrasus_ea_change", 
                                        "Charonia_ea_change", 
                                        "Murex.Fasciolaria_ea_change", 
                                        "Lambis_ea_change", 
                                        "Anadara_ea_change", 
                                        "Pinctada.Isognomon.Atrina.Pinna_ea_change", 
                                        "Tridacna_ea_change", 
                                        "Tripneustes_ea_change", 
                                        "Aristeidae_ea_change", 
                                        "Scylla_ea_change", 
                                        "Palunirus_ea_change"))

tableofmeansabundance <- as.data.frame(colMeans(abundanceaverages, na.rm=TRUE))



rownames(tableofmeansabundance)[rownames(tableofmeansabundance) == "Mullidae_ea_change"] <- "Mullidae"
rownames(tableofmeansabundance)[rownames(tableofmeansabundance) == "Pomacentridae_ea_change"] <- "Pomacentridae"
rownames(tableofmeansabundance)[rownames(tableofmeansabundance) == "Acanthuridae_ea_change"] <- "Acanthuridae"
rownames(tableofmeansabundance)[rownames(tableofmeansabundance) == "Scaridae_ea_change"] <- "Scaridae"
rownames(tableofmeansabundance)[rownames(tableofmeansabundance) == "Scorpaenidae_ea_change"] <- "Scorpaenidae"
rownames(tableofmeansabundance)[rownames(tableofmeansabundance) == "Siganidae_ea_change"] <- "Siganidae"
rownames(tableofmeansabundance)[rownames(tableofmeansabundance) == "Sphyraenidae_ea_change"] <- "Sphyraenidae"
rownames(tableofmeansabundance)[rownames(tableofmeansabundance) == "Labridae_ea_change"] <- "Labridae"
rownames(tableofmeansabundance)[rownames(tableofmeansabundance) == "Priacanthidae_ea_change"] <- "Priacanthidae"
rownames(tableofmeansabundance)[rownames(tableofmeansabundance) == "Chaetodontidae_ea_change"] <- "Chaetodontidae"
rownames(tableofmeansabundance)[rownames(tableofmeansabundance) == "Fistularidae_ea_change"] <- "Fistularidae"
rownames(tableofmeansabundance)[rownames(tableofmeansabundance) == "Gerreidae_ea_change"] <- "Gerreidae"
rownames(tableofmeansabundance)[rownames(tableofmeansabundance) == "Haemulidae_ea_change"] <- "Haemulidae"
rownames(tableofmeansabundance)[rownames(tableofmeansabundance) == "Lethrinidae_ea_change"] <- "Lethrinidae"
rownames(tableofmeansabundance)[rownames(tableofmeansabundance) == "Ostraciidae_ea_change"] <- "Ostraciidae"
rownames(tableofmeansabundance)[rownames(tableofmeansabundance) == "Pomacanthidae_ea_change"] <- "Pomacanthidae"
rownames(tableofmeansabundance)[rownames(tableofmeansabundance) == "Rajidae_ea_change"] <- "Rajidae"
rownames(tableofmeansabundance)[rownames(tableofmeansabundance) == "Rhinobatidae_ea_change"] <- "Rhinobatidae"
rownames(tableofmeansabundance)[rownames(tableofmeansabundance) == "Soleidae_ea_change"] <- "Soleidae"
rownames(tableofmeansabundance)[rownames(tableofmeansabundance) == "Tetraodontidae_ea_change"] <- "Tetraodontidae"
rownames(tableofmeansabundance)[rownames(tableofmeansabundance) == "Carangidae_ea_change"] <- "Carangidae"
rownames(tableofmeansabundance)[rownames(tableofmeansabundance) == "Balistidae_ea_change"] <- "Balistidae"
rownames(tableofmeansabundance)[rownames(tableofmeansabundance) == "Zanclidae_ea_change"] <- "Zanclidae"
rownames(tableofmeansabundance)[rownames(tableofmeansabundance) == "Octopus_ea_change"] <- "Octopoda"
rownames(tableofmeansabundance)[rownames(tableofmeansabundance) == "Loligo_ea_change"] <- "Loligo"
rownames(tableofmeansabundance)[rownames(tableofmeansabundance) == "Sepia_ea_change"] <- "Sepia"
rownames(tableofmeansabundance)[rownames(tableofmeansabundance) == "Cypraea.Conus_ea_change"] <- "Cypraeidae/Conus"
rownames(tableofmeansabundance)[rownames(tableofmeansabundance) == "Pyrasus_ea_change"] <- "Pyrasus"
rownames(tableofmeansabundance)[rownames(tableofmeansabundance) == "Charonia_ea_change"] <- "Charonia"
rownames(tableofmeansabundance)[rownames(tableofmeansabundance) == "Murex.Fasciolaria_ea_change"] <- "Muricidae"
rownames(tableofmeansabundance)[rownames(tableofmeansabundance) == "Lambis_ea_change"] <- "Lambis"
rownames(tableofmeansabundance)[rownames(tableofmeansabundance) == "Anadara_ea_change"] <- "Anadara"
rownames(tableofmeansabundance)[rownames(tableofmeansabundance) == "Pinctada.Isognomon.Atrina.Pinna_ea_change"] <- "Pteriidae/Isognomonidae"
rownames(tableofmeansabundance)[rownames(tableofmeansabundance) == "Tridacna_ea_change"] <- "Tridacna"
rownames(tableofmeansabundance)[rownames(tableofmeansabundance) == "Tripneustes_ea_change"] <- "Tripneustes"
rownames(tableofmeansabundance)[rownames(tableofmeansabundance) == "Aristeidae_ea_change"] <- "Aristeidae"
rownames(tableofmeansabundance)[rownames(tableofmeansabundance) == "Scylla_ea_change"] <- "Scylla"
rownames(tableofmeansabundance)[rownames(tableofmeansabundance) == "Palunirus_ea_change"] <- "Palunirus"

colnames(tableofmeansabundance)[colnames(tableofmeansabundance) == "colMeans(abundanceaverages, na.rm = TRUE)"] <- "Average"

tableofmeansabundance$row_names <- row.names(tableofmeansabundance) 

write.csv(tableofmeansabundance, "AverageAbundance.csv", row.names=FALSE)


jpeg("AverageAbundance", width =4000, height = 2000, res=300)
ggplot(tableofmeansabundance, aes(x=fct_inorder(row_names), y=Average), )+
  geom_bar(position = "dodge", stat = "identity", width = 0.80, fill="forestgreen") +
  ggtitle("Average Perceived Change in Ecological Abundance over 10 Years by Marine Resource Group") +
  labs(x="Marine Resource Group", y="Perceived Change") +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_text(aes(y = Average + 0.1* sign(Average), label = round(Average, digits = 2)), position = position_dodge(width=0.75), vjust=1, check_overlap = TRUE, size=2.75) + 
  theme(plot.title = element_text(hjust = 0.5, margin=margin(b=20)))+
  theme(axis.title.x = element_text(margin=margin(t=20)), axis.title.y = element_text(margin=margin(r=20)))
dev.off()









harvestaverages <- select(averages, c("Mullidae_harv_change", 
                                      "Pomacentridae_harv_change",
                                      "Acanthuridae_harv_change", 
                                      "Scaridae_harv_change", 
                                      "Scorpaenidae_harv_change", 
                                      "Siganidae_harv_change", 
                                      "Sphyraenidae_harv_change", 
                                      "Labridae_harv_change", 
                                      "Priacanthidae_harv_change", 
                                      "Chaetodontidae_harv_change", 
                                      "Fistularidae_harv_change", 
                                      "Gerreidae_harv_change", 
                                      "Haemulidae_harv_change", 
                                      "Lethrinidae_harv_change", 
                                      "Ostraciidae_harv_change", 
                                      "Pomacanthidae_harv_change", 
                                      "Rajidae_harv_change", 
                                      "Rhinobatidae_harv_change", 
                                      "Soleidae_harv_change", 
                                      "Tetraodontidae_harv_change", 
                                      "Carangidae_harv_change", 
                                      "Balistidae_harv_change", 
                                      "Zanclidae_harv_change", 
                                      "Octopus_harv_change", 
                                      "Loligo_harv_change", 
                                      "Sepia_harv_change", 
                                      "Cypraea.Conus_harv_change", 
                                      "Pyrasus_harv_change", 
                                      "Charonia_harv_change", 
                                      "Murex.Fasciolaria_harv_change", 
                                      "Lambis_harv_change", 
                                      "Anadara_harv_change", 
                                      "Pinctada.Isognomon.Atrina.Pinna_harv_change", 
                                      "Tridacna_harv_change", 
                                      "Tripneustes_harv_change", 
                                      "Aristeidae_harv_change", 
                                      "Scylla_harv_change", 
                                      "Palunirus_harv_change"))

tableofmeansharvest <- as.data.frame(colMeans(harvestaverages, na.rm=TRUE))



rownames(tableofmeansharvest)[rownames(tableofmeansharvest) == "Mullidae_harv_change"] <- "Mullidae"
rownames(tableofmeansharvest)[rownames(tableofmeansharvest) == "Pomacentridae_harv_change"] <- "Pomacentridae"
rownames(tableofmeansharvest)[rownames(tableofmeansharvest) == "Acanthuridae_harv_change"] <- "Acanthuridae"
rownames(tableofmeansharvest)[rownames(tableofmeansharvest) == "Scaridae_harv_change"] <- "Scaridae"
rownames(tableofmeansharvest)[rownames(tableofmeansharvest) == "Scorpaenidae_harv_change"] <- "Scorpaenidae"
rownames(tableofmeansharvest)[rownames(tableofmeansharvest) == "Siganidae_harv_change"] <- "Siganidae"
rownames(tableofmeansharvest)[rownames(tableofmeansharvest) == "Sphyraenidae_harv_change"] <- "Sphyraenidae"
rownames(tableofmeansharvest)[rownames(tableofmeansharvest) == "Labridae_harv_change"] <- "Labridae"
rownames(tableofmeansharvest)[rownames(tableofmeansharvest) == "Priacanthidae_harv_change"] <- "Priacanthidae"
rownames(tableofmeansharvest)[rownames(tableofmeansharvest) == "Chaetodontidae_harv_change"] <- "Chaetodontidae"
rownames(tableofmeansharvest)[rownames(tableofmeansharvest) == "Fistularidae_harv_change"] <- "Fistularidae"
rownames(tableofmeansharvest)[rownames(tableofmeansharvest) == "Gerreidae_harv_change"] <- "Gerreidae"
rownames(tableofmeansharvest)[rownames(tableofmeansharvest) == "Haemulidae_harv_change"] <- "Haemulidae"
rownames(tableofmeansharvest)[rownames(tableofmeansharvest) == "Lethrinidae_harv_change"] <- "Lethrinidae"
rownames(tableofmeansharvest)[rownames(tableofmeansharvest) == "Ostraciidae_harv_change"] <- "Ostraciidae"
rownames(tableofmeansharvest)[rownames(tableofmeansharvest) == "Pomacanthidae_harv_change"] <- "Pomacanthidae"
rownames(tableofmeansharvest)[rownames(tableofmeansharvest) == "Rajidae_harv_change"] <- "Rajidae"
rownames(tableofmeansharvest)[rownames(tableofmeansharvest) == "Rhinobatidae_harv_change"] <- "Rhinobatidae"
rownames(tableofmeansharvest)[rownames(tableofmeansharvest) == "Soleidae_harv_change"] <- "Soleidae"
rownames(tableofmeansharvest)[rownames(tableofmeansharvest) == "Tetraodontidae_harv_change"] <- "Tetraodontidae"
rownames(tableofmeansharvest)[rownames(tableofmeansharvest) == "Carangidae_harv_change"] <- "Carangidae"
rownames(tableofmeansharvest)[rownames(tableofmeansharvest) == "Balistidae_harv_change"] <- "Balistidae"
rownames(tableofmeansharvest)[rownames(tableofmeansharvest) == "Zanclidae_harv_change"] <- "Zanclidae"
rownames(tableofmeansharvest)[rownames(tableofmeansharvest) == "Octopus_harv_change"] <- "Octopoda"
rownames(tableofmeansharvest)[rownames(tableofmeansharvest) == "Loligo_harv_change"] <- "Loligo"
rownames(tableofmeansharvest)[rownames(tableofmeansharvest) == "Sepia_harv_change"] <- "Sepia"
rownames(tableofmeansharvest)[rownames(tableofmeansharvest) == "Cypraea.Conus_harv_change"] <- "Cypraeidae/Conus"
rownames(tableofmeansharvest)[rownames(tableofmeansharvest) == "Pyrasus_harv_change"] <- "Pyrasus"
rownames(tableofmeansharvest)[rownames(tableofmeansharvest) == "Charonia_harv_change"] <- "Charonia"
rownames(tableofmeansharvest)[rownames(tableofmeansharvest) == "Murex.Fasciolaria_harv_change"] <- "Muricidae"
rownames(tableofmeansharvest)[rownames(tableofmeansharvest) == "Lambis_harv_change"] <- "Lambis"
rownames(tableofmeansharvest)[rownames(tableofmeansharvest) == "Anadara_harv_change"] <- "Anadara"
rownames(tableofmeansharvest)[rownames(tableofmeansharvest) == "Pinctada.Isognomon.Atrina.Pinna_harv_change"] <- "Pteriidae/Isognomonidae"
rownames(tableofmeansharvest)[rownames(tableofmeansharvest) == "Tridacna_harv_change"] <- "Tridacna"
rownames(tableofmeansharvest)[rownames(tableofmeansharvest) == "Tripneustes_harv_change"] <- "Tripneustes"
rownames(tableofmeansharvest)[rownames(tableofmeansharvest) == "Aristeidae_harv_change"] <- "Aristeidae"
rownames(tableofmeansharvest)[rownames(tableofmeansharvest) == "Scylla_harv_change"] <- "Scylla"
rownames(tableofmeansharvest)[rownames(tableofmeansharvest) == "Palunirus_harv_change"] <- "Palunirus"

colnames(tableofmeansharvest)[colnames(tableofmeansharvest) == "colMeans(harvestaverages, na.rm = TRUE)"] <- "Average"

tableofmeansharvest$row_names <- row.names(tableofmeansharvest) 

write.csv(tableofmeansharvest, "AverageHarvest.csv", row.names=FALSE)


jpeg("AverageHarvest", width =4000, height = 2000, res=300)
ggplot(tableofmeansharvest, aes(x=fct_inorder(row_names), y=Average), )+
  geom_bar(position = "dodge", stat = "identity", width = 0.80, fill="blue3") +
  ggtitle("Average Perceived Change in Harvest over 10 Years by Marine Resource Group") +
  labs(x="Marine Resource Group", y="Perceived Change") +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_text(aes(y = Average + 0.1* sign(Average), label = round(Average, digits = 2)), position = position_dodge(width=0.75), vjust=1, check_overlap = TRUE, size=2.75) + 
  theme(plot.title = element_text(hjust = 0.5, margin=margin(b=20)))+
  theme(axis.title.x = element_text(margin=margin(t=20)), axis.title.y = element_text(margin=margin(r=20)))
dev.off()






##Randomeffects Consumption

Mullidaerancons<-as.data.frame(ranef(testMullidaecons, probs = c(0.1,0.9)))
Mullidaerancons$Village<-row.names(Mullidaerancons)
Mullidaerancons$group<- "Mullidae"
colnames(Mullidaerancons)<-c("estimate","se","conf.low","conf.high","Village", "group")

Pomacentridaerancons<-as.data.frame(ranef(testPomacentridaecons, probs = c(0.1,0.9)))
Pomacentridaerancons$Village<-row.names(Pomacentridaerancons)
Pomacentridaerancons$group<- "Pomacentridae"
colnames(Pomacentridaerancons)<-c("estimate","se","conf.low","conf.high","Village", "group")

Acanthuridaerancons<-as.data.frame(ranef(testAcanthuridaecons, probs = c(0.1,0.9)))
Acanthuridaerancons$Village<-row.names(Acanthuridaerancons)
Acanthuridaerancons$group<- "Acanthuridae"
colnames(Acanthuridaerancons)<-c("estimate","se","conf.low","conf.high","Village", "group")

Scaridaerancons<-as.data.frame(ranef(testScaridaecons, probs = c(0.1,0.9)))
Scaridaerancons$Village<-row.names(Scaridaerancons)
Scaridaerancons$group<- "Scaridae"
colnames(Scaridaerancons)<-c("estimate","se","conf.low","conf.high","Village", "group")


Scorpaenidaerancons<-as.data.frame(ranef(testScorpaenidaecons, probs = c(0.1,0.9)))
Scorpaenidaerancons$Village<-row.names(Scorpaenidaerancons)
Scorpaenidaerancons$group<- "Scorpaenidae"
colnames(Scorpaenidaerancons)<-c("estimate","se","conf.low","conf.high","Village", "group")

Siganidaerancons<-as.data.frame(ranef(testSiganidaecons, probs = c(0.1,0.9)))
Siganidaerancons$Village<-row.names(Siganidaerancons)
Siganidaerancons$group<- "Siganidae"
colnames(Siganidaerancons)<-c("estimate","se","conf.low","conf.high","Village", "group")


Sphyraenidaerancons<-as.data.frame(ranef(testSphyraenidaecons, probs = c(0.1,0.9)))
Sphyraenidaerancons$Village<-row.names(Sphyraenidaerancons)
Sphyraenidaerancons$group<- "Sphyraenidae"
colnames(Sphyraenidaerancons)<-c("estimate","se","conf.low","conf.high","Village", "group")

Labridaerancons<-as.data.frame(ranef(testLabridaecons, probs = c(0.1,0.9)))
Labridaerancons$Village<-row.names(Labridaerancons)
Labridaerancons$group<- "Labridae"
colnames(Labridaerancons)<-c("estimate","se","conf.low","conf.high","Village", "group")


Priacanthidaerancons<-as.data.frame(ranef(testPriacanthidaecons, probs = c(0.1,0.9)))
Priacanthidaerancons$Village<-row.names(Priacanthidaerancons)
Priacanthidaerancons$group<- "Priacanthidae"
colnames(Priacanthidaerancons)<-c("estimate","se","conf.low","conf.high","Village", "group")

Chaetodontidaerancons<-as.data.frame(ranef(testChaetodontidaecons, probs = c(0.1,0.9)))
Chaetodontidaerancons$Village<-row.names(Chaetodontidaerancons)
Chaetodontidaerancons$group<- "Chaetodontidae"
colnames(Chaetodontidaerancons)<-c("estimate","se","conf.low","conf.high","Village", "group")

Fistularidaerancons<-as.data.frame(ranef(testFistularidaecons, probs = c(0.1,0.9)))
Fistularidaerancons$Village<-row.names(Fistularidaerancons)
Fistularidaerancons$group<- "Fistularidae"
colnames(Fistularidaerancons)<-c("estimate","se","conf.low","conf.high","Village", "group")

Gerreidaerancons<-as.data.frame(ranef(testGerreidaecons, probs = c(0.1,0.9)))
Gerreidaerancons$Village<-row.names(Gerreidaerancons)
Gerreidaerancons$group<- "Gerreidae"
colnames(Gerreidaerancons)<-c("estimate","se","conf.low","conf.high","Village", "group")



Haemulidaerancons<-as.data.frame(ranef(testHaemulidaecons, probs = c(0.1,0.9)))
Haemulidaerancons$Village<-row.names(Haemulidaerancons)
Haemulidaerancons$group<- "Haemulidae"
colnames(Haemulidaerancons)<-c("estimate","se","conf.low","conf.high","Village", "group")


Lethrinidaerancons<-as.data.frame(ranef(testLethrinidaecons, probs = c(0.1,0.9)))
Lethrinidaerancons$Village<-row.names(Lethrinidaerancons)
Lethrinidaerancons$group<- "Lethrinidae"
colnames(Lethrinidaerancons)<-c("estimate","se","conf.low","conf.high","Village", "group")



Ostraciidaerancons<-as.data.frame(ranef(testOstraciidaecons, probs = c(0.1,0.9)))
Ostraciidaerancons$Village<-row.names(Ostraciidaerancons)
Ostraciidaerancons$group<- "Ostraciidae"
colnames(Ostraciidaerancons)<-c("estimate","se","conf.low","conf.high","Village", "group")



Pomacanthidaerancons<-as.data.frame(ranef(testPomacanthidaecons, probs = c(0.1,0.9)))
Pomacanthidaerancons$Village<-row.names(Pomacanthidaerancons)
Pomacanthidaerancons$group<- "Pomacanthidae"
colnames(Pomacanthidaerancons)<-c("estimate","se","conf.low","conf.high","Village", "group")


Rajidaerancons<-as.data.frame(ranef(testRajidaecons, probs = c(0.1,0.9)))
Rajidaerancons$Village<-row.names(Rajidaerancons)
Rajidaerancons$group<- "Rajidae"
colnames(Rajidaerancons)<-c("estimate","se","conf.low","conf.high","Village", "group")



Rhinobatidaerancons<-as.data.frame(ranef(testRhinobatidaecons, probs = c(0.1,0.9)))
Rhinobatidaerancons$Village<-row.names(Rhinobatidaerancons)
Rhinobatidaerancons$group<- "Rhinobatidae"
colnames(Rhinobatidaerancons)<-c("estimate","se","conf.low","conf.high","Village", "group")


Soleidaerancons<-as.data.frame(ranef(testSoleidaecons, probs = c(0.1,0.9)))
Soleidaerancons$Village<-row.names(Soleidaerancons)
Soleidaerancons$group<- "Soleidae"
colnames(Soleidaerancons)<-c("estimate","se","conf.low","conf.high","Village", "group")


Tetraodontidaerancons<-as.data.frame(ranef(testTetraodontidaecons, probs = c(0.1,0.9)))
Tetraodontidaerancons$Village<-row.names(Tetraodontidaerancons)
Tetraodontidaerancons$group<- "Tetraodontidae"
colnames(Tetraodontidaerancons)<-c("estimate","se","conf.low","conf.high","Village", "group")



Carangidaerancons<-as.data.frame(ranef(testCarangidaecons, probs = c(0.1,0.9)))
Carangidaerancons$Village<-row.names(Carangidaerancons)
Carangidaerancons$group<- "Carangidae"
colnames(Carangidaerancons)<-c("estimate","se","conf.low","conf.high","Village", "group")


Balistidaerancons<-as.data.frame(ranef(testBalistidaecons, probs = c(0.1,0.9)))
Balistidaerancons$Village<-row.names(Balistidaerancons)
Balistidaerancons$group<- "Balistidae"
colnames(Balistidaerancons)<-c("estimate","se","conf.low","conf.high","Village", "group")



Zanclidaerancons<-as.data.frame(ranef(testZanclidaecons, probs = c(0.1,0.9)))
Zanclidaerancons$Village<-row.names(Zanclidaerancons)
Zanclidaerancons$group<- "Zanclidae"
colnames(Zanclidaerancons)<-c("estimate","se","conf.low","conf.high","Village", "group")




Octopusrancons<-as.data.frame(ranef(testOctopuscons, probs = c(0.1,0.9)))
Octopusrancons$Village<-row.names(Octopusrancons)
Octopusrancons$group<- "Octopus"
colnames(Octopusrancons)<-c("estimate","se","conf.low","conf.high","Village", "group")




Loligorancons<-as.data.frame(ranef(testLoligocons, probs = c(0.1,0.9)))
Loligorancons$Village<-row.names(Loligorancons)
Loligorancons$group<- "Loligo"
colnames(Loligorancons)<-c("estimate","se","conf.low","conf.high","Village", "group")




Sepiarancons<-as.data.frame(ranef(testSepiacons, probs = c(0.1,0.9)))
Sepiarancons$Village<-row.names(Sepiarancons)
Sepiarancons$group<- "Sepia"
colnames(Sepiarancons)<-c("estimate","se","conf.low","conf.high","Village", "group")



Cypraea.Conusrancons<-as.data.frame(ranef(testCypraea.Conuscons, probs = c(0.1,0.9)))
Cypraea.Conusrancons$Village<-row.names(Cypraea.Conusrancons)
Cypraea.Conusrancons$group<- "Cypraea.Conus"
colnames(Cypraea.Conusrancons)<-c("estimate","se","conf.low","conf.high","Village", "group")



Pyrasusrancons<-as.data.frame(ranef(testPyrasuscons, probs = c(0.1,0.9)))
Pyrasusrancons$Village<-row.names(Pyrasusrancons)
Pyrasusrancons$group<- "Pyrasus"
colnames(Pyrasusrancons)<-c("estimate","se","conf.low","conf.high","Village", "group")



Charoniarancons<-as.data.frame(ranef(testCharoniacons, probs = c(0.1,0.9)))
Charoniarancons$Village<-row.names(Charoniarancons)
Charoniarancons$group<- "Charonia"
colnames(Charoniarancons)<-c("estimate","se","conf.low","conf.high","Village", "group")


Murex.Fasciolariarancons<-as.data.frame(ranef(testMurex.Fasciolariacons, probs = c(0.1,0.9)))
Murex.Fasciolariarancons$Village<-row.names(Murex.Fasciolariarancons)
Murex.Fasciolariarancons$group<- "Murex.Fasciolaria"
colnames(Murex.Fasciolariarancons)<-c("estimate","se","conf.low","conf.high","Village", "group")



Lambisrancons<-as.data.frame(ranef(testLambiscons, probs = c(0.1,0.9)))
Lambisrancons$Village<-row.names(Lambisrancons)
Lambisrancons$group<- "Lambis"
colnames(Lambisrancons)<-c("estimate","se","conf.low","conf.high","Village", "group")


Anadararancons<-as.data.frame(ranef(testAnadaracons, probs = c(0.1,0.9)))
Anadararancons$Village<-row.names(Anadararancons)
Anadararancons$group<- "Anadara"
colnames(Anadararancons)<-c("estimate","se","conf.low","conf.high","Village", "group")



Pinctada.Isognomon.Atrina.Pinnarancons<-as.data.frame(ranef(testPinctada.Isognomon.Atrina.Pinnacons, probs = c(0.1,0.9)))
Pinctada.Isognomon.Atrina.Pinnarancons$Village<-row.names(Pinctada.Isognomon.Atrina.Pinnarancons)
Pinctada.Isognomon.Atrina.Pinnarancons$group<- "Pinctada.Isognomon.Atrina.Pinna"
colnames(Pinctada.Isognomon.Atrina.Pinnarancons)<-c("estimate","se","conf.low","conf.high","Village", "group")



Tridacnarancons<-as.data.frame(ranef(testTridacnacons, probs = c(0.1,0.9)))
Tridacnarancons$Village<-row.names(Tridacnarancons)
Tridacnarancons$group<- "Tridacna"
colnames(Tridacnarancons)<-c("estimate","se","conf.low","conf.high","Village", "group")


Tripneustesrancons<-as.data.frame(ranef(testTripneustescons, probs = c(0.1,0.9)))
Tripneustesrancons$Village<-row.names(Tripneustesrancons)
Tripneustesrancons$group<- "Tripneustes"
colnames(Tripneustesrancons)<-c("estimate","se","conf.low","conf.high","Village", "group")


Aristeidaerancons<-as.data.frame(ranef(testAristeidaecons, probs = c(0.1,0.9)))
Aristeidaerancons$Village<-row.names(Aristeidaerancons)
Aristeidaerancons$group<- "Aristeidae"
colnames(Aristeidaerancons)<-c("estimate","se","conf.low","conf.high","Village", "group")



Scyllarancons<-as.data.frame(ranef(testScyllacons, probs = c(0.1,0.9)))
Scyllarancons$Village<-row.names(Scyllarancons)
Scyllarancons$group<- "Scylla"
colnames(Scyllarancons)<-c("estimate","se","conf.low","conf.high","Village", "group")



Palunirusrancons<-as.data.frame(ranef(testPaluniruscons, probs = c(0.1,0.9)))
Palunirusrancons$Village<-row.names(Palunirusrancons)
Palunirusrancons$group<- "Palunirus"
colnames(Palunirusrancons)<-c("estimate","se","conf.low","conf.high","Village", "group")




Randomeffectsizecons <- rbind(Mullidaerancons,
                              Pomacentridaerancons,
                              Acanthuridaerancons,
                              Scaridaerancons,
                              Scorpaenidaerancons,
                              Siganidaerancons,
                              Sphyraenidaerancons,
                              Labridaerancons,
                              Priacanthidaerancons,
                              Chaetodontidaerancons,
                              Fistularidaerancons,
                              Gerreidaerancons,
                              Haemulidaerancons,
                              Lethrinidaerancons,
                              Ostraciidaerancons,
                              Pomacanthidaerancons,
                              Rajidaerancons,
                              Rhinobatidaerancons,
                              Soleidaerancons,
                              Tetraodontidaerancons,
                              Carangidaerancons,
                              Balistidaerancons,
                              Zanclidaerancons,
                              Octopusrancons,
                              Loligorancons,
                              Sepiarancons,
                              Cypraea.Conusrancons,
                              Pyrasusrancons,
                              Charoniarancons,
                              Murex.Fasciolariarancons,
                              Lambisrancons,
                              Anadararancons,
                              Pinctada.Isognomon.Atrina.Pinnarancons,
                              Tridacnarancons,
                              Tripneustesrancons,
                              Aristeidaerancons,
                              Scyllarancons,
                              Palunirusrancons)


Randomeffectsizecons[Randomeffectsizecons == "Cypraea.Conus"] <- "Cypraeidae/Conus"
Randomeffectsizecons[Randomeffectsizecons == "Pinctada.Isognomon.Atrina.Pinna"] <- "Pteriidae/Isognomonidae"
Randomeffectsizecons[Randomeffectsizecons == "Murex.Fasciolaria"] <- "Muricidae"
Randomeffectsizecons[Randomeffectsizecons == "Octopus"] <- "Octopoda"

write.csv(Randomeffectsizecons, "ConsEffectSizesRANDOM.csv", row.names=FALSE)

Randomeffectsizecons <- read.csv("ConsEffectSizesRANDOM.csv")



jpeg("VillageEffectSizeConsumption", width =2000, height = 5000, res=300)
ggplot(Randomeffectsizecons,aes(x=estimate, fill=Village, y=fct_rev(fct_inorder(group)),xmin=conf.low,xmax=conf.high))+geom_errorbar(position=position_dodge(width=1))+geom_point(pch=21,size=3, position=position_dodge(width =1))+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  geom_smooth() +
  labs(title = str_wrap("Effect of Village on Perception of Change in Consumption of each Marine Resource Group", 45))+
  labs(x="Effect Size by Village", y="Marine Resource Group")+
  theme(plot.title = element_text(hjust = 0.5, margin=margin(b=20)))+
  theme(axis.title.x = element_text(margin=margin(t=20)), axis.title.y = element_text(margin=margin(r=20)))
dev.off()









##Randomeffects eaumption

Mullidaeranea<-as.data.frame(ranef(testMullidaeea, probs = c(0.1,0.9)))
Mullidaeranea$Village<-row.names(Mullidaeranea)
Mullidaeranea$group<- "Mullidae"
colnames(Mullidaeranea)<-c("estimate","se","conf.low","conf.high","Village", "group")

Pomacentridaeranea<-as.data.frame(ranef(testPomacentridaeea, probs = c(0.1,0.9)))
Pomacentridaeranea$Village<-row.names(Pomacentridaeranea)
Pomacentridaeranea$group<- "Pomacentridae"
colnames(Pomacentridaeranea)<-c("estimate","se","conf.low","conf.high","Village", "group")

Acanthuridaeranea<-as.data.frame(ranef(testAcanthuridaeea, probs = c(0.1,0.9)))
Acanthuridaeranea$Village<-row.names(Acanthuridaeranea)
Acanthuridaeranea$group<- "Acanthuridae"
colnames(Acanthuridaeranea)<-c("estimate","se","conf.low","conf.high","Village", "group")

Scaridaeranea<-as.data.frame(ranef(testScaridaeea, probs = c(0.1,0.9)))
Scaridaeranea$Village<-row.names(Scaridaeranea)
Scaridaeranea$group<- "Scaridae"
colnames(Scaridaeranea)<-c("estimate","se","conf.low","conf.high","Village", "group")


Scorpaenidaeranea<-as.data.frame(ranef(testScorpaenidaeea, probs = c(0.1,0.9)))
Scorpaenidaeranea$Village<-row.names(Scorpaenidaeranea)
Scorpaenidaeranea$group<- "Scorpaenidae"
colnames(Scorpaenidaeranea)<-c("estimate","se","conf.low","conf.high","Village", "group")

Siganidaeranea<-as.data.frame(ranef(testSiganidaeea, probs = c(0.1,0.9)))
Siganidaeranea$Village<-row.names(Siganidaeranea)
Siganidaeranea$group<- "Siganidae"
colnames(Siganidaeranea)<-c("estimate","se","conf.low","conf.high","Village", "group")


Sphyraenidaeranea<-as.data.frame(ranef(testSphyraenidaeea, probs = c(0.1,0.9)))
Sphyraenidaeranea$Village<-row.names(Sphyraenidaeranea)
Sphyraenidaeranea$group<- "Sphyraenidae"
colnames(Sphyraenidaeranea)<-c("estimate","se","conf.low","conf.high","Village", "group")

Labridaeranea<-as.data.frame(ranef(testLabridaeea, probs = c(0.1,0.9)))
Labridaeranea$Village<-row.names(Labridaeranea)
Labridaeranea$group<- "Labridae"
colnames(Labridaeranea)<-c("estimate","se","conf.low","conf.high","Village", "group")


Priacanthidaeranea<-as.data.frame(ranef(testPriacanthidaeea, probs = c(0.1,0.9)))
Priacanthidaeranea$Village<-row.names(Priacanthidaeranea)
Priacanthidaeranea$group<- "Priacanthidae"
colnames(Priacanthidaeranea)<-c("estimate","se","conf.low","conf.high","Village", "group")

Chaetodontidaeranea<-as.data.frame(ranef(testChaetodontidaeea, probs = c(0.1,0.9)))
Chaetodontidaeranea$Village<-row.names(Chaetodontidaeranea)
Chaetodontidaeranea$group<- "Chaetodontidae"
colnames(Chaetodontidaeranea)<-c("estimate","se","conf.low","conf.high","Village", "group")

Fistularidaeranea<-as.data.frame(ranef(testFistularidaeea, probs = c(0.1,0.9)))
Fistularidaeranea$Village<-row.names(Fistularidaeranea)
Fistularidaeranea$group<- "Fistularidae"
colnames(Fistularidaeranea)<-c("estimate","se","conf.low","conf.high","Village", "group")

Gerreidaeranea<-as.data.frame(ranef(testGerreidaeea, probs = c(0.1,0.9)))
Gerreidaeranea$Village<-row.names(Gerreidaeranea)
Gerreidaeranea$group<- "Gerreidae"
colnames(Gerreidaeranea)<-c("estimate","se","conf.low","conf.high","Village", "group")



Haemulidaeranea<-as.data.frame(ranef(testHaemulidaeea, probs = c(0.1,0.9)))
Haemulidaeranea$Village<-row.names(Haemulidaeranea)
Haemulidaeranea$group<- "Haemulidae"
colnames(Haemulidaeranea)<-c("estimate","se","conf.low","conf.high","Village", "group")


Lethrinidaeranea<-as.data.frame(ranef(testLethrinidaeea, probs = c(0.1,0.9)))
Lethrinidaeranea$Village<-row.names(Lethrinidaeranea)
Lethrinidaeranea$group<- "Lethrinidae"
colnames(Lethrinidaeranea)<-c("estimate","se","conf.low","conf.high","Village", "group")



Ostraciidaeranea<-as.data.frame(ranef(testOstraciidaeea, probs = c(0.1,0.9)))
Ostraciidaeranea$Village<-row.names(Ostraciidaeranea)
Ostraciidaeranea$group<- "Ostraciidae"
colnames(Ostraciidaeranea)<-c("estimate","se","conf.low","conf.high","Village", "group")



Pomacanthidaeranea<-as.data.frame(ranef(testPomacanthidaeea, probs = c(0.1,0.9)))
Pomacanthidaeranea$Village<-row.names(Pomacanthidaeranea)
Pomacanthidaeranea$group<- "Pomacanthidae"
colnames(Pomacanthidaeranea)<-c("estimate","se","conf.low","conf.high","Village", "group")


Rajidaeranea<-as.data.frame(ranef(testRajidaeea, probs = c(0.1,0.9)))
Rajidaeranea$Village<-row.names(Rajidaeranea)
Rajidaeranea$group<- "Rajidae"
colnames(Rajidaeranea)<-c("estimate","se","conf.low","conf.high","Village", "group")



Rhinobatidaeranea<-as.data.frame(ranef(testRhinobatidaeea, probs = c(0.1,0.9)))
Rhinobatidaeranea$Village<-row.names(Rhinobatidaeranea)
Rhinobatidaeranea$group<- "Rhinobatidae"
colnames(Rhinobatidaeranea)<-c("estimate","se","conf.low","conf.high","Village", "group")


Soleidaeranea<-as.data.frame(ranef(testSoleidaeea, probs = c(0.1,0.9)))
Soleidaeranea$Village<-row.names(Soleidaeranea)
Soleidaeranea$group<- "Soleidae"
colnames(Soleidaeranea)<-c("estimate","se","conf.low","conf.high","Village", "group")


Tetraodontidaeranea<-as.data.frame(ranef(testTetraodontidaeea, probs = c(0.1,0.9)))
Tetraodontidaeranea$Village<-row.names(Tetraodontidaeranea)
Tetraodontidaeranea$group<- "Tetraodontidae"
colnames(Tetraodontidaeranea)<-c("estimate","se","conf.low","conf.high","Village", "group")



Carangidaeranea<-as.data.frame(ranef(testCarangidaeea, probs = c(0.1,0.9)))
Carangidaeranea$Village<-row.names(Carangidaeranea)
Carangidaeranea$group<- "Carangidae"
colnames(Carangidaeranea)<-c("estimate","se","conf.low","conf.high","Village", "group")


Balistidaeranea<-as.data.frame(ranef(testBalistidaeea, probs = c(0.1,0.9)))
Balistidaeranea$Village<-row.names(Balistidaeranea)
Balistidaeranea$group<- "Balistidae"
colnames(Balistidaeranea)<-c("estimate","se","conf.low","conf.high","Village", "group")



Zanclidaeranea<-as.data.frame(ranef(testZanclidaeea, probs = c(0.1,0.9)))
Zanclidaeranea$Village<-row.names(Zanclidaeranea)
Zanclidaeranea$group<- "Zanclidae"
colnames(Zanclidaeranea)<-c("estimate","se","conf.low","conf.high","Village", "group")




Octopusranea<-as.data.frame(ranef(testOctopusea, probs = c(0.1,0.9)))
Octopusranea$Village<-row.names(Octopusranea)
Octopusranea$group<- "Octopus"
colnames(Octopusranea)<-c("estimate","se","conf.low","conf.high","Village", "group")




Loligoranea<-as.data.frame(ranef(testLoligoea, probs = c(0.1,0.9)))
Loligoranea$Village<-row.names(Loligoranea)
Loligoranea$group<- "Loligo"
colnames(Loligoranea)<-c("estimate","se","conf.low","conf.high","Village", "group")




Sepiaranea<-as.data.frame(ranef(testSepiaea, probs = c(0.1,0.9)))
Sepiaranea$Village<-row.names(Sepiaranea)
Sepiaranea$group<- "Sepia"
colnames(Sepiaranea)<-c("estimate","se","conf.low","conf.high","Village", "group")



Cypraea.Conusranea<-as.data.frame(ranef(testCypraea.Conusea, probs = c(0.1,0.9)))
Cypraea.Conusranea$Village<-row.names(Cypraea.Conusranea)
Cypraea.Conusranea$group<- "Cypraea.Conus"
colnames(Cypraea.Conusranea)<-c("estimate","se","conf.low","conf.high","Village", "group")



Pyrasusranea<-as.data.frame(ranef(testPyrasusea, probs = c(0.1,0.9)))
Pyrasusranea$Village<-row.names(Pyrasusranea)
Pyrasusranea$group<- "Pyrasus"
colnames(Pyrasusranea)<-c("estimate","se","conf.low","conf.high","Village", "group")



Charoniaranea<-as.data.frame(ranef(testCharoniaea, probs = c(0.1,0.9)))
Charoniaranea$Village<-row.names(Charoniaranea)
Charoniaranea$group<- "Charonia"
colnames(Charoniaranea)<-c("estimate","se","conf.low","conf.high","Village", "group")


Murex.Fasciolariaranea<-as.data.frame(ranef(testMurex.Fasciolariaea, probs = c(0.1,0.9)))
Murex.Fasciolariaranea$Village<-row.names(Murex.Fasciolariaranea)
Murex.Fasciolariaranea$group<- "Murex.Fasciolaria"
colnames(Murex.Fasciolariaranea)<-c("estimate","se","conf.low","conf.high","Village", "group")



Lambisranea<-as.data.frame(ranef(testLambisea, probs = c(0.1,0.9)))
Lambisranea$Village<-row.names(Lambisranea)
Lambisranea$group<- "Lambis"
colnames(Lambisranea)<-c("estimate","se","conf.low","conf.high","Village", "group")


Anadararanea<-as.data.frame(ranef(testAnadaraea, probs = c(0.1,0.9)))
Anadararanea$Village<-row.names(Anadararanea)
Anadararanea$group<- "Anadara"
colnames(Anadararanea)<-c("estimate","se","conf.low","conf.high","Village", "group")



Pinctada.Isognomon.Atrina.Pinnaranea<-as.data.frame(ranef(testPinctada.Isognomon.Atrina.Pinnaea, probs = c(0.1,0.9)))
Pinctada.Isognomon.Atrina.Pinnaranea$Village<-row.names(Pinctada.Isognomon.Atrina.Pinnaranea)
Pinctada.Isognomon.Atrina.Pinnaranea$group<- "Pinctada.Isognomon.Atrina.Pinna"
colnames(Pinctada.Isognomon.Atrina.Pinnaranea)<-c("estimate","se","conf.low","conf.high","Village", "group")



Tridacnaranea<-as.data.frame(ranef(testTridacnaea, probs = c(0.1,0.9)))
Tridacnaranea$Village<-row.names(Tridacnaranea)
Tridacnaranea$group<- "Tridacna"
colnames(Tridacnaranea)<-c("estimate","se","conf.low","conf.high","Village", "group")


Tripneustesranea<-as.data.frame(ranef(testTripneustesea, probs = c(0.1,0.9)))
Tripneustesranea$Village<-row.names(Tripneustesranea)
Tripneustesranea$group<- "Tripneustes"
colnames(Tripneustesranea)<-c("estimate","se","conf.low","conf.high","Village", "group")


Aristeidaeranea<-as.data.frame(ranef(testAristeidaeea, probs = c(0.1,0.9)))
Aristeidaeranea$Village<-row.names(Aristeidaeranea)
Aristeidaeranea$group<- "Aristeidae"
colnames(Aristeidaeranea)<-c("estimate","se","conf.low","conf.high","Village", "group")



Scyllaranea<-as.data.frame(ranef(testScyllaea, probs = c(0.1,0.9)))
Scyllaranea$Village<-row.names(Scyllaranea)
Scyllaranea$group<- "Scylla"
colnames(Scyllaranea)<-c("estimate","se","conf.low","conf.high","Village", "group")



Palunirusranea<-as.data.frame(ranef(testPalunirusea, probs = c(0.1,0.9)))
Palunirusranea$Village<-row.names(Palunirusranea)
Palunirusranea$group<- "Palunirus"
colnames(Palunirusranea)<-c("estimate","se","conf.low","conf.high","Village", "group")




Randomeffectsizeea <- rbind(Mullidaeranea,
                            Pomacentridaeranea,
                            Acanthuridaeranea,
                            Scaridaeranea,
                            Scorpaenidaeranea,
                            Siganidaeranea,
                            Sphyraenidaeranea,
                            Labridaeranea,
                            Priacanthidaeranea,
                            Chaetodontidaeranea,
                            Fistularidaeranea,
                            Gerreidaeranea,
                            Haemulidaeranea,
                            Lethrinidaeranea,
                            Ostraciidaeranea,
                            Pomacanthidaeranea,
                            Rajidaeranea,
                            Rhinobatidaeranea,
                            Soleidaeranea,
                            Tetraodontidaeranea,
                            Carangidaeranea,
                            Balistidaeranea,
                            Zanclidaeranea,
                            Octopusranea,
                            Loligoranea,
                            Sepiaranea,
                            Cypraea.Conusranea,
                            Pyrasusranea,
                            Charoniaranea,
                            Murex.Fasciolariaranea,
                            Lambisranea,
                            Anadararanea,
                            Pinctada.Isognomon.Atrina.Pinnaranea,
                            Tridacnaranea,
                            Tripneustesranea,
                            Aristeidaeranea,
                            Scyllaranea,
                            Palunirusranea)


Randomeffectsizeea[Randomeffectsizeea == "Cypraea.Conus"] <- "Cypraeidae/Conus"
Randomeffectsizeea[Randomeffectsizeea == "Pinctada.Isognomon.Atrina.Pinna"] <- "Pteriidae/Isognomonidae"
Randomeffectsizeea[Randomeffectsizeea == "Murex.Fasciolaria"] <- "Muricidae"
Randomeffectsizeea[Randomeffectsizeea == "Octopus"] <- "Octopoda"

write.csv(Randomeffectsizeea, "eaEffectSizesRANDOM.csv", row.names=FALSE)

Randomeffectsizeea <- read.csv("eaEffectSizesRANDOM.csv")



jpeg("VillageEffectSizeEA", width =2000, height = 5000, res=300)
ggplot(Randomeffectsizeea,aes(x=estimate, fill=Village, y=fct_rev(fct_inorder(group)),xmin=conf.low,xmax=conf.high))+geom_errorbar(position=position_dodge(width=1))+geom_point(pch=21,size=3, position=position_dodge(width =1))+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  geom_smooth() +
  labs(title = str_wrap("Effect of Village on Perception of Change in Ecological Abundance of each Marine Resource Group", 45))+
  labs(x="Effect Size by Village", y="Marine Resource Group")+
  theme(plot.title = element_text(hjust = 0.5, margin=margin(b=20)))+
  theme(axis.title.x = element_text(margin=margin(t=20)), axis.title.y = element_text(margin=margin(r=20)))
dev.off()





##Randomeffects harvumption

Mullidaeranharv<-as.data.frame(ranef(testMullidaeharv, probs = c(0.1,0.9)))
Mullidaeranharv$Village<-row.names(Mullidaeranharv)
Mullidaeranharv$group<- "Mullidae"
colnames(Mullidaeranharv)<-c("estimate","se","conf.low","conf.high","Village", "group")

Pomacentridaeranharv<-as.data.frame(ranef(testPomacentridaeharv, probs = c(0.1,0.9)))
Pomacentridaeranharv$Village<-row.names(Pomacentridaeranharv)
Pomacentridaeranharv$group<- "Pomacentridae"
colnames(Pomacentridaeranharv)<-c("estimate","se","conf.low","conf.high","Village", "group")

Acanthuridaeranharv<-as.data.frame(ranef(testAcanthuridaeharv, probs = c(0.1,0.9)))
Acanthuridaeranharv$Village<-row.names(Acanthuridaeranharv)
Acanthuridaeranharv$group<- "Acanthuridae"
colnames(Acanthuridaeranharv)<-c("estimate","se","conf.low","conf.high","Village", "group")

Scaridaeranharv<-as.data.frame(ranef(testScaridaeharv, probs = c(0.1,0.9)))
Scaridaeranharv$Village<-row.names(Scaridaeranharv)
Scaridaeranharv$group<- "Scaridae"
colnames(Scaridaeranharv)<-c("estimate","se","conf.low","conf.high","Village", "group")


Scorpaenidaeranharv<-as.data.frame(ranef(testScorpaenidaeharv, probs = c(0.1,0.9)))
Scorpaenidaeranharv$Village<-row.names(Scorpaenidaeranharv)
Scorpaenidaeranharv$group<- "Scorpaenidae"
colnames(Scorpaenidaeranharv)<-c("estimate","se","conf.low","conf.high","Village", "group")

Siganidaeranharv<-as.data.frame(ranef(testSiganidaeharv, probs = c(0.1,0.9)))
Siganidaeranharv$Village<-row.names(Siganidaeranharv)
Siganidaeranharv$group<- "Siganidae"
colnames(Siganidaeranharv)<-c("estimate","se","conf.low","conf.high","Village", "group")


Sphyraenidaeranharv<-as.data.frame(ranef(testSphyraenidaeharv, probs = c(0.1,0.9)))
Sphyraenidaeranharv$Village<-row.names(Sphyraenidaeranharv)
Sphyraenidaeranharv$group<- "Sphyraenidae"
colnames(Sphyraenidaeranharv)<-c("estimate","se","conf.low","conf.high","Village", "group")

Labridaeranharv<-as.data.frame(ranef(testLabridaeharv, probs = c(0.1,0.9)))
Labridaeranharv$Village<-row.names(Labridaeranharv)
Labridaeranharv$group<- "Labridae"
colnames(Labridaeranharv)<-c("estimate","se","conf.low","conf.high","Village", "group")


Priacanthidaeranharv<-as.data.frame(ranef(testPriacanthidaeharv, probs = c(0.1,0.9)))
Priacanthidaeranharv$Village<-row.names(Priacanthidaeranharv)
Priacanthidaeranharv$group<- "Priacanthidae"
colnames(Priacanthidaeranharv)<-c("estimate","se","conf.low","conf.high","Village", "group")

Chaetodontidaeranharv<-as.data.frame(ranef(testChaetodontidaeharv, probs = c(0.1,0.9)))
Chaetodontidaeranharv$Village<-row.names(Chaetodontidaeranharv)
Chaetodontidaeranharv$group<- "Chaetodontidae"
colnames(Chaetodontidaeranharv)<-c("estimate","se","conf.low","conf.high","Village", "group")

Fistularidaeranharv<-as.data.frame(ranef(testFistularidaeharv, probs = c(0.1,0.9)))
Fistularidaeranharv$Village<-row.names(Fistularidaeranharv)
Fistularidaeranharv$group<- "Fistularidae"
colnames(Fistularidaeranharv)<-c("estimate","se","conf.low","conf.high","Village", "group")

Gerreidaeranharv<-as.data.frame(ranef(testGerreidaeharv, probs = c(0.1,0.9)))
Gerreidaeranharv$Village<-row.names(Gerreidaeranharv)
Gerreidaeranharv$group<- "Gerreidae"
colnames(Gerreidaeranharv)<-c("estimate","se","conf.low","conf.high","Village", "group")



Haemulidaeranharv<-as.data.frame(ranef(testHaemulidaeharv, probs = c(0.1,0.9)))
Haemulidaeranharv$Village<-row.names(Haemulidaeranharv)
Haemulidaeranharv$group<- "Haemulidae"
colnames(Haemulidaeranharv)<-c("estimate","se","conf.low","conf.high","Village", "group")


Lethrinidaeranharv<-as.data.frame(ranef(testLethrinidaeharv, probs = c(0.1,0.9)))
Lethrinidaeranharv$Village<-row.names(Lethrinidaeranharv)
Lethrinidaeranharv$group<- "Lethrinidae"
colnames(Lethrinidaeranharv)<-c("estimate","se","conf.low","conf.high","Village", "group")



Ostraciidaeranharv<-as.data.frame(ranef(testOstraciidaeharv, probs = c(0.1,0.9)))
Ostraciidaeranharv$Village<-row.names(Ostraciidaeranharv)
Ostraciidaeranharv$group<- "Ostraciidae"
colnames(Ostraciidaeranharv)<-c("estimate","se","conf.low","conf.high","Village", "group")



Pomacanthidaeranharv<-as.data.frame(ranef(testPomacanthidaeharv, probs = c(0.1,0.9)))
Pomacanthidaeranharv$Village<-row.names(Pomacanthidaeranharv)
Pomacanthidaeranharv$group<- "Pomacanthidae"
colnames(Pomacanthidaeranharv)<-c("estimate","se","conf.low","conf.high","Village", "group")


Rajidaeranharv<-as.data.frame(ranef(testRajidaeharv, probs = c(0.1,0.9)))
Rajidaeranharv$Village<-row.names(Rajidaeranharv)
Rajidaeranharv$group<- "Rajidae"
colnames(Rajidaeranharv)<-c("estimate","se","conf.low","conf.high","Village", "group")



Rhinobatidaeranharv<-as.data.frame(ranef(testRhinobatidaeharv, probs = c(0.1,0.9)))
Rhinobatidaeranharv$Village<-row.names(Rhinobatidaeranharv)
Rhinobatidaeranharv$group<- "Rhinobatidae"
colnames(Rhinobatidaeranharv)<-c("estimate","se","conf.low","conf.high","Village", "group")


Soleidaeranharv<-as.data.frame(ranef(testSoleidaeharv, probs = c(0.1,0.9)))
Soleidaeranharv$Village<-row.names(Soleidaeranharv)
Soleidaeranharv$group<- "Soleidae"
colnames(Soleidaeranharv)<-c("estimate","se","conf.low","conf.high","Village", "group")


Tetraodontidaeranharv<-as.data.frame(ranef(testTetraodontidaeharv, probs = c(0.1,0.9)))
Tetraodontidaeranharv$Village<-row.names(Tetraodontidaeranharv)
Tetraodontidaeranharv$group<- "Tetraodontidae"
colnames(Tetraodontidaeranharv)<-c("estimate","se","conf.low","conf.high","Village", "group")



Carangidaeranharv<-as.data.frame(ranef(testCarangidaeharv, probs = c(0.1,0.9)))
Carangidaeranharv$Village<-row.names(Carangidaeranharv)
Carangidaeranharv$group<- "Carangidae"
colnames(Carangidaeranharv)<-c("estimate","se","conf.low","conf.high","Village", "group")


Balistidaeranharv<-as.data.frame(ranef(testBalistidaeharv, probs = c(0.1,0.9)))
Balistidaeranharv$Village<-row.names(Balistidaeranharv)
Balistidaeranharv$group<- "Balistidae"
colnames(Balistidaeranharv)<-c("estimate","se","conf.low","conf.high","Village", "group")



Zanclidaeranharv<-as.data.frame(ranef(testZanclidaeharv, probs = c(0.1,0.9)))
Zanclidaeranharv$Village<-row.names(Zanclidaeranharv)
Zanclidaeranharv$group<- "Zanclidae"
colnames(Zanclidaeranharv)<-c("estimate","se","conf.low","conf.high","Village", "group")




Octopusranharv<-as.data.frame(ranef(testOctopusharv, probs = c(0.1,0.9)))
Octopusranharv$Village<-row.names(Octopusranharv)
Octopusranharv$group<- "Octopus"
colnames(Octopusranharv)<-c("estimate","se","conf.low","conf.high","Village", "group")




Loligoranharv<-as.data.frame(ranef(testLoligoharv, probs = c(0.1,0.9)))
Loligoranharv$Village<-row.names(Loligoranharv)
Loligoranharv$group<- "Loligo"
colnames(Loligoranharv)<-c("estimate","se","conf.low","conf.high","Village", "group")




Sepiaranharv<-as.data.frame(ranef(testSepiaharv, probs = c(0.1,0.9)))
Sepiaranharv$Village<-row.names(Sepiaranharv)
Sepiaranharv$group<- "Sepia"
colnames(Sepiaranharv)<-c("estimate","se","conf.low","conf.high","Village", "group")



Cypraea.Conusranharv<-as.data.frame(ranef(testCypraea.Conusharv, probs = c(0.1,0.9)))
Cypraea.Conusranharv$Village<-row.names(Cypraea.Conusranharv)
Cypraea.Conusranharv$group<- "Cypraea.Conus"
colnames(Cypraea.Conusranharv)<-c("estimate","se","conf.low","conf.high","Village", "group")



Pyrasusranharv<-as.data.frame(ranef(testPyrasusharv, probs = c(0.1,0.9)))
Pyrasusranharv$Village<-row.names(Pyrasusranharv)
Pyrasusranharv$group<- "Pyrasus"
colnames(Pyrasusranharv)<-c("estimate","se","conf.low","conf.high","Village", "group")



Charoniaranharv<-as.data.frame(ranef(testCharoniaharv, probs = c(0.1,0.9)))
Charoniaranharv$Village<-row.names(Charoniaranharv)
Charoniaranharv$group<- "Charonia"
colnames(Charoniaranharv)<-c("estimate","se","conf.low","conf.high","Village", "group")


Murex.Fasciolariaranharv<-as.data.frame(ranef(testMurex.Fasciolariaharv, probs = c(0.1,0.9)))
Murex.Fasciolariaranharv$Village<-row.names(Murex.Fasciolariaranharv)
Murex.Fasciolariaranharv$group<- "Murex.Fasciolaria"
colnames(Murex.Fasciolariaranharv)<-c("estimate","se","conf.low","conf.high","Village", "group")



Lambisranharv<-as.data.frame(ranef(testLambisharv, probs = c(0.1,0.9)))
Lambisranharv$Village<-row.names(Lambisranharv)
Lambisranharv$group<- "Lambis"
colnames(Lambisranharv)<-c("estimate","se","conf.low","conf.high","Village", "group")


Anadararanharv<-as.data.frame(ranef(testAnadaraharv, probs = c(0.1,0.9)))
Anadararanharv$Village<-row.names(Anadararanharv)
Anadararanharv$group<- "Anadara"
colnames(Anadararanharv)<-c("estimate","se","conf.low","conf.high","Village", "group")



Pinctada.Isognomon.Atrina.Pinnaranharv<-as.data.frame(ranef(testPinctada.Isognomon.Atrina.Pinnaharv, probs = c(0.1,0.9)))
Pinctada.Isognomon.Atrina.Pinnaranharv$Village<-row.names(Pinctada.Isognomon.Atrina.Pinnaranharv)
Pinctada.Isognomon.Atrina.Pinnaranharv$group<- "Pinctada.Isognomon.Atrina.Pinna"
colnames(Pinctada.Isognomon.Atrina.Pinnaranharv)<-c("estimate","se","conf.low","conf.high","Village", "group")



Tridacnaranharv<-as.data.frame(ranef(testTridacnaharv, probs = c(0.1,0.9)))
Tridacnaranharv$Village<-row.names(Tridacnaranharv)
Tridacnaranharv$group<- "Tridacna"
colnames(Tridacnaranharv)<-c("estimate","se","conf.low","conf.high","Village", "group")


Tripneustesranharv<-as.data.frame(ranef(testTripneustesharv, probs = c(0.1,0.9)))
Tripneustesranharv$Village<-row.names(Tripneustesranharv)
Tripneustesranharv$group<- "Tripneustes"
colnames(Tripneustesranharv)<-c("estimate","se","conf.low","conf.high","Village", "group")


Aristeidaeranharv<-as.data.frame(ranef(testAristeidaeharv, probs = c(0.1,0.9)))
Aristeidaeranharv$Village<-row.names(Aristeidaeranharv)
Aristeidaeranharv$group<- "Aristeidae"
colnames(Aristeidaeranharv)<-c("estimate","se","conf.low","conf.high","Village", "group")



Scyllaranharv<-as.data.frame(ranef(testScyllaharv, probs = c(0.1,0.9)))
Scyllaranharv$Village<-row.names(Scyllaranharv)
Scyllaranharv$group<- "Scylla"
colnames(Scyllaranharv)<-c("estimate","se","conf.low","conf.high","Village", "group")



Palunirusranharv<-as.data.frame(ranef(testPalunirusharv, probs = c(0.1,0.9)))
Palunirusranharv$Village<-row.names(Palunirusranharv)
Palunirusranharv$group<- "Palunirus"
colnames(Palunirusranharv)<-c("estimate","se","conf.low","conf.high","Village", "group")




Randomeffectsizeharv <- rbind(Mullidaeranharv,
                              Pomacentridaeranharv,
                              Acanthuridaeranharv,
                              Scaridaeranharv,
                              Scorpaenidaeranharv,
                              Siganidaeranharv,
                              Sphyraenidaeranharv,
                              Labridaeranharv,
                              Priacanthidaeranharv,
                              Chaetodontidaeranharv,
                              Fistularidaeranharv,
                              Gerreidaeranharv,
                              Haemulidaeranharv,
                              Lethrinidaeranharv,
                              Ostraciidaeranharv,
                              Pomacanthidaeranharv,
                              Rajidaeranharv,
                              Rhinobatidaeranharv,
                              Soleidaeranharv,
                              Tetraodontidaeranharv,
                              Carangidaeranharv,
                              Balistidaeranharv,
                              Zanclidaeranharv,
                              Octopusranharv,
                              Loligoranharv,
                              Sepiaranharv,
                              Cypraea.Conusranharv,
                              Pyrasusranharv,
                              Charoniaranharv,
                              Murex.Fasciolariaranharv,
                              Lambisranharv,
                              Anadararanharv,
                              Pinctada.Isognomon.Atrina.Pinnaranharv,
                              Tridacnaranharv,
                              Tripneustesranharv,
                              Aristeidaeranharv,
                              Scyllaranharv,
                              Palunirusranharv)


Randomeffectsizeharv[Randomeffectsizeharv == "Cypraea.Conus"] <- "Cypraeidae/Conus"
Randomeffectsizeharv[Randomeffectsizeharv == "Pinctada.Isognomon.Atrina.Pinna"] <- "Pteriidae/Isognomonidae"
Randomeffectsizeharv[Randomeffectsizeharv == "Murex.Fasciolaria"] <- "Muricidae"
Randomeffectsizeharv[Randomeffectsizeharv == "Octopus"] <- "Octopoda"

write.csv(Randomeffectsizeharv, "harvEffectSizesRANDOM.csv", row.names=FALSE)

Randomeffectsizeharv <- read.csv("harvEffectSizesRANDOM.csv")



jpeg("VillageEffectSizeHarvest", width =2000, height = 5000, res=300)
ggplot(Randomeffectsizeharv,aes(x=estimate, fill=Village, y=fct_rev(fct_inorder(group)),xmin=conf.low,xmax=conf.high))+geom_errorbar(position=position_dodge(width=1))+geom_point(pch=21,size=3, position=position_dodge(width =1))+geom_vline(xintercept = 0,lty=2)+theme_classic()+
  geom_smooth() +
  labs(title = str_wrap("Effect of Village on Perception of Change in Harvest of each Marine Resource Group", 45))+
  labs(x="Effect Size by Village", y="Marine Resource Group")+
  theme(plot.title = element_text(hjust = 0.5, margin=margin(b=20)))+
  theme(axis.title.x = element_text(margin=margin(t=20)), axis.title.y = element_text(margin=margin(r=20)))
dev.off()


