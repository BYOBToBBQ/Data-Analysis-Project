##install packages ggplot2,car,GGally, knitr

##stargazer for table
library(ggplot2)
library(car)
library(knitr)
library(GGally)
library(leaps)


prostate <- read.table("http://sma.epfl.ch/~lbelzile/math341/prostate.dat", header = TRUE)
##Exploratory data analysis
str(prostate)
apply(prostate, 2, range)
print(GGally::ggpairs(prostate[,-c(5,7)]), progress = FALSE)
##drop only one with gleason 8
prostate <- prostate[c(-37),]
attach(prostate)
##lpsa has linear relationship with lcavol, (less)lweight, lcp, rests seems uncorrelated
##collinearity between lcavol-lcp 
##no need for transformations
##default value observed before taking log need to clean
gleason <- factor(prostate$gleason, ordered=FALSE)
svi <- factor(svi, ordered=FALSE)

##Fit full model
full_model <- lm(lpsa ~ lcavol + lweight + age + lbph + svi + lcp + gleason + pgg45)
##do diagnostics
studentized<-rstudent(full_model)

##, mar = c(5, 5, 1.5, 0.5)
par(mfrow = c(3, 3))
plot(x=lcavol, y=studentized,xlab="Log cancer volumne",ylab="Residuals",cex.lab=1.4)
plot(x=lweight, y=studentized,xlab="Log prostate weight",ylab="Residuals",cex.lab=1.4)
plot(x=age, y=studentized,xlab="Age of patient",ylab="Residuals",cex.lab=1.4)
plot(x=lbph, y=studentized,xlab="Log benign prostatic hyperplasia",ylab="Residuals",cex.lab=1.4)
plot(x=svi, y=studentized,xlab="Seminal vesicle invasion",ylab="Residuals",cex.lab=1.4)
plot(x=lcp, y=studentized,xlab="Log capsular penetration",ylab="Residuals",cex.lab=1.4)
plot(x=gleason, y=studentized,xlab="Gleason score",ylab="Residuals",cex.lab=1.4)
plot(x=pgg45, y=studentized,xlab="Percentages of Gleason score 4 or 5",ylab="Residuals",cex.lab=1.4)


par(mfrow = c(2, 2), mar = c(5, 5, 1.5, 0.5))
bl <- scales::alpha("black", 0.5) #semi-transparent black
n <- nrow(prostate)
#Student Q-Q plot
qqPlot(full_model, simulate = 1999, envelope = TRUE,
       ylab = "Externally studentized residuals", 
       xlab = "Theoretical student quantiles", 
       pch = 20, col = bl)
#Residuals vs fitted values
residualPlot(full_model, type = "rstudent", quadratic = FALSE, 
             pch = 20, ylab = "Externally studentized residuals")
#Cook distance hat values
plot(cooks.distance(full_model), col = bl, pch = 20, ylab = "Cook's distances")
abline(h = 8/(n-2*length(coef(full_model))), col = 2)
influencePlot(full_model)
##colinearity
vif(full_model)
M <- model.matrix(full_model)
kappa(M, exact = TRUE, norm='2')

##gleason progress pas de manière linéaire
summary_fm <- summary(full_model)
knitr::kable(coef(summary_fm), digits = 3)
##have intercept table, interpret t values
##strong evidence against the hypothesis that lcavol, lweight, svi souldn't be included
##pgg45 close to 0

##Automatic Model Building

##perform backward selection
##drop multiple term each step to avoid repeated testing
drop1(lm(lpsa ~ lcavol + lweight + age + lbph + svi + lcp + gleason + pgg45), test = "F")
drop1(lm(lpsa ~ lcavol + lweight + age + lbph + svi + lcp), test = "F")
drop1(lm(lpsa ~ lcavol + lweight + lbph + svi), test = "F")
drop1(lm(lpsa ~ lcavol + lweight + svi), test = "F")

backward_selection_model<- lm(lpsa ~ lcavol + lweight + svi)

##fit all possible models
leaps=regsubsets(lpsa ~ lcavol + lweight + age + lbph + svi + lcp + gleason + pgg45,data=prostate, nbest=8)
par(mfrow=c(1,1))
plot(leaps, scale = "bic", main = "BIC")
##end up with same model

##add_interactions with lcavol to answer question 1
interaction1_model<- lm(lpsa ~ lcavol + lweight + svi + lcavol:lweight)
interaction2_model<- lm(lpsa ~ lcavol + lweight + svi + lcavol:svi)
##interaction3_model<- lm(lpsa ~ lcavol + lweight + svi + lcavol:svi + lcavol:lweight)

##perform model selection on models built

AIC(full_model)
AIC(backward_selection_model)
AIC(interaction1_model)
AIC(interaction2_model)
##AIC(interaction3_model)

BIC(full_model)
BIC(backward_selection_model)
BIC(interaction1_model)
BIC(interaction2_model)
##BIC(interaction3_model)

summary(full_model)$adj.r.squared
summary(backward_selection_model)$adj.r.squared
summary(interaction1_model)$adj.r.squared
summary(interaction2_model)$adj.r.squared
##summary(interaction3_model)$adj.r.squared

vif(full_model)
vif(backward_selection_model)
vif(interaction1_model)
vif(interaction2_model)
##vif(interaction3_model)

##select backward_selection_model