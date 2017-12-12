##install packages ggplot2,car,GGally, knitr
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
prostate <- prostate[c(-37),]
attach(prostate)
gleason <- factor(gleason, ordered=FALSE)
svi <- factor(svi, ordered=FALSE)

##Fit full model
full_model <- lm(lpsa ~ lcavol + lweight + age + lbph + svi + lcp + gleason + pgg45)

##Diagnostics
studentized<-rstudent(full_model)
par(mfrow = c(3, 3))
plot(x=lcavol, y=studentized,xlab="Log cancer volumne",ylab="Standardized residuals",cex.lab=1.4)
plot(x=lweight, y=studentized,xlab="Log prostate weight",ylab="Standardized residuals",cex.lab=1.4)
plot(x=age, y=studentized,xlab="Age of patient",ylab="Standardized residuals",cex.lab=1.4)
plot(x=lbph, y=studentized,xlab="Log benign prostatic hyperplasia",ylab="Standardized residuals",cex.lab=1.4)
plot(x=svi, y=studentized,xlab="Seminal vesicle invasion",ylab="Standardized residuals",cex.lab=1.4)
plot(x=lcp, y=studentized,xlab="Log capsular penetration",ylab="Standardized residuals",cex.lab=1.4)
plot(x=gleason, y=studentized,xlab="Gleason score",ylab="Standardized residuals",cex.lab=1.4)
plot(x=pgg45, y=studentized,xlab="Pct of Gleason score 4 or 5",ylab="Standardized residuals",cex.lab=1.4)

par(mfrow = c(2, 2), mar = c(5, 5, 1.5, 0.5))
bl <- scales::alpha("black", 0.5) #semi-transparent black
n <- nrow(prostate)

#Student Q-Q plot
qqPlot(full_model, simulate = 1999, envelope = TRUE,
       ylab = "Standardized residuals", 
       xlab = "Theoretical student quantiles", 
       pch = 20, col = bl, main="(a) QQPlot")

#Residuals vs fitted values
residualPlot(full_model, type = "rstudent", quadratic = FALSE, 
             pch = 20, ylab = "Standardized residuals", main="(b) Standardized residuals vs fitted values")

#Cook distance hat values
plot(cooks.distance(full_model), col = bl, pch = 20, ylab = "Cook's distances", main="(c) Cook's distance")
abline(h = 8/(n-2*length(coef(full_model))), col = 2)
influencePlot(full_model, ylab = "Standardized residuals", main="(d) Hat-Values")

##Look at influential observations
prostate[38,] ##outlier value
prostate[40,] ##leverage point
prostate[46,] ##leverage point

##Colinearity
vif(full_model)
M <- model.matrix(full_model)
kappa(M, exact = TRUE, norm='2')

##Automatic Model Building

##perform backward selection
##drop multiple term each step to preserve power
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

##add_interactions with lcavol to try to answer question 1
interaction1_model<- lm(lpsa ~ lcavol + lweight + svi + lcavol:lweight)
interaction2_model<- lm(lpsa ~ lcavol + lweight + svi + lcavol:svi)

##perform model selection on models built

AIC(full_model)
AIC(backward_selection_model)
AIC(interaction1_model)
AIC(interaction2_model)

BIC(full_model)
BIC(backward_selection_model)
BIC(interaction1_model)
BIC(interaction2_model)

summary(full_model)$adj.r.squared
summary(backward_selection_model)$adj.r.squared
summary(interaction1_model)$adj.r.squared
summary(interaction2_model)$adj.r.squared

vif(full_model)
vif(backward_selection_model)
vif(interaction1_model)

##select backward_selection_model
final_model <- backward_selection_model

summary_final <- summary(final_model)
knitr::kable(coef(summary_final), digits = 3)
vif(final_model)

par(mfrow = c(1, 1))
car::avPlot(lm(lcavol ~ lpsa + lweight + age + lbph + svi + lcp + gleason + pgg45), variable = "lpsa", ellipse = TRUE, col.lines = NULL)
