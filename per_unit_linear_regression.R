# Library
library("raveio")
library("ggplot2")
library("lme4")
library("MASS")


############ Linear regression ############
############ Variance Explained ############

for (i in 1:15) {

unit = i

############## Load & format data ###############
print(paste("Iteration ",as.character(i)))

220407_14650_datasetComplete1_allFlights_9

#dat <- read_mat(paste("GLM/datasets/training/220407_14650_dataset100msBatCentric_allFlights_",as.character(unit),".mat",sep=""))
#colnames(dat$table_) <- c("n","Bx","By","Bat_to_Mx","Bat_to_My","Bat_to_Kx","Bat_to_Ky","Bat_to_BOx","Bat_to_BOy","C","Pre","Dur","Post")
#dat <- read_mat(paste("GLM/datasets/training/220407_14650_dataset100ms_allFlights_",as.character(unit),".mat",sep=""))
#colnames(dat$table_) <- c("n","Bx","By","BOx","BOy","Kx","Ky","Mx","My","C","P","D","PO")
#dat <- read_mat(paste("GLM/datasets/training/220407_14650_datasetDistance_allFlights_",as.character(unit),".mat",sep=""))
#colnames(dat$table_) <- c("n","B2M","B2K","C","PDP")
#dat <- read_mat(paste("GLM/datasets/training/220407_14650_datasetDistanceOnlyPre_allFlights_",as.character(unit),".mat",sep=""))
#colnames(dat$table_) <- c("n","B2M","B2K","C")
dat <- read_mat(paste("GLM/datasets/training/220407_14650_datasetComplete1_allFlights_",as.character(unit),".mat",sep=""))
colnames(dat$table_) <- c("n","B2M","B2K","B2OB","D2Goal","D2Home","WhoL","WhoT","Goal_Factor","C","DFPF","PDP")

dat_frame <- data.frame(dat$table_)
rm(dat)
gc()

# Change the pre, during and post to factors
dat_frame$P <- as.factor(dat_frame$P)
dat_frame$D <- as.factor(dat_frame$D)
dat_frame$P <- as.factor(dat_frame$P)

dat_frame$Goal_Factor = as.factor(dat_frame$Goal_Factor)
dat_frame$PDP = as.factor(dat_frame$PDP)




########## Fit the GLM model ##########
m_glm <- glm(n ~ 1 + B2M + B2K + B2OB + D2Goal + D2Home + WhoL + WhoT + C + DFPF,family=gaussian(link='identity'),data=dat_frame)
m_glm <- glm(n ~ 1 + B2M + B2K + B2OB + C + DFPF,family=gaussian(link='identity'),data=dat_frame)
m_glm_summary <- summary(m_glm)

# Plot the residuals to check model fit
res_glm <- resid(m_glm)
plot(fitted(m_glm), res_glm)
abline(0,0)

# Make qq plot to check residuals follow gaussian distribution
qqnorm(res_glm)
qqline(res_glm)

# Check for under or over dispersion
E2 <- resid(m_glm,type='pearson')
N <-nrow(dat_frame)
pl <- length(coef(m_glm))
mglm_disp <- sum(E2^2)/(N-pl)

# Not really sure how to interpret variance explained here.



########## Fit the LM model ##########
m_lm <- lm(n ~ 1 + B2M + B2K + C + B2M*C + B2K*C, data=dat_frame)
m_lm <- lm(n ~ 1 + B2M + B2K + B2OB + D2Goal + D2Home + WhoL + WhoT + Goal_Factor + C + DFPF + PDP, data=dat_frame)

colnames(dat$table_) <- c("n","B2M","B2K","B2OB","D2Goal","D2Home","WhoL","WhoT","Goal_Factor","C","DFPF","PDP")


anova_lm <- anova(m_lm)
lm_ss <- anova_lm$"Sum Sq"
print(cbind(anova_lm,PctExp=lm_ss/sum(lm_ss)*100))
variance_explained <- cbind(anova_lm,PctExp=lm_ss/sum(lm_ss)*100)

# Plot the residuals to check model fit
res_lm <- resid(m_lm)
plot(fitted(m_lm), res_lm)
abline(0,0)

# Make qq plot to check residuals follow gaussian distribution
qqnorm(res_lm)
qqline(res_lm)



########## Fit Ordinal Logistic Regression model ##########
library("MASS")
dat_frame$n = factor(dat_frame$n)

o_lr = polr(n ~ 1 + B2M + B2K + C + PDP, data = dat_frame, Hess = TRUE)
summary(o_lr)

# Perform chi squared test
chi_olr <- 1-pchisq(deviance(o_lr),df.residual(o_lr))





########## Fit GLM Gamma model ##########
library("ggplot2")

glm_gamma <- glm(n ~ B2M + B2K + C + PDP, family=Gamma,data=dat_frame)
res_gamma <- resid(glm_gamma)
plot(fitted(glm_gamma), res_gamma)
abline(0,0)

qqnorm(res_gamma)
qqline(res_gamma)






########## Fit GLM Poisson model ##########

glm_pois <- glm(n ~ B2M + B2K + C + B2K*C + B2M*C, family=poisson,data=dat_frame)
res_gamma <- resid(glm_gamma)
plot(fitted(glm_gamma), res_gamma)
abline(0,0)

qqnorm(res_gamma)
qqline(res_gamma)

# FOr the poisson and zero inflated stuff change 1 to zeros in n
dat_frame$n = dat_frame$n-1

# Pois
pois <- glm(n ~ B2K + B2M + C + B2K*C + B2M*C, family = 'poisson',data=dat_frame)
# Check for under or over dispersion
E2 <- resid(pois,type='pearson')
N <-nrow(dat_frame)
pl <- length(coef(pois))
pois_disp <- sum(E2^2)/(N-pl)

# Negative binomial
negbin <- glm.nb(formula = n ~ B2K + B2M + C + B2K*C + B2M*C, data=dat_frame)
# Check for under or over dispersion
E2 <- resid(negbin,type='pearson')
N <-nrow(dat_frame)
pl <- length(coef(negbin))
negbin_disp <- sum(E2^2)/(N-pl)

# Zip Poisson
zippois <- zeroinfl(n ~ B2K + B2M + B2OB + C + PDP,dist = 'poisson',data=dat_frame)
summary(zippois)
E2 <- resid(zippois,type='pearson')
N <-nrow(dat_frame)
pl <- length(coef(zippois))
zippois_disp <- sum(E2^2)/(N-pl)

# Zip neg bin
ZNB <- zeroinfl(n ~ B2K + B2M + B2OB + B2K* C + DFPF + PDP, dist = 'negbin',data=dat_frame)
summary(ZNB)
E2 <- resid(ZNB,type='pearson')
N <-nrow(dat_frame)
pl <- length(coef(ZNB))
ZNB_disp <- sum(E2^2)/(N-pl)








########## Multiple Mixed Regression ##########
m_mmr <- lmer(n ~ 1 + (1|C) + (1 | Bx) + (1 | BOx) + (1|Kx) + (1|Mx) + (1|P) + (1 | D), dat_frame)

m_mmr_summary <- summary(m_mmr)
write.table(m_mmr_summary$coefficients, "GLM/datasets/results/220407/m_mmr",sep="\t")
anova_results <- anova(m_mmr, test="Chisq")
write.table(anova_results$Pr, "GLM/datasets/results/220407/m_mmr_anova",sep="\t")





########### Write results ###########
write.table(variance_explained, paste("GLM/datasets/results/220407/m_lm_variance_explained_",as.character(unit),sep=""),sep="\t")

rm(dat_frame)
rm(m_lr)
rm(m_lr_summary)
rm(anova_results)
gc()

}
