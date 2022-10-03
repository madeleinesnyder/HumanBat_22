# Library
library("lme4")
library("raveio")

############### Repeated Measures Logistic Regression ################
####### Run 10x each with different data subsets ################
for (unit in 1:15) {
for (i in 1:10) {

############## Load & format data #################
print(paste("Iteration ",as.character(i)))

dat <- read_mat(paste("GLM/datasets/training/220407_14650_dataset3_",as.character(unit),".mat",sep=""))
colnames(dat$table1) <- c("n","B","BO","K","M","C","P","D","PO")
dat_frame <- data.frame(dat$table1)
rm(dat)
gc()

# Randomly choose subset of data
train <- dat_frame[sample(nrow(dat_frame), round(nrow(dat_frame)/4)), ]
rm(dat_frame)
gc()

########## Fit the Logistic GLM model ##########
m_logitc <- glm(n ~ 1 + B + BO + M + K + C + P + D,family=binomial(link='logit'),data=train)


########## Look at goodness of fit ##########
deviance_diff <- m_logitc$null.deviance-m_logitc$deviance
# IF this value is small, shitty model


########## Look at Log Odds Ratio for each coefficient ##########
# predict gives the predicted value in terms of logits
plot.dat <- data.frame(prob = train$n,
                       B = train$B,
                       fit = predict(m_logitc, train))
#convert those logit values to probabilities
plot.dat$fit_prob <- exp(plot.dat$fit)/(1+exp(plot.dat$fit))



exp(coef(m_logitc))

m_logitc_summary <- summary(m_logitc)
#write.table(m_logitc_summary$coefficients, paste("GLM/datasets/results/220407/unit_",as.character(unit),"_m_logitc_",as.character(i)),sep="\t")
#anova_results <- anova(m_logitc, test="Chisq")
#write.table(anova_results$Pr, paste("GLM/datasets/results/220407/unit_",as.character(unit),"_m_logitc_anova_",as.character(i)),sep="\t")

rm(train)
rm(m_logitc)
rm(m_logitc_summary)
rm(anova_results)
gc()
}
}
