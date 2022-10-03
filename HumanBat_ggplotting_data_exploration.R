######### Data Exploration ##########

######### Load Libraries #############
library("lme4")
library("raveio")
library("ggplot2")

# Load in dataframe (Only Pre, -1s and takeoff)
dat <- read_mat(paste("GLM/datasets/training/220407_14650_datasetDistanceOnlyPre_allFlights_",as.character(unit),".mat",sep=""))
colnames(dat$table_) <- c("n","B2M","B2K","C")

# Load in dataframe (Only During)
#dat <- read_mat(paste("GLM/datasets/training/220407_14650_datasetDistanceOnlyDur_allFlights_",as.character(unit),".mat",sep=""))
#colnames(dat$table_) <- c("n","B2M","B2K","C")

# Make into a dataframe
dat_frame <- data.frame(dat$table_)
rm(dat)
gc()

######### Plot variable distributions #############
ggplot(data=dat_frame) + geom_bar(mapping = aes(x = D2Home))
ggplot(data = dat_frame) +geom_histogram(mapping = aes(x = D2Home), binwidth = 50)

with(dat_frame, plot(n ~ B2K))

plot(dat_frame$C,resid(m_lm),ylab="Residuals" )

# Plot a Poisson
Pois <- rpoisson(1000,shape=1,scale=1)
