rm(list=ls()) # clear workspace
setwd("C:/Users/.../")

library(haven)
library(tidyverse)
library(stargazer)
library(readr)
library(ggplot2)
library(sandwich)
library(lmtest)
library(AER)#overdispersion test
library(margins)
library(MASS) #negative binomial
library(pscl) #Hurdle
library(msme)#Pearson chi-square goodness-of-fit test
library(censReg)
library(VGAM)
library(estimatr)


# Import data
tradedata <- read_csv("Assignment/Trade Data.csv",show_col_types = FALSE)
attach(tradedata)
view(tradedata)

summary(tradedata)

##############
#   Zeros    #
##############

trade_zero<-mean(!trade)*100 #% of trade made up of zeros
cat("Zeros make up ",round(trade_zero),"% of the observations for the trade flow variable 
where its minimum value is",min(trade), ", its maximum value is", max(trade), 
    ", and its mean is", mean(trade))

# The distribution has a long right tail.
# Trade has a the max value is 101000000 yet more than half of the values are under 172129.5.
# This makes sense given that the proportion of zeros is 47.6%.


###############################################################################
#                       OLS Linear Regression: trade                          #
###############################################################################
trade_OLS<-lm(trade ~ logGDP_im + logGDP_ex + logGDPC_im + logGDPC_ex +
                logdist + comborder + comlang +
                colonial + comFTA + landlocked_ex + landlocked_im +  
                remoteness_im + remoteness_ex + openness
              ,data = tradedata)
summary(trade_OLS)
stargazer(trade_OLS,type="text")

#ROBUST standard errors
coeftest(trade_OLS, vcov = vcovHC(trade_OLS, type = "HC1"))

##############
#  Problem   #
##############

# The regression has negative values
tradedata$yhat_trade_OLS <- predict(trade_OLS, tradedata, type = "response")
summary(tradedata[, c("trade","yhat_trade_OLS")])

###############################################################################
#                       MLE Poisson Regression: trade                         #
###############################################################################

poisson_model <- glm(trade ~ logGDP_im + logGDP_ex + logGDPC_im + logGDPC_ex +
                       logdist + comborder + comlang +
                       colonial + comFTA + landlocked_ex + landlocked_im +  
                       remoteness_im + remoteness_ex + openness,
                     data = tradedata, family = poisson(link = "log"))
summary(poisson_model)
stargazer(poisson_model, type="text")

# The residual for deviation is 1.7404e+09  on 18345  degrees of freedom (), given
# that the rule of thumb is ratio of deviance to degrees of freedom should be 1
# indicates overdispersion.
# [(1,740,390,150/18345)= 94870.5369]


#ROBUST standard errors
robust_poisson_model<-coeftest(poisson_model, vcov. = vcovHC(poisson_model, type = "HC1"))
robust_poisson_model

stargazer(poisson_model,robust_poisson_model,type="text")

#####################
#  OVERDISPERSION   #
#####################

dispersiontest(poisson_model)
# the p-value is significant therefore we can conclude: overdispersion. 

#####################
#  Goodness of Fit  #
#####################

#The squared correlation coefficient is another way to assess goodness of fit.
# R-squared for Poisson

# Predict the number of counts (exp(x'beta))
tradedata$yPhat <- predict(poisson_model, type = "response")

summary(tradedata[c("yPhat", "trade")])

correlation <- cor(tradedata$trade, tradedata$yPhat)

cat("Squared correlation between y and yhat = ", correlation^2, "\n")

# The squared correlation coefficient is 0.86, it means that 86% of the
# variation in trade can be explained by variation in the independent variables

###########################################
# Pearson chi-square goodness-of-fit test #
###########################################

P__disp(poisson_model)

# pearson.chi2   dispersion 
# 3141793875.4     171261.6 

# Pearson chi-square goodness-of-fit test measures of how well the model fits observations.

# With a large pearson.chi2 statistic of 3141793875.4 it is clear
# there is a big the difference between the observations and the expectations.

################################
# ZERO PREDICTIONS VS OBSERVED #
################################

# yPhat is the expected mean count for each observation
# rounded sum of predict probability of zero count using expected mean counts
zero_predict_P <- round(sum(dpois(x = 0, lambda = tradedata$yPhat)))

# observed number of 0's
zero_actual_P <- sum(tradedata$trade < 1)

cat("The model predicts ",zero_predict_P, " zero counts but we observe", zero_actual_P, 
    "zero counts \nThis indicates some severe underfitting of zero counts")

###############################################################################
#                     MLE NEGATIVE BINOMIAL Regression: trade                 #
###############################################################################

# Standard negative binomial (NB2) with default SEs
nb_model<- glm.nb(trade ~ logGDP_im + logGDP_ex + logGDPC_im + logGDPC_ex + logdist + comborder + comlang +
                    colonial + comFTA + landlocked_ex + landlocked_im +  
                    remoteness_im + remoteness_ex + openness, data = tradedata)


summary(nb_model)
stargazer(nb_model,type="text")

#####################
#  Goodness of Fit  #
#####################

# COMPARABLE: squared correlation between fitted and actual counts
# "For the exponential model, the R-squared is computed as the squared correlation 
# coefficient between \textit{trade} and $\widehat{\textit{trade}}{=exp(x_{i}\widehat{\beta})}$" (w2 p653)

tradedata$yNBhat <- predict(nb_model, type = "response") # Predict expected NEAB counts

summary(tradedata[c("yNBhat", "trade")])

correlation_NB <- cor(tradedata$trade, tradedata$yNBhat)
squared_correlation_NB <- correlation_NB^2

cat("Squared correlation between y and yhat = ", squared_correlation_NB, "\n")

# the result is 0.4889337 which is not at all similar to the Poisson model's 
# 0.8616187, however this is likely indicative of the Poisson model's overdispersion.

###########################################
# Pearson chi-square goodness-of-fit test #
###########################################

P__disp(nb_model)

# pearson.chi2   dispersion 
# 95357.459377     5.198008

# Pearson chi-square goodness-of-fit test measures of how well the model fits observations.

# With a large pearson.chi2 statistic of 95357.459377 this model appears a better fit.

################################
# ZERO PREDICTIONS VS OBSERVED #
################################

# yNBhat is the expected mean count for each observation
# rounded sum of predict probability of zero count using expected mean counts
zero_predict_NB <- round(sum(dpois(x = 0, lambda = tradedata$yNBhat)))

# observed number of 0's
zero_actual_NB <- sum(tradedata$trade < 1)

cat("The model predicts ",zero_predict_NB, " zero counts but we observe", zero_actual_NB, 
    "zero counts \nThis indicates some underfitting of zero counts \nthough, not as bads as the Poisson's prediction of ",zero_predict_P)

#################################
#  Poisson & NEGATIVE BINOMIAL  #
#          Correlation          #
#################################

# Predict the linear index (X'B)
tradedata$XBP_hat <- predict(poisson_model, type = "link")
tradedata$XBNBhat <- predict(nb_model, type="link")

# Correlation between Poisson and Negative Binomial linear predictors
correlation_XB <- cor(tradedata$XBP_hat, tradedata$XBNBhat)

# Correlation between Poisson and Negative Binomial expected counts
correlation_y <- cor(tradedata$yPhat, tradedata$yNBhat)

cat("Correlation between Poisson and NegBin linear predictors: ", correlation_XB, "\n")
#Correlation between Poisson and NegBin linear predictors:  0.9876102 

cat("Correlation between Poisson and NegBin expected counts: ", correlation_y, "\n")
#Correlation between Poisson and NegBin expected counts:  0.790757 


# Using base R (remove unwanted columns)
tradedata$yPhat <- NULL
tradedata$XBP_hat <- NULL
tradedata$yNBhat <- NULL
tradedata$XBNBhat <- NULL

###############################################################################
#                 MARGINAL EFFECTS: POISSON & NEGATIVE BINOMIAL               #
###############################################################################

################
#     AME:     #
#   POISSON    #
################

# Compute average marginal effects for log distance
P_AMElogdist <- margins(poisson_model, var="logdist")

# Print the results
summary(P_AMElogdist)

#####################
#       AME:        #
# NEGATIVE BINOMIAL #
#####################

# Compute average marginal effects for log distance
NB_AMElogdist <- margins(nb_model, var="logdist")

# Print the results
summary(NB_AMElogdist)

#############################################
#       MARGINAL EFFECTS AT THE MEAN        #
#############################################

MEMtradedata <- data.frame(logGDP_im = rep(mean(tradedata$logGDP_im),1),
                           logGDP_ex = rep(mean(tradedata$logGDP_ex),1),
                           logGDPC_im = rep(mean(tradedata$logGDPC_im),1),
                           logGDPC_ex = rep(mean(tradedata$logGDPC_ex),1),
                           logdist = rep(mean(tradedata$logdist),1),
                           comborder = rep(mean(tradedata$comborder),1),
                           comlang = rep(mean(tradedata$comlang),1),
                           colonial = rep(mean(tradedata$colonial),1),
                           comFTA = rep(mean(tradedata$comFTA),1),
                           landlocked_ex = rep(mean(tradedata$landlocked_ex),1),
                           landlocked_im = rep(mean(tradedata$landlocked_im),1),
                           remoteness_im = rep(mean(tradedata$remoteness_im),1),
                           remoteness_ex = rep(mean(tradedata$remoteness_ex),1),
                           openness = rep(mean(tradedata$openness),1))

################
#     MEM:     #
#   POISSON    #
################

MEM_P <- margins(poisson_model, data=MEMtradedata, var="logdist")
summary(MEM_P)

#####################
#       MEM:        #
# NEGATIVE BINOMIAL #
#####################

MEM_NB <- margins(nb_model, data=MEMtradedata, var="logdist")
summary(MEM_NB)


###########
#   NAs   #
###########

logtrade_na<-mean(!trade, na.rm=TRUE)*100 #% of trade made up of NA
logtrade_na

###############################################################################
#                     OLS Linear Regression: log(trade)                       #
###############################################################################
logtrade_OLS<-lm(logtrade ~ logGDP_im + logGDP_ex + logGDPC_im + logGDPC_ex +
                   logdist + comborder + comlang +
                   colonial + comFTA + landlocked_ex + landlocked_im +  
                   remoteness_im + remoteness_ex + openness, data = tradedata)
summary(logtrade_OLS)
stargazer(logtrade_OLS,type="text")


##############
#  Problem   #
##############

# We have lost 18,360 - 9,613 observations

tradedata$yhat_logtrade_OLS <- predict(logtrade_OLS, tradedata, type = "response")

summary(tradedata[, c("logtrade","yhat_logtrade_OLS")])

ggplot(tradedata, aes(x = yhat_logtrade_OLS)) +
  geom_point(aes(y=logtrade), color="blue") +
  geom_point(aes(y=yhat_logtrade_OLS), color="red") +
  labs(title = "Scatter Plot of logtrade vs. logtrade_hat")

###############################################################################
#                     HURDLE: NEGATIVE BINOMIAL Regression: trade             #
###############################################################################

# Fit zero-truncated Negative Binomial
hurdle_NB_model <- hurdle(formula = trade ~ logGDP_im + logGDP_ex + logGDPC_im + logGDPC_ex + logdist + comborder + comlang + colonial + comFTA + landlocked_ex + landlocked_im + remoteness_im + remoteness_ex + openness|
                            logGDP_im + logGDP_ex + logGDPC_im + logGDPC_ex + logdist + comborder + comlang + colonial + comFTA + landlocked_ex + landlocked_im + remoteness_im + remoteness_ex + openness, 
                          data=tradedata, dist = "negbin", zero.dist = "binomial")
summary(hurdle_NB_model)
stargazer(hurdle_NB_model, type="text")

#####################
#   DATA Subset:    #
# Positive Outcomes #
#####################

tradedata_positive <- subset(tradedata, trade > 0)

############################
#         AME:             #
# HURDLE NEGATIVE BINOMIAL #
############################

# Saving coefficients for the respective parts
binary_logit_NB <- coef(hurdle_NB_model, model = "zero")
Zero_truncated_NegBin <- coef(hurdle_NB_model, model = "count")

#linear content of the logit step
tradedata$binary_logit_NBXB <- with(tradedata, binary_logit_NB["logGDP_im"]*logGDP_im + 
                                      binary_logit_NB["logGDP_ex"]*logGDP_ex + 
                                      binary_logit_NB["logGDPC_im"]*logGDPC_im + 
                                      binary_logit_NB["logGDPC_ex"]*logGDPC_ex + 
                                      binary_logit_NB["logdist"]*logdist + 
                                      binary_logit_NB["comborder"]*comborder + 
                                      binary_logit_NB["comlang"]*comlang + 
                                      binary_logit_NB["colonial"]*colonial + 
                                      binary_logit_NB["comFTA"]*comFTA + 
                                      binary_logit_NB["landlocked_ex"]*landlocked_ex + 
                                      binary_logit_NB["landlocked_im"]*landlocked_im + 
                                      binary_logit_NB["remoteness_im"]*remoteness_im + 
                                      binary_logit_NB["remoteness_ex"]*remoteness_ex + 
                                      binary_logit_NB["openness"]*openness + 
                                      binary_logit_NB["(Intercept)"])

#linear content of the poisson step
tradedata_positive$Zero_truncated_NegBinXB <- with(tradedata_positive, Zero_truncated_NegBin["logGDP_im"]*logGDP_im + 
                                                     Zero_truncated_NegBin["logGDP_ex"]*logGDP_ex + 
                                                     Zero_truncated_NegBin["logGDPC_im"]*logGDPC_im + 
                                                     Zero_truncated_NegBin["logGDPC_ex"]*logGDPC_ex + 
                                                     Zero_truncated_NegBin["logdist"]*logdist + 
                                                     Zero_truncated_NegBin["comborder"]*comborder + 
                                                     Zero_truncated_NegBin["comlang"]*comlang + 
                                                     Zero_truncated_NegBin["colonial"]*colonial + 
                                                     Zero_truncated_NegBin["comFTA"]*comFTA + 
                                                     Zero_truncated_NegBin["landlocked_ex"]*landlocked_ex + 
                                                     Zero_truncated_NegBin["landlocked_im"]*landlocked_im + 
                                                     Zero_truncated_NegBin["remoteness_im"]*remoteness_im + 
                                                     Zero_truncated_NegBin["remoteness_ex"]*remoteness_ex + 
                                                     Zero_truncated_NegBin["openness"]*openness + 
                                                     Zero_truncated_NegBin["(Intercept)"])

# exp[XB]/(1+exp[XB])
tradedata$NBT1 <- exp(tradedata$binary_logit_NBXB) / (1 + exp(tradedata$binary_logit_NBXB))

# (exp[XB]/(1+exp[XB]))^2
tradedata$NBT1d <- exp(tradedata$binary_logit_NBXB) / (1 + exp(tradedata$binary_logit_NBXB))^2

# exp[truncatedXB]/(1-exp(-exp[truncatedXB]))
nbt2 <- exp(tradedata_positive$Zero_truncated_NegBinXB) / (1 - exp(-exp(tradedata_positive$Zero_truncated_NegBinXB)))

# 1 - (exp(-exp[truncatedXB])*exp[truncatedXB])/(1-exp(-exp[truncatedXB]))
nbt3 <- 1 - (exp(-exp(tradedata_positive$Zero_truncated_NegBinXB)) * exp(tradedata_positive$Zero_truncated_NegBinXB)) / (1 - exp(-exp(tradedata_positive$Zero_truncated_NegBinXB)))

#PROBABILITY OF TRADE TAKING ANY VALUE GIVEN THAT TRADE>0 AND GIVEN X
# Ensure to only calculate T2 and T3 where y > 0
tradedata$NBT2 <- ifelse(tradedata$trade > 0, nbt2, NA)

tradedata$NBT3 <- ifelse(tradedata$trade > 0, nbt3, NA)


tradedata$Hurdle_NB_AMElogdist <- with(tradedata, NBT1d * NBT2 * binary_logit_NB["logdist"] + NBT1 * NBT2 * NBT3 * Zero_truncated_NegBin["logdist"])
mean(tradedata$Hurdle_NB_AMElogdist,var="logdist", na.rm = TRUE)#AME mean

############################
#           MEM:           #
# HURDLE NEGATIVE BINOMIAL #
############################

MEMhurdledata <- data.frame(logGDP_im = rep(mean(tradedata$logGDP_im),1),
                            logGDP_ex = rep(mean(tradedata$logGDP_ex),1),
                            logGDPC_im = rep(mean(tradedata$logGDPC_im),1),
                            logGDPC_ex = rep(mean(tradedata$logGDPC_ex),1),
                            logdist = rep(mean(tradedata$logdist),1),
                            comborder = rep(mean(tradedata$comborder),1),
                            comlang = rep(mean(tradedata$comlang),1),
                            colonial = rep(mean(tradedata$colonial),1),
                            comFTA = rep(mean(tradedata$comFTA),1),
                            landlocked_ex = rep(mean(tradedata$landlocked_ex),1),
                            landlocked_im = rep(mean(tradedata$landlocked_im),1),
                            remoteness_im = rep(mean(tradedata$remoteness_im),1),
                            remoteness_ex = rep(mean(tradedata$remoteness_ex),1),
                            openness = rep(mean(tradedata$openness),1),
                            trade = rep(floor(mean(tradedata$trade)),1))  

#linear content of the logit step
MEMhurdledata$binary_logit_NBXB <- with(MEMhurdledata, binary_logit_NB["logGDP_im"]*logGDP_im + 
                                          binary_logit_NB["logGDP_ex"]*logGDP_ex + 
                                          binary_logit_NB["logGDPC_im"]*logGDPC_im + 
                                          binary_logit_NB["logGDPC_ex"]*logGDPC_ex + 
                                          binary_logit_NB["logdist"]*logdist + 
                                          binary_logit_NB["comborder"]*comborder + 
                                          binary_logit_NB["comlang"]*comlang + 
                                          binary_logit_NB["colonial"]*colonial + 
                                          binary_logit_NB["comFTA"]*comFTA + 
                                          binary_logit_NB["landlocked_ex"]*landlocked_ex + 
                                          binary_logit_NB["landlocked_im"]*landlocked_im + 
                                          binary_logit_NB["remoteness_im"]*remoteness_im + 
                                          binary_logit_NB["remoteness_ex"]*remoteness_ex + 
                                          binary_logit_NB["openness"]*openness + 
                                          binary_logit_NB["(Intercept)"])

#linear content of the poisson step
MEMhurdledata$Zero_truncated_NegBinXB <- with(MEMhurdledata, Zero_truncated_NegBin["logGDP_im"]*logGDP_im + 
                                                Zero_truncated_NegBin["logGDP_ex"]*logGDP_ex + 
                                                Zero_truncated_NegBin["logGDPC_im"]*logGDPC_im + 
                                                Zero_truncated_NegBin["logGDPC_ex"]*logGDPC_ex + 
                                                Zero_truncated_NegBin["logdist"]*logdist + 
                                                Zero_truncated_NegBin["comborder"]*comborder + 
                                                Zero_truncated_NegBin["comlang"]*comlang + 
                                                Zero_truncated_NegBin["colonial"]*colonial + 
                                                Zero_truncated_NegBin["comFTA"]*comFTA + 
                                                Zero_truncated_NegBin["landlocked_ex"]*landlocked_ex + 
                                                Zero_truncated_NegBin["landlocked_im"]*landlocked_im + 
                                                Zero_truncated_NegBin["remoteness_im"]*remoteness_im + 
                                                Zero_truncated_NegBin["remoteness_ex"]*remoteness_ex + 
                                                Zero_truncated_NegBin["openness"]*openness + 
                                                Zero_truncated_NegBin["(Intercept)"])

# exp[XB]/(1+exp[XB])
MEMhurdledata$NBT1 <- exp(MEMhurdledata$binary_logit_NBXB) / (1 + exp(MEMhurdledata$binary_logit_NBXB))

# (exp[XB]/(1+exp[XB]))^2
MEMhurdledata$NBT1d <- exp(MEMhurdledata$binary_logit_NBXB) / (1 + exp(MEMhurdledata$binary_logit_NBXB))^2

# exp[truncated_XB]/(1-exp(-exp[truncated_XB]))
MEMnbt2 <- exp(MEMhurdledata$Zero_truncated_NegBinXB) / (1 - exp(-exp(MEMhurdledata$Zero_truncated_NegBinXB)))

# 1 - (exp(-exp[truncated_XB])*exp[truncated_XB])/(1-exp(-exp[truncated_XB]))
MEMnbt3 <- 1 - (exp(-exp(MEMhurdledata$Zero_truncated_NegBinXB)) * exp(MEMhurdledata$Zero_truncated_NegBinXB)) / (1 - exp(-exp(MEMhurdledata$Zero_truncated_NegBinXB)))

# PROBABILITY OF TRADE TAKING ANY VALUE GIVEN THAT TRADE>0 AND GIVEN X
# Ensure to only calculate T2 and T3 where y > 0

MEMhurdledata$NBT2 <- ifelse(MEMhurdledata$trade > 0, MEMnbt2, NA)
MEMhurdledata$NBT3 <- ifelse(MEMhurdledata$trade > 0, MEMnbt3, NA)

# logdist
MEMhurdledata$Hurdle_NB_MEMlogdist <- with(MEMhurdledata, NBT1d * NBT2 * binary_logit_NB["logdist"] + NBT1 * NBT2 * NBT3 * Zero_truncated_NegBin["logdist"])

#MEM
mean(MEMhurdledata$Hurdle_NB_MEMlogdist, na.rm = TRUE)

##############################################################################
#                        Zero-inflated: NEGATIVE BINOMIAL                    #
##############################################################################

Zero_Inflated_NBmodel <- zeroinfl(trade ~ logGDP_im + logGDP_ex + logGDPC_im + logGDPC_ex + logdist + comborder + comlang + colonial + comFTA + landlocked_ex + landlocked_im + remoteness_im + remoteness_ex + openness|
                                    logGDP_im + logGDP_ex + logGDPC_im + logGDPC_ex + logdist + comborder + comlang + colonial + comFTA + landlocked_ex + landlocked_im + remoteness_im + remoteness_ex + openness, 
                                  data=tradedata, dist = "negbin")
summary(Zero_Inflated_NBmodel)
stargazer(Zero_Inflated_NBmodel,type="text")

####################################
#              AME:                #
# Zero-inflated: NEGATIVE BINOMIAL #
####################################

marginsNB_count <- margins(Zero_Inflated_NBmodel, variables = "logdist", atmeans = TRUE)
summary(marginsNB_count)

#########################################
#                MEM:                   #
#    Zero-inflated: NEGATIVE BINOMIAL   # 
#########################################

MEMZIdata <- data.frame(logGDP_im = rep(mean(tradedata$logGDP_im),1),
                        logGDP_ex = rep(mean(tradedata$logGDP_ex),1),
                        logGDPC_im = rep(mean(tradedata$logGDPC_im),1),
                        logGDPC_ex = rep(mean(tradedata$logGDPC_ex),1),
                        logdist = rep(mean(tradedata$logdist),1),
                        comborder = rep(mean(tradedata$comborder),1),
                        comlang = rep(mean(tradedata$comlang),1),
                        colonial = rep(mean(tradedata$colonial),1),
                        comFTA = rep(mean(tradedata$comFTA),1),
                        landlocked_ex = rep(mean(tradedata$landlocked_ex),1),
                        landlocked_im = rep(mean(tradedata$landlocked_im),1),
                        remoteness_im = rep(mean(tradedata$remoteness_im),1),
                        remoteness_ex = rep(mean(tradedata$remoteness_ex),1),
                        openness = rep(mean(tradedata$openness),1),
                        trade = rep(floor(mean(tradedata$trade)),1))  


count_coefs_ziNB <- coef(Zero_Inflated_NBmodel, model = "count")  # Count model coefficients

inflate_coefs_ziNB <- coef(Zero_Inflated_NBmodel, model = "zero")  # Inflation model coefficients


MEMZIdata$NB_MEMinflXB_0 <- with(MEMZIdata,inflate_coefs_ziNB["(Intercept)"] +
                                   inflate_coefs_ziNB["logGDP_im"]*logGDP_im + 
                                   inflate_coefs_ziNB["logGDP_ex"]*logGDP_ex + 
                                   inflate_coefs_ziNB["logGDPC_im"]*logGDPC_im + 
                                   inflate_coefs_ziNB["logGDPC_ex"]*logGDPC_ex + 
                                   0* inflate_coefs_ziNB["logdist"]+
                                   inflate_coefs_ziNB["comborder"]*comborder + 
                                   inflate_coefs_ziNB["comlang"]*comlang + 
                                   inflate_coefs_ziNB["colonial"]*colonial + 
                                   inflate_coefs_ziNB["comFTA"]*comFTA + 
                                   inflate_coefs_ziNB["landlocked_ex"]*landlocked_ex + 
                                   inflate_coefs_ziNB["landlocked_im"]*landlocked_im + 
                                   inflate_coefs_ziNB["remoteness_im"]*remoteness_im + 
                                   inflate_coefs_ziNB["remoteness_ex"]*remoteness_ex + 
                                   inflate_coefs_ziNB["openness"]*openness)


MEMZIdata$NB_MEMXB_0 <- with(MEMZIdata, count_coefs_ziNB["(Intercept)"] +
                               count_coefs_ziNB["logGDP_im"]*logGDP_im + 
                               count_coefs_ziNB["logGDP_ex"]*logGDP_ex + 
                               count_coefs_ziNB["logGDPC_im"]*logGDPC_im + 
                               count_coefs_ziNB["logGDPC_ex"]*logGDPC_ex + 
                               0* count_coefs_ziNB["logdist"]+
                               count_coefs_ziNB["comborder"]*comborder + 
                               count_coefs_ziNB["comlang"]*comlang + 
                               count_coefs_ziNB["colonial"]*colonial + 
                               count_coefs_ziNB["comFTA"]*comFTA + 
                               count_coefs_ziNB["landlocked_ex"]*landlocked_ex + 
                               count_coefs_ziNB["landlocked_im"]*landlocked_im + 
                               count_coefs_ziNB["remoteness_im"]*remoteness_im + 
                               count_coefs_ziNB["remoteness_ex"]*remoteness_ex + 
                               count_coefs_ziNB["openness"]*openness)


MEMZIdata$NB_ziMEMlogdist <- with(MEMZIdata, (exp(NB_MEMXB_0 + count_coefs_ziNB["logdist"]) / (1 + exp(NB_MEMinflXB_0 + count_coefs_ziNB["logdist"]))) - 
                                    (exp(NB_MEMXB_0) / (1 + exp(NB_MEMinflXB_0))))

mean(MEMZIdata$NB_ziMEMlogdist)  #Zero-inflated NEGATIVE BINOMIAL

##################
# CENSORED TRADE #
##################

tradedata$trade_cen = tradedata$trade
tradedata$trade_cen <- tradedata$trade_cen+1
tradedata$log_trade_cen <- log(tradedata$trade_cen)

summary(tradedata$log_trade_cen) #Mean:4.42 
var(tradedata$log_trade_cen) #Variance: 23.32249
hist(tradedata$log_trade_cen)

###############################################################################
#                   OLS Linear Regression: log(trade + 1)                     #
###############################################################################

logtrade_cen_OLS<-lm(log_trade_cen ~ logGDP_im + logGDP_ex + logGDPC_im + logGDPC_ex +
                       logdist + comborder + comlang +
                       colonial + comFTA + landlocked_ex + landlocked_im +  
                       remoteness_im + remoteness_ex + openness, data = tradedata)
summary(logtrade_cen_OLS)
stargazer(logtrade_cen_OLS,type="text")

###############################################################################
#                  CENSORED Regression Model: log(trade + 1)                  #
###############################################################################

# Tobit regression for wage_cens
Censored_log <- censReg(log_trade_cen ~ logGDP_im + logGDP_ex + logGDPC_im + logGDPC_ex +
                          logdist + comborder + comlang +
                          colonial + comFTA + landlocked_ex + landlocked_im +  
                          remoteness_im + remoteness_ex + openness, left = 0, data = tradedata)
summary(Censored_log)
stargazer(Censored_log, type="text")

#########################################################
#       MARGINAL EFFECTS: Tobit aka Censored            #
#########################################################

C <- vglm(log_trade_cen ~ logGDP_im + logGDP_ex + logGDPC_im + logGDPC_ex +
            logdist + comborder + comlang +
            colonial + comFTA + landlocked_ex + landlocked_im +  
            remoteness_im + remoteness_ex + openness, data = tradedata, tobit(Lower = 0, Upper = Inf, type.f = "cens"), trace = TRUE)
summary(C)

#################
#     AME:      #
#   Censored    #
#################

XB_C <- predict(C, type = "response")    # Linear predictor
sigma_C <- exp(coef(C)["(Intercept):2"]) # Extract sigma (standard error) from model
meanXB_C <- mean(XB_C) # Mean of xb1
mPHIC <- pnorm(meanXB_C / sigma_C) # Normal CDF 

# Marginal Effects at the Mean (MEM)
AME_C <- data.frame(logGDP_im= mPHIC * coef(C)["logGDP_im"],1,
                    logGDP_ex=mPHIC * coef(C)["logGDP_ex"],1, 
                    logGDPC_im=mPHIC * coef(C)["logGDPC_im"],1,
                    logGDPC_ex=mPHIC * coef(C)["logGDPC_ex"],1, 
                    logdist=mPHIC * coef(C)["logdist"],1, 
                    comborder=mPHIC * coef(C)["comborder"],1,
                    comlang=mPHIC * coef(C)["comlang"],1,
                    colonial=mPHIC * coef(C)["colonial"],1, 
                    comFTA=mPHIC * coef(C)["comFTA"],1, 
                    landlocked_ex=mPHIC * coef(C)["landlocked_ex"],1,
                    landlocked_im=mPHIC * coef(C)["landlocked_im"],1,
                    remoteness_im=mPHIC * coef(C)["remoteness_im"],1,
                    remoteness_ex=mPHIC * coef(C)["remoteness_ex"],1,
                    openness=mPHIC * coef(C)["openness"],1)
mean(AME_C$logdist)

#################
#     MEM:      #
#   Censored    #
#################

MxB_C <- coef(C)["(Intercept):1"]+coef(C)["logGDP_im"]*mean(tradedata$logGDP_im)+coef(C)["logGDP_ex"]*mean(tradedata$logGDP_ex)+
  coef(C)["logGDPC_im"]*mean(tradedata$logGDPC_im)+coef(C)["logGDPC_ex"]*mean(tradedata$logGDPC_ex)+coef(C)["logdist"]*mean(tradedata$logdist)+
  coef(C)["comborder"]*mean(tradedata$comborder)+coef(C)["comlang"]*mean(tradedata$comlang)+
  coef(C)["colonial"]*mean(tradedata$colonial)+coef(C)["comFTA"]*mean(tradedata$comFTA)+
  coef(C)["landlocked_ex"]*mean(tradedata$landlocked_ex)+coef(C)["landlocked_im"]*mean(tradedata$landlocked_im)+
  coef(C)["remoteness_im"]*mean(tradedata$remoteness_im)+coef(C)["remoteness_ex"]*mean(tradedata$remoteness_ex)+
  coef(C)["openness"]*mean(tradedata$openness)

# Compute necessary components for marginal effects calculation
Msigma_C <- exp(coef(C)["(Intercept):2"]) # Extract sigma (standard error) from model
MeansxB_C <- mean(MxB_C) # Mean of xb
mPHIC_2 <- pnorm(MeansxB_C / Msigma_C) # Normal CDF of meanxb1/sigma1


MEM_C  <- data.frame(logGDP_im= mPHIC_2 * coef(C)["logGDP_im"],1,
                     logGDP_ex=mPHIC_2 * coef(C)["logGDP_ex"],1, 
                     logGDPC_im=mPHIC_2 * coef(C)["logGDPC_im"],1,
                     logGDPC_ex=mPHIC_2 * coef(C)["logGDPC_ex"],1, 
                     logdist=mPHIC_2 * coef(C)["logdist"],1, 
                     comborder=mPHIC_2 * coef(C)["comborder"],1,
                     comlang=mPHIC_2 * coef(C)["comlang"],1,
                     colonial=mPHIC_2 * coef(C)["colonial"],1, 
                     comFTA=mPHIC_2 * coef(C)["comFTA"],1, 
                     landlocked_ex=mPHIC_2 * coef(C)["landlocked_ex"],1,
                     landlocked_im=mPHIC_2 * coef(C)["landlocked_im"],1,
                     remoteness_im=mPHIC_2 * coef(C)["remoteness_im"],1,
                     remoteness_ex=mPHIC_2 * coef(C)["remoteness_ex"],1,
                     openness=mPHIC_2 * coef(C)["openness"],1)

mean(MEM_C$logdist)

detach(tradedata)

