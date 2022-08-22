rm(list = ls())

library("AER")
library("ggplot2")
library("sandwich")
library("stargazer")
library("survival")
library("plm")
library("margins")
library("readxl")
library("xts")
library("dynlm")
library("zoo")
library("urca")
library("strucchange")
library("orcutt")
library("fGarch")
library("quantmod")
library("tidyr")

setwd("~/Econometrics/Fama")

# Load CSV file into R
ff_data <- read.csv(file = "FF_Data_Annual.csv")
ff_data$Date <- format(as.yearmon(as.character(ff_data$Date), "%Y%m"), "%Y-%m")

##Creating xts
ff_data_xts <- xts(ff_data$Mkt.RF,order.by = as.yearmon(ff_data$Date))["2010::2018"]


##Estimating AR1
ff_data_ar1 <- lm(ff_data_xts ~ lag(ff_data_xts))

png("Summary_Results_of_Fama_French_Weekly_AR1.png")
summary(ff_data_ar1)
graphics.off()

coeftest(ff_data_ar1)
stargazer(ff_data_ar1,
          digits = 3,
          header = F,
          type = "html",
          title = "Fama French Risk Premium",
          out = "ff_data_ar1.doc")

##Using Bayes Information Criterion
BIC <- function(model) {
  ssr <- sum(model$residuals^2)
  t <- length(model$residuals)
  p <- length(model$coef) - 1
  return(
    round(c("p" = p, "BIC" = log(ssr/t) + (p + 1)*log(t)/t), 
          4)
  )
}

for (p in 1:5) {
  print(BIC(lm(ff_data_xts ~ lag(xts(ff_data_xts), 1:p))))
}##highest BIC(1:21) = 5.1658

##Testing for Violations of Key Time Series Assumptions
##Coefficients Biased Towards Zero?
beta_1_hat = ff_data_ar1$coefficients[2]
t_statistic = coeftest(ff_data_ar1)[2,3]
cat("Mean of beta_1_hat:", mean(beta_1_hat), "\n")
expected_value_analytical = 1-5.3/T
cat("Analytical Expected Value of beta_1_hat:", expected_value_analytical, "\n")

#testing for unit root in AR1
##H_0: beta_1 = 1: unit root process
##H_1: beta_1 < 1: stationarity present
AR1_DickeyFuller_lm <- lm(diff(ff_data_xts) ~ lag(ff_data_xts, 1))
coeftest(AR1_DickeyFuller_lm)
stargazer(AR1_DickeyFuller_lm,
          digits = 3,
          header = F,
          type = "html",
          title = "Fama French Weekly AR(1) First Differenced",
          out = "AR1_DickeyFuller_Test.doc")

unit_root_test = ur.df(ff_data_xts, lags = 0, type = 'drift')
DickeyFuller_teststat_automatic = unit_root_test@teststat[1]
DickeyFuller_critical_values = unit_root_test@cval[1]

png("Unit_Root_FF_Data.png")
summary(unit_root_test)
graphics.off()

DickeyFuller_critical_values < DickeyFuller_teststat_automatic

##Critical Values to test:
print("Dickey-Fuller Intercept-Only Critical Values:")
print(DickeyFuller_critical_values)
print(DickeyFuller_teststat_automatic)
##accept null hypothesis

##Taking first differences
ff_data_fdif <- lm(diff(ff_data_xts) ~ lag(ff_data_xts))
png("Coeftest_First_Dif_Fama.png")
coeftest(ff_data_fdif)
graphics.off()

for (p in 1:5) {
  print(BIC(lm(diff(ff_data_xts) ~ lag(xts(ff_data_xts), 1:p))))
}


ff_data_fdiff_1 <- lm(diff(ff_data_xts) ~ lag(ff_data_xts, 1))
coeftest(ff_data_fdiff_1)
summary(ff_data_fdiff_1)

##Accepting that there may be breaks using QLR Test
##QLR Test
num_periods <- length(ff_data_xts)
tau_0 = round(0.15*num_periods,digits=0)
tau_1 = round((1-0.15)*num_periods, digits = 0)
num_tests <- tau_1 - tau_0 + 1
tau <- seq(tau_0, tau_1)

D <- 1*(time(ff_data_xts) > time(ff_data_xts)[tau[1]]) 
Chow_test_model = lm(ff_data_xts ~ lag(ff_data_xts,1) + D + (D*lag(ff_data_xts,1)))
png("Chow_Test_Model_QLR.png")
coeftest(Chow_test_model)
graphics.off()

##set up the array of test statistics, equal in length to the number of
##tests we run
chow_test_statistics <- array(num_tests)
for (i in 1:num_tests) {
  D <- 1*(time(ff_data_xts) > time(ff_data_xts)[tau[i]])  
  Chow_test_model = lm(ff_data_xts ~ lag(ff_data_xts,1) + D + (D*lag(ff_data_xts,1)))
  chow_test = linearHypothesis(Chow_test_model, c("D=0", "lag(ff_data_xts, 1):D"), test="F", white.adjust = FALSE)
  chow_test_statistics[i] = chow_test$F[2]
}
data.frame("Level" = tau, "F-stat" = chow_test_statistics)

##pick the maximum Chow test stat
##This is the test statistic you need to compare against the appropriate critical values.
QLR_test_stat = max(chow_test_statistics)
cat("QLR Test Statistic: ", QLR_test_stat, "\n")
# This line tells us at WHICH tau we get the highest Chow test stat
tau_hat <- tau[which.max(chow_test_statistics)]
# Now, let's see which observation the likely break is at:
ff_data_xts[tau_hat]

# Now, we can see what happens AT the period where our QLR test thinks the break 
# was most likely by running the Chow test model and interpreting specific 
# coefficients.
D <- 1*(time(ff_data_xts) > time(ff_data_xts)[tau_hat])  
Chow_test_model = lm(ff_data_xts ~ lag(ff_data_xts,1) + D + (D*lag(ff_data_xts,1)))
png("Chow_Test_QLR_Coeft.png")
coeftest(Chow_test_model)
graphics.off()

##plot our data with a line at
##break period to see if the results look like what we see on the plot
png("QLR_Break_2009.png")
plot(as.zoo(ff_data_xts))
abline(v=as.yearmon("2009-02"))
graphics.off()
##Is coeficient on Yt-1 significant

##Municipal Solid Waste
msw_data <- read.csv(file = "Total_MSW_1960_2018.csv")
msw_data$Date <- format(as.yearmon(as.character(msw_data$Date), "%Y%m"), "%Y-%m")
##Creating xts
msw_data_xts <- xts(msw_data$Total.MSW,order.by = as.yearmon(msw_data$Date))["2010::2018"]

##Granger Test
ff_data_granger_22 <- lm(diff(ff_data_xts) ~ lag(ff_data_xts) + lag(msw_data_xts))
coeftest(ff_data_granger_22)

ff_granger_fd <- lm(diff(ff_data_xts) ~ lag(ff_data_xts) + diff(lag(msw_data_xts)))
coeftest(ff_granger_fd)
jpeg("granger-test")
linearHypothesis(ff_data_granger_22, 
                 c("lag(ff_data_xts)=0", 
                   "lag(msw_data_xts)=0"), 
                 test = "F", white.adjust = FALSE)
graphics.off()

##================================================================================================
##GLS Estimate of Dynamic Multipliers
##Commencing with Distributed Lag Model: DLM
## r= 3
ff_dlm <- lm(ff_data_xts ~ lag(ff_data_xts) + msw_data_xts + lag(msw_data_xts, 1) + lag(msw_data_xts, 2))
coeftest(ff_dlm)
##Assuming model exhibits exogeneity as lags of Yt and Xt
ols_uhat <- ff_dlm$residuals
v <- (msw_data_xts - mean(msw_data_xts))*ols_uhat
plot(as.zoo(ols_uhat), main = "DLM Residuals", ylab = "u_t", xlab = "t")
acf(ols_uhat, main = "Residual Autocorrelations")

delta_0_hat <- ff_dlm$coefficients[3]
delta_1_hat <- ff_dlm$coefficients[4]
phi_1_hat <- ff_dlm$coefficients[2] ##Confirm validity of selections

##Computing 2 Dynamic Mulitpliers
c("hat_beta_1" = delta_0_hat, "hat_beta_2" = delta_1_hat + (phi_1_hat*delta_0_hat))

##Feasible GLS
Residuals <- lm(ols_uhat ~ lag(ols_uhat))
coeftest(Residuals)
phi1_hat <- Residuals$coefficients[2]
cat("Phi_1_hat from manual Cochrane-Orcutt:", phi1_hat, "\n")

##==================================================================================================
##==================================================================================================


##Testing for Cointegration
co_ff_msw <- ff_data_xts - msw_data_xts
##Plotting ff data and msw data
plot(as.zoo(ff_data_xts),
     plot.type = "single",
     lty = 2,
     lwd = 2,
     col = "blue",
     xlab = "Date",
     ylab = "Risk Premium",
     ylim = c(0, 150),
     main = "Riskn Premia")
lines(as.zoo(msw_data_xts),
      col = "orange",
      lwd = 2,
      xlab = "Date",
      ylab = "Increase",
      main = "MSW Waste")
lines(as.zoo(co_ff_msw),
      col = "purple",
      lwd = 2,
      xlab = "Date",
      ylab = "ff data - msw data",
      main = "Spread for ff and msw")

##Cointegration test 1: Theta is known
z = co_ff_msw
ADFTest_ff_msw_spread <- ur.ers(co_ff_msw, 
                                        lag.max = 1, type = "DF-GLS",
                                        model = "constant")
jpeg("Cointegration With KNown Theta.jpeg")
summary(ADFTest_ff_msw_spread)
graphics.off()


##Cointegration Test: Estimating Theta
Cointegration_First_Step <- lm(ff_data_xts ~ msw_data_xts)
coeftest(Cointegration_First_Step)
theta_hat = Cointegration_First_Step$coefficients[2]

##Using theta hat to generate z
z_hat = ff_data_xts - (theta_hat*msw_data_xts)

##Testing for staionarity
ADF_Z_hat <- ur.df(z_hat, lags = 1, selectlags = "BIC", type = "none")
jpeg("Augmented DF.jpeg")
summary(ADF_Z_hat)
graphics.off()

##DF-GLS test on z
GLS_DF_test_z_hat <- ur.ers(z_hat, lag.max = 1, type = "DF-GLS", model = "constant")
jpeg("GLS_DF.jpeg")
summary(GLS_DF_test_z_hat)
graphics.off()