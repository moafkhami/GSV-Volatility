#check for required packages
if (!require(xts, tseries, lmtest, fGarch)){ 
  install.packages('xts', 'tseries', 'lmtest', 'fGarch')}
if (!require(dplyr, sandwich)) install.packages('dplyr', 'sandwich')
#importing daily prices
df <- read.csv("Prices.csv")
df$Date <- as.Date( as.character(df$Date), "%d-%b-%y")
# library(dplyr) 
# df <- df %>% #alternative to line 10
#   filter (Date > as.Date("2004-01-03") & Date < as.Date("2016-07-23"))
df <- subset(df, as.Date("2004-01-01") < Date & Date < as.Date("2016-07-24"))

#returns
returns <- diff(as.matrix(log(df[,-1])))
returns <- data.frame(df$Date[-1], returns)
colnames(returns)[1] <- c("Date")

#weekly volatilites
library(xts)
colSd <- function (x, na.rm=FALSE) {apply(X=x, MARGIN=2,
                                          FUN=sd, na.rm=na.rm)}
returns <- as.xts(returns[, 2:ncol(returns)], order.by = as.Date(returns$Date))
volatilities <- apply.weekly(returns, colSd)
volatilities <- data.frame(Date=index(volatilities),
                           coredata(volatilities))

#merge with trends data
trends <- read.csv("trends.csv")
df <- cbind(volatilities, trends[,-1])

#ADF unit-root test
library(tseries)
suppressWarnings(apply(volatilities[2:ncol(volatilities)], 2, adf.test))

#Granger Causality
library(lmtest)
GrangerTable <- matrix(NA, nrow = ncol(trends)-1, ncol = ncol(volatilities)-1)
for (i in 1:nrow(GrangerTable)){
  for (j in 1:ncol(GrangerTable)){
    pval = grangertest(df[,j+1] ~ df[,i+ncol(GrangerTable)+1], order=2, data=df)$Pr[2]
    GrangerTable[i,j] = ifelse(pval < 0.05, signif(pval, digits = 2), "--")
    #GrangerTable[i,j] = grangertest(df[,j+1] ~ df[,i+ncol(GrangerTable)+1], order=2, data=df)$Pr[2]
  }
}
colnames(GrangerTable) <- colnames(volatilities)[-1]
rownames(GrangerTable) <- colnames(trends)[-1]
# library(xtable)
# xtable(GrangerTable)

#Garch estimates
library(fGarch)
WeeklyRet <- apply.weekly(returns, mean)
#GARCH (1,1) model
GarchEst <- matrix(NA, nrow = 2*ncol(WeeklyRet), ncol = 7)
#Residual and cond. var. extraction
Params <- matrix(NA, nrow = nrow(WeeklyRet), ncol = 2*ncol(WeeklyRet))
for (i in 1:ncol(WeeklyRet)){
  fitted <- garchFit(formula = ~garch(1,1), data = WeeklyRet[,i], cond.dist="std")
  GarchEst[2*i-1, 1] = colnames(WeeklyRet)[i]
  GarchEst[2*i-1, 2:5] = fitted@fit$coef[1:4]
  GarchEst[2*i, 2:5] = paste("(",fitted@fit$tval[1:4], ")")
  GarchEst[2*i-1, 6] = fitted@fit$llh
  GarchEst[2*i-1, 7] = fitted@fit$ics[1]
  Params[, 2*i-1] = log(residuals(fitted)^2)
  Params[, 2*i] = fitted@h.t
  #Extracting residuals and conditional variance
}
GarchEst[is.na(GarchEst)] <- ""
colnames(GarchEst) <- c("Commodity", "mu", "omega", "alpha", "beta", "logL", "AIC" )
#xtable(GarchEst)
colnames(Params) <- c(1:ncol(Params))
for (i in 1:ncol(WeeklyRet)){
  colnames(Params)[2*i-1] <- paste0("s_", colnames(WeeklyRet)[i], sep='')
  colnames(Params)[2*i] <- paste0("v_", colnames(WeeklyRet)[i], sep='')
}

df <- cbind (Params[-1,], trends[-ncol(trends),]) #using lagged trends
#df <- cbind(Params, trends)
library(dplyr)
df <- df %>%
  select(Date, everything())

#F-test and tables 4-5-6
library(lmtest)
library(sandwich)
onelist <- list()
twolist <- list()
threelist <- list()
k1 = k2 = k3 = 1
for (i in 1:ncol(WeeklyRet)){
  z1 <- lm(df[,2*i] ~ df[,2*i+1], data = df) #original regression
  sumz1 <- summary(z1)
  err1 <- coeftest(z1, vcov = vcovHC(z1, type="HC1")) #White Errors
  for (j in 2:ncol(trends)){
    if (GrangerTable[j - 1, i]!= "--"){ #Granger causality condition
      z2 <- lm(df[,2*i] ~ df[,2*i + 1] + df[,2*ncol(WeeklyRet) + j], data = df)
      sumz2 <- summary(z2) #new equation
      err2 <- coeftest(z2, vcov = vcovHC(z2, type="HC1"))
      if (anova(z1, z2)$Pr[2] <= 0.05 && err2[3,4] <= 0.05
          && sumz2$adj.r.squared > 1.1* sumz1$adj.r.squared){ #enhancement requirements
        dat1 = c(colnames(WeeklyRet)[i], colnames(trends[j]),
                signif(sumz2$coefficients[1:3,1], 3),
                signif(sumz2$coefficients[1:3,3], 3),
                signif(sumz2$adj.r.squared, 3)) 
        onelist[[k1]] = dat1 #append
        k1 = k1 + 1
        for (l in (j+1):ncol(trends)){ #two keywords
          z3 <- lm(df[,2*i] ~ df[,2*i+1] +
                     df[,2*ncol(WeeklyRet) + j] + df[,2*ncol(WeeklyRet) + l] ,
                   data = df)
          sumz3 <- summary(z3)
          err3 <- coeftest(z3, vcov = vcovHC(z3, type="HC1"))
          if (anova(z2, z3)$Pr[2] <= 0.05 && err3[4,4] <= 0.05
              && sumz3$adj.r.squared > 1.1* sumz2$adj.r.squared){ 
            dat2 = c(colnames(WeeklyRet)[i], colnames(trends[j]),
                     colnames(trends[l]),
                    signif(sumz3$coefficients[1:4,1], 3),
                    signif(sumz3$coefficients[1:4,3], 3),
                    signif(sumz3$adj.r.squared, 3)) 
            twolist[[k2]] = dat2
            k2 = k2 + 1
            if (l < ncol(trends)){
              for (m in (l+1):ncol(trends)){
                z4 <- lm(df[,2*i] ~ df[,2*i+1] +
                           df[,2*ncol(WeeklyRet) + j] + df[,2*ncol(WeeklyRet) + l] +
                           df[,2*ncol(WeeklyRet) + m], data = df)
                sumz4 <- summary(z4)
                err4 <- coeftest(z4, vcov = vcovHC(z4, type="HC1"))
                if (anova(z3, z4)$Pr[2] <= 0.05 && err4[5,4] <= 0.05
                    && sumz4$adj.r.squared > 1.1* sumz3$adj.r.squared){
                  imp=(100*(sumz4$adj.r.squared-sumz1$adj.r.squared))/sumz1$adj.r.squared
                  dat3 = c(colnames(WeeklyRet)[i], colnames(trends[j]),
                           colnames(trends[l]),
                           colnames(trends[m]),
                           signif(sumz4$coefficients[1:5,1], 3),
                           signif(sumz4$coefficients[1:5,3], 3),
                           signif(sumz4$adj.r.squared, 3),
                           #imp)
                           round(imp, 3)) 
                  threelist[[k3]] = dat3
                  k3 = k3 + 1
                }
            }
            }
          }
        }
      }
    }
  }
}
OneKey <- data.frame(matrix(unlist(onelist), nrow=length(onelist), byrow=T))
colnames(OneKey) <- c("Commodity", "Term", "beta_0", "beta_1", "k_1",
                      "t_beta_0", "t_beta_1", "t_k_1", "Adj.R.Sq")
#OneKey <- OneKey[order(OneKey$Adj.R.Sq,decreasing = TRUE),]
TwoKey <- data.frame(matrix(unlist(twolist), nrow=length(twolist), byrow=T))
colnames(TwoKey) <- c("Commodity", "Term1", "Term2", "beta_0", "beta_1", "k_1", "k_2",
                      "t_beta_0", "t_beta_1", "t_k_1", "t_k_2", "Adj.R.Sq")
#TwoKey <- TwoKey[order(TwoKey$Adj.R.Sq,decreasing = TRUE),]
ThreeKey <- data.frame(matrix(unlist(threelist), nrow=length(threelist), byrow=T))
colnames(ThreeKey) <- c("Commodity", "Term1", "Term2", "Term3", "beta_0", "beta_1",
                      "k_1", "k_2", "k_3", "t_beta_0", "t_beta_1", "t_k_1", "t_k_2",
                      "t_k_3", "Adj.R.Sq", "improvement")
ThreeKey <- ThreeKey[order(ThreeKey$Adj.R.Sq, decreasing = TRUE),]
#rm(onelist, twolist, threelist)
