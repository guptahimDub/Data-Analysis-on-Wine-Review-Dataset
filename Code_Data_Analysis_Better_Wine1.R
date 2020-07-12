library(sqldf)
library(readxl)
df<-read_excel("winemag-data-130k-v2.xlsx")

#View(data_all)
#remove na rows
#df<-droplevels(na.omit(data_all))

#Data Handling
dataframe <- df[ which(df$price =='15'),]
dataframe<-sqldf("select country, price, points,variety from dataframe")

dfa <- dataframe[ which(dataframe$country== "Chile"),]
dfa <- dfa[ which(dfa$variety== "Chardonnay"),]

dfb <- dataframe[ which(dataframe$country=='South Africa'),]
dfb <- dfb[ which(dfb$variety== "Sauvignon Blanc"),]

new_df = rbind(dfa,dfb) #combine both dataframes
new_df <- droplevels(na.omit(new_df)) #drop the records having NA values

library(ggplot2)
ggplot(new_df) + geom_boxplot(aes(variety, points, fill = variety)) + 
  geom_jitter(aes(variety, points, shape = new_df$variety))

findQuantile1 <- new_df[ which(new_df$variety =='Chardonnay'),]
quantile(findQuantile1$points)

findQuantile2 <- new_df[ which(new_df$variety =='Sauvignon Blanc'),]
quantile(findQuantile2$points)


tapply(new_df$points, new_df$variety, mean)
tapply(new_df$points, new_df$variety, median)
tapply(new_df$points, new_df$variety, sd)

c = tapply(new_df$points, new_df$variety, mean)

#Difference in mean in observed data between South Africa and Chile
DiffOfMean = c['Sauvignon Blanc'] - c['Chardonnay']
print('Difference in mean in observed data between South Africa and Chile:')
print(DiffOfMean)

#t-test sample:
t.test(points ~ variety, data=new_df, var.equal = TRUE)


#Bayesian Model
compare_points_gibbs <- function(y, ind, mu0 = 50, tau0 = 1/400, del0 = 0, gamma0 = 1/400,
                                 a0 = 1, b0 = 50, maxiter = 10000)
{
  print('Value of IND')
  print(ind)
  y1 <- y[ind == 'Chardonnay']
  y2 <- y[ind == 'Sauvignon Blanc']
  n1 <- length(y1)
  n2 <- length(y2)
  ##### starting values
  mu <- (mean(y1) + mean(y2)) / 2
  del <- (mean(y1) - mean(y2)) / 2
  mat_store <- matrix(0, nrow = maxiter, ncol = 3)
  #####
  ##### Gibbs sampler
  an <- a0 + (n1 + n2)/2
  for(s in 1 : maxiter)
  {
    ##update tau
    bn <- b0 + 0.5 * (sum((y1 - mu - del) ^ 2) + sum((y2 - mu + del) ^ 2))
    tau <- rgamma(1, an, bn)
    ##
    ##update mu
    taun <- tau0 + tau * (n1 + n2)
    mun <- (tau0 * mu0 + tau * (sum(y1 - del) + sum(y2 + del))) / taun
    #print('MUN')
    #print(mun)
    mu <- rnorm(1, mun, sqrt(1/taun))
    ##
    ##update del
    gamman <- tau0 + tau*(n1 + n2)
    deln <- ( del0 * tau0 + tau * (sum(y1 - mu) - sum(y2 - mu))) / gamman
    del<-rnorm(1, deln, sqrt(1/gamman))
    ##
    ## store parameter values
    mat_store[s, ] <- c(mu, del, tau)
    #print(mat_store[s, ])
  }
  colnames(mat_store) <- c("mu", "del", "tau")
  #print(mat_store)
  return(mat_store)
}


library(MCMCpack)

##Calling function
fit <- compare_points_gibbs(new_df$points, as.factor(new_df$variety))

library(MCMCpack)
plot(as.mcmc(fit))

raftery.diag(as.mcmc(fit))

apply(fit, 2, mean)

apply(fit, 2, sd)


mean(1/sqrt(fit[, 3]))

sd(1/sqrt(fit[, 3]))

y1_sim <- rnorm(10000, fit[, 1] + fit[, 2], sd = 1/sqrt(fit[, 3]))
y2_sim <- rnorm(10000, fit[, 1] - fit[, 2], sd = 1/sqrt(fit[, 3]))

ggplot(data.frame(y_sim_diff = y1_sim - y2_sim)) + stat_bin(aes(y_sim_diff))

#How much better is the Wine:

#difference  = mean(y2_sim) - mean(y1_sim)

print('Probability that the Sauvignon Blanc will be better')
mean(y2_sim > y1_sim)
difference  = mean(y2_sim - y1_sim)
print(difference)

ggplot(data.frame(y1_sim, y2_sim)) + geom_point(aes(y1_sim, y2_sim), alpha = 0.3) +
  geom_abline(slope = 1, intercept = 0)
