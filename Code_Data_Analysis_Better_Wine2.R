library(sqldf)
library(readxl)
library(dplyr)
dfa <- read_excel("winemag-data-130k-v2.xlsx")
df <- sqldf("select country, price, points,variety,region_1 from dfa")
#View(df)

df1 <- df[ which(df$country== "Italy"),]
#dim(df1)


df1 <- df1[ which(df1$price < 20),]

length(which(is.na(df1)))
df1 <- na.omit((df1))
dim(df1)
df1 <- unique(df1)
#View(df1)
#################################
#Finding the dataframe having the region_1 values/reviews greater than 3
df1<-df1 %>% 
  group_by(region_1) %>% 
  filter(n()>3)
#View(df1)

####################################
#Find the average of points on observed data (price < 20 and region occurrence > 3 and country = Italy)
averagePoints <- mean(df1$points)

#####################################
#Comparing Multiple Countries w.r.t points
df1$region_1<- factor(df1$region_1)
nlevels(df1$region_1)

#Plot
ggplot(df1) + geom_boxplot(aes(x = reorder(region_1, points, median), points,
                               fill = reorder(region_1, points, median)), show.legend=FALSE)

#Distribution
ggplot(df1, aes(x = reorder(region_1, region_1, length))) + stat_count()

ggplot(df1, aes(points)) + stat_bin(binwidth = 1)

              ggplot(data.frame(size = tapply(df1$points, df1$region_1, length),
                  mean_score = tapply(df1$points, df1$region_1, mean)), aes(size, mean_score)) +
  geom_point() + xlab("Region 1 sample size") + ylab("Mean Score") +
  ggtitle("Mean Rating versus Sample size of Region")

#str(df1)

#Comparing using the Gibbs sampling:
compare_m_gibbs <- function(y, ind, mu0 = 50, tau0 = 1/400,
                            a0 = 1, b0 = 50, alpha0 =1, beta0 = 50, maxiter = 2)
{
  ### weakly informative priors
  a = ind
  b = as.numeric(ind)
  a0 <- 1/2 ; b0 <- 50 ## tau_w hyperparameters
  alpha0 <-1/2 ; beta0 <- 50 ## tau_b hyperparameters
  mu0<-50 ; tau0 <- 1/25
  #f= 1
  #f1 = 1
  # print('Print IND')
  # print(ind)
  #f=ind
  ####
  ### starting values
  #print(ind)
  m <- nlevels(ind)
  #print('Value of M Start')
  #print(m)
  #print('Value of M End')
  ybar <- theta <- tapply(y, ind, mean)
  #print('tapply(y, ind, var)')
  #print(mean(1/tapply(y, ind, var)))
  tau_w <- 1/mean(tapply(y, ind, var)) ##within group precision
  #print('tau_w')
  #print(tau_w)
  mu <- mean(theta)
  #print('mu')
  #print(mu)
  tau_b <-var(theta) ##between group precision
  #print('tau_b')
  #print(tau_b)
  n_m <- tapply(y, ind, length)
  #print('n_m')
  #print(n_m)
  alphan <- alpha0 + sum(n_m)/2
  #print('alphan')
  #print(alphan)
  ###
  ### setup MCMC
  theta_mat <- matrix(0, nrow=maxiter, ncol=m)
  mat_store <- matrix(0, nrow=maxiter, ncol=3)
  ###
  ### MCMC algorithm
  for(s in 1:maxiter)
  {
    # sample new values of the thetas
    for(j in 1:m)
    {
      #print('n_m[j] * tau_w  ')
      #print(n_m[j] * tau_w)
      #print(n_m[j])
      #print('tau_w')
      #print(tau_w)
      taun <- n_m[j] * tau_w + tau_b
      thetan <- (ybar[j] * n_m[j] * tau_w + mu * tau_b) / taun
      #print('rnorm(1, thetan, 1/sqrt(taun))')
      #print(rnorm(1, thetan, 1/sqrt(taun)))
      theta[j]<-rnorm(1, thetan, 1/sqrt(taun))
    }
    #sample new value of tau_w
    ss <- 0
    for(j in 1:m){
      j <- toString(j)
      ind <- as.numeric(ind)
      ss <- ss + sum((y[ind == j] - theta[as.numeric(j)])^2)
    }
    betan <- beta0 + ss/2
    tau_w <- rgamma(1, alphan, betan)
    #sample a new value of mu
    taum <- m * tau_b + tau0
    mum <- (mean(theta) * m * tau_b + mu0 * tau0) / taum
    mu <- rnorm(1, mum, 1/ sqrt(taum))
    # sample a new value of tau_b
    am <- a0 + m/2
    bm <- b0 + sum((theta - mu)^2) / 2
    tau_b <- rgamma(1, am, bm)
    #store results
    theta_mat[s,] <- theta
    mat_store[s, ] <- c(mu, tau_w, tau_b)
  }
  #print(mat_store)
  colnames(mat_store) <- c("mu", "tau_w", "tau_b")
  return(list(params = mat_store, theta = theta_mat, a=a, b=b))
}


fit2 <- compare_m_gibbs(df1$points, df1$region_1)
df3 <- data.frame(fit2$a,fit2$b)
df3 <- unique(df3)
df3 <- df3[order(df3$fit2.b, decreasing = FALSE),]

apply(fit2$params, 2, mean)
apply(fit2$params, 2, sd)
mean(1/sqrt(fit2$params[, 2]))
sd(1/sqrt(fit2$params[, 2]))
mean(1/sqrt(fit2$params[, 3]))
sd(1/sqrt(fit2$params[, 3]))

theta_hat <- apply(fit2$theta, 2, mean) ## get basic posterior summary

names(theta_hat) <- df3$fit2.a
sorted <- sort(theta_hat, decreasing = TRUE) ## which region's wine rated best and worst?
index<-  sorted > averagePoints # TRUE or FALSE
cat("Regions having ratings better than average value:",names(sorted[index]),sep=",",fill=TRUE)


print("Regions having ratings better than average value:")
print(names(sorted[index]))
#write(paste(as.character(names(sorted[index])), collapse=","),"note.txt",append="TRUE")
cat("Regions having ratings better than average value:",names(sorted[index]),sep=",",fill=TRUE)
      
View(names(sorted[index]))

theta_ci <- apply(fit2$theta, 2, quantile, prob = c(0.025, .975)) ## upper/lower bounds for thetas
df_error <- data.frame(lower = theta_ci[1, ], upper = theta_ci[2, ], mean = theta_hat,
                       school = factor(1:150))
ggplot(df_error, aes(x = reorder(school, mean), mean)) + geom_errorbar(aes(ymin = lower, ymax = upper))

## reformat samples for ggplot
theta_df <- data.frame(samples = as.numeric(fit2$theta),
                       region_1 = rep(1:ncol(fit2$theta), each = nrow(fit2$theta)))
ggplot(theta_df) + geom_boxplot(aes(x = reorder(region_1, samples, median), samples,
                                    fill = reorder(region_1, samples, median)), show.legend=FALSE)

ggplot(data.frame(size = tapply(df1$points, df1$region_1, length), theta_hat = theta_hat),
       aes(size, theta_hat)) + geom_point()

ggplot(data.frame(ybar = tapply(df1$points, df1$region_1, mean), theta_hat = theta_hat),
       aes(ybar, theta_hat)) + geom_point()


