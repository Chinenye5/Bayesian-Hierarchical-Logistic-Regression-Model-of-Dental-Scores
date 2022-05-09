library(mistr); library(readxl); library(tidyr);
library(tidyverse) ;library(standardize) ;
library(coda) ;library(lme4) ;library(MCMCvis);
library(bayesplot) ;library(R2jags); 
library(R2OpenBUGS) ;library(plyr) ;library(pander);
library(Matrix) ;library(lattice);
library("readxl");

BDA1 <- read_excel("C:\\Users\\Data.xlsx")

BDA1$Gums_d_m[BDA1$Gums_d_m=="Healthy"] <- "1"
BDA1$Gums_d_m[BDA1$Gums_d_m=="Not healthy"] <- "0"

BDA1$training[BDA1$training=="video"] <- "1"
BDA1$training[BDA1$training=="standard"] <- "0"

BDA1$Gender[BDA1$Gender=="female"] <- "1"
BDA1$Gender[BDA1$Gender=="male"] <- "0"

BDA1$Assessment_ease[BDA1$Assessment_ease=="Difficult"] <- "1"
BDA1$Assessment_ease[BDA1$Assessment_ease=="Easy"] <- "0"

BDA1$education[BDA1$education=="Yes"] <- "1"
BDA1$education[BDA1$education=="No"] <- "0"

BDA1$Gums[BDA1$Gums=="Healthy"] <- "1"
BDA1$Gums[BDA1$Gums=="Not healthy"] <- "0"

colnames(BDA1)<-c("FAC","SUB","SEX","AGE","EASY",
                  "TRAIN","EDU","SEX_ND","GUM_ND","GUM_D")




v=table(BDA1$SUB)

measurement=unlist(lapply(v, function(i) 1:i))

BDA1$Measurement = unlist(lapply(v, function(i) 1:i))

#p= pivot_wider(data=BDA1,id_cols =c("FAC","SUB","SEX","AGE","EASY",
#          "SEX_ND","Measurement"),  
#values_from =C("EDU","TRAIN","GUM_D"), 
#names_from =Measurement )






####### Part 1#######

BDA2H<- BDA1[BDA1$GUM_D==1, ]

BDA2H$GUM_ND=as.numeric(BDA2H$GUM_ND)
BDA2H$TRAIN=as.numeric(BDA2H$TRAIN)

nrep=table(BDA2H$SUB)

library(data.table)

s=dcast(setDT(BDA2H) , FAC+SUB+SEX+AGE+EASY ~ Measurement, 
        value.var = c("TRAIN", "EDU","SEX_ND", "GUM_ND"
        ))

nfac =table(BDA2H$FAC)

glm(GUM_ND~TRAIN, 
    data=BDA2H,family = binomial(link = logit))



model1<-"model{

for(i in 1:Nobs) {			
  for(j in 1:nrep[i]){		
    y[i, j] ~ dbern(pi[i,j])	
    logit(pi[i, j]) = FAC_random[i] + SUB[i,j] 
    SUB[i,j] ~ dnorm(beta[1] + beta[2]*TRAIN[i,j] , 
    sigma_sub_inv)
    
     
  }		
  
  FAC_random[i] ~ dnorm(0, sigma_fac_inv)		
  		
  
}			

#Hyperpriors for SUB Random Effect

sigma.sub ~ dunif(0, 20)
sigma_sub_inv <- 1/pow(sigma.sub ,2)

#Hyperpriors for FAC Random Effect

sigma.fac ~ dunif(0, 20)
sigma_fac_inv <- 1/pow(sigma.fac ,2)


#Priors for Fixed Effects

for (p in 1:2){

beta[p] ~ dnorm(0,0.0001)
 }



}"



library(rjags)

library(R2jags)

dat <- list(y=s[,c("GUM_ND_1","GUM_ND_2", "GUM_ND_3", "GUM_ND_4")]
            ,Nobs=nrow(s), nrep=nrep,
            TRAIN=s[,c("TRAIN_1" ,
                       "TRAIN_2" , "TRAIN_3",  "TRAIN_4")])

inits = list(
  list(beta=c(2.0502,-0.3638),sigma.sub=0.5,sigma.fac=1.5),
  list(beta=c(2.0502,-0.3638)*2,sigma.sub=0.5*2,sigma.fac=1.5*2),
  list(beta=c(2.0502,-0.3638)/2,sigma.sub=0.5/2, sigma.fac=1.5/2)
)


model <- jags(model.file = textConnection
              (model1),data = dat,n.chains=3, inits = inits,
              parameters.to.save=c("beta","sigma.sub",
                                   "sigma.fac"), n.burnin = 40000, 
              n.iter = 50000, n.thin = 1, DIC=T)



model_11 <- jags(model.file = textConnection
                 (model1),data = dat,n.chains=3, inits = inits,
                 parameters.to.save=c("SUB","FAC_random"), n.burnin = 40000, 
                 n.iter = 50000, n.thin = 1, DIC=T)

model_11$BUGSoutput$DIC

model_11$BUGSoutput$pD


mmm6=model_11$BUGSoutput$summary[, c("mean", "2.5%",
                                     "97.5%")] 

hist(mmm6[str_detect(row.names(mmm6),"SUB")],
     xlab = "Subjects", main=
       "Distribution plot for Random effect Subject")

hist(mmm6[str_detect(row.names(mmm6),"FAC_random")],
     xlab = "Facilities", main=
       "Distribution plot for Random effect Facilities")

o=mmm6[str_detect(row.names(mmm6),"SUB"),]
f=mmm6[str_detect(row.names(mmm6),"FAC_random"),]


subject=o[,1]
facilities=f[,1]


pp=data.frame(subject)
rr=data.frame(facilities)


qqnorm(o)
qqline(o, col = "red", lwd = 2)


qqnorm(n)
qqline(n, col = "red", lwd = 2)




ggplot(pp, aes(x=subject))+
  geom_histogram(color="black", fill="#800000")+
  theme_bw()


ggplot(rr, aes(x=facilities))+
  geom_histogram(color="black", fill="#808000")+
  theme_bw()


qplot(r, data = rr, geom = "histogram", bins=30,
      main = "Histogram for Random effect Facilities",
      fill=I("blue"))+
  theme_bw()

qplot(p, data = pp, geom = "histogram", bins=30,
      main = "Histogram for Random effect Subjects",
      fill=I("pink") ,
      xlab ="Subjects")+
  theme_bw()

ggplot(pp, aes(sample=subject))+stat_qq() + theme_bw()
ggplot(rr, aes(sample=facilities))+stat_qq() + theme_bw()

ggplot(rr) + geom_qq(aes(sample = facilities))+
  geom_abline(intercept = -18, slope = 6,
              color = "red", size = 1.3, alpha = 0.5)



#t-distribution for the random effacets


model16<-"model{

for(i in 1:Nobs) {			
  for(j in 1:nrep[i]){		
    y[i, j] ~ dbern(pi[i,j])	
    logit(pi[i, j]) = FAC_random[i] + SUB[i,j] 
    SUB[i,j] ~ dt(beta[1] + beta[2]*TRAIN[i,j],
    sigma_sub_inv, df_sub)
    
  }		
  
  FAC_random[i] ~ dt(0, sigma_fac_inv, df_fac)		
  

}			

#Hyperpriors for SUB Random Effect

sigma.sub ~ dunif(0, 20)
sigma_sub_inv <- 1/pow(sigma.sub ,2)
df_sub ~ dgamma(0.1,0.1)

#Hyperpriors for FAC Random Effect

sigma.fac ~ dunif(0, 20)
sigma_fac_inv <- 1/pow(sigma.fac ,2)
df_fac ~ dgamma(0.1,0.1)

#Priors for Fixed Effects

for (p in 1:2){

beta[p] ~ dnorm(0,0.0001)
 }



}"



dat <- list(y=s[,c("GUM_ND_1","GUM_ND_2", "GUM_ND_3", "GUM_ND_4")]
            ,Nobs=nrow(s), nrep=nrep,
            TRAIN=s[,c("TRAIN_1" ,
                       "TRAIN_2" , "TRAIN_3",  "TRAIN_4")])

inits = list(
  list(beta=c(2.0502,-0.3638),sigma.sub=0.5,sigma.fac=1.5),
  list(beta=c(2.0502,-0.3638)*2,sigma.sub=0.5*2,sigma.fac=1.5*2),
  list(beta=c(2.0502,-0.3638)/2,sigma.sub=0.5/2, sigma.fac=1.5/2)
)



model_16 <- jags(model.file = textConnection
                 (model16),data = dat,n.chains=3, inits = inits,
                 parameters.to.save=c("SUB","FAC_random"), n.burnin = 40000, 
                 n.iter = 50000, n.thin = 1, DIC=T)


model_16$BUGSoutput$DIC
model_16$BUGSoutput$pD




summary(model)


model$BUGSoutput$DIC

model$BUGSoutput$pD

###Checking convergence of randomly 
#selected random effects

library(mcmcplots)
library(ggmcmc)

model_mcmc= as.mcmc(model)
model_reffcts_ggs <- ggs(model_mcmc)

summary(model_mcmc)

#posterior means across the 3 chains
psmean_reffects <- model_reffcts_ggs %>% group_by(Parameter)  %>% summarise(Mean=mean(value))

###Histogram of Random Effects
hist(psmean_reffects$Mean[-101],xlab="b",main="Histogram of Posterior Means of Normal Random Intercepts",probability=T,col="cyan4")




mmm2=model$BUGSoutput$summary[, c("mean", "2.5%",
                                  "97.5%")] 

hist(mmm2[str_detect(row.names(mmm2),"SUB")],
     xlab = "Subject", main=
       "Distribution plot for Random effect Subject")

hist(mmm2[str_detect(row.names(mmm2),"FAC_random")],
     xlab = "Facilities", main=
       "Distribution plot for Random effect Facilities")

m=mmm2[str_detect(row.names(mmm2),"SUB"),]
n=mmm2[str_detect(row.names(mmm2),"FAC_random"),]


h=m[,1]
g=n[,1]


hh=data.frame(h)
mn=data.frame(g)

qqnorm(m)
qqline(m, col = "red", lwd = 2)


qqnorm(n)
qqline(n, col = "red", lwd = 2)


qplot(sample= h, data=hh)
qplot(sample= g, data=mn)

qplot(g, data = mn, geom = "histogram", bins=30,
      main = "Histogram for Random effect Facilities",
      fill=I("blue"))+
  theme_bw()

qplot(h, data = hh, geom = "histogram", bins=30,
      main = "Histogram for Random effect Subjects",
      fill=I("pink") ,
      xlab ="Subjects")+
  theme_bw()

ggplot(hh, aes(sample=h))+stat_qq() + theme_bw()
ggplot(mn, aes(sample=g))+stat_qq() + theme_bw()

ggplot(hh) + geom_qq(aes(sample = h))+
  geom_abline(intercept = 20, slope = 7,
              color = "red", size = 1.3, alpha = 0.5)

#p = ggplot(mn, aes(x=g))
#hh=p+geom_bar



model_mcmc= as.mcmc(model)
hist(model_mcmc)

ms=ggs(md)

ggs_traceplot(ms)

ggs_autocorrelation(ms)+theme_classic()


gelman.plot(md)


R2jags::traceplot(model, mfrow = c(2, 2), 
                  varname = c("beta","sigma.sub", "sigma.fac"))








## PPC

model12<-"model{

for(i in 1:Nobs) {			
  for(j in 1:nrep[i]){		
    y[i, j] ~ dbern(pi[i,j])
    y.rep[i,j] ~ dbern(pi.rep[i,j])
    logit(pi[i, j]) =  FAC_random[i] + SUB[i,j]
    SUB[i,j] ~ dnorm(beta[1] + beta[2]*TRAIN[i,j],
    sigma_sub_inv)
    
    #distribution of future observation of the probability
    
    logit(pi.rep[i,j])=FAC_random.rep[i] + SUB.rep[i,j]
    SUB.rep[i,j] ~ dnorm(beta[1] + beta[2]*TRAIN[i,j],
    sigma_sub_inv)
    
    #distribution of future observation y (pi.rep)
    
        #pi.rep[i,j] ~ dbern(pi[i,j])
  }		
  
  FAC_random[i] ~ dnorm(0, sigma_fac_inv)		
  FAC_random.rep[i] ~ dnorm(0, sigma_fac_inv)		
  
}			

#Hyperpriors for SUB Random Effect

sigma.sub ~ dunif(0, 20)
sigma_sub_inv <- 1/pow(sigma.sub ,2)

#Hyperpriors for FAC Random Effect

sigma.fac ~ dunif(0, 20)
sigma_fac_inv <- 1/pow(sigma.fac ,2)


#Priors for Fixed Effects

for (p in 1:2){

beta[p] ~ dnorm(0,0.0001)
 }


for(i in 1:Nobs){
#nrep[i]
## minimum and maximum across rows
  		
  a[i]=min(pi[i,1:nrep[i]])
  b[i]=min(pi.rep[i,1:nrep[i]])
  
  c[i]=max(pi[i,1:nrep[i]])
  d[i]=max(pi.rep[i,1:nrep[i]])
  
  #Calculate the average probability per patient
  
  pi_mean[i]= mean(pi[i,1:nrep[i]])
  pi.rep_mean[i]=mean(pi.rep[i,1:nrep[i]])
  

  
}

##PPCs for whole dataset

pi_min=min(a[])
pi_max=max(c[])
pi.rep_min = min(b[])
pi.rep_max= max(d[])
  
m1=mean(pi_mean[])
m2=mean(pi.rep_mean[])

### test for kolmogorov

rank_pi=sort(pi_mean[])
rank_pi.rep=sort(pi.rep_mean[])
for(k in 1:Nobs){
F_pi[k] = phi(rank_pi[k])
F_pi.rep[k] = phi(rank_pi.rep[k])
F_diff[k]=max(F_pi[k] - (k-1)/Nobs, k/Nobs-
F_pi[k])
F_diff.rep[k]=max(F_pi.rep[k]-(k-1)/Nobs,
k/Nobs-F_pi.rep[k])
}
ks_pi=max(F_diff)
ks_pi.rep=max(F_diff.rep)

##PPCs

        #test1 min and max
        #y_min <- min(a[])
        #y.rep_min <- min(b[])

        #y_max <- max(c[])
        #y.rep_max <- max(d[])
        
        
        #tests
        #t_min <- step(y.rep_min - y_min)
        #t_max <- step(y.rep_max - y_max)
        

##tests2 Skewness and Kurtosis
        
        #for(i in 1:Nobs) {			
  #for(j in 1:nrep2[i]){
            #skewness
    #sky[i,j] <- pow((y[i,j] - pi[i,j])/(pi[i,j]*(1-pi[i,j])),
    #3)
  #sky.rep[i,j] <- pow((y.rep2[i,j] - pi[i,j])/(pi[i,j]*(1-pi[i,j]))
         #   ,3)
  #kurtosis
# kky[i,j] <- pow((y[i,j] - pi[i,j])/(pi[i,j]*(1-pi[i,j]))
            #,4)
# kky.rep[i,j] <- pow((y.rep2[i,j] - pi[i,j])/(pi[i,j]*(1-pi[i,j]))
# ,4)


 ## }


#sky2[i]=mean(sky[i,1:nrep2[i]])
  #sky.rep2[i]=mean(sky.rep[i,1:nrep2[i]])
  #kky2[i]=mean(kky[i,1:nrep2[i]])
  #kky.rep2[i]=mean(kky.rep[i,1:nrep2[i]])
#}

#skew_y <- mean(sky2[])
        #skew_y.rep <- mean(sky.rep2[])
        #skew_test <- step(skew_y.rep-skew_y)
      
        #kurt_y <- mean(kky2[])-3
        #kurt_y.rep <- mean(kky.rep2[])-3
        #kurt_test <- step(kurt_y.rep - kurt_y)


###tests
      ppc_test[1] <- step(pi.rep_min - pi_min)
      ppc_test[2] <- step(pi.rep_max - pi_max)
      ppc_test[3] <- step(ks_pi.rep - ks_pi)
      ppc_test[4] <- step(m2 - m1)
      
      
###PPC measures
      ppc_measure[1] <- pi_min
      ppc_measure[2] <- pi.rep_min
      ppc_measure[3] <- pi_max
      ppc_measure[4] <- pi.rep_max
      ppc_measure[5] <- ks_pi
      ppc_measure[6] <- ks_pi.rep
      ppc_measure[7] <- m1
      ppc_measure[8] <- m2 

}"




dat <- list(y=s[,c("GUM_ND_1","GUM_ND_2", "GUM_ND_3", "GUM_ND_4")]
            ,Nobs=nrow(s), nrep=nrep,
            TRAIN=s[,c("TRAIN_1" ,
                       "TRAIN_2" , "TRAIN_3",  "TRAIN_4")])


inits = list(
  list(beta=c(2.0502,-0.3638),sigma.sub=0.5,sigma.fac=1.5),
  list(beta=c(2.0502,-0.3638)*2,sigma.sub=0.5*2,sigma.fac=1.5*2),
  list(beta=c(2.0502,-0.3638)/2,sigma.sub=0.5/2, sigma.fac=1.5/2)
)

model_14 <- jags(model.file = textConnection
                 (model12),data = dat,n.chains=3, inits = inits,
                 parameters.to.save=c("ppc_test",
                                      "ppc_measure"), 
                 n.burnin = 40000, 
                 n.iter = 50000, n.thin=10, DIC=T)

mmm2=model_14$BUGSoutput$summary[, c("mean", "2.5%",
                                     "97.5%")]

a = mmm2[str_detect(row.names(mmm2),"pi_mean"),]
b = mmm2[str_detect(row.names(mmm2),"pi.rep_mean"),]
v=a[,1]
w=b[,1]
sv=data.frame(v)
sw=data.frame(w)
data.frame(v=a[,1], w=b[,1])
ab=data.frame(pi=a[,1], pi.rep=b[,1])


library(ggplot2)

p= ggplot(ab, aes(x=pi))
p + geom_histogram(bins=15, color ="black",
                   fill="gray") +
  geom_histogram(aes(x=pi.rep), bins=15, color="black",
                 fill="green") + theme_bw()

ab=data.frame(pi_posterior=a[,1],
              pi_posterior_predictive=b[,1])

data.frame(pi_posterior=a[,1],
           pi_posterior_predictive=b[,1]) %>%
  pivot_longer(cols=everything(),
               names_to = "Distribution", values_to = "Posterior")
ggplot(ab, aes(x=Posterior, fill =Distribution))+
  geom_histogram(bins = 6, color= "black") +
  theme_bw()  +
  theme(axis.title = element_text(size =15),
        legend.position = "top")




save(model_14, file = "model_14s.RData")
save.image("model_14.RData")



##convert to MCMC

model_14_ppc_mcmc= as.mcmc(model_14)

summary(model_14_ppc_mcmc)

### Checking convergence of the posterior 
#discrepancies(normal random effects okay)

model_14_ppc_ggs = ggs(model_14_ppc_mcmc)
p_D <- model_14_ppc_ggs %>% group_by(Parameter) %>% 
  summarise(Mean=mean(value))



#as a ggs object
#ms2=ggs(md2)

#traceplot
#ggs_traceplot(ms2)

#autocorrelation plot
#autocorr.plot(md2)

#BGR diagnostic version 2
#gelman.plot(md2)


#R2jags::traceplot(model_2, mfrow = c(2, 2), 
#varname = c("beta","sigma.sub",
# "sigma.fac"))

#autocorr.plot(md2, c("beta"))




##### using t-distribution for the random effects


model22<-"model{

for(i in 1:Nobs) {			
  for(j in 1:nrep[i]){		
    y[i, j] ~ dbern(pi[i,j])
    y.rep[i,j] ~ dbern(pi.rep[i,j])
    logit(pi[i, j]) =  FAC_random[i] + SUB[i,j]
    SUB[i,j] ~ dt(beta[1] + beta[2]*TRAIN[i,j],
    sigma_sub_inv, df_sub)
    
    
    #distribution of future observation of the probability
    
    logit(pi.rep[i,j])=FAC_random.rep[i] + SUB.rep[i,j]
    SUB.rep[i,j] ~ dt(beta[1] + beta[2]*TRAIN[i,j],
    sigma_sub_inv, df_sub)
    
  }		
  
  FAC_random[i] ~ dt(0, sigma_fac_inv, df_fac)		
  FAC_random.rep[i] ~ dt(0, sigma_fac_inv, df_fac)

}			

#Hyperpriors for SUB Random Effect

sigma.sub ~ dunif(0, 20)
sigma_sub_inv <- 1/pow(sigma.sub ,2)
df_sub ~ dgamma(0.001,0.001)

#Hyperpriors for FAC Random Effect

sigma.fac ~ dunif(0, 20)
sigma_fac_inv <- 1/pow(sigma.fac ,2)
df_fac ~ dgamma(0.001,0.001)

#Priors for Fixed Effects

for (p in 1:2){

beta[p] ~ dnorm(0,0.0001)
}
 
for(i in 1:Nobs){

e[i]=min(SUB[i,1:nrep[i]])
f[i]=min(SUB.rep[i,1:nrep[i]])
  
g[i]=max(SUB[i,1:nrep[i]])
h[i]=max(SUB.rep[i,1:nrep[i]]) 


SUB_mean[i]= mean(SUB[i,1:nrep[i]])
SUB.rep_mean[i]=mean(SUB.rep[i,1:nrep[i]])
}



#FAC_random_min=min(FAC_random[]) 
#FAC_random.rep_min=min(FAC_random.rep[])
#FAC_random_max=max(FAC_random[])
#FAC_random.rep_max=max(FAC_random.rep[])

SUB_min=min(e[])
SUB.rep_min=min(f[])
SUB_max=max(g[])
SUB.rep_max=max(h[])



#for(i in 1:Nobs){
#nrep[i]
## minimum and maximum across rows
  		
#a[i]=min(pi[i,1:nrep[i]])
#b[i]=min(pi.rep[i,1:nrep[i]])
  
#c[i]=max(pi[i,1:nrep[i]])
#d[i]=max(pi.rep[i,1:nrep[i]])
  
#Calculate the average probability per patient
  
#pi_mean[i]= mean(pi[i,1:nrep[i]])
#pi.rep_mean[i]=mean(pi.rep[i,1:nrep[i]])

#}

##PPCs for whole dataset

#pi_min=min(a[])
#pi_max=max(c[])
#pi.rep_min = min(b[])
#pi.rep_max= max(d[])
  
#m1=mean(pi_mean[])
#m2=mean(pi.rep_mean[])



### test for kolmogorov

#rank_FAC_random=sort(FAC_random[])
#rank_FAC_random.rep=sort(FAC_random.rep[])
#for(k in 1:Nobs){
#F_FAC_random[k] = phi(rank_FAC_random[k])
#F_FAC_random.rep[k] = phi(rank_FAC_random.rep[k])
#F_diff1[k]=max(F_FAC_random[k] - (k-1)/Nobs, k/Nobs-
#F_FAC_random[k])
#F_diff1.rep[k]=max(F_FAC_random.rep[k]-(k-1)/Nobs,
#k/Nobs-F_FAC_random.rep[k])
#}
#ks_FAC_random=max(F_diff1)
#ks_FAC_random.rep=max(F_diff1.rep)


rank_SUB=sort(SUB_mean[])
rank_SUB.rep=sort(SUB.rep_mean[])
for(k in 1:Nobs){
F_SUB[k] = phi(rank_SUB[k])
F_SUB.rep[k] = phi(rank_SUB.rep[k])
F_diff[k]=max(F_SUB[k] - (k-1)/Nobs, k/Nobs-
F_SUB[k])
F_diff.rep[k]=max(F_SUB.rep[k]-(k-1)/Nobs,
k/Nobs-F_SUB.rep[k])
}
ks_SUB=max(F_diff)
ks_SUB.rep=max(F_diff.rep)


###tests

ppc_test[1] <- step(FAC_random.rep_min - FAC_random_min)
ppc_test[2] <- step(FAC_random.rep_max - FAC_random_max)
ppc_test[3] <- step(SUB.rep_min - SUB_min)
ppc_test[4] <- step(SUB.rep_max - SUB.rep_max)
ppc_test[5] <- step(ks_FAC_random.rep - ks_FAC_random)
ppc_test[6] <- step(ks_SUB.rep - ks_SUB)


#ppc_test[4] <- step(m2 - m1)
      
      
###PPC measures

#ppc_measure[01] <- FAC_random_min
#ppc_measure[02] <- FAC_random.rep_min
#ppc_measure[03] <- FAC_random_max
#ppc_measure[04] <- FAC_random.rep_max
ppc_measure[05] <-SUB_min
ppc_measure[06] <- SUB.rep_min
ppc_measure[07] <- SUB_max
ppc_measure[08] <- SUB.rep_max
#ppc_measure[09] <- ks_FAC_random
#ppc_measure[10] <- ks_FAC_random.rep
ppc_measure[11] <- ks_SUB
ppc_measure[12] <- ks_SUB.rep

#ppc_measure[7] <- m1
#ppc_measure[8] <- m2 

}"



dat <- list(y=s[,c("GUM_ND_1","GUM_ND_2", "GUM_ND_3", "GUM_ND_4")]
            ,Nobs=nrow(s), nrep=nrep,
            TRAIN=s[,c("TRAIN_1" ,
                       "TRAIN_2" , "TRAIN_3",  "TRAIN_4")])


inits = list(
  list(beta=c(2.0502,-0.3638),sigma.sub=0.5,sigma.fac=1.5),
  list(beta=c(2.0502,-0.3638)*2,sigma.sub=0.5*2,sigma.fac=1.5*2),
  list(beta=c(2.0502,-0.3638)/2,sigma.sub=0.5/2, sigma.fac=1.5/2)
)

model_22 <- jags(model.file = textConnection
                 (model22),data = dat,n.chains=3, inits = inits,
                 parameters.to.save=c("ppc_test",
                                      "ppc_measure"), 
                 n.burnin = 40000, 
                 n.iter = 50000, n.thin=10, DIC=T)




model_22_t_mcmc= as.mcmc(model_22)
summary(model_22_t_mcmc)



ms2=ggs(md2)


summary(md2)

ms2$DIC 

model_22$BUGSoutput$DIC

model_22$BUGSoutput$pD





##### using normal distribution for the random effects


model23<-"model{

for(i in 1:Nobs) {			
  for(j in 1:nrep[i]){		
    y[i, j] ~ dbern(pi[i,j])
    y.rep[i,j] ~ dbern(pi.rep[i,j])
    logit(pi[i, j]) =  FAC_random[i] + SUB[i,j]
    SUB[i,j] ~ dnorm(beta[1] + beta[2]*TRAIN[i,j],
    sigma_sub_inv)
    
    #distribution of future observation of the probability
    
    logit(pi.rep[i,j])=FAC_random.rep[i] + SUB.rep[i,j]
    SUB.rep[i,j] ~ dnorm(beta[1] + beta[2]*TRAIN[i,j],
    sigma_sub_inv)
    
  }		
  
  FAC_random[i] ~ dnorm(0, sigma_fac_inv)		
  FAC_random.rep[i] ~ dnorm(0, sigma_fac_inv)
}			

#Hyperpriors for SUB Random Effect

sigma.sub ~ dunif(0, 20)
sigma_sub_inv <- 1/pow(sigma.sub ,2)

#Hyperpriors for FAC Random Effect

sigma.fac ~ dunif(0, 20)
sigma_fac_inv <- 1/pow(sigma.fac ,2)



#Priors for Fixed Effects

for (p in 1:2){

beta[p] ~ dnorm(0,0.0001)
}
 
for(i in 1:Nobs){

e[i]=min(SUB[i,1:nrep[i]])
f[i]=min(SUB.rep[i,1:nrep[i]])
  
g[i]=max(SUB[i,1:nrep[i]])
h[i]=max(SUB.rep[i,1:nrep[i]]) 


SUB_mean[i]= mean(SUB[i,1:nrep[i]])
SUB.rep_mean[i]=mean(SUB.rep[i,1:nrep[i]])
}



FAC_random_min=min(FAC_random[]) 
FAC_random.rep_min=min(FAC_random.rep[])
FAC_random_max=max(FAC_random[])
FAC_random.rep_max=max(FAC_random.rep[])

SUB_min=min(e[])
SUB.rep_min=min(f[])
SUB_max=max(g[])
SUB.rep_max=max(h[])



#for(i in 1:Nobs){
#nrep[i]
## minimum and maximum across rows
  		
#a[i]=min(pi[i,1:nrep[i]])
#b[i]=min(pi.rep[i,1:nrep[i]])
  
#c[i]=max(pi[i,1:nrep[i]])
#d[i]=max(pi.rep[i,1:nrep[i]])
  
#Calculate the average probability per patient
  
#pi_mean[i]= mean(pi[i,1:nrep[i]])
#pi.rep_mean[i]=mean(pi.rep[i,1:nrep[i]])

#}

##PPCs for whole dataset

#pi_min=min(a[])
#pi_max=max(c[])
#pi.rep_min = min(b[])
#pi.rep_max= max(d[])
  
#m1=mean(pi_mean[])
#m2=mean(pi.rep_mean[])



### test for kolmogorov

rank_FAC_random=sort(FAC_random[])
rank_FAC_random.rep=sort(FAC_random.rep[])
for(k in 1:Nobs){
F_FAC_random[k] = phi(rank_FAC_random[k])
F_FAC_random.rep[k] = phi(rank_FAC_random.rep[k])
F_diff1[k]=max(F_FAC_random[k] - (k-1)/Nobs, k/Nobs-
F_FAC_random[k])
F_diff1.rep[k]=max(F_FAC_random.rep[k]-(k-1)/Nobs,
k/Nobs-F_FAC_random.rep[k])
}
ks_FAC_random=max(F_diff1)
ks_FAC_random.rep=max(F_diff1.rep)


rank_SUB=sort(SUB_mean[])
rank_SUB.rep=sort(SUB.rep_mean[])
for(k in 1:Nobs){
F_SUB[k] = phi(rank_SUB[k])
F_SUB.rep[k] = phi(rank_SUB.rep[k])
F_diff[k]=max(F_SUB[k] - (k-1)/Nobs, k/Nobs-
F_SUB[k])
F_diff.rep[k]=max(F_SUB.rep[k]-(k-1)/Nobs,
k/Nobs-F_SUB.rep[k])
}
ks_SUB=max(F_diff)
ks_SUB.rep=max(F_diff.rep)


###tests

ppc_test[1] <- step(FAC_random.rep_min - FAC_random_min)
ppc_test[2] <- step(FAC_random.rep_max - FAC_random_max)
ppc_test[3] <- step(SUB.rep_min - SUB_min)
ppc_test[4] <- step(SUB.rep_max - SUB.rep_max)
ppc_test[5] <- step(ks_FAC_random.rep - ks_FAC_random)
ppc_test[6] <- step(ks_SUB.rep - ks_SUB)


#ppc_test[4] <- step(m2 - m1)
      
      
###PPC measures

ppc_measure[1] <- FAC_random_min
ppc_measure[2] <- FAC_random.rep_min
ppc_measure[3] <- FAC_random_max
ppc_measure[4] <- FAC_random.rep_max
ppc_measure[5] <-SUB_min
ppc_measure[6] <- SUB.rep_min
ppc_measure[7] <- SUB_max
ppc_measure[8] <- SUB.rep_max
ppc_measure[9] <- ks_FAC_random
ppc_measure[10] <- ks_FAC_random.rep
ppc_measure[11] <- ks_SUB
ppc_measure[12] <- ks_SUB.rep

#ppc_measure[7] <- m1
#ppc_measure[8] <- m2 

}"



dat <- list(y=s[,c("GUM_ND_1","GUM_ND_2", "GUM_ND_3", "GUM_ND_4")]
            ,Nobs=nrow(s), nrep=nrep,
            TRAIN=s[,c("TRAIN_1" ,
                       "TRAIN_2" , "TRAIN_3",  "TRAIN_4")])


inits = list(
  list(beta=c(2.0502,-0.3638),sigma.sub=0.5,sigma.fac=1.5),
  list(beta=c(2.0502,-0.3638)*2,sigma.sub=0.5*2,sigma.fac=1.5*2),
  list(beta=c(2.0502,-0.3638)/2,sigma.sub=0.5/2, sigma.fac=1.5/2)
)

model_23 <- jags(model.file = textConnection
                 (model23),data = dat,n.chains=3, inits = inits,
                 parameters.to.save=c("ppc_test",
                                      "ppc_measure"), 
                 n.burnin = 40000, 
                 n.iter = 50000, n.thin=10, DIC=T)









### Model without random effect for care facility


model23<-"model{

for(i in 1:Nobs) {			
  for(j in 1:nrep[i]){		
    y[i, j] ~ dbern(pi[i,j])	
    logit(pi[i, j]) =   SUB[i,j]
    SUB[i,j] ~ dnorm(beta[1] + beta[2]*TRAIN[i,j],
    sigma_sub_inv)
    
  }		
  
  
}			

#Hyperpriors for SUB Random Effect

sigma.sub ~ dunif(0, 20)
sigma_sub_inv <- 1/pow(sigma.sub ,2)

#Priors for Fixed Effects

for (p in 1:2){

beta[p] ~ dnorm(0,0.0001)
 }

}"



dat <- list(y=s[,c("GUM_ND_1","GUM_ND_2", "GUM_ND_3", "GUM_ND_4")]
            ,Nobs=nrow(s), nrep=nrep,
            TRAIN=s[,c("TRAIN_1" ,
                       "TRAIN_2" , "TRAIN_3",  "TRAIN_4")])

inits = list(
  list(beta=c(2.0502,-0.3638),sigma.sub=0.5,sigma.fac=1.5),
  list(beta=c(2.0502,-0.3638)*2,sigma.sub=0.5*2,sigma.fac=1.5*2),
  list(beta=c(2.0502,-0.3638)/2,sigma.sub=0.5/2, sigma.fac=1.5/2)
)


model_23 <- jags(model.file = textConnection
                 (model23),data = dat,n.chains=3, inits = inits,
                 parameters.to.save=c("beta","sigma.sub"),
                 n.burnin = 40000, n.iter = 50000, n.thin=1)


model_23$BUGSoutput$DIC

model_23$BUGSoutput$pD


#as an mcmc object
md2= as.mcmc(model_2)




#as a ggs object
ms2=ggs(md2)

#traceplot
ggs_traceplot(ms2)

#autocorrelation plot
autocorr.plot(md2)

#BGR diagnostic version 2
gelman.plot(md2)


R2jags::traceplot(model_2, mfrow = c(2, 2), 
                  varname = c("beta","sigma.sub",
                              "sigma.fac"))


#autocorr.plot(md2, c("beta"))











