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

BDA1$gender_care[BDA1$gender_care=="Female"] <- "1"
BDA1$gender_care[BDA1$gender_care=="Male"] <- "0"

colnames(BDA1)<-c("FAC","SUB","SEX","AGE","EASY",
                  "TRAIN","EDU","SEX_ND","GUM_ND","GUM_D")


v=table(BDA1$SUB)



##### Part 2 #####



measurement=unlist(lapply(v, function(i) 1:i))

BDA1$Measurement = unlist(lapply(v, function(i) 1:i))


BDA2NH<- BDA1[BDA1$GUM_D == 0, ]


BDA2NH$GUM_ND=as.numeric(BDA2NH$GUM_ND)
BDA2NH$TRAIN=as.numeric(BDA2NH$TRAIN)
BDA2NH$SEX=as.numeric(BDA2NH$SEX)
BDA2NH$AGE=as.numeric(BDA2NH$AGE)
BDA2NH$EASY=as.numeric(BDA2NH$EASY)
BDA2NH$EDU=as.numeric(BDA2NH$EDU)
BDA2NH$SEX_ND=as.numeric(BDA2NH$SEX_ND)


nrep2=table(BDA2NH$SUB)

library(data.table)

s2=dcast(setDT(BDA2NH) , FAC+SUB+SEX+AGE+EASY ~ Measurement, 
         value.var = c("TRAIN", "EDU","SEX_ND", "GUM_ND"
         ))

nfac2 =table(BDA2NH$FAC)


glm(GUM_ND~SEX+AGE+EASY+TRAIN+ EDU+SEX_ND, 
    data=BDA2NH,family = binomial(link = logit))




model2<-"model{

for(i in 1:Nobs) {			
  for(j in 1:nrep2[i]){		
    y[i, j] ~ dbern(pi[i,j])	
    logit(pi[i, j]) =  FAC_random[i] + SUB[i,j]
    SUB[i,j] ~ dnorm(beta[1] + beta[2]*SEX[i] + beta[3]*AGE[i]
    + beta[4]*EASY[i]+ beta[5]*TRAIN[i,j] + 
    beta[6]*EDU[i,j]+ beta[7]*SEX_ND[i,j],
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

for (p in 1:7){

beta[p] ~ dnorm(0,0.0001)
 }

}"



library(rjags)

library(R2jags)

dat <- list(y=s2[,c("GUM_ND_1","GUM_ND_2", "GUM_ND_3", "GUM_ND_4")]
            ,Nobs=nrow(s2), nrep2=nrep2,
            TRAIN=s2[,c("TRAIN_1" ,
                        "TRAIN_2" , "TRAIN_3",  "TRAIN_4")],
            SEX=s2$SEX, AGE=s2$AGE, EASY=s2$EASY, EDU=s2[,
                                                         c("EDU_1","EDU_2","EDU_3" ,"EDU_4")], 
            SEX_ND=s2[, c("SEX_ND_1", "SEX_ND_2", "SEX_ND_3", "SEX_ND_4")])

inits = list(
  list(beta=c(1.28486,0.69710,0.01199,0.17843,-0.48455,
              -0.47600,1.83693),
       sigma.sub=0.5,sigma.fac=1.5),
  list(beta=c(1.28486,0.69710,0.01199,0.17843,-0.48455,
              -0.47600,1.83693)*2,
       sigma.sub=0.5*2,sigma.fac=1.5*2),
  list(beta=c(1.28486,0.69710,0.01199,0.17843,-0.48455,
              -0.47600,1.83693)/2,
       sigma.sub=0.5/2, sigma.fac=1.5/2)
)


model_2 <- jags(model.file = textConnection
                (model2),data = dat,n.chains=3, inits = inits,
                parameters.to.save=c("beta","sigma.sub",
                                     "sigma.fac"), n.burnin = 40000, 
                n.iter = 50000, n.thin=1)


library(mcmcplots)
library(ggmcmc)

#as an mcmc object
model_2_mcmc= as.mcmc(model_2)
model_2_ggs <- ggs(model_2_mcmc)



mmm4=model_2$BUGSoutput$summary[, c("mean", "2.5%",
                                    "97.5%")] 
model_2$BUGSoutput$DIC
model_2$BUGSoutput$pD


summary(model_2_mcmc)




model_16 <- jags(model.file = textConnection
                 (model2),data = dat,n.chains=3, inits = inits,
                 parameters.to.save=c("SUB","FAC_random"), n.burnin = 40000, 
                 n.iter = 50000, n.thin=1)


#as an mcmc object
model_random_mcmc= as.mcmc(model_2_random)
model_reffcts_ggs <- ggs(model_random_mcmc)


#posterior means across the 3 chains

psmean_random <- model_reffcts_ggs %>% 
  group_by(Parameter)  %>% summarise(Mean=mean(value))

###Histogram of Random Effects
hist(psmean_random$Mean[-101],xlab="b",
     main="Histogram of Posterior Means of Normal Random Intercepts",
     probability=T,col="cyan4")





mmm2=model_2_random$BUGSoutput$summary[, c("mean", "2.5%",
                                           "97.5%")] 

hist(mmm2[str_detect(row.names(mmm2),"SUB")],
     xlab = "Subject", main=
       "Distribution plot for Random effect Subject")

hist(mmm2[str_detect(row.names(mmm2),"FAC_random")],
     xlab = "Facilities", main=
       "Distribution plot for Random effect Facilities")

m=mmm2[str_detect(row.names(mmm2),"SUB"),]
n=mmm2[str_detect(row.names(mmm2),"FAC_random"),]


subject=m[,1]
facilities=n[,1]


hh=data.frame(subject)
mn=data.frame(facilities)

qqnorm(m)
qqline(m, col = "red", lwd = 2)


qqnorm(n)
qqline(n, col = "red", lwd = 2)


qplot(sample= subject, data=hh)
qplot(sample= facilities, data=mn)

qplot(subject, data = hh, geom = "histogram", bins=30,
      main = "Histogram for Random effect Facilities",
      fill=I("blue"))+
  theme_bw()

qplot(facilities, data = mn, geom = "histogram", bins=30,
      main = "Histogram for Random effect Subjects",
      fill=I("pink") ,
      xlab ="Subjects")+
  theme_bw()

ggplot(hh, aes(x=subject))+
  geom_histogram(color="darkblue", fill="lightblue")+
  theme_bw()


ggplot(mn, aes(x=facilities))+
  geom_histogram(color="darkblue", fill="pink")+
  theme_bw()




ggplot(hh, aes(sample=subject))+stat_qq()+ theme_bw()
ggplot(mn, aes(sample=facilities))+stat_qq() + theme_bw()


ggplot(hh) + geom_qq(aes(sample = subject))+
  geom_abline(intercept = 10, slope = 10,
              color = "red", size = 1.3, alpha = 0.5)+ theme_bw()


ggplot(mn) + geom_qq(aes(sample = facilities))+
  geom_abline(intercept = -09, slope = 13,
              color = "red", size = 1.3, alpha = 0.5)+ theme_bw()



library(ggplot2)

#as a ggs object
ms2=ggs(md2)

#traceplot
ggs_traceplot(ms2)+theme_classic()

ggs_autocorrelation(ms2)+theme_classic()

#BGR diagnostic version 2
gelman.plot(md2)


R2jags::traceplot(model_2, mfrow = c(2, 2), 
                  varname = c("beta","sigma.sub",
                              "sigma.fac"))


#autocorr.plot(md2, c("beta"))



#without random effect for facilities


model14<-"model{

for(i in 1:Nobs) {			
  for(j in 1:nrep2[i]){		
    y[i, j] ~ dbern(pi[i,j])	
    logit(pi[i, j]) =   SUB[i,j]
    SUB[i,j] ~ dnorm(beta[1] + beta[2]*SEX[i] + beta[3]*AGE[i]
    + beta[4]*EASY[i]+ beta[5]*TRAIN[i,j] + 
    beta[6]*EDU[i,j]+ beta[7]*SEX_ND[i,j],
    sigma_sub_inv)
    
  }		
  
  #FAC_random[i] ~ dnorm(0, sigma_fac_inv)		
  		
  
}			

#Hyperpriors for SUB Random Effect

sigma.sub ~ dunif(0, 20)
sigma_sub_inv <- 1/pow(sigma.sub ,2)

#Hyperpriors for FAC Random Effect

#sigma.fac ~ dunif(0, 20)
#sigma_fac_inv <- 1/pow(sigma.fac ,2)


#Priors for Fixed Effects

for (p in 1:7){

beta[p] ~ dnorm(0,0.0001)
 }

}"



dat <- list(y=s2[,c("GUM_ND_1","GUM_ND_2", "GUM_ND_3", "GUM_ND_4")]
            ,Nobs=nrow(s2), nrep2=nrep2,
            TRAIN=s2[,c("TRAIN_1" ,
                        "TRAIN_2" , "TRAIN_3",  "TRAIN_4")],
            SEX=s2$SEX, AGE=s2$AGE, EASY=s2$EASY, EDU=s2[,
                                                         c("EDU_1","EDU_2","EDU_3" ,"EDU_4")], 
            SEX_ND=s2[, c("SEX_ND_1", "SEX_ND_2", "SEX_ND_3", "SEX_ND_4")])

inits = list(
  list(beta=c(1.28486,0.69710,0.01199,0.17843,-0.48455,
              -0.47600,1.83693),
       sigma.sub=0.5,sigma.fac=1.5),
  list(beta=c(1.28486,0.69710,0.01199,0.17843,-0.48455,
              -0.47600,1.83693)*2,
       sigma.sub=0.5*2,sigma.fac=1.5*2),
  list(beta=c(1.28486,0.69710,0.01199,0.17843,-0.48455,
              -0.47600,1.83693)/2,
       sigma.sub=0.5/2, sigma.fac=1.5/2)
)


model_14 <- jags(model.file = textConnection
                 (model14),data = dat,n.chains=3, inits = inits,
                 parameters.to.save=c("beta","sigma.sub"), n.burnin = 40000, 
                 n.iter = 50000, n.thin=1)




#as an mcmc object
model_2_mcmc= as.mcmc(model_2)
model_2_ggs <- ggs(model_2_mcmc)



mmm5=model_14$BUGSoutput$summary[, c("mean", "2.5%",
                                     "97.5%")] 
model_14$BUGSoutput$DIC
model_14 $BUGSoutput$pD


summary(model_2_mcmc)




#t-distribution for the random effects


model13<-"model{

for(i in 1:Nobs) {			
  for(j in 1:nrep2[i]){		
    y[i, j] ~ dbern(pi[i,j])	
    logit(pi[i, j]) =  FAC_random[i] + SUB[i,j]
    SUB[i,j] ~ dt(beta[1] + beta[2]*SEX[i] + beta[3]*AGE[i]
    + beta[4]*EASY[i]+ beta[5]*TRAIN[i,j] + 
    beta[6]*EDU[i,j]+ beta[7]*SEX_ND[i,j],
    sigma_sub_inv, df_sub)
    
  }		
  
  FAC_random[i] ~ dt(0, sigma_fac_inv, df_fac)		
  		
  
}			

#Hyperpriors for SUB Random Effect

sigma.sub ~ dunif(0, 20)
sigma_sub_inv <- 1/pow(sigma.sub ,2)
df_sub ~ dgamma(0.01,0.01)


#Hyperpriors for FAC Random Effect

sigma.fac ~ dunif(0, 20)
sigma_fac_inv <- 1/pow(sigma.fac ,2)
df_fac ~ dgamma(0.01,0.01)


#Priors for Fixed Effects

for (p in 1:7){

beta[p] ~ dnorm(0,0.0001)
 }

}"




dat <- list(y=s2[,c("GUM_ND_1","GUM_ND_2", "GUM_ND_3", "GUM_ND_4")]
            ,Nobs=nrow(s2), nrep2=nrep2,
            TRAIN=s2[,c("TRAIN_1" ,
                        "TRAIN_2" , "TRAIN_3",  "TRAIN_4")],
            SEX=s2$SEX, AGE=s2$AGE, EASY=s2$EASY, EDU=s2[,
                                                         c("EDU_1","EDU_2","EDU_3" ,"EDU_4")], 
            SEX_ND=s2[, c("SEX_ND_1", "SEX_ND_2", "SEX_ND_3", "SEX_ND_4")])

inits = list(
  list(beta=c(1.28486,0.69710,0.01199,0.17843,-0.48455,
              -0.47600,1.83693),
       sigma.sub=0.5,sigma.fac=1.5),
  list(beta=c(1.28486,0.69710,0.01199,0.17843,-0.48455,
              -0.47600,1.83693)*2,
       sigma.sub=0.5*2,sigma.fac=1.5*2),
  list(beta=c(1.28486,0.69710,0.01199,0.17843,-0.48455,
              -0.47600,1.83693)/2,
       sigma.sub=0.5/2, sigma.fac=1.5/2)
)


model_13_random <- jags(model.file = textConnection
                        (model13),data = dat,n.chains=3, inits = inits,
                        parameters.to.save=c("SUB","FAC_random"), n.burnin = 40000, 
                        n.iter = 50000, n.thin=1, DIC=T)



#as an mcmc object

model13_random_mcmc= as.mcmc(model_13_random)
model13_reffcts_ggs <- ggs(model13_random_mcmc)




model_13_random$BUGSoutput$DIC
model_13_random$BUGSoutput$pD


mmm3=model_13_random$BUGSoutput$summary[, c("mean", "2.5%",
                                            "97.5%")] 

hist(mmm3[str_detect(row.names(mmm3),"SUB")],
     xlab = "Subject", main=
       "Distribution plot for Random effect Subject", col="#800000" )


hist(mmm3[str_detect(row.names(mmm3),"FAC_random")],
     xlab = "Facilities", main=
       "Distribution plot for Random effect Facilities", col="#808000")


j=mmm3[str_detect(row.names(mmm3),"SUB"),]
k=mmm3[str_detect(row.names(mmm3),"FAC_random"),]


subject=j[,1]
facilities=k[,1]


jj=data.frame(subject)
kk=data.frame(facilities)

qqnorm(j)
qqline(j, col = "red", lwd = 2)


qqnorm(k)
qqline(k, col = "red", lwd = 2)


qplot(sample= subject, data=hh)
qplot(sample= facilities, data=mn)

qplot(subject, data = jj, geom = "histogram", bins=30,
      main = "Histogram for Random effect Facilities",
      fill=I("blue"))+
  theme_bw()

qplot(facilities, data = kk, geom = "histogram", bins=30,
      main = "Histogram for Random effect Subjects",
      fill=I("pink") ,
      xlab ="Subjects")+
  theme_bw()


ggplot(jj, aes(x=subject))+
  geom_histogram(color="darkblue", fill="#800000")+
  theme_bw()


ggplot(kk, aes(x=facilities))+
  geom_histogram(color="darkblue", fill="#808000")+
  theme_bw()



ggplot(jj, aes(sample=subject))+stat_qq()+ theme_bw()
ggplot(kk, aes(sample=facilities))+stat_qq() + theme_bw()


ggplot(jj) + geom_qq(aes(sample = subject))+
  geom_abline(intercept = -300, slope =100 ,
              color = "red", size = 1.3, alpha = 0.5)+ theme_bw()


ggplot(kk) + geom_qq(aes(sample = facilities))+
  geom_abline(intercept = 10, slope = 10,
              color = "red", size = 1.3, alpha = 0.5)+ theme_bw()






## PPC

model21<-"model{

for(i in 1:Nobs) {			
  for(j in 1:nrep2[i]){		
    y[i, j] ~ dbern(pi[i,j])
y.rep[i,j] ~ dbern(pi.rep[i,j])	
    logit(pi[i, j]) =  FAC_random[i] + SUB[i,j]
    SUB[i,j] ~ dnorm(beta[1] + beta[2]*SEX[i] + beta[3]*AGE[i]
    + beta[4]*EASY[i]+ beta[5]*TRAIN[i,j] + 
    beta[6]*EDU[i,j]+ beta[7]*SEX_ND[i,j],
    sigma_sub_inv)
    

#distribution of future observation of the probability

logit(pi.rep[i,j])=FAC_random.rep[i] + SUB.rep[i,j]
SUB.rep[i,j] ~ dnorm(beta[1] + beta[2]*SEX[i] + beta[3]*AGE[i]
    + beta[4]*EASY[i]+ beta[5]*TRAIN[i,j] + 
    beta[6]*EDU[i,j]+ beta[7]*SEX_ND[i,j],
    sigma_sub_inv)

    
  #distribution of future observation y (y.rep)
    
        #y.rep2[i,j] ~ dbern(pi[i,j])
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

for (p in 1:7){

beta[p] ~ dnorm(0,0.0001)
 }


for(i in 1:Nobs){
#nrep[i]
## minimum and maximum across rows
  		
  a[i]=min(pi[i,1:nrep2[i]])
  b[i]=min(pi.rep[i,1:nrep2[i]])
  
  c[i]=max(pi[i,1:nrep2[i]])
  d[i]=max(pi.rep[i,1:nrep2[i]])

#Calculate the average probability per patient
  
  pi_mean[i]= mean(pi[i,1:nrep2[i]])
  pi.rep_mean[i]=mean(pi.rep[i,1:nrep2[i]])
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

##PPCs
        #test1 min and max
        #y_min <- min(a[])
        #y.rep_min <- min(b[])

        #y_max <- max(c[])
        #y.rep_max <- max(d[])
        
        
        #tests
        #t_min <- step(y.rep_min - y_min)
       # t_max <- step(y.rep_max - y_max)
        

##tests2 Skewness and Kurtosis
        
        #for(i in 1:Nobs) {			
  #for(j in 1:nrep2[i]){
            #skewness
    #sky[i,j] <- pow((y[i,j] - pi[i,j])/(pi[i,j]*(1-pi[i,j])),
    #3)
  #sky.rep[i,j] <- pow((y.rep2[i,j] - pi[i,j])/(pi[i,j]*(1-pi[i,j]))
          #  ,3)
  #kurtosis
 #kky[i,j] <- pow((y[i,j] - pi[i,j])/(pi[i,j]*(1-pi[i,j]))
           # ,4)
 #kky.rep[i,j] <- pow((y.rep2[i,j] - pi[i,j])/(pi[i,j]*(1-pi[i,j]))
 #,4)


 # }


#sky2[i]=mean(sky[i,1:nrep2[i]])
 # sky.rep2[i]=mean(sky.rep[i,1:nrep2[i]])
 # kky2[i]=mean(kky[i,1:nrep2[i]])
 # kky.rep2[i]=mean(kky.rep[i,1:nrep2[i]])
#}

#skew_y <- mean(sky2[])
       # skew_y.rep <- mean(sky.rep2[])
        #skew_test <- step(skew_y.rep-skew_y)
      
        #kurt_y <- mean(kky2[])-3
        #kurt_y.rep <- mean(kky.rep2[])-3
       # kurt_test <- step(kurt_y.rep - kurt_y)


###tests
     # ppc_test[1] <- t_min
     # ppc_test[2] <- t_max
     # ppc_test[3] <- skew_test
     # ppc_test[4] <- kurt_test
      
      
###PPC measures
     # ppc_measure[1] <- y_min
      #ppc_measure[2] <- y.rep_min
      #ppc_measure[3] <- y_max
     # ppc_measure[4] <- y.rep_max
     # ppc_measure[5] <- skew_y
     # ppc_measure[6] <- skew_y.rep
     # ppc_measure[7] <- kurt_y
     # ppc_measure[8] <- kurt_y.rep 

}"




dat <- list(y=s2[,c("GUM_ND_1","GUM_ND_2", "GUM_ND_3", "GUM_ND_4")]
            ,Nobs=nrow(s2), nrep2=nrep2,
            TRAIN=s2[,c("TRAIN_1" ,
                        "TRAIN_2" , "TRAIN_3",  "TRAIN_4")],
            SEX=s2$SEX, AGE=s2$AGE, EASY=s2$EASY, EDU=s2[,
                                                         c("EDU_1","EDU_2","EDU_3" ,"EDU_4")], 
            SEX_ND=s2[, c("SEX_ND_1", "SEX_ND_2", "SEX_ND_3", "SEX_ND_4")])

inits = list(
  list(beta=c(1.28486,0.69710,0.01199,0.17843,-0.48455,
              -0.47600,1.83693),
       sigma.sub=0.5,sigma.fac=1.5),
  list(beta=c(1.28486,0.69710,0.01199,0.17843,-0.48455,
              -0.47600,1.83693)*2,
       sigma.sub=0.5*2,sigma.fac=1.5*2),
  list(beta=c(1.28486,0.69710,0.01199,0.17843,-0.48455,
              -0.47600,1.83693)/2,
       sigma.sub=0.5/2, sigma.fac=1.5/2)
)


model_21 <- jags(model.file = textConnection
                 (model21),data = dat,n.chains=3, inits = inits,
                 parameters.to.save=c("ppc_test","ppc_measure"), 
                 n.burnin = 40000, 
                 n.iter = 50000, n.thin=1)



##convert to MCMC

model_21_ppc_mcmc= as.mcmc(model_21)

summary(model_21_ppc_mcmc)


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



##### using t-distribution for the random effects


model22<-"model{

for(i in 1:Nobs) {			
  for(j in 1:nrep2[i]){		
    y[i, j] ~ dbern(pi[i,j])	
y.rep[i,j] ~ dbern(pi.rep[i,j])
    logit(pi[i, j]) =  FAC_random[i] + SUB[i,j]
    SUB[i,j] ~ dt(beta[1] + beta[2]*SEX[i] + beta[3]*AGE[i]
    + beta[4]*EASY[i]+ beta[5]*TRAIN[i,j] + 
    beta[6]*EDU[i,j]+ beta[7]*SEX_ND[i,j],
    sigma_sub_inv, df_sub)
    


#distribution of future observation of the probability
    
    logit(pi.rep[i,j])=FAC_random.rep[i] + SUB.rep[i,j]

SUB.rep[i,j] ~ dt(beta[1] + beta[2]*SEX[i] + beta[3]*AGE[i]
    + beta[4]*EASY[i]+ beta[5]*TRAIN[i,j] + 
    beta[6]*EDU[i,j]+ beta[7]*SEX_ND[i,j],
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

for (p in 1:7){

beta[p] ~ dnorm(0,0.0001)
 }


for(i in 1:Nobs){

e[i]=min(SUB[i,1:nrep2[i]])
f[i]=min(SUB.rep[i,1:nrep2[i]])
  
g[i]=max(SUB[i,1:nrep2[i]])
h[i]=max(SUB.rep[i,1:nrep2[i]]) 


SUB_mean[i]= mean(SUB[i,1:nrep2[i]])
SUB.rep_mean[i]=mean(SUB.rep[i,1:nrep2[i]])
}



SUB_min=min(e[])
SUB.rep_min=min(f[])
SUB_max=max(g[])
SUB.rep_max=max(h[])

#FAC_random_min=min(FAC_random[]) 
#FAC_random.rep_min=min(FAC_random.rep[])
#FAC_random_max=max(FAC_random[])
#FAC_random.rep_max=max(FAC_random.rep[])


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

#ppc_test[1] <- step(FAC_random.rep_min - FAC_random_min)
#ppc_test[2] <- step(FAC_random.rep_max - FAC_random_max)
ppc_test[3]<-(SUB.rep_min - SUB_min)
ppc_test[4]<-(SUB.rep_max - SUB_max)
#ppc_test[5] <- step(ks_FAC_random.rep - ks_FAC_random)
ppc_test[6] <- step(ks_SUB.rep - ks_SUB)


###PPC measures

#ppc_measure[01] <- FAC_random_min
#ppc_measure[02] <- FAC_random.rep_min
#ppc_measure[03] <- FAC_random_max
#ppc_measure[04] <- FAC_random.rep_max
ppc_measure[05]<-SUB_min
ppc_measure[06] <- SUB.rep_min
ppc_measure[07] <- SUB_max
ppc_measure[08] <- SUB.rep_max
#ppc_measure[09] <- ks_FAC_random
#ppc_measure[10] <- ks_FAC_random.rep
ppc_measure[11] <- ks_SUB
ppc_measure[12] <- ks_SUB.rep

}"




dat <- list(y=s2[,c("GUM_ND_1","GUM_ND_2", "GUM_ND_3", "GUM_ND_4")]
            ,Nobs=nrow(s2), nrep2=nrep2,
            TRAIN=s2[,c("TRAIN_1" ,
                        "TRAIN_2" , "TRAIN_3",  "TRAIN_4")],
            SEX=s2$SEX, AGE=s2$AGE, EASY=s2$EASY, EDU=s2[,
                                                         c("EDU_1","EDU_2","EDU_3" ,"EDU_4")], 
            SEX_ND=s2[, c("SEX_ND_1", "SEX_ND_2", "SEX_ND_3", "SEX_ND_4")])

inits = list(
  list(beta=c(1.28486,0.69710,0.01199,0.17843,-0.48455,
              -0.47600,1.83693),
       sigma.sub=0.5,sigma.fac=1.5),
  list(beta=c(1.28486,0.69710,0.01199,0.17843,-0.48455,
              -0.47600,1.83693)*2,
       sigma.sub=0.5*2,sigma.fac=1.5*2),
  list(beta=c(1.28486,0.69710,0.01199,0.17843,-0.48455,
              -0.47600,1.83693)/2,
       sigma.sub=0.5/2, sigma.fac=1.5/2)
)


model_22 <- jags(model.file = textConnection
                 (model22),data = dat,n.chains=3, inits = inits,
                 parameters.to.save=c("ppc_test","ppc_measure" ),
                 n.burnin = 40000, 
                 n.iter = 50000, n.thin=1)





model_22_ppc_mcmc= as.mcmc(model_22)

summary(model_22_ppc_mcmc)





ms2=ggs(md2)


summary(md2)

ms2$DIC 

model_22$BUGSoutput$DIC

model_22$BUGSoutput$pD


### Model without random effect for care facility


model23<-"model{

for(i in 1:Nobs) {			
  for(j in 1:nrep2[i]){		
    y[i, j] ~ dbern(pi[i,j])	
    logit(pi[i, j]) =   SUB[i,j]
    SUB[i,j] ~ dnorm(beta[1] + beta[2]*SEX[i] + beta[3]*AGE[i]
    + beta[4]*EASY[i]+ beta[5]*TRAIN[i,j] + 
    beta[6]*EDU[i,j]+ beta[7]*SEX_ND[i,j],
    sigma_sub_inv)
    
  }		
  
  		
  		
  
}			

#Hyperpriors for SUB Random Effect

sigma.sub ~ dunif(0, 20)
sigma_sub_inv <- 1/pow(sigma.sub ,2)

#Priors for Fixed Effects

for (p in 1:7){

beta[p] ~ dnorm(0,0.0001)
 }

}"




dat <- list(y=s2[,c("GUM_ND_1","GUM_ND_2", "GUM_ND_3", "GUM_ND_4")]
            ,Nobs=nrow(s2), nrep2=nrep2,
            TRAIN=s2[,c("TRAIN_1" ,
                        "TRAIN_2" , "TRAIN_3",  "TRAIN_4")],
            SEX=s2$SEX, AGE=s2$AGE, EASY=s2$EASY, EDU=s2[,
                                                         c("EDU_1","EDU_2","EDU_3" ,"EDU_4")], 
            SEX_ND=s2[, c("SEX_ND_1", "SEX_ND_2", "SEX_ND_3", "SEX_ND_4")])

inits = list(
  list(beta=c(1.28486,0.69710,0.01199,0.17843,-0.48455,
              -0.47600,1.83693),
       sigma.sub=0.5,sigma.fac=1.5),
  list(beta=c(1.28486,0.69710,0.01199,0.17843,-0.48455,
              -0.47600,1.83693)*2,
       sigma.sub=0.5*2,sigma.fac=1.5*2),
  list(beta=c(1.28486,0.69710,0.01199,0.17843,-0.48455,
              -0.47600,1.83693)/2,
       sigma.sub=0.5/2, sigma.fac=1.5/2)
)


model_23 <- jags(model.file = textConnection
                 (model23),data = dat,n.chains=3, inits = inits,
                 parameters.to.save=c("beta","sigma.sub",
                                      "sigma.fac"), n.burnin = 5000, 
                 n.iter = 10000)


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










