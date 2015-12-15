### power calculation
#install.packages('pwr')
library(pwr)
pwr.2p.test(power=0.8,sig.level=0.05,h=-0.05,alternative='less')


#### Meta-analysis for endemic countries symptom probability

alphas = c(33,16,8,7,161,29)
betas = c(23,35,22,48,404,239)

expected = digamma(alphas) - digamma(alphas + betas)
vars = trigamma(alphas) - trigamma(alphas + betas)

#install.packages('metafor')
library(metafor)
mod = rma.uni(yi=expected,vi=vars)

exp(c(mod$b,mod$ci.lb,mod$ci.ub))

for (i in 1:6){
  print(qbeta(c(0.025,0.975),alphas[i],betas[i]))
}
alphas/(alphas+betas)

set.seed(1)
sigma = exp(rnorm(10000,mod$b,mod$se))

prevE = prevI = matrix(NA,10000,6)
inc = c(0.5,1,2,5,10,1.8)
for (i in 1:6){
  prevE[,i] = (inc[i]/1000)*((1-sigma)/sigma)*(1.55/365.25)
  prevI[,i] = (inc[i]/1000)*((1-sigma)/sigma)*(5.09/365.25)
}

quantile(prevE[,5],c(0.5,0.025,0.975))
quantile(prevI[,5],c(0.5,0.025,0.975))
quantile((prevE+prevI)[,5],c(0.5,0.025,0.975))*100000

#p = mean((prevE+prevI)[,5])

EE = EI = II = matrix(NA,10000,6)
delta = 1/1.55; gammaA = 1/5.09
for (i in 1:6){
  EE[,i] = prevE[,i]*(1-pexp(2,delta))
  EI[,i] = prevE[,i]*(delta*(1-exp(-2*gammaA))-gammaA*(1-exp(-2*delta)))/(delta-gammaA)
  II[,i] = prevI[,i]*(1-pexp(2,gammaA))
}

p = EE + EI + II
quantile(p[,5]*100000,c(0.5,0.025,0.975))


i = 6
quantile(((dbinom(3,454,p[,i]))/(dbinom(1,454,p[,i])+dbinom(2,454,p[,i])+dbinom(3,454,p[,i])))*100,c(0.5,0.025,0.975))

i = 6
quantile(dbinom(3,454,p[,i])*100,c(0.5,0.025,0.975))

quantile(((1-dbinom(0,454,p[,i])))*100,c(0.5,0.025,0.975))

for (i in 1:5){
  print(quantile(p[,i],c(0.5,0.025,0.975))*100000)
  #print(quantile((prevE[,i]+prevI[,i])*100000,c(0.5,0.025,0.975)))
}


#save(p,file='p.Rdata')

#############################
### Sensitivity specificity
#############################

set.seed(1)
sens = rbeta(10000,65,3)
spec = rbeta(10000,61,5)
quantile(spec,c(0.5,0.025,0.975))
pRDT = matrix(NA,10000,6)
for (i in 1:6){
  pRDT[,i] = spec*prevE[,i]*(exp(-2*delta)+(delta*(1-exp(-2*gammaA))-gammaA*(1-exp(-2*delta)))/(delta-gammaA)) + (1-sens)*prevI[,i]*exp(-2*gammaA)
}
for (i in 1:5){
  print(quantile(pRDT[,i],c(0.5,0.025,0.975))*100000)
}


quantile((1-dbinom(0,454,pRDT[,j]))*100,c(0.5,0.025,0.975))

out = matrix(NA,6,3)
for (j in 1:6){
  out[j,] = quantile((1 - (1-dbinom(0,454,pRDT[,j]))/(1-dbinom(0,454,p[,j])))*100,c(0.5,0.025,0.975))
}
out

###########################
### Abx time of departure
###########################
library(fitdistrplus)

set.seed(1)
abx.durpars = qmedist(c(2.4,2.4,2.4,2.4,rep(2.74,94),3.07,3.07,3.07,3.07),distr='norm',probs=c(0.025,0.975))$estimate ##
gammaA = 0.1964637
abx.rskpars = qmedist(c(rep(log(0.18),4),rep(log(0.34),96),rep(log(0.66),4)),dist='norm',probs=c(0.025,.975))$estimate ### RR of shedding in prophylaxed vs. unproph exposed
gammaAbx = 1/(1/gammaA - rnorm(10000,abx.durpars[1],abx.durpars[2]))
nuAbx = rlnorm(10000,abx.rskpars[1],abx.rskpars[2])

pTOD = matrix(NA,1e4,6)
for (i in 1:6){
  pTOD[,i] = nuAbx*prevE[,i]*(exp(-2*delta) + (delta*(1-exp(-2*gammaAbx)) - gammaAbx*(1-exp(-2*delta)))/(delta-gammaAbx)) + prevI[,i]*exp(-2*gammaAbx)
}
for (i in 1:5){
  print(quantile(pTOD[,i],c(0.5,0.025,0.975))*100000)
}

j = 1
quantile((1-dbinom(0,454,pTOD[,j]))*100,c(0.5,0.025,0.975))

out = matrix(NA,6,3)
for (j in 1:6){
  out[j,] = quantile((1 - (1-dbinom(0,454,pTOD[,j]))/(1-dbinom(0,454,p[,j])))*100,c(0.5,0.025,0.975))
}
out


###########################
### Abx early initiated
###########################

prevEAbx = prevIAbx = matrix(NA,10000,6)
inc = c(0.5,1,2,5,10,1.8)
for (i in 1:6){
  prevEAbx[,i] = nuAbx*(inc[i]/1000)*((1-sigma)/sigma)*((1/delta)/365.25)
  prevIAbx[,i] = nuAbx*(inc[i]/1000)*((1-sigma)/sigma)*((1/gammaAbx)/365.25)
}
for (i in 1:5){
  print(quantile((prevEAbx+prevIAbx)[,i],c(0.5,0.025,0.975))*100000)
}


pEI = matrix(NA,1e4,6)
for (i in 1:6){
  pEI[,i] = nuAbx*prevEAbx[,i]*(exp(-2*delta) + (delta*(1-exp(-2*gammaAbx)) - gammaAbx*(1-exp(-2*delta)))/(delta-gammaAbx)) + prevIAbx[,i]*exp(-2*gammaAbx)
}
for (i in 1:5){
  print(quantile(pEI[,i],c(0.5,0.025,0.975))*100000)
}

quantile(pEI[,5],c(0.5,0.025,0.975))*100000
j = 6
quantile((1-dbinom(0,454,pEI[,j]))*100,c(0.5,0.025,0.975))

out = matrix(NA,6,3)
for (j in 1:6){
  out[j,] = quantile((1 - (1-dbinom(0,454,pEI[,j]))/(1-dbinom(0,454,p[,j])))*100,c(0.5,0.025,0.975))
}
out

###########################
### Abx sensitivity
###########################
set.seed(1)
pars.abx = c(2.7350040,0.1709258)
gammaA = 0.1964637
gammaAbxR = 1/(1/gammaA - 1.5*rnorm(5e3,pars.abx[1],pars.abx[2]))
nuAbxR = 0.5*rlnorm(10000,abx.rskpars[1],abx.rskpars[2])

pTODR = matrix(NA,1e4,5)
for (i in 1:5){
  pTODR[,i] = nuAbxR*prevE[,i]*(exp(-2*delta) + (delta*(1-exp(-2*gammaAbxR)) - gammaAbxR*(1-exp(-2*delta)))/(delta-gammaAbxR)) + prevI[,i]*exp(-2*gammaAbxR)
}

prevEAbxR = prevIAbxR = matrix(NA,10000,5)
for (i in 1:5){
  prevEAbxR[,i] = nuAbxR*(inc[i]/1000)*((1-sigma)/sigma)*((1/delta)/365.25)
  prevIAbxR[,i] = nuAbxR*(inc[i]/1000)*((1-sigma)/sigma)*((1/gammaAbxR)/365.25)
}

pEIR = matrix(NA,1e4,5)
for (i in 1:5){
  pEIR[,i] = nuAbxR*prevEAbxR[,i]*(exp(-2*delta) + (delta*(1-exp(-2*gammaAbxR)) - gammaAbxR*(1-exp(-2*delta)))/(delta-gammaAbxR)) + prevIAbxR[,i]*exp(-2*gammaAbxR)
}

########################
### Vax importation prob
########################


prevVaxE = prevVaxI = matrix(NA,10000,6)
inc = c(0.5,1,2,5,10,1.8)
for (i in 1:6){
  prevVaxE[,i] = (inc[i]/1000)*((1-sigma + sigma*nuVaxDis)/sigma)*(1.55/365.25)
  prevVaxI[,i] = (inc[i]/1000)*((1-sigma + sigma*nuVaxDis)/sigma)*(5.09/365.25)
}

prevAbxE = prevAbxI = matrix(NA,10000,6)
for (i in 1:6){
  prevAbxE[,i] = nuAbx*(inc[i]/1000)*((1-sigma + sigma*nuAbxDis)/sigma)*(1.55/365.25)
  prevAbxI[,i] = nuAbx*(inc[i]/1000)*((1-sigma + sigma*nuAbxDis)/sigma)*(5.09/365.25)
}

quantile(prevAbxE[,j]/prevE[,j],c(0.5,0.025,0.975))

########################
### Case probabilities
########################

### This is the workhorse section for calculating main outcomes (case probabilities and
### effectiveness estimates) prior to adjustment for pre-deployment probabilities.

### The format is as follows: the "input1a, 2a, and 3a" vectors
### store the date of the first cholera case (9999 if transmission is exhausted) under the
### intervention scenario. The "input 1b, 2b and 3b" vectors store the date of the first
### case under the relevant "base-case" scenario.


set.seed(101)
probsA = probsB = matrix(NA,10000,3)
input1a = vax.chol1sim[,1]; input1b = chol1sim[,1]
input2a = vax.chol2sim[,1]; input2b = chol2sim[,1]
input3a = vax.chol3sim[,1]; input3b = chol3sim[,1]

### Case probability under intervention (A) and base-case (B)
for (i in 1:10000){
  ind = sample(1:5000,5000,replace=T)
  probsA[i,] = c(mean(input1a[ind]!=9999),mean(input2a[ind]!=9999),mean(input3a[ind]!=9999))
  probsB[i,] = c(mean(input1b[ind]!=9999),mean(input2b[ind]!=9999),mean(input3b[ind]!=9999))
}
### Effectiveness
for (j in 1:3){
  print(100*quantile(probsA[,j],c(0.5,0.025,0.975)))
  print(quantile((1-probsA[,j]/probsB[,j])*100,c(0.5,0.025,0.975)))
}


############################
### Direct RDT sensitivity
############################


truepos = c(65,66,189)
falseneg = c(2,6,14)

trueneg = c(24,102,209)
falsepos = c(8,38,213)

#### sens
alpha = truepos; alpha.beta = truepos + falseneg
truepos+falseneg
E = digamma(alpha) - digamma(alpha.beta)
V = trigamma(alpha) - trigamma(alpha.beta)

library(metafor)
mod = rma.uni(E,V,method='REML')
summary(mod)

set.seed(1)
out = rnorm(10000,mod$b,mod$se)
hist(exp(out))

sensDIR = exp(out)
sensDIR[sensDIR>1] = 1

alpha = trueneg; alpha.beta = trueneg + falsepos

E = digamma(alpha) - digamma(alpha.beta)
V = trigamma(alpha) - trigamma(alpha.beta)

quantile(exp(rnorm(10000,-1.31,sqrt(0.08))),c(0.5,0.025,0.975))

mod = rma.uni(E,V,method='REML')
set.seed(10)
out = rnorm(10000,mod$b,mod$se)

specDIR = spec# exp(out)
specDIR[specDIR>1] = 1

#sensDIR = 0.9; specDIR = 0.9
#sensDIR = 0.99; specDIR = 0.99

pRDTdir = matrix(NA,10000,5)
for (i in 1:5){
  pRDTdir[,i] = specDIR*prevE[,i]*(exp(-2*delta)+(delta*(1-exp(-2*gammaA))-gammaA*(1-exp(-2*delta)))/(delta-gammaA)) + (1-sensDIR)*prevI[,i]*exp(-2*gammaA)
}

j = 2
quantile((1-dbinom(0,454,pRDTdir[,j]))*100,c(0.5,0.025,0.975))

out = matrix(NA,5,3)
for (j in 1:5){
  out[j,] = quantile((1 - (1-dbinom(0,454,pRDTdir[,j]))/(1-dbinom(0,454,p[,j])))*100,c(0.5,0.025,0.975))
}
out

#################################
### Specificity/False positives
#################################
specDIR = 0.5
rdtspec = (1 - specDIR)*(1-prevI)
cispec = array(NA,dim=c(10000,5,3))
for (j in 1:5) for (i in 1:10000){
  cispec[i,j,] = qbinom(c(0.5,0.025,0.975),454,rdtspec[i,j])
}
mean(cispec[,1,1]); mean(cispec[,1,2]); mean(cispec[,1,3])
mean(cispec[,5,1]); mean(cispec[,5,2]); mean(cispec[,5,3])

fpos = array(NA,dim=c(10000,6))
for (i in 1:length(spec)) for (j in 1:6){
  fpos[i,j] = rbinom(1,454,(1-spec[i])*(1-prevI[i,j]-prevE[i,j]))
}

quantile(fpos[,1],c(0.5,0.025,0.975))/454

########################
## Considering vaccine protection against infection (sensitivity analysis)
########################

pvax5 = 0.95*p
pvax10 = 0.9*p
pvax25 = 0.75*p
pvax50 = 0.5*p

########################
### Total effects
########################

input1 = chol1sim[,1]; input2 = chol2sim[,1]; input3 = chol3sim[,1]
temp = c(); out = matrix(NA,6,3)


totalprob.fn = function(input1,input2,input3,prev){
  temp = c(); out = matrix(NA,6,3)
  for (j in 1:6){
    set.seed(101)
    for (i in 1:10000){
      ind = sample(1:5000,5000,replace=T)
      probs = c(mean(input1[ind]!=9999),mean(input2[ind]!=9999),mean(input3[ind]!=9999))
      pre = dbinom(1:3,454,prev[i,j])
      temp[i] = pre %*% probs
    }
    out[j,] = quantile(temp,c(0.5,0.025,0.975))
  }
  print(out*100)
}
totalprob.fn(input1=vax.chol1sim[,1],input2=vax.chol2sim[,1],input3=vax.chol3sim[,1],prev=pvax50)

totaleff.fn = function(input1a,input2a,input3a,input1b,input2b,input3b,preva,prevb){
  tempA = tempB = c(); out = matrix(NA,6,3)
  for (j in 1:6){
    set.seed(101)
    for (i in 1:10000){
      ind = sample(1:5000,5000,replace=T)
      probsA = c(mean(input1a[ind]!=9999),mean(input2a[ind]!=9999),mean(input3a[ind]!=9999))
      probsB = c(mean(input1b[ind]!=9999),mean(input2b[ind]!=9999),mean(input3b[ind]!=9999))
      preA = dbinom(1:3,454,preva[i,j])
      preB = dbinom(1:3,454,prevb[i,j])
      tempA[i] = preA %*% probsA
      tempB[i] = preB %*% probsB
    }
    out[j,] = quantile(1 - tempA/tempB,c(0.5,0.025,0.975))
  }
  print(out*100)
}

totaleff.fn(input1a=vax.chol1sim[,1],input2a=vax.chol2sim[,1],input3a=vax.chol3sim[,1],preva=pvax50,
            input1b=chol1sim[,1],input2b=chol2sim[,1],input3b=chol3sim[,1],prevb=p)

mean(chol1sim[chol1sim[,1]!=9999,1]<8)

#################################
### Protection to peacekeepers
#################################


set.seed(1)
### probability of disease | infection, baseline
piBase = rbeta(10000,6,10)

### probability of disease | infection, cipro
piClin = rbeta(10000,1,14)

nuAbxDis = piClin/piBase ### conditional prob for disease given infection
quantile(nuAbxDis,c(0.5,0.025,0.975))
quantile(1 - nuAbxDis,c(0.5,0.025,0.975))
quantile(1 - nuAbx*nuAbxDis,c(0.5,0.025,0.975))
gammaD = 0.3012048

bp = (inc/1000)*(((1/delta)/365.25) + (1/gammaD)/365.25) ### baseline probabilities

### Individual efficacy
int = (((1/delta)/365.25)*nuAbx*nuAbxDis*1.1 + (1/gammaD)/365.25)
base = (((1/delta)/365.25) + (1/gammaD)/365.25)
quantile(1 - int/base,c(0.5,0.025,0.975))*100

int = ((1/delta + 1/gammaD)/365.25)*nuAbx*nuAbxDis
base = ((1/delta + 1/gammaD)/365.25)
quantile(1 - int/base,c(0.5,0.025,0.975))*100

j = 5
p.tod = ((inc[j]/1000)*(((1/delta)/365.25)*(1.25*nuAbx*1.25*nuAbxDis) + (1/gammaD)/365.25))
quantile(p.tod,c(0.5,0.025,0.975))*100
quantile(1 - p.tod/bp[j],c(0.5,0.025,0.975))*100

j = 5
p.ei = (inc[j]/1000)*((1/delta + 1/gammaD)/365.25)*(nuAbxDis*1.25*1.25)
quantile(p.ei,c(0.5,0.025,0.975))*100
quantile(1-p.ei/bp[j],c(0.5,0.025,0.975))*100

set.seed(1)
piBaseV = rbeta(10000,7+6,0+2)
piClinV = rbeta(10000,3,6)
DisVax = piClinV/piBaseV
quantile(1-DisVax,c(0.5,0.025,0.975))

piBaseVDI = rbeta(10000,7+6,0+2)
piClinVDI = rbeta(10000,4+3,6+4)
DisVaxInfect = piClinVDI/piBaseVDI
quantile(1-piClinVDI/piBaseVDI,c(0.5,0.025,0.975))

quantile(p[,5],c(0.5,0.025,0.975))*100000

quantile(1 - nuAbx*nuAbxDis*DisVaxInfect,c(0.5,0.025,0.975))
quantile(prevI[,6],c(0.5,0.025,0.975))*100000
quantile((10.0/1000)*((1-sigma)/sigma)*((1/delta)/365.25),c(0.5,0.025,0.975))*100000

#piBaseInfectV = rbeta(10000,)
#piClinInfectV

#piBaseInfectV = rbeta(10000,7+,0)
#piClinInfectV = rbeta(10000,10,1)

#InfectV = c(rep(1,10+7),rep(0,1+2))
#InfectCtl = c(rep(1,7+8))

DisGInfectV = c(rep(1,4+3),rep(0,6+4))
DisGInfectCtl = c(rep(1,7+6),rep(0,0+2))
#DisV = c(rep(1,4+3),rep(0,7+6))
#DisCtl = c(rep(1,))

#set.seed(1)
#pVaxInfect = pVaxDisGInfect = c()
#for (i in 1:10000){
#  pVaxInfect[i] = mean(sample(InfectV,length(InfectV),replace=T))
#  pVaxDisGInfect[i] = mean(sample(DisGInfectV,length(DisGInfectV),replace=T))/mean(sample(DisGInfectCtl,length(DisGInfectCtl),replace=T))
# # pVaxDis[i] = mean(sample())
#}

quantile(1 - pVaxDisGInfect,c(0.5,0.025,0.975))

nuVaxDis = piClinV/piBaseV
quantile(1 - nuVaxDis,c(0.5,0.025,0.975))

mean(chol1sim[,1][chol1sim[,1]!=9999]<=4)

j = 3
p.vax = (inc[j]/1000)*((1/delta + 1/gammaD)/365.25)*nuVaxDis
quantile(p.vax,c(0.5,0.025,0.975))*100
quantile(1-p.vax/bp[j],c(0.5,0.025,0.975))*100

j = 3
p.todV = (inc[j]/1000)*nuVaxDis*((nuAbx*(1/delta) + 1/gammaD)/365.25)
quantile(p.todV,c(0.5,0.025,0.975))*100
quantile(1-p.todV/bp[j],c(0.5,0.025,0.975))*100

j = 3
p.eiV = (inc[j]/1000)*nuAbx*nuAbxDis*nuVaxDis*((1/delta + 1/gammaD)/365.25)
quantile(p.eiV,c(0.5,0.025,0.975))*100
quantile(1-p.eiV/bp[j],c(0.5,0.025,0.975))*100

p.tod = p.ei = p.vax = p.todV = p.eiV = matrix(NA,10000,5)
for (i in 1:5){
  p.tod[,i] = (inc[i]/1000)*(((1/delta)/365.25)*nuAbx*nuAbxDis*1.25 + (1/gammaD)/365.25)
  p.ei[,i] = (inc[i]/1000)*((1/delta + 1/gammaD)/365.25)*nuAbx*nuAbxDis*1.25
  #p.vax[,i] = (inc[i]/1000)*((1/delta + 1/gammaD)/365.25)*nuVaxDis
  #p.todV[,i] = (inc[i]/1000)*nuVaxDis*((nuAbx*nuAbxDis*(1/delta) + 1/gammaD)/365.25)
  #p.eiV[,i] = (inc[i]/1000)*nuAbx*nuAbxDis*nuVaxDis*((1/delta + 1/gammaD)/365.25)
}

#################################
### Protection to battalions
#################################


#set.seed = 1
#prob.base = matrix(NA,10000,5)
#for (i in 1:10000) for (j in 1:5){
#  prob.base[i,j] = mean(rbinom(10000,454,bp[j]))
#}

j = 1
1 - dbinom(0,454,p.tod[,1])
1 - dbinom(0,454,prob.base[,1])
j = 1


prob.base[,1]
dbinom(0,454,prob.base[,1])
dbinom(0,454,bp[1])
j = 4
quantile(1-dbinom(0,454,p.tod[,j]),c(0.5,0.025,0.975))*100
quantile(1 - (1-dbinom(0,454,p.tod[,j]))/(1-dbinom(0,454,bp[j])),c(0.5,0.025,0.975))*100

quantile(1-dbinom(0,454,p.ei[,j]),c(0.5,0.025,0.975))*100
quantile(1 - (1-dbinom(0,454,p.ei[,j]))/(1-dbinom(0,454,bp[j])),c(0.5,0.025,0.975))*100

quantile(1-dbinom(0,454,p.vax[,j]),c(0.5,0.025,0.975))*100
quantile(1 - (1-dbinom(0,454,p.vax[,j]))/(1-dbinom(0,454,bp[j])),c(0.5,0.025,0.975))*100

quantile(1-dbinom(0,454,p.todV[,j]),c(0.5,0.025,0.975))*100
quantile(1 - (1-dbinom(0,454,p.todV[,j]))/(1-dbinom(0,454,bp[j])),c(0.5,0.025,0.975))*100

quantile(1-dbinom(0,454,p.eiV[,j]),c(0.5,0.025,0.975))*100
quantile(1 - (1-dbinom(0,454,p.eiV[,j]))/(1-dbinom(0,454,bp[j])),c(0.5,0.025,0.975))*100


#############################
### R0
############################

pars = c(1/1.55,1/5.09,1/3.32,1/1.77,100,1,0.1,0.025,9923243,879644,9043599,1,1/30,8.58)
parnames = c('delta','gammaA','gammaD','gammaC','eta','xi','kappa','zeta','N','Nart','Noth','mu1','mu2','nu',
             'pi','betaW','betaL','k','omega','popA')
RW1 = RH1 = R01 = RA1 = c()
r = 1 + log10(nu)

set.seed(202)
nsims=5000
for (l in 1:nsims){
  par = c(pars,state1[sample(2001:25000,1),,sample(1:3,1)])
  for (z in 1:length(parnames)){
    assign(parnames[z],par[z])
  }
  RW1[l] = betaW*(gammaC*gammaD*(1-zeta)*kappa*(mu1+eta*mu2)+betaW*(mu1+eta*eta*mu2)*(gammaA*(gammaD*nu*gammaC)-gammaC*gammaD*(1-zeta)))*xi*popA*Nart/(gammaA*gammaC*gammaD*(1-zeta)*kappa*kappa*mu1*mu2)
  RH1[l] = betaL*((1-pi)/gammaA + pi*(gammaA*r + gammaD)/(gammaC*gammaD*(1-zeta)))
  R01[l] = (RW1[l]*popA*Nart + RH1[l]*N)/N
  RA1[l] = RW1[l]*popA + RH1[l]
}

RW2 = RH2 = R02 = RA2 = c()
for (l in 1:nsims){
  par = c(pars,state2[sample(2001:25000,1),,sample(1:3,1)])
  for (z in 1:length(parnames)){
    assign(parnames[z],par[z])
  }
  RW2[l] = betaW*(gammaC*gammaD*(1-zeta)*kappa*(mu1+eta*mu2)+betaW*(mu1+eta*eta*mu2)*(gammaA*(gammaD*nu*gammaC)-gammaC*gammaD*(1-zeta)))*xi*popA*Nart/(gammaA*gammaC*gammaD*(1-zeta)*kappa*kappa*mu1*mu2)
  RH2[l] = betaL*((1-pi)/gammaA + pi*(gammaA*r + gammaD)/(gammaC*gammaD*(1-zeta)))
  R02[l] = (RW2[l]*popA*Nart + RH2[l]*N)/N
  RA2[l] = RW2[l]*popA + RH2[l]
}


RW3 = RH3 = R03 = RA3 = c()
for (l in 1:nsims){
  par = c(pars,state3[sample(2001:25000,1),,sample(1:3,1)])
  for (z in 1:length(parnames)){
    assign(parnames[z],par[z])
  }
  RW3[l] = betaW*(gammaC*gammaD*(1-zeta)*kappa*(mu1+eta*mu2)+betaW*(mu1+eta*eta*mu2)*(gammaA*(gammaD*nu*gammaC)-gammaC*gammaD*(1-zeta)))*xi*popA*Nart/(gammaA*gammaC*gammaD*(1-zeta)*kappa*kappa*mu1*mu2)
  RH3[l] = betaL*((1-pi)/gammaA + pi*(gammaA*r + gammaD)/(gammaC*gammaD*(1-zeta)))
  R03[l] = (RW3[l]*popA*Nart + RH3[l]*N)/N
  RA3[l] = RW3[l]*popA + RH3[l]
}

for (j in 1:6){
  for (i in 1:length(R01)){
    out[i] = c(RH1[i],RH2[i],RH3[i]) %*% dbinom(1:3,454,p[i,j])/sum(dbinom(1:3,454,p[i,j]))
  }
  print(quantile(out,c(0.5,0.025,0.975)))
}

quantile(RA3,c(0.5,0.025,0.975))


################################
### Plotting
################################


dist1 = density(state.weakprior[10001:50000,1,]); dist2 = density(state1[10001:50000,1,])

layout(mat=matrix(c(1:4),nrow=2,ncol=2,byrow=T),widths=c(1.25,1),heights=c(1,1.35))

par(mar=c(1,5,2,1))
plot(y=dbeta(seq(0.1,0.4,by=0.001),26.3,81.2),x=seq(0.1,0.4,by=0.001),xlim=c(0.1,0.4),type='l',col='blue',lwd=2,bty='n',yaxt='n',xaxt='n',ylab=NA,xlab=NA,main=NA)
par(new=T); plot(y=dist1$y/sum(dist1$y),x=dist1$x,xlim=c(0.1,0.4),type='l',lwd=4,col=rgb(1,0,0,alpha=0.5),bty='n',yaxt='n',ylab='Density',xlab=NA)
mtext(expression(paste(pi,' ~ ','Beta(26.3, 81.2)',sep='')),3,adj=0)
axis(2,labels=F)

par(mar=c(1,1,2,1))
plot(y=dbeta(seq(0.1,0.4,by=0.001),263,812),x=seq(0.1,0.4,by=0.001),xlim=c(0.1,0.4),type='l',col='blue',lwd=2,bty='n',yaxt='n',xaxt='n',ylab=NA,xlab=NA,main=NA)
par(new=T); plot(y=dist2$y,x=dist2$x,xlim=c(0.1,0.4),type='l',lwd=4,col=rgb(1,0,0,alpha=0.5),bty='n',ylab=NA,xlab=NA,yaxt='n')
mtext(expression(paste(pi,' ~ ','Beta(263, 812)',sep='')),3,adj=0)
axis(2,labels=F)

par(mar=c(5,5,2,1))
plot((state.weakprior[10001:50000,1,]),state.weakprior[10001:50000,3,],cex=0.01,ylim=c(0.14,0.4),xlim=c(0.1,0.4),bty='n',ylab=expression(hat(beta)[L]),xlab=expression(hat(pi)))
par(mar=c(5,1,2,1))
plot((state1[10001:50000,1,]),state1[10001:50000,3,],ylim=c(0.14,0.4),xlim=c(0.1,0.4),cex=0.01,bty='n',ylab=NA,xlab=expression(hat(pi)),yaxt='n')
axis(2,labels=F)


####################
### MCMC trace plot
###################
h = 6
w = var(state3[10001:50000,h,1])+var(state3[10001:50000,h,2])+var(state3[10001:50000,h,3])
n = 40000
m = 3
b = (n/(m-1))*((mean(state3[10001:50000,h,1]) - mean(state3[10001:50000,h,]))^2
               +(mean(state3[10001:50000,h,2]) - mean(state3[10001:50000,h,]))^2
               +(mean(state3[10001:50000,h,3]) - mean(state3[10001:50000,h,]))^2)
theta = (1-1/n)*w + (1/n)*b

sqrt(theta/w)

statep = state1; statep[,2,] = log(state1[,2,]); statep[,4,] = log(state1[,4,])
library(coda)
?gelman.diag

tryobj = as.mcmc(statep[10001:50000,,])


par(mfrow=c(2,3))
for (i in 1:6){
  plot(statep[,i,1],type='l',col='red',ylim=range(statep[,i,]))
  lines(statep[,i,2],col='blue')
  lines(statep[,i,3],col='green')
}
plot(acf(state1[10001:50000,4,1]))
for (i in 1:6){
  distr = list()
  for (j in 1:3){
    distr[[j]] = density(statep[10001:50000,i,j])
  }
  plot(y=log(distr[[1]]$y),x=distr[[1]]$x,col='red',type='l',ylim=c(0,2*max(log(distr[[1]]$y))))
  lines(y=log(distr[[2]]$y),x=distr[[2]]$x,col='blue')
  lines(y=log(distr[[3]]$y),x=distr[[3]]$x,col='green')
}

#####################
### First case plot
#####################

w = matrix(NA,10000,3)
for (i in 1:10000){
  w[i,] = dbinom(c(1,2,3),454,p[i,6])/sum(dbinom(c(1,2,3),454,p[i,6]))
}

first = c()
for (i in 1:100000){
  ind = sample(1:10000,1)
  weight = sample(1:3,1,prob=w[ind,])
  if (weight==1){
    first[i] = sample(chol1sim[chol1sim[,1]!=9999,1],1)
  } else{
    if (weight==2){
      first[i] = sample(chol2sim[chol2sim[,1]!=9999,1],1)
    } else{
      if (weight==3){
        first[i] = sample(chol3sim[chol3sim[,1]!=9999,1],1)
      }
    }
  }
}

df = density(first)

ivers.case = data.frame(y=(df$y[1:which.min(abs(df$x-4))]),
                        x=(df$x[1:which.min(abs(df$x-4))]))

cuba.case = data.frame(y=(df$y[1:which.min(abs(df$x-6))]),
                       x=(df$x[1:which.min(abs(df$x-6))]))

piarroux.case = data.frame(y=(df$y[1:which.min(abs(df$x-8))]),
                           x=(df$x[1:which.min(abs(df$x-8))]))
mean(first<6)
#layout(matrix(1:2,nrow=2),heights=c(1,1.35))
#par(mar=c(2,5,1,1))
par(mar=c(5,1,1,1))
plot(y=df$y[df$x<=22],x=df$x[df$x<=22],type='l',lwd=0,bty='n',yaxt='n',ylab='Probability density',main=NA,xaxt='n',xlab='October 2010')
polygon(y=c(ivers.case$y,0,0),x=c(ivers.case$x,max(ivers.case$x),0),col=rgb(1,0,0,alpha=0.6),border=NA)
polygon(y=c(cuba.case$y,0,0),x=c(cuba.case$x,max(cuba.case$x),0),col=rgb(0.8,0,0.2,alpha=0.4),border=NA)
polygon(y=c(piarroux.case$y,0,0),x=c(piarroux.case$x,max(piarroux.case$x),0),col=rgb(0.5,0,0.5,alpha=0.25),border=NA)
polygon(y=c(df$y[df$x<=22],0,0),x=c(df$x[df$x<=22],0,0),col=rgb(0,0,1,alpha=0.15),border=NA)
lines(y=c(0,max(ivers.case$y)),x=rep(max(ivers.case$x),2),lty='dotted',lwd=2)
lines(y=c(0,cuba.case$y[which.max(cuba.case$x)]),x=rep(max(cuba.case$x),2),lty='dotted',lwd=2)
lines(y=c(0,piarroux.case$y[which.max(piarroux.case$x)]),x=rep(max(piarroux.case$x),2),lty='dotted',lwd=2)
lines(y=df$y[df$x<=22],x=df$x[df$x<=22],col='dark blue',lwd=4)
axis(1,at=0:22,labels=F)
axis(1,at=c(0.5,3.5,5.5,7.5,21.5),tick=F,labels=c('9','12','14','16','31'))
text(x=c(2,4.5,6.5)+0.5,y=rep(0.025,3),labels=c('A','B','C'),font=2)
#axis(2,labels=F)
?metafor
par(mar=c(5,5,1,1))
plot(y=cumsum(df$y[df$x<=22])/sum(df$y[df$x<=22]),x=df$x[df$x<=22],type='l',lwd=0,bty='n',yaxt='n',ylab='Cumulative probability',main=NA,xaxt='n',xlab='October 2010')
polygon(y=c(cumsum(ivers.case$y),0,0)/sum(df$y[df$x<=22]),x=c(ivers.case$x,max(ivers.case$x),0),col=rgb(1,0,0,alpha=0.6),border=NA)
polygon(y=c(cumsum(piarroux.case$y),rep(cumsum(piarroux.case$y[which.min(abs(piarroux.case$x-4))]),2))/sum(df$y[df$x<=22]),x=c(piarroux.case$x,max(piarroux.case$x),0),col=rgb(1,0.2,0,alpha=0.25),border=NA)
#polygon(y=c(cumsum(df$y[df$x<=22]),0,0)/sum(df$y[df$x<=22]),x=c(df$x[df$x<=22],max(df$x[df$x<=22]),0),col=rgb(0,0,1,alpha=0.25),border=NA)
lines(y=c(0,sum(ivers.case$y))/sum(df$y[df$x<=22]),x=rep(max(ivers.case$x),2),lty='dotted',lwd=2)
lines(y=c(0,sum(piarroux.case$y))/sum(df$y[df$x<=22]),x=rep(max(piarroux.case$x),2),lty='dotted',lwd=2)
lines(y=cumsum(df$y[df$x<=22])/sum(df$y[df$x<=22]),x=df$x[df$x<=22],col='dark blue',lwd=4)
axis(1,at=0:22,labels=F)
axis(1,at=c(0.5,3.5,7.5,21.5),tick=F,labels=c('9','12','16','31'))
