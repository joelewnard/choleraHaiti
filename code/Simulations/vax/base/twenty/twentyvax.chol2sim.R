
varnames = c('PK',
             'S1a','S1b','S2',
             'Ea1a','Ea1b','Ea2',
             'Es1a','Es1b','Es2',
             'Ia1a','Ia1b','Ia2',
             'Is1a','Is1b','Is2',
             'Y1a','Y1b','Y2',
             'Ya1a','Ya1b','Ya2',
             'Ic1a','Ic1b','Ic2',
             'Ra1a','Ra1b','Ra2',
             'Rs1a','Rs1b','Rs2',
             'B1','B2'
)

parnames = c('delta','gammaA','gammaD','gammaC','eta','xi','kappa','zeta','N','Nart','Noth','mu1','mu2','nu',
             'pi','betaW','betaL','k','omega','popA')
nu=8.58
vaxd = 20*(9.9*10^5)/(5.1*10^7)

######## INFECTIONS: -1 Suscept, +1 latent
#### infections by Water (Asx and Sx)
infectW.A = rep(0,33); infectW.A[which(varnames=='S1a')] = -1; infectW.A[which(varnames=='Ea1a')] = 1
infectW.S = rep(0,33); infectW.S[which(varnames=='S1a')] = -1; infectW.S[which(varnames=='Es1a')] = 1

#### infections by Local (Asx and Sx, pops 1a, 1b, 2)
infectL.A.1a = rep(0,33); infectL.A.1a[which(varnames=='S1a')] = -1; infectL.A.1a[which(varnames=='Ea1a')] = 1
infectL.S.1a = rep(0,33); infectL.S.1a[which(varnames=='S1a')] = -1; infectL.S.1a[which(varnames=='Es1a')] = 1
infectL.A.1b = rep(0,33); infectL.A.1b[which(varnames=='S1b')] = -1; infectL.A.1b[which(varnames=='Ea1b')] = 1
infectL.S.1b = rep(0,33); infectL.S.1b[which(varnames=='S1b')] = -1; infectL.S.1b[which(varnames=='Es1b')] = 1
infectL.A.2 = rep(0,33); infectL.A.2[which(varnames=='S2')] = -1; infectL.A.2[which(varnames=='Ea2')] = 1
infectL.S.2 = rep(0,33); infectL.S.2[which(varnames=='S2')] = -1; infectL.S.2[which(varnames=='Es2')] = 1

#### migration
mig.S.1a = rep(0,33); mig.S.1a[which(varnames=='S1a')] = -1; mig.S.1a[which(varnames=='S1b')] = 1
mig.S.1b = rep(0,33); mig.S.1b[which(varnames=='S1b')] = -1; mig.S.1b[which(varnames=='S1a')] = 1
mig.Ea.1a = rep(0,33); mig.Ea.1a[which(varnames=='Ea1a')] = -1; mig.Ea.1a[which(varnames=='Ea1b')] = 1
mig.Ea.1b = rep(0,33); mig.Ea.1b[which(varnames=='Ea1b')] = -1; mig.Ea.1b[which(varnames=='Ea1a')] = 1
mig.Es.1a = rep(0,33); mig.Es.1a[which(varnames=='Es1a')] = -1; mig.Es.1a[which(varnames=='Es1b')] = 1
mig.Es.1b = rep(0,33); mig.Es.1b[which(varnames=='Es1b')] = -1; mig.Es.1b[which(varnames=='Es1a')] = 1
mig.Ia.1a = rep(0,33); mig.Ia.1a[which(varnames=='Ia1a')] = -1; mig.Ia.1a[which(varnames=='Ia1b')] = 1
mig.Ia.1b = rep(0,33); mig.Ia.1b[which(varnames=='Ia1b')] = -1; mig.Ia.1b[which(varnames=='Ia1a')] = 1
mig.Ic.1a = rep(0,33); mig.Ic.1a[which(varnames=='Ic1a')] = -1; mig.Ic.1a[which(varnames=='Ic1b')] = 1
mig.Ic.1b = rep(0,33); mig.Ic.1b[which(varnames=='Ic1b')] = -1; mig.Ic.1b[which(varnames=='Ic1a')] = 1

######## PROGRESSION: -1 latent, +1 infectious; +1 obs if severe
prog.A.1a = rep(0,33); prog.A.1a[which(varnames=='Ea1a')] = -1; prog.A.1a[which(varnames=='Ia1a')] = 1; prog.A.1a[which(varnames=='Ya1a')] = 1
prog.S.1a = rep(0,33); prog.S.1a[which(varnames=='Es1a')] = -1; prog.S.1a[which(varnames=='Is1a')] = 1; prog.S.1a[which(varnames=='Y1a')] = 1
prog.A.1b = rep(0,33); prog.A.1b[which(varnames=='Ea1b')] = -1; prog.A.1b[which(varnames=='Ia1b')] = 1; prog.A.1b[which(varnames=='Ya1b')] = 1
prog.S.1b = rep(0,33); prog.S.1b[which(varnames=='Es1b')] = -1; prog.S.1b[which(varnames=='Is1b')] = 1; prog.S.1b[which(varnames=='Y1b')] = 1
prog.A.2 = rep(0,33); prog.A.2[which(varnames=='Ea2')] = -1; prog.A.2[which(varnames=='Ia2')] = 1; prog.A.2[which(varnames=='Ya2')] = 1
prog.S.2 = rep(0,33); prog.S.2[which(varnames=='Es2')] = -1; prog.S.2[which(varnames=='Is2')] = 1; prog.S.2[which(varnames=='Y2')] = 1

######## CONVALESCENCE: -1 severe, +1 convalesc
conval.1a = rep(0,33); conval.1a[which(varnames=='Is1a')] = -1; conval.1a[which(varnames=='Ic1a')] = 1
conval.1b = rep(0,33); conval.1b[which(varnames=='Is1b')] = -1; conval.1b[which(varnames=='Ic1b')] = 1
conval.2 = rep(0,33); conval.2[which(varnames=='Is2')] = -1; conval.2[which(varnames=='Ic2')] = 1

######## MORTALITY: -1 severe
mortal.1a = rep(0,33); mortal.1a[which(varnames=='Is1a')] = -1
mortal.1b = rep(0,33); mortal.1b[which(varnames=='Is1b')] = -1
mortal.2 = rep(0,33); mortal.2[which(varnames=='Is2')] = -1

######## RECOVERY: -1 infectious (Asx and Cx); +1 recovered (Asx and Sx)
recov.A.1a = rep(0,33); recov.A.1a[which(varnames=='Ia1a')] = -1; recov.A.1a[which(varnames=='Ra1a')] = 1
recov.S.1a = rep(0,33); recov.S.1a[which(varnames=='Ic1a')] = -1; recov.S.1a[which(varnames=='Rs1a')] = 1
recov.A.1b = rep(0,33); recov.A.1b[which(varnames=='Ia1b')] = -1; recov.A.1b[which(varnames=='Ra1b')] = 1
recov.S.1b = rep(0,33); recov.S.1b[which(varnames=='Ic1b')] = -1; recov.S.1b[which(varnames=='Rs1b')] = 1
recov.A.2 = rep(0,33); recov.A.2[which(varnames=='Ia2')] = -1; recov.A.2[which(varnames=='Ra2')] = 1
recov.S.2 = rep(0,33); recov.S.2[which(varnames=='Ic2')] = -1; recov.S.2[which(varnames=='Rs2')] = 1

######## PK RECOVERIES: -1 PK
recov.PK = rep(0,33); recov.PK[which(varnames=='PK')] = -1

######## SHED: +1 B1
shed.PK = rep(0,33); shed.PK[which(varnames=='B1')] = vaxd
shed.A = rep(0,33); shed.A[which(varnames=='B1')] = 1
shed.S = rep(0,33); shed.S[which(varnames=='B1')] = nu ### shed nu vibrios relative to Asx and Conval
shed.C = rep(0,33); shed.C[which(varnames=='B1')] = 1

trans = rep(0,33); trans[which(varnames=='B1')] = -vaxd; trans[which(varnames=='B2')] = vaxd
die = rep(0,33); die[which(varnames=='B2')] = -vaxd

transitions = matrix(c(
  infectW.A,
  infectW.S,
  infectL.A.1a,
  infectL.S.1a,
  infectL.A.1b,
  infectL.S.1b,
  infectL.A.2,
  infectL.S.2,
  mig.S.1a,
  mig.S.1b,
  mig.Ea.1a,
  mig.Ea.1b,
  mig.Es.1a,
  mig.Es.1b,
  mig.Ia.1a,
  mig.Ia.1b,
  mig.Ic.1a,
  mig.Ic.1b,
  prog.A.1a,
  prog.S.1a,
  prog.A.1b,
  prog.S.1b,
  prog.A.2,
  prog.S.2,
  conval.1a,
  conval.1b,
  conval.2,
  mortal.1a,
  mortal.1b,
  mortal.2,
  recov.A.1a,
  recov.S.1a,
  recov.A.1b,
  recov.S.1b,
  recov.A.2,
  recov.S.2,
  recov.PK,
  shed.PK,
  shed.A,
  shed.S,
  shed.C,
  trans,
  die),byrow=T,ncol=33)

onestep = function(x,pars){
  
  for (z in 1:length(x)){
    x[z] = max(0,x[z])
  }
  
  for (i in 1:33){
    assign(varnames[i],x[i+1])
  }

  for (h in 1:length(parnames)){
    assign(parnames[h],pars[h])
  }
  r = 1+log10(nu)
  
  lambdaW = betaW*(eta*B1+B2)/(betaW*(eta*B1+B2) + kappa)
  lambdaL = k*log(1 + betaL*(Ia1a + Ia1b + Ia2 + r*(Is1a + Is1b + Is2) + Ic1a + Ic1b + Ic2)/k)/(Nart+Noth)
  
  rates = c(
    infectW.A = lambdaW*(1-lambdaW)*S1a,
    infectW.S = (lambdaW^2)*S1a,
    infectL.A.1a = (1-pi)*lambdaL*S1a,
    infectL.S.1a = pi*lambdaL*S1a,
    infectL.A.1b = (1-pi)*lambdaL*S1b,
    infectL.S.1b = pi*lambdaL*S1b,
    infectL.A.2 = (1-pi)*lambdaL*S2,
    infectL.S.2 = pi*lambdaL*S2,
    
    mig.S.1a = 0,#omega*(1-popA)*S1a,
    mig.S.1b = 0,#omega*popA*S1b,
    mig.Ea.1a = omega*(1-popA)*Ea1a,
    mig.Ea.1b = omega*popA*Ea1b,
    mig.Es.1a = omega*(1-popA)*Es1a,
    mig.Es.1b = omega*popA*Es1b,
    mig.Ia.1a = omega*(1-popA)*Ia1a,
    mig.Ia.1b = omega*popA*Ia1b,
    mig.Ic.1a = omega*(1-popA)*Ic1a,
    mig.Ic.1b = omega*popA*Ic1b,
    
    prog.A.1a = delta*Ea1a,
    prog.S.1a = delta*Es1a,
    prog.A.1b = delta*Ea1b,
    prog.S.1b = delta*Es1b,
    prog.A.2 = delta*Ea2,
    prog.S.2 = delta*Es2,
    
    conval.1a = gammaD*Is1a, ### (1 - zeta)*gammaD/(1-zeta) ### probability of recovery times exit rate
    conval.1b = gammaD*Is1b,
    conval.2 = gammaD*Is2,
    
    mortal.1a = zeta*gammaD*Is1a/(1-zeta), ## zeta*(gammaD/(1-zeta)) ### probability of death times exit rate
    mortal.1b = zeta*gammaD*Is1b/(1-zeta),
    mortal.2 = zeta*gammaD*Is2/(1-zeta),
    
    recov.A.1a = gammaA*Ia1a,
    recov.S.1a = gammaC*Ic1a,
    recov.A.1b = gammaA*Ia1b,
    recov.S.1b = gammaC*Ic1b,
    recov.A.2 = gammaA*Ia2,
    recov.S.2 = gammaC*Ic2,
    
    recov.PK = gammaA*PK,
    
    shed.PK = PK,
    shed.A = Ia1a,
    shed.S = Is1a, 
    shed.C = Ic1a,
    trans = mu1*B1/vaxd,
    die = mu2*B2/vaxd
  )
  if ((Ia1a==0)&(Ia1b==0)&(Ea1a==0)&(Ea1b==0)&(Is1a==0)&(Is1b==0)&(Es1a==0)&(Es1b==0)&(Ea2==0)&(Es2==0)&(Ia2==0)&(Is2==0)){
    rates['mig.S.1a'] = rates['mig.S.1b'] = 0
    rates['mig.Ea.1a'] = rates['mig.Ea.1b'] = 0
    rates['mig.Es.1a'] = rates['mig.Es.1b'] = 0
    rates['mig.Ia.1a'] = rates['mig.Ia.1b'] = 0
    rates['mig.Ic.1a'] = rates['mig.Ic.1b'] = 0
  }
  
  tot.rate = sum(rates)
  tau = rexp(n=1,rate=tot.rate)
  if (is.na(tau)){
    return('no transmission')
  } else{
    event = sample(1:length(rates),size=1,prob=rates/tot.rate)
    return(x+c(tau,transitions[event,]))
  }  
}


simul.fn = function(x,params,maxstep,tmax){
  names(x) = c('time',varnames)
  j = 0

  while (j<=maxstep){
    if (j>1){
      if ((x['time']>tmax)&(y['PK']==0)&(y['Ea1a']==0)&(y['Ia1a']==0)&(y['Es1a']==0)){
        return(c(9999,j))
      }
    }
    j = j+1
    y = onestep(x,params)
    
    names(y) = c('time',varnames)
    
    if (y[1]=='no transmission'){
      return(c(9999,j))
    }
    if (sum(c(y['PK']==0),(y['B1']==0),(y['B2']==0),
            (y['Ea1a']==0),(y['Es1a']==0),(y['Ia1a']==0),(y['Is1a']==0),(y['Ic1a']==0),
            (y['Ea1b']==0),(y['Es1b']==0),(y['Ia1b']==0),(y['Is1b']==0),(y['Ic1b']==0),
            (y['Ea2']==0),(y['Es2']==0),(y['Ia2']==0),(y['Is2']==0),(y['Ic2']==0))==18){
      return(c(9999,j))
    }
    if (sum(c(y['Y1a'],y['Y1b'],y['Y2'])>0)){
      return(c(y['time'],j))
    }
    x = y
    names(x) = c('time',varnames)
  }
  return(c(9999,j))
}




set.seed(20102)
load(file='~/chol2.mcmc.Rdata')
pars = c(1/1.55,1/5.09,1/3.32,1/1.77,100,1,0.1,0.025,9923243,879644,9043599,1,1/30,8.58)

nsims = 5e3
out = matrix(NA,nsims,2)

for (l in 1:nsims){
  par = c(pars,state2[sample(2001:25000,1),,sample(1:3,1)])
  for (z in 1:length(parnames)){
    assign(parnames[z],par[z])
  }
  init = c(2,popA*Nart,(1-popA)*Nart,Noth,rep(0,29))
  out[l,] = simul.fn(x=c(0,init),params=par,maxstep=5e4,tmax=100)
  if ((l/1e1)==ceiling(l/1e1)){
    print(l)
  }
}

twentyvax.chol2sim = out

save(twentyvax.chol2sim,file='twentyvax.chol2sim.Rdata')

