
dat1 = c(rep(0,7),as.matrix(read.csv('piarrouxAall.csv',header=F))) + c(rep(0,7),as.matrix(read.csv('piarrouxBall.csv',header=F)))
dat2 = c(rep(0,7),as.matrix(read.csv('piarrouxCall.csv',header=T)))

library(deSolve)


varnames = c('PK',
             'S1a','S1b','S2',
             'Ea1a','Ea1b','Ea2',
             'Es1a','Es1b','Es2',
             'Ia1a','Ia1b','Ia2',
             'Is1a','Is1b','Is2',
             'Y1a','Y1b','Y2',
             'Ic1a','Ic1b','Ic2',
             'Ra1a','Ra1b','Ra2',
             'Rs1a','Rs1b','Rs2',
             'B1','B2'
)

parnames = c('delta','gammaA','gammaD','gammaC','eta','xi','kappa','zeta','N','Nart','Noth','mu1','mu2','nu',
             'pi','betaW','betaL','k','omega','popA')


SIR.fn = function(t,y,parms){
  
  for (i in 1:length(varnames)){
    assign(varnames[i],y[i])
  }
  
  for (i in 1:length(parnames)){
    assign(parnames[i],parms[i])
  }
  
  r = 1 + log(nu,10)
  
  lambdaL = k*log(1 + betaL*(Ia1a + Ia1b + Ia2 + r*(Is1a + Is1b + Is2) + Ic1a + Ic1b + Ic2)/k)/(Nart+Noth)
  
  lambdaW = betaW*(eta*B1 + B2)/(betaW*(eta*B1 + B2) + kappa)
  piW = betaW*(eta*B1 + B2)/(betaW*(eta*B1 + B2) + kappa)
  
  pops = c(popA*Nart,(1-popA)*Nart,Noth)   

  dS1a = -(lambdaL + lambdaW)*S1a + omega*(popA*S1b - (1-popA)*S1a) 
  dEa1a = ((1-pi)*lambdaL + (1-piW)*lambdaW)*S1a - delta*Ea1a + omega*(popA*Ea1b - (1-popA)*Ea1a)
  dEs1a = (pi*lambdaL + lambdaW*piW)*S1a - delta*Es1a + omega*(popA*Es1b - (1-popA)*Es1a)
  dIa1a = delta*Ea1a - gammaA*Ia1a + omega*(popA*Ia1b - (1-popA)*Ia1a)
  dIs1a = delta*Es1a - (gammaD/(1-zeta))*Is1a
  dIc1a = gammaD*Is1a - gammaC*Ic1a + omega*(popA*Ic1b - (1-popA)*Ic1a)
  dRa1a = gammaA*Ia1a 
  dRs1a = gammaC*Ic1a
  dY1a = delta*Es1a
  
  dPK = -gammaA*PK
  dB1 = xi*(PK + Ia1a + nu*Is1a + Ic1a) - mu1*B1
  dB2 = mu1*B1 - mu2*B2

  dS1b = -lambdaL*S1b + omega*((1-popA)*S1a - popA*S1b)
  dEa1b = (1-pi)*lambdaL*S1b - delta*Ea1b + omega*((1-popA)*Ea1a - popA*Ea1b)
  dEs1b = pi*lambdaL*S1b - delta*Es1b + omega*((1-popA)*Es1a - popA*Es1b)
  dIa1b = delta*Ea1b - gammaA*Ia1b + omega*((1-popA)*Ia1a - popA*Ia1b)
  dIs1b = delta*Es1b - (gammaD/(1-zeta))*Is1b 
  dIc1b = gammaD*Is1b - gammaC*Ic1b + omega*((1-popA)*Ic1a - popA*Ic1b)
  dRa1b = gammaA*Ia1b 
  dRs1b = gammaC*Ic1b 
  dY1b = delta*Es1b
  
  dS2 = -lambdaL*S2 + omega*Ra2
  dEa2 = (1-pi)*lambdaL*S2 - delta*Ea2 
  dEs2 = pi*lambdaL*S2 - delta*Es2 
  dIa2 = delta*Ea2 - gammaA*Ia2 
  dIs2 = delta*Es2 - (gammaD/(1-zeta))*Is2 
  dIc2 = gammaD*Is2 - gammaC*Ic2 
  dRa2 = gammaA*Ia2 - omega*Ra2
  dRs2 = gammaC*Ic2
  dY2 = delta*Es2
  
  return(list(c(dPK,
                dS1a,dS1b,dS2,
                dEa1a,dEa1b,dEa2,
                dEs1a,dEs1b,dEs2,
                dIa1a,dIa1b,dIa2,
                dIs1a,dIs1b,dIs2,
                dY1a,dY1b,dY2,
                dIc1a,dIc1b,dIc2,
                dRa1a,dRa1b,dRa2,
                dRs1a,dRs1b,dRs2,
                dB1,dB2)))
}


pars = list()
pars[[1]] = c(1/1.55,1/5.09,1/3.32,1/1.77,100,1,0.1,0.025,9923243,879644,9043599,1,1/30,8.58,
         #### now fitted
         0.23, ## pi
         8e-8, ## betaW
         0.35, ## betaL
         5e2, ## k
         1e-2, ## omega
         0.035) ## popA

pars[[2]] = c(1/1.55,1/5.09,1/3.32,1/1.77,100,1,0.1,0.025,9923243,879644,9043599,1,1/30,8.58,
              #### now fitted
              0.2, ## pi
              1.3e-7, ## betaW
              0.4, ## betaL
              7e2, ## k
              5e-3, ## omega
              0.02) ## popA

pars[[3]] = c(1/1.55,1/5.09,1/3.32,1/1.77,100,1,0.1,0.025,9923243,879644,9043599,1,1/30,8.58,
              #### now fitted
              0.25, ## pi
              2.3e-7, ## betaW
              0.45, ## betaL
              4e2, ## k
              4e-2, ## omega
              0.03) ## popA

names(pars[[1]]) = parnames

set.seed(1)

nit = 5e4
state = array(NA,dim=c(nit,6,3))
state[1,,1] = pars[[1]][15:20]; state[1,,2] = pars[[2]][15:20]; state[1,,3] = pars[[3]][15:20]

post.state = rep(-1e10,3)
RL = RW = R0 = c()
for (i in 2:nit){

  for (c in 1:3){
    prior = -Inf
    
    while(prior==(-Inf)){
      
      step = vector(length=6)
      step[1:5] = c(
        rnorm(1,0,5e-2),
        rnorm(1,0,1e-9),
        0,
        rnorm(1,0,1e1),
        rnorm(1,0,1e-3)
      )
      step[3] = runif(1,-3,1)*step[1]
      #   step[4] = runif(1,-1,2)*step[1]
      #  step[4] = (1e3)*runif(1,-2,1)*step[5]
      step[6] = (1e6)*runif(1,-3,1)*step[2]
      step[sample(1:6,3)] = 0
      
      cand = state[i-1,,c] + step
      
      par = c(pars[[1]][1:14],cand)
      
      prior = sum(c(
        dbeta(par[15],263/100,812/100,log=T),
        dunif(par[16],0,1e-3,log=T),
        dunif(par[17],0.1,0.6,log=T),
         dunif(par[18],0,1e6,log=T),
        dunif(par[19],1e-4,1e-1,log=T),
       # dunif(par[19],1/(7*8.3),1/(0.3*7),log=T),
        #dgamma(1/par[19],mu*theta,mu/(sigma^2),log=T), ### gamma dist from King estimate (2-path model)
        dunif(par[20],1e-3,1e-1,log=T)
      ))
    }
    init = rep(0,length(varnames)); names(init)=varnames
    names(par) = parnames
    init = rep(0,length(varnames)); names(init)=varnames
    init['PK'] = 1;
    init['S1a'] = par['popA']*par['Nart'];
    init['S1b'] = (1-par['popA'])*par['Nart'];
    init['S2'] = par['Noth']
    
    out = ode(y=init,times=0:27,func=SIR.fn,parms=par,method='rk4')
    case1 = diff(out[,'Y1a']) + diff(out[,'Y1b']); case2 = diff(out[,'Y2']) 
    
    post.cand = prior + sum(dpois(round(dat1[1:27]),case1,log=T)) + sum(dpois(round(dat2[1:27]),case2,log=T))
      
    ratio = exp(post.cand - post.state[c])
    if (runif(1)<ratio){
      state[i,,c] = cand
      post.state[c] = post.cand
    } else{
      state[i,,c] = state[i-1,,c]
    }
  }
  
#  if ((i/1000)==round(i/1000)){
#    par(mfrow=c(4,2))
#    for (h in 1:6){
#      plot(state[1:(i-1),h,1],ylim=c(min(state[1:(i-1),h,]),max(state[1:(i-1),h,])),type='l',col='red')
#      lines(state[1:(i-1),h,2],col='blue'); lines(state[1:(i-1),h,3],col='green')
#    }
#    plot(case1,type='l',col='red',ylim=c(0,3000)); points(dat1[1:27]); plot(case2,type='l',col='red',ylim=c(0,600)); points(dat2[1:27])
#  }  
}
state1 = state
save(state1,file='chol1.mcmc.senspi.Rdata')




