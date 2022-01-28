#set.seed(191322)
#set.seed(1913268)
#SEIR on network
#setwd("~/Dropbox/Newyork_University/TB_Project")# setting work directory
remove (list = objects())

#library(statnet)
library(deSolve)
library(Matrix)
# class(cm_C) = "dsCMatrix"
library(fastnet)
#save.image('myEnvironment_SensitivityTest.RData')
#load('myEnvironment.RData')# #R.version.string #Verify installed R version


#sessionInfo() #prints latest R version in use and other functions/libraries installed 
# require(devtools)
# install_version("fastnet", version = "0.1.6", repos="http://cran.us.r-project.org")
# saveRDS(cm,file=<file name>)



# BarabasiAlbert=function(N, K){
#   CM=matrix(0, ncol=N, nrow=N)
#   CM[1,2]=1
#   CM[2,1]=1
#   for(i in 3:N){
#     probs=apply(CM, 1, sum)
#     link=unique(sample(c(1:N)[-i],
#                        size=min(c(K, i-1)), prob=probs[-i]))
#     CM[i, link]=CM[link, i]=1
#   }
#   class(CM)="cm"
#   return(CM)
# }

BarabasiAlbert=function(N, K){
  CM = Matrix(data = 0,nrow=N,ncol=N,sparse = T)
  edges = to.edgelist(net.barabasi.albert(N,K))
  CM[edges] = 1
  CM[cbind(edges[,2],edges[,1])] = 1
  return(CM)
}

rxnom <- function(probabilities, indicatorVector) {
  probsVec = runif(length(indicatorVector))
  sump = 0.0
  lastp = 0.0
  outMat = matrix(nrow=length(indicatorVector), ncol = length(probabilities))
  i = 1
  for (p in probabilities) {
    lastp = sump
    sump = sump + p
    thisInd = ((probsVec > lastp) & (probsVec <= sump))*indicatorVector
    outMat[,i] = thisInd
    i = i + 1
  }
  outMat
  
}

rStateFromSusceptable <- function(p_expose_vec,p_immune_scalar, p_die_scalar, indicatorVector) {
  probsVec = runif(length(indicatorVector))
  sump = 0.0
  lastp = 0.0
  sus2exp = as.numeric(probsVec <= p_expose_vec)*indicatorVector
  sus2immune = as.numeric((probsVec > p_expose_vec) & (probsVec <= (p_expose_vec + p_immune_scalar)))*indicatorVector
  sus2die = as.numeric((probsVec > (p_expose_vec + p_immune_scalar)) & (probsVec <= (p_expose_vec + p_immune_scalar + p_die_scalar)))*indicatorVector
  # if (any(sus2die == sus2die)) {
  #   breakpoint = 1
  # }
  list(exposed = sus2exp, immune = sus2immune, dead = sus2die)
}

# From V state
rStateFromVaccinated <- function(p_exposeVac_vec, p_dieVac_scalar, indicatorVector) {
  probsVec = runif(length(indicatorVector))
  sump = 0.0
  lastp = 0.0
  Vac2exp = as.numeric(probsVec <= p_exposeVac_vec)*indicatorVector
  Vac2die = as.numeric((probsVec > p_exposeVac_vec) & (probsVec <= (p_exposeVac_vec + p_dieVac_scalar)))*indicatorVector
  list(exposed_Vac = Vac2exp, dead_Vac = Vac2die)
}


##################################################################################


NetworkSELI <- function(CM, tau,lambda,mu,delta, sigma, gamma, omega,p, tau_vac,
                        lambda_vac,lambda_vacbirth,delta_tx, alpha) {  
  #generate SEIR epidemic on a CM-network  
  #CM = contact matrix 
  #tau = probability of infection across an edge 
  #delta = probability of exposed across an edge
  #gamma = probability of removal per time step 
  #set.seed(3)
  N=dim(CM)[1] 
  S=matrix(rep(1,N),nrow=N,ncol=1) #First susceptibles
  E=matrix(rep(0,N),nrow=N,ncol=1) #Early Exposed
  L=matrix(rep(0,N),nrow=N,ncol=1) #Later Exposed
  I=matrix(rep(0,N),nrow=N,ncol=1) #First infected
  V=matrix(rep(0,N),nrow=N,ncol=1) #First immunized
  I1=sample(1:N, size=10)           #Pick first random infected  
  I[I1,1]=1  
  S[I1,1]=0  
  t=1 
  daysToRun <- 700
  currentTau = tau
  currentDelta = delta
  currentLambda = lambda
  #currentKappa = kappa
 # currentZeta = zeta
  # sum(I[,t])
  #while ((sum(I[, ncol(I)], na.rm = T)>0) | (t==1)) 
  while (t < daysToRun){
    #n <- t - 1
    n = t
    #print(t)
    # if (t > 100){currentTau <- tau_vac
    # currentDelta <- delta_tx
    #   }
    
    if (t >=400){currentTau <- tau_vac
    currentDelta <- delta_tx
    currentLambda <- lambda_vac
    #currentKappa <- kappa_vac
    }
    
    if (t >=401){
      currentLambda <- lambda_vacbirth
    }
    
    infneigh=CM %*% I[, n]
    #print(infneigh)
    pinf = 1-(1-currentTau)^infneigh
    if (any(is.na(pinf))) {
      print("pinf or pinf_2 have NA(s)")
      break
    }
    
    infneigh=CM %*% I[, n]
    #print(infneigh)
    pinf_V = 1-(1-alpha)^infneigh
    if (any(is.na(pinf_V))) {
      print("pinf_V have NA(s)")
      break
    }
    
    dependentNextS = rStateFromSusceptable(pinf, currentLambda, mu,S[,n])#rxnom(c(pinf,mu,1-pinf-mu),S[, n])
    newE =  dependentNextS[["exposed"]]
    newDS = dependentNextS[["dead"]]
    newIM = dependentNextS[["immune"]]
    
    #vaccination using vaccine with X% efficacy 
    dependentNextV = rStateFromVaccinated(pinf_V, mu ,V[, n])
    newVE = dependentNextV[["exposed_Vac"]] #new exposed from vaccinated state
    newDV = dependentNextV[["dead_Vac"]] #new natural death from vaccinated state
    
    dependentNextE = rxnom(c((sigma*(1-p)),(sigma*p),mu,
                             1-(sigma*(1-p))-(sigma*p)-mu),E[, n])
    newEL = dependentNextE[,1]#rbinom(N, E[,n], (1-p))
    newEI = dependentNextE[,2]#rbinom(N, E[,n], p)
    newDE =  dependentNextE[,3]#rbinom(N, E[,n], mu)
    
    
    dependentNextL = rxnom(c(gamma,mu,1-gamma-mu),L[, n])
    newI = dependentNextL[,1]#rbinom(N, L[, n], pinf)
    newDL =  dependentNextL[,2]#rbinom(N, L[,n], mu)
    
    
    dependentNextI = rxnom(c(omega,currentDelta,1-omega-currentDelta),I[, n])
    newDI =  dependentNextI[,1]#rbinom(N, I[,n], omega)
    newR =   dependentNextI[,2]#rbinom(N,I[,n], delta)
    
    # #vaccination using vaccine with 100% efficacy
    # dependentNextV = rxnom(c(mu),V[, n])#100% vaccine efficacy
    # newDV = dependentNextV[,1]
    
    #vaccination using vaccine with X% efficacy 
    # dependentNextV = rxnom(c((pinf_V*alpha),mu,1-(pinf_V*alpha)-mu),V[, n])
    # newVE = dependentNextV[,1] #new exposed from vaccinated state
    # newDV = dependentNextV[,2] #new natural death from vaccinated state
    
    
    nextS = S[, n] - newE - newIM - newDS + newDS  + newR + ((1-currentLambda)*(newDE + newDL+ newDI + newDV))
    nextE = E[, n] + newE + newVE - newEL - newEI -newDE 
    nextL = L[, n] + newEL - newI - newDL
    nextI = I[, n]  + newI + newEI - newDI - newR
    nextV = V[, n] + newIM - newVE - newDV + (currentLambda*(newDE+newDL+newDI+newDV))
    if (any(nextS <0))
    {
      breakpoint = 1
    }
    
    S=cbind(S, nextS)
    E=cbind(E, nextE)
    L=cbind(L, nextL)
    I=cbind(I, nextI)
    V=cbind(V, nextV)
    
    t = t+1 
    
    # print(ncol(I))
    # if (ncol(I) == 100) {
    #   break
    # }
  }  
  res=list(S=S,E=E,L=L,I=I, V=V)  
  class(res)="netSELI"  
  return(res)  
} 

#######################################################################
#Without any intervention
set.seed(2020)
taus = c(0.048, 0.028, 0.014, 0.0096) #rate of getting infected 
Ks =   c(6,     10,    20,    30) #number of neighbors
N = 100 #Population size
n_sims = 1 #number of simulations to run
n_data_points = n_sims * length(Ks) #number of data points (n_sims data point per model)
all_sims_w = data.frame(k=numeric(n_data_points),n_inf=numeric(n_data_points)) # number of infected per model per simulation #
for (i_sim in 1:n_sims) {
  for (iK in 1:length(Ks)) {
    k = Ks[iK]
    tau = taus[iK]
    cm_c = BarabasiAlbert(N,k)
    cepW = NetworkSELI(CM = cm_c, tau = tau, tau_vac = tau,
                       mu = 1-exp(-0.0008*1), delta = 1-exp(-0.02*1),lambda = 1-exp(-0*1),
                       lambda_vac = 1-exp(-0.6*1),lambda_vacbirth = 1-exp(-0*1),sigma = 1-exp(-0.2232*1),
                       gamma = 1-exp(-0.2232*1), omega = 1-exp(-0.69315*1), p = 1-exp(-0.07796*1), 
                       delta_tx = 1-exp(-0.02*1))
    result = apply(cepW$I, 2, sum, na.rm = T)[[700]]
    result_index = ( ( i_sim - 1 ) * length(Ks) ) + iK
    all_sims_w$k[result_index] = k
    all_sims_w$n_inf[result_index] = result
    print(result_index)
  }
}

#######################################################################
#Vaccine intervention
##loading network with ten thousands nodes and 10 neighbors
cm_X = readRDS('C:/Users/milalm01/Dropbox/Newyork_University/TB_Project/cm_X.rds');
##loading network with hundred thousands nodes and 10 neighbors
#cm_C_10 = readRDS('C:/Users/milalm01/Dropbox/Newyork_University/TB_Project/cm_C_10.rds');
#mean(apply(cm_C_10,2,sum)) #computes average number of neighbors per node
#mean(rowSums(cm_C_10))#works for large sparse matrices
#min(rowSums(cm_C_10))# minimum number of neighbors
#max(rowSums(cm_C_10))#max number of neighbors
cepV_85 = NetworkSELI(CM = cm_C_10, tau = 0.028, tau_vac = 0.028,alpha = 0.0042,
                   mu = 1-exp(-0.0008*1), delta = 1-exp(-0.02*1),lambda = 1-exp(-0*1),
                   lambda_vac = 1-exp(-0.6*1),lambda_vacbirth = 1-exp(-0.6*1*0.0012),
                   sigma = 1-exp(-0.2*1),gamma = 1-exp(-0.2*1), omega = 1-exp(-0.5*1),
                   p = 1-exp(-0.075*1),delta_tx = 1-exp(-0.02*1))
cepV = cepV_85
smV = summary(cepV)
inf_V=ifelse(apply(cepV$I,2,sum)>0,2,1)
#nwt_V=network(cm_X, directed=FALSE)
#plot(nwt_V, vertex.col=inf_V)
top_y_val_V <- max(apply(cepV$S,2,sum)) * 1.1
top_x_val_V <- ncol(cepV$S) * 1.1


set.seed(2020)
vac_intervention <- function(lambda_vac, lambda_vacbirth){
  taus = c(0.048, 0.028, 0.014, 0.0096) #rate of getting infected 
  #alphas = c(0.01632, 0.00952, 0.00476, 0.003264)# vaccine with 66% efficacy
  alphas = c(0.0024, 0.0014, 0.0007, 0.00048) # vaccine with 95% efficacy
  Ks =   c(6, 10, 20, 30) #number of neighbors
  N = 100000
  n_sims = 2 #number of simulations to run
  n_data_points = n_sims * length(Ks) #number of data points (n_sims data point per model)
  all_sims_vac = data.frame(k=numeric(n_data_points),n_inf=numeric(n_data_points), p_size=numeric(n_data_points)) # number of infected per model per simulation #
  for (i_sim in 1:n_sims) {
    for (iK in 1:length(Ks)) {
      # random seed here = everything will be same
      k = Ks[iK]
      tau = taus[iK]
      alpha = alphas[iK]
      cm_c = BarabasiAlbert(N,k)
      cepV = NetworkSELI(CM = cm_c, tau = tau, tau_vac = tau,alpha = alpha,
                         mu = 1-exp(-0.0008*1), delta = 1-exp(-0.02*1),lambda = 1-exp(-0*1),
                         lambda_vac = lambda_vac,lambda_vacbirth = lambda_vacbirth,sigma = 1-exp(-0.2*1),
                         gamma = 1-exp(-0.2*1), omega = 1-exp(-0.5*1), p = 1-exp(-0.075*1), delta_tx = 1-exp(-0.02*1))
      result_index = ( ( i_sim - 1 ) * length(Ks) ) + iK
      all_sims_vac$k[result_index] = k
      all_sims_vac$n_inf_700[result_index] = apply(cepV$I, 2, sum, na.rm = T)[[700]]
      all_sims_vac$n_inf_390[result_index] = apply(cepV$I, 2, sum, na.rm = T)[[390]]
      all_sims_vac$n_inf_400[result_index] = apply(cepV$I, 2, sum, na.rm = T)[[400]]
      all_sims_vac$n_inf_401[result_index] = apply(cepV$I, 2, sum, na.rm = T)[[401]]
      
      all_sims_vac$p_size[result_index] = lambda_vac
      print(result_index)
    }
  }
return(all_sims_vac)
}
all_sims_result_vac <- list()
lambda_vac_r = c(0.6, 0.92, 1.9, 3.0)
#lambda_vac_r = c(3.0)
lambda_vacs = c(1 - exp(-lambda_vac_r*1))
lambda_vacbirths = c(1 - exp(-lambda_vac_r*1*0.0012))
for (i in 1:4){
all_sims_result_vac[[i]] <- vac_intervention(lambda_vacs[i], lambda_vacbirths[i])
}

#saveRDS(all_sims_result_vac, file="all_sims_result_vac_P_sizes_0.45_0.6_0.85.rds")


#######################################################################
#Treatment intervention
cepT_4 = NetworkSELI(CM = cm_C_10, tau = 0.028, tau_vac = 0.028,alpha = 0,
                   mu = 1-exp(-0.0008*1), delta = 1-exp(-0.02*1),lambda = 1-exp(-0*1),
                   lambda_vac = 1-exp(-0*1),lambda_vacbirth = 1-exp(-0*1),sigma = 1-exp(-0.2*1),
                   gamma = 1-exp(-0.2*1), omega = 1-exp(-0.5*1), p = 1-exp(-0.075*1),
                   delta_tx = 0.4)

cepT = cepT_4
smT = summary(cepT)
inf_T=ifelse(apply(cepT$I,2,sum)>0,2,1)
#nwt_V=network(cm_X, directed=FALSE)
#plot(nwt_V, vertex.col=inf_V)
top_y_val_T <- max(apply(cepT$S,2,sum)) * 1.1
top_x_val_T <- ncol(cepT$S) * 1.1



set.seed(2020)
treatment_intervention <- function(delta_tx){
  taus = c(0.048, 0.028, 0.014, 0.0096) #rate of getting infected 
  Ks =   c(6, 10, 20, 30) #number of neighbors
  N = 100000
  n_sims = 2 #number of simulations to run
  n_data_points = n_sims * length(Ks) #number of data points (n_sims data point per model)
  all_sims_tx = data.frame(k=numeric(n_data_points),n_inf=numeric(n_data_points), rec_rate=numeric(n_data_points)) # number of infected per model per simulation #
  for (i_sim in 1:n_sims) {
    for (iK in 1:length(Ks)) {
      # random seed here = everything will be same
      k = Ks[iK]
      tau = taus[iK]
      cm_c = BarabasiAlbert(N,k)
      cepT = NetworkSELI(CM = cm_c, tau = tau, tau_vac = tau,alpha = 0,
                         mu = 1-exp(-0.0008*1), delta = 1-exp(-0.02*1),lambda = 1-exp(-0*1),
                         lambda_vac = 1-exp(-0*1),lambda_vacbirth = 1-exp(-0*1),sigma = 1-exp(-0.22*1),
                         gamma = 1-exp(-0.22*1), omega = 1-exp(-0.69315*1), p = 1-exp(-0.07796*1), 
                         delta_tx = delta_tx)
      result = apply(cepT$I, 2, sum, na.rm = T)[[700]]
      result_index = ( ( i_sim - 1 ) * length(Ks) ) + iK
      all_sims_tx$k[result_index] = k
      all_sims_tx$n_inf[result_index] = result
      all_sims_tx$p_size[result_index] = delta_tx
      print(result_index)
    }
  }
  return(all_sims_tx)
}
all_sims_result_tx <- list()
#delta_txs = c(0.2,0.3,0.4,0.5)
delta_txs = c(0.4,0.5)
for (i in 1:4){
  all_sims_result_tx[[i]] <- treatment_intervention(delta_txs[i])
}



###################################################################
#All to all mixing model
seir_ode<-function(t,Y,par){
  S<-Y[1]
  E<-Y[2]
  L<-Y[3]
  I<-Y[4]
  V<-Y[5]
  
  beta<-par[1]
  sigma<-par[2]
  gamma<-par[3]
  mu<-par[4]
  p<-par[5]
  omega<-par[6]
  delta<-par[7]
  lambda<-par[8]
  
  if (t >= 400){
    beta <- par[9]
    delta <- par[10]
    lambda <-par[11]
  }
  
  if (t >=401){
    lambda <-par[12]
  } 
  
  dYdt<-vector(length=4)
  dYdt[1]=-beta*I*S - lambda*S - mu*S + mu*S + delta*I + (1-lambda)*(mu*E + mu*L + omega*I + mu*V) # Susceptible
  dYdt[2]=beta*I*S - (1-p)*sigma*E - p*sigma*E - mu*E #Early latent
  dYdt[3]=(1-p)*sigma*E - gamma*L - mu*L #Later latent
  dYdt[4]=p*sigma*E + gamma*L - omega*I - delta*I #Active infected
  dYdt[5]=lambda*S - mu*V + lambda*(mu*E + mu*L + omega*I + mu*V)#Vaccinated
  return(list(dYdt))
}

# Set parameter values
beta<-0.00000767
#beta<-0.00000689 # for N = 100k, 6, 10, 20 and 30 nodes when Tau is pertubed
sigma<-0.2
gamma<-0.2
mu<-0.0008
p<-0.075
#omega<-0.54
omega<-0.5
delta<-0.02
lambda<-0

time<-seq(1,700)

#Initial values for sub-populations
X<-99990       # susceptible hosts
W<-0           # Early latent exposed
Y<-0           # later latent exposed
Z<-10          # Active TB
U<-0           # immune

init<-c(S<-X, E<-W, L<-Y, I<-Z, V<-U)


##############################################################################
#Interventions
#1) Without any Treatment
par<-c(beta,sigma,gamma,mu,p,omega,delta,lambda,0.00000689,0.02,0,0)
# Solve system using lsoda
sol_W<-lsoda(init,time,seir_ode,par)
sol_W=as.data.frame(sol_W) 
head(round(sol_W, 4))

# #Plot
# par(mfrow=c(1,1))
# plot(time,sol_W[,2],type="l",col="black",ylab="Population size",xlab="Time (days)",
#      ylim=c(0,100000) ,xlim=c(0,800))
# lines(time,sol_W[,3],col="orange")
# lines(time,sol_W[,4],col="brown") 
# lines(time,sol_W[,5],col="red")
# lines(time,sol_W[,6],col="blue")
# legend("topright", legend=c("S","E","L","I","V"), lty=c(1,1,1,1,1),
#        col=c("black", "orange", "brown","red","blue"), bty = "n")

###############################################################
#With drug Treatment improving recovery rate by 20 times 
par<-c(beta,sigma,gamma,mu,p,omega,delta,lambda,0.00000767,0.5,0,0)
sol_T<-lsoda(init,time,seir_ode,par)
sol_T=as.data.frame(sol_T) 
head(round(sol_T, 4))

# #Plot
# par(mfrow=c(2,2))
# plot(time,sol_T[,2],type="l",col="black",ylab="Population size", 
#      ylim=c(0,100000) ,xlim=c(0,800))
# lines(time,sol_T[,3],col="orange")
# lines(time,sol_T[,4],col="brown") 
# lines(time,sol_T[,5],col="red")
# lines(time,sol_T[,6],col="blue")
# legend("topright", legend=c("S","E","L","I","V"), lty=c(1,1,1,1,1),
#       col=c("black", "orange", "brown","red","blue"))

###############################################################
# With Vaccine intervention with 100% efficacy moving one third of susceptibles to immune state
par<-c(beta,sigma,gamma,mu,p,omega,delta,lambda,0.00000767,0.02,3.0,(0.0012*3.0))

#With Vaccine intervention with 66% efficacy decreasing rate of exposure by one third
par<-c(beta,sigma,gamma,mu,p,omega,delta,lambda,0.0000026,0.02,0,0)

# Solve system using lsoda
sol_V<-lsoda(init,time,seir_ode,par)
sol_V=as.data.frame(sol_V) 
head(round(sol_V, 4))

# #Plot
# #par(mfrow=c(2,2))
# plot(time,sol_V[,2],type="l",col="black",lty = 2, ylab="Population size",xlab="Time (days)",
#      ylim=c(0,11000),xlim=c(0,330))
# lines(time,sol_V[,3],col="orange", lty = 2)
# lines(time,sol_V[,4],col="brown", lty = 2) 
# lines(time,sol_V[,5],col="red", lty = 2)
# legend("topright", legend=c("SV","EV","LV","IV"), lty=c(2,2,2,2),
#        col=c("black", "orange", "brown","red"), bty = "n")

#################################################################
# With Vaccine intervention making 60% of the population immune and tx improve recovery 
# by 20 times
par<-c(beta,sigma,gamma,mu,p,omega,delta,lambda,0.00000766,0.4,1.9,(0.0012*1.9))
# Solve system using lsoda
sol_VT<-lsoda(init,time,seir_ode,par)
sol_VT=as.data.frame(sol_VT) 
head(round(sol_VT, 4))

# #Plot
# #par(mfrow=c(2,2))
# plot(time,sol_VT[,2],type="l",col="black",ylab="Numbers", ylim=c(0,13000)
#      ,xlim=c(0,330))
# lines(time,sol_VT[,3],col="yellow")
# lines(time,sol_VT[,4],col="orange") 
# lines(time,sol_VT[,5],col="red")
# #legend("topright", legend=c("S","E","L","I"), lty=c(1,1,1,1),
# #      col=c("black", "yellow", "orange","red"))

####################################################################

#Figure all graphs including vaccinated 
par(mfrow=c(2,2))
plot(apply(cepW$S,2, sum),type="l",col="black", xlab="Time (days)", ylab="Population size",
     xlim = c(0,top_x_val_W), ylim=c(0,top_y_val_W))
lines(apply(cepW$E, 2, sum),type="l", col="orange")
lines(apply(cepW$L, 2, sum), type="l", col="brown")
lines(apply(cepW$I, 2, sum), type="l", col="red")
lines(apply(cepW$V, 2, sum), type="l", col="blue")
lines(time,sol_W[,2],type="l", col="black", lty = 2)
lines(time,sol_W[,3],type="l", col="orange", lty = 2)
lines(time,sol_W[,4],type="l", col="brown", lty = 2)
lines(time,sol_W[,5],type="l", col="red", lty = 2)
lines(time,sol_W[,6],type="l", col="blue", lty =2)
#legend("topright", legend=c("SV_NB","SV", "Ig_NB", "Ig"), lty=c(1,2,1,2),
# col=c("black", "black", "blue", "blue"), bty = "n")
legend("bottomright", legend=c("SV_NB","SV"), lty=c(1,2),
       col=c("black", "black"), bty = "n")

par(mfrow=c(1,2))
plot(apply(cepV$S,2, sum),type="l",col="black", xlab="Time (days)", ylab="Population size",
     xlim = c(0,top_x_val_V), ylim=c(0,top_y_val_V))
lines(apply(cepV$E, 2, sum),type="l", col="orange")
lines(apply(cepV$L, 2, sum), type="l", col="brown")
lines(apply(cepV$I, 2, sum), type="l", col="red")
lines(apply(cepV$V, 2, sum), type="l", col="blue")
lines(time,sol_V[,2],type="l", col="black", lty = 2)
lines(time,sol_V[,3],type="l", col="orange", lty = 2)
lines(time,sol_V[,4],type="l", col="brown", lty = 2)
lines(time,sol_V[,5],type="l", col="red", lty = 2)
lines(time,sol_V[,6],type="l", col="blue", lty =2)
#legend("topright", legend=c("SV_NB","SV", "Ig_NB", "Ig"), lty=c(1,2,1,2),
# col=c("black", "black", "blue", "blue"), bty = "n")
legend("bottomright", legend=c("SV_NB","SV"), lty=c(1,2),
       col=c("black", "black"), bty = "n")
legend("topright", legend=c("100 efficacy", "95% efficacy","85% efficacy","66% efficacy"), lty=c(1,1,1,1),
        col=c("black", "brown","orange","blue", ), bty = "n")

par(mfrow=c(1,2))
plot(apply(cepT$S,2, sum),type="l",col="black", xlab="Time (days)", ylab="Population size",
     xlim = c(0,top_x_val_T), ylim=c(0,top_y_val_T))
lines(apply(cepT$E, 2, sum),type="l", col="orange")
lines(apply(cepT$L, 2, sum), type="l", col="brown")
lines(apply(cepT$I, 2, sum), type="l", col="red")
lines(time,sol_T[,2],type="l", col="black", lty = 2)
lines(time,sol_T[,3],type="l", col="orange", lty = 2)
lines(time,sol_T[,4],type="l", col="brown", lty = 2)
lines(time,sol_T[,5],type="l", col="red", lty = 2)
#legend("topright", legend=c("SV_NB","SV", "Ig_NB", "Ig"), lty=c(1,2,1,2),
# col=c("black", "black", "blue", "blue"), bty = "n")



#####################################################################
#Fig_1
par(mfrow=c(2,2))
plot(apply(cepW$S,2, sum),type="l", xlab="Time (days)", ylab="Population size",
     xlim = c(0,top_x_val_W), ylim=c(0,top_y_val_W))
#lines(apply(cepV$S, 2, sum, na.rm = T), type="l", col="black")
lines(time,sol_W[,2],type="l", col="black", lty = 2)
#lines(time,sol_V[,2],type="l", col="blue", lty =2)
#legend("topright", legend=c("SV_NB","SV", "Ig_NB", "Ig"), lty=c(1,2,1,2),
#col=c("black", "black", "blue", "blue"), bty = "n")
legend("topright", legend=c("SW_NB","SW"), lty=c(1,2),
       col=c("black", "black"), bty = "n")

plot(apply(cepW$E,2, sum),type="l",col="orange", xlab="Time (days)", ylab="Population size",
     xlim = c(0,top_x_val_W), ylim=c(0,top_y_val_W))
lines(time,sol_W[,3],type="l", col="orange", lty = 2)
legend("topright", legend=c("EW_NB","EW"), lty=c(1,2),
       col=c("orange", "orange"), bty = "n")

plot(apply(cepW$L,2, sum),type="l",col="brown", xlab="Time (days)", ylab="Population size",
     xlim = c(0,top_x_val_W), ylim=c(0,top_y_val_W))
lines(time,sol_W[,4],type="l", col="brown", lty = 2)
legend("topright", legend=c("LW_NB","LW"), lty=c(1,2),
       col=c("brown", "brown"), bty = "n")

plot(apply(cepW$I,2, sum),type="l",col="red", xlab="Time (days)", ylab="Population size",
     xlim = c(0,top_x_val_W), ylim=c(0,top_y_val_W))
lines(time,sol_W[,5],type="l", col="red", lty = 2)
legend("topright", legend=c("IW_NB","IW"), lty=c(1,2),
       col=c("red", "red"), bty = "n")

##################################################################################
#Fig_2 & 3
par(mfrow=c(2,2))
plot(apply(cepV$S,2, sum),type="l",col="black", xlab="Time (days)", ylab="Population size",
     xlim = c(0,top_x_val_V), ylim=c(0,top_y_val_V))
#lines(apply(cepV$V, 2, sum, na.rm = T), type="l", col="blue")
lines(time,sol_V[,2],type="l", col="black", lty = 2)
#lines(time,sol_V[,6],type="l", col="blue", lty =2)
#legend("topright", legend=c("SV_NB","SV", "Ig_NB", "Ig"), lty=c(1,2,1,2),
# col=c("black", "black", "blue", "blue"), bty = "n")
legend("bottomright", legend=c("SV_NB","SV"), lty=c(1,2),
       col=c("black", "black"), bty = "n")
#legend("topright", legend=c("Ig_NB", "Ig"), lty=c(1,2),
#col=c("blue", "blue"), bty = "n")

plot (apply(cepV$E,2, sum),type="l", xlab="Time (days)", col = 'orange', ylab="Population size",
      xlim=c(0,top_x_val_V), ylim=c(0, top_y_val_V)) 
#lines(apply(cepW$E,2,sum, na.rm = T), type="l",col= 'orange', lty = 2)
lines(time,sol_V[,3],type="l", col="orange", lty = 2)
#lines(time,sol_W[,3],type="l", col="green", lty =2)
legend("topright", legend=c("EV_NB", "EV"), lty=c(1,2),
       col=c("orange", "orange"), bty = "n")
plot (apply(cepV$L,2, sum),type="l", xlab="Time (days)", col = 'brown', ylab="Population size",
      xlim=c(0,top_x_val_V), ylim=c(0, top_y_val_V)) 
#lines(apply(cepW$L,2,sum, na.rm = T), type="l",col= 'orange', lty = 2)
lines(time,sol_V[,4],type="l", col="brown", lty = 2)
#lines(time,sol_W[,4],type="l", col="purple", lty =2)
legend("topright", legend=c("LV_NB", "LV"), lty=c(1,2),
       col=c("brown", "brown"), bty = "n")
plot (apply(cepV$I,2, sum),type="l", xlab="Time (days)", col = 'red', ylab="Population size",
      xlim=c(0,top_x_val_V), ylim=c(0,top_y_val_V))
#lines(apply(cepW$I,2,sum, na.rm = T), type="l",col= 'red', lty = 2)
lines(time,sol_V[,5],type="l", col="red", lty = 2)
#lines(time,sol_W[,5],type="l", col="pink", lty =2)
legend("topright", legend=c("IV_NB", "IV"), lty=c(1,2),
       col=c("red", "red"),bty = "n")

####################################################################
#Fig_4
par(mfrow=c(1,2))
plot(apply(cepT$S,2, sum),type="l",col="black", xlab="Time (days)", ylab="Population size",
     xlim = c(0,top_x_val_T), ylim=c(0,top_y_val_T))
#lines(apply(cepV$V, 2, sum, na.rm = T), type="l", col="blue")
lines(time,sol_T[,2],type="l", col="black", lty = 2)
#lines(time,sol_V[,6],type="l", col="blue", lty =2)
#legend("topright", legend=c("SV_NB","SV", "Ig_NB", "Ig"), lty=c(1,2,1,2),
# col=c("black", "black", "blue", "blue"), bty = "n")
legend("bottomright", legend=c("ST_NB","ST"), lty=c(1,2),
       col=c("black", "black"), bty = "n")
#legend("topright", legend=c("Ig_NB", "Ig"), lty=c(1,2),
#col=c("blue", "blue"), bty = "n")

plot (apply(cepT$E,2, sum),type="l", xlab="Time (days)", col = 'orange', ylab="Population size",
      xlim=c(0,top_x_val_T), ylim=c(0, top_y_val_T)) 
#lines(apply(cepW$E,2,sum, na.rm = T), type="l",col= 'orange', lty = 2)
lines(time,sol_T[,3],type="l", col="orange", lty = 2)
#lines(time,sol_W[,3],type="l", col="green", lty =2)
legend("topright", legend=c("ET_NB", "ET"), lty=c(1,2),
       col=c("orange", "orange"), bty = "n")
plot (apply(cepT$L,2, sum),type="l", xlab="Time (days)", col = 'brown', ylab="Population size",
      xlim=c(0,top_x_val_T), ylim=c(0, top_y_val_T)) 
#lines(apply(cepW$L,2,sum, na.rm = T), type="l",col= 'orange', lty = 2)
lines(time,sol_T[,4],type="l", col="brown", lty = 2)
#lines(time,sol_W[,4],type="l", col="purple", lty =2)
legend("topright", legend=c("LT_NB", "LT"), lty=c(1,2),
       col=c("brown", "brown"), bty = "n")
plot (apply(cepT$I,2, sum),type="l", xlab="Time (days)", col = 'red', ylab="Population size",
      xlim=c(0,top_x_val_T), ylim=c(0,top_y_val_T))
#lines(apply(cepW$I,2,sum, na.rm = T), type="l",col= 'red', lty = 2)
lines(time,sol_T[,5],type="l", col="red", lty = 2)
#lines(time,sol_W[,5],type="l", col="pink", lty =2)
legend("topright", legend=c("IT_NB", "IT"), lty=c(1,2),
       col=c("red", "red"),bty = "n")

######################################################################
#Fig_5
par(mfrow=c(2,2))
plot(apply(cepVT$S,2, sum),type="l",col="black", xlab="Time (days)", ylab="Population size",
     xlim = c(0,top_x_val_VT), ylim=c(0,top_y_val_VT))
#lines(apply(cepV$V, 2, sum, na.rm = T), type="l", col="blue")
lines(time,sol_VT[,2],type="l", col="black", lty = 2)
#lines(time,sol_V[,6],type="l", col="blue", lty =2)
#legend("topright", legend=c("SV_NB","SV", "Ig_NB", "Ig"), lty=c(1,2,1,2),
# col=c("black", "black", "blue", "blue"), bty = "n")
legend("topright", legend=c("SVT_NB","SVT"), lty=c(1,2),
       col=c("black", "black"), bty = "n")
#legend("topright", legend=c("Ig_NB", "Ig"), lty=c(1,2),
#col=c("blue", "blue"), bty = "n")

plot (apply(cepVT$E,2, sum),type="l", xlab="Time (days)", col = 'orange', ylab="Population size",
      xlim=c(0,top_x_val_VT), ylim=c(0, top_y_val_VT)) 
#lines(apply(cepW$E,2,sum, na.rm = T), type="l",col= 'orange', lty = 2)
lines(time,sol_VT[,3],type="l", col="orange", lty = 2)
#lines(time,sol_W[,3],type="l", col="green", lty =2)
legend("topright", legend=c("EVT_NB", "EVT"), lty=c(1,2),
       col=c("orange", "orange"), bty = "n")
plot (apply(cepVT$L,2, sum),type="l", xlab="Time (days)", col = 'brown', ylab="Population size",
      xlim=c(0,top_x_val_VT), ylim=c(0, top_y_val_VT)) 
#lines(apply(cepW$L,2,sum, na.rm = T), type="l",col= 'orange', lty = 2)
lines(time,sol_VT[,4],type="l", col="brown", lty = 2)
#lines(time,sol_W[,4],type="l", col="purple", lty =2)
legend("topright", legend=c("LVT_NB", "LVT"), lty=c(1,2),
       col=c("brown", "brown"), bty = "n")
plot (apply(cepVT$I,2, sum),type="l", xlab="Time (days)", col = 'red', ylab="Population size",
      xlim=c(0,top_x_val_VT), ylim=c(0,top_y_val_VT))
#lines(apply(cepW$I,2,sum, na.rm = T), type="l",col= 'red', lty = 2)
lines(time,sol_VT[,5],type="l", col="red", lty = 2)
#lines(time,sol_W[,5],type="l", col="pink", lty =2)
legend("topright", legend=c("IVT_NB", "IVT"), lty=c(1,2),
       col=c("red", "red"),bty = "n")

###############################################################
