##################################################################################################
################################## CODE BAYES : AIR ##############################################
##################################################################################################

# Gibbs sampler car dim 3 
## log random walk symétrique autour de 0 donc proba d'acceptation = min(1, g(X¨*)X*/g(Xt))Xt

#### Données
"alpha" <- 4.48        
"beta" <- 0.76         
"sigma2" <- 81.14      
"J" <- 3               
"y" <-
  c(21, 20, 15)
"n" <-
  c(48, 34, 21)
"Z" <-
  c(10, 30, 50)


#### Modèle 
calcul.pj=function(theta1,theta2,xj){
  e=exp(theta1+theta2*xj)
  return(e/(1+e))
}


air <- function(Z,n,y,J,sigma2,beta,alpha, nchain = 10^5, prop_sd = c(2,0.01,0.01)){
  #initialisation des Xj
  chain.x= matrix(NA,nchain+1,J)
  init.x=alpha+beta*Z
  acc.rates.x=rep(0,J)
  chain.x[1,] <- init.x
  prop.sd.x=prop_sd[1]
  
  #initialisation des theta
  chain.theta <- matrix(NA, nchain + 1, 2)
  acc.rates.theta <- rep(0, 2)
  init.theta <- c( 0 , 0 ) #=> pj=1/2
  chain.theta[1,] <- init.theta
  prop.sd.theta1=prop_sd[2]
  prop.sd.theta2=prop_sd[3]
  
  for (iter in 1:nchain){
    current.x = chain.x[iter,]
    curent.theta = chain.theta[iter,]
    theta1=curent.theta[1]
    theta2=curent.theta[2]
    chain.x[iter+1,]
    
    #mise à jour de Xj
    for (j in 1:J){ #ici j=1,2,3
      X=current.x[j]
      
      mu=alpha+beta*Z[j] 
      prop=rnorm(1,X, prop.sd.x)
      pj=calcul.pj(theta1,theta2,X)
      pj.prop=calcul.pj(theta1,theta2,prop)
      
      top=dnorm(prop,mu,sqrt(sigma2),log=TRUE)+ dbinom(y[j],n[j],prob=pj.prop, log = TRUE)
      bottom=dnorm(X,mu,sqrt(sigma2),log=TRUE)+ dbinom(y[j],n[j],prob=pj, log = TRUE)
      
      acc_prob.x <- min(1,exp(top - bottom))
      if (runif(1) < acc_prob.x){
        current.x[j] <- prop
        acc.rates.x[j] <- acc.rates.x[j] + 1
        
      }
    }
    
    chain.x[iter+1,]=current.x
    
    #Mise à jour de theta1
    prop=rnorm(1,theta1,prop.sd.theta1) #K = gaussien => acceptation min(1,g(X*)/g(X))
    p=calcul.pj(theta1,theta2,current.x)
    p.prop=calcul.pj(prop,theta2,current.x)
    
    top=dnorm(prop,0,1/sqrt(0.001),log=TRUE) + sum(dbinom(y, n, prob=p.prop, log = TRUE))
    bottom=dnorm(theta1,0,1/sqrt(0.001),log=TRUE) + sum(dbinom(y, n, prob=p, log = TRUE))
    
    acc_proba.theta1=min(1,exp(top-bottom))
    if (runif(1) < acc_proba.theta1){
      curent.theta[1] <- prop
      acc.rates.theta[1] <- acc.rates.theta[1] + 1
    }
    
    chain.theta[iter+1,1]=curent.theta[1]
    
    #Mise à jour de theta2
    prop=rnorm(1,theta2,prop.sd.theta2) #K = gaussien => acceptation min(1,g(X*)/g(X))
    p=calcul.pj(theta1,theta2,current.x)
    p.prop=calcul.pj(theta1,prop,current.x)
    
    top=dnorm(prop,0,1/sqrt(0.001),log=TRUE)+ sum( dbinom(y ,n ,prob=p.prop, log = TRUE))
    bottom=dnorm(theta2,0,1/sqrt(0.001),log=TRUE)+ sum( dbinom(y, n, prob=p, log = TRUE))
    
    acc_proba.theta2=min(1,exp(top-bottom))
    
    if (runif(1) < acc_proba.theta2){
      curent.theta[2] <- prop
      acc.rates.theta[2] <- acc.rates.theta[2] + 1
    }
    
    chain.theta[iter+1,2]=curent.theta[2]
  }
  
  return(list(chain.x = chain.x, acc.rates.x = acc.rates.x / nchain, chain.theta=chain.theta,acc.rates.theta=acc.rates.theta/nchain))
}

out<-air(Z,n,y,J,sigma2,beta,alpha)


#on enlève la période de chauffe
out$chain.x<-out$chain.x[-(1:1000),]
out$chain.theta<-out$chain.theta[-(1:1000),]

#on vérifie nos résultats par rapport au tableau donné
mean(out$chain.theta[,1])  
sd(out$chain.theta[,1])   
mean(out$chain.theta[,2])  
sd(out$chain.theta[,2])
mean(out$chain.x[,1]) 
sd(out$chain.x[,1]) 
mean(out$chain.x[,2]) 
sd(out$chain.x[,2])
mean(out$chain.x[,3]) 
sd(out$chain.x[,3])

apply(out$chain.x,2,quantile, prob=c(0.025, 0.975))
apply(out$chain.theta,2,quantile, prob=c(0.025, 0.975))

#on regarde les acc.rates
out$acc.rates.x[1]
out$acc.rates.x[2]
out$acc.rates.x[3]
out$acc.rates.theta[1]
out$acc.rates.theta[2]

############################ Analyse

# chaines de Markov
par(mfrow = c(1, 3), mar = c(4, 5, 0.5, 0.5))
ylabs <- c(expression(x[0]), expression(x[1]), expression(x[2]))
for (j in 1:3){
  plot(out$chain.x[,j], type = "l", ylab = ylabs[j])
} 

par(mfrow = c(1, 2), mar = c(4, 5, 0.5, 0.5))
ylabs <- c(expression(theta[1]), expression(theta[2]))
for (j in 1:2)
  plot(out$chain.theta[,j], type = "l", ylab = ylabs[j])

par(mfrow = c(1, 3), mar = c(4, 5, 0.5, 0.5))
xlabs <- c(expression(x[0]), expression(x[1]), expression(x[2]))
for (j in 1:3)
  plot(density(out$chain.x[,j]), type = "l", xlab = xlabs[j], main = "")

plot(out$chain.theta) #thetas corrélés


x_hat <- apply(out$chain.x, 2, median)
theta_hat <- apply(out$chain.theta, 2, median)
proba <- as.numeric(exp(theta_hat[1]+theta_hat[2]*x_hat) / (1 + exp(theta_hat[1]+theta_hat[2]*x_hat)))
# Y ~ Bin(proba, n) 
rbind(predicted = proba*n, observed = y)

# Pour chaque categorie et chaque etat de la chaine, nous simulons
#le nombre de yes: malades respiratoires

simu <- matrix(NA, nrow(out$chain.x), 3)
for (t in 1:nrow(out$chain.x)){
  x <- out$chain.x[t,]
  theta <- out$chain.theta[t,]
  proba <- as.numeric(exp(theta[1]+theta[2]*x) / (1 + exp(theta[1]+theta[2]*x)))
  simu[t,] <- rbinom(3, n, proba)
}

par(mfrow = c(1, 3))
xlabs <- c("Nb de malades en catégorie 1",
           "Nb de malades en catégorie 2",
           "Nb de malades en catégorie 3")
for (j in 1:3){
  plot(table(simu[,j]) / nrow(simu),
       xlab = xlabs[j], xlim=c(0,45),ylim=c(0,0.15),ylab = "Probabilite")
}
