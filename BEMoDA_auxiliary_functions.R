# BEMoDA shiny v1.1 - BiowaivEr aid for Model Dependent-Independent Approach application for in-vitro dissolution profile comparison
# 
# Model Dependent-Independent Approach application for in-vitro dissolution profile comparison as proposed by Sathe et al. in 1996
# (Sathe PM, Tsong Y, Shah VP. In-vitro dissolution profile comparison: statistics and analysis, model dependent approach. Pharm Res. 1996 Dec;13(12):1799-803).
# 
# Copyright (C) 2017 Jakub Szlęk, Aleksander Mendyk
# 
# Authors: 
# Jakub Szlęk, Aleksander Mendyk
# 
# Affiliation: 
# Jagiellonian University Medical College,
# Faculty of Pharmacy,
# Department of Pharmaceucial Technology and Biopharmaceutics,
# Medyczna 9 st.,
# 30-688 Kraków
# Poland
# 
# Bugs, issues, please e-mail to maintainer
# Jakub Szlęk: j.szlek@uj.edu.pl
# 
# Copyright (C) 2017 Jakub Szlęk, Aleksander Mendyk
# 
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General 
# Public License as published by the Free Software Foundation, either version 3 of the License,
# or (at your option) any later version.
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with this program.
# If not, see <http://www.gnu.org/licenses/>.
# # 
# File contains auxiliary functions used by the script:
# BEMoDA_shiny_Dep.R
# 


#### Time points extraction from data table

extract.timepoints <- function(dataframe){

    time <- gsub("X","",colnames(dataframe))
    time <- gsub(".min","",time)
    time <- as.numeric(time)

    return(time)
}


#### To find out the points on the borded of critical region (CR)

draw.ellipse <- function(x, y){
#
# Estimate the parameters.
#
X <- as.data.frame(cbind(One=1, x, y, xy=x*y, x2=x^2, y2=y^2))
fit <- lm(One ~ . - 1, X)
beta.hat <- coef(fit)

#
# Plot the estimate, the point, and the original axes.
#
evaluate <- function(x, y, beta) {
  if (missing(y)) {
    y <- x[, 2]; x <- x[, 1]
  }
    as.vector(cbind(x, y, x*y, x^2, y^2) %*% beta - 1)
}
e.x <- diff(range(x)) / 40
e.y <- diff(range(y)) / 40
n.x <- 100
n.y <- 60
u <- seq(min(x)-e.x, max(x)+e.x, length.out=n.x)
v <- seq(min(y)-e.y, max(y)+e.y, length.out=n.y)
z <- matrix(evaluate(as.matrix(expand.grid(u, v)), beta=beta.hat), n.x)

return(list(u=u,v=v,z=z, beta.hat=beta.hat))
}

#### Get ellipse params from ellipse equation
# The equation of ellipse is:
# A*x^2 + B*xy + C*y^2 + D*x + E*y - 1 = 0

get.ellipse.params <- function(params){
  
  # params should be provided in form that
  # > sapply(myParams, class)
  # 
  # x         y        xy        x2        y2 
  # "numeric" "numeric" "numeric" "numeric" "numeric"
  
  A <- params[["x2"]]
  B <- params[["xy"]]
  C <- params[["y2"]]
  D <- params[["x"]]
  E <- params[["y"]]
  F <- -1
  
  M0 <- matrix(c(F,D/2,E/2, D/2, A, B/2, E/2, B/2, C), nrow=3, byrow=TRUE)
  M <- matrix(c(A,B/2,B/2,C), nrow=2)
  lambda <- eigen(M)$values
  
  # assuming abs(lambda[1] - A) < abs(lambda[1] - C), if not, swap lambda[1] and lambda[2] in the following equations:
  if(abs(lambda[1] - A) < abs(lambda[1] - C)){
    
    a <- sqrt(-det(M0)/(det(M)*lambda[1]))  
    b <- sqrt(-det(M0)/(det(M)*lambda[2]))
    
    } else if(abs(lambda[1] - A) >= abs(lambda[1] - C)){
    
    a <- sqrt(-det(M0)/(det(M)*lambda[2]))  
    b <- sqrt(-det(M0)/(det(M)*lambda[1]))
    
  }
  
  xc <- (B*E-2*C*D)/(4*A*C-B^2)
  yc <- (B*D-2*A*E)/(4*A*C-B^2)
  phi <- (pi/2 - atan((A-C)/B))/2
  
  return(list(major.axis.length=a,minor.axis.length=b,x.center=xc,y.center=yc,major.angle=phi))
}


#### Function to search the par[1] and par[2] via nloptr function for critical region (CR)

cr.ellipse <- function(par,S.cr,k.val,mean.diff,p,nr,level){

    f.val <- qf(level,p,2*nr-p-1)
    
    
# The part of an below equation 
# 
#   left.eq <- S[1,1]*(par[1]-mean.diff[1])^2 + S[1,2]*(par[1]-mean.diff[1])*S[2,1]*(par[2]-mean.diff[2]) + S[2,2]*(par[2]-mean.diff[2])^2    
# 
#  can be substituted by: 
# 
    left.eq <- t(par-mean.diff) %*% solve(S.cr) %*% (par-mean.diff)
   
    right.eq <- f.val / k.val # k.val <- ((nt+nr-p-1)/((nt+nr-2)*p))*((nt*nr)/(nt+nr)) 

    res <- mean((left.eq-right.eq)^2)
    
    if (is.na(res)){
        res<-Inf
    }
    
    return(res)
}


#### Function to search the par[1] and par[2] via nloptr function for similarity region (SR)

sr.ellipse <- function(par,S, k.val, p, nr, level){

    alpha <- par[1]
    beta <- par[2]
    f.val <- qf(level,p,2*nr-p-1)
    
# The part of an below equation 
# 
#   left.eq <- S[1,1]*(par[1]-mean.diff[1])^2 + S[1,2]*(par[1]-mean.diff[1])*S[2,1]*(par[2]-mean.diff[2]) + S[2,2]*(par[2]-mean.diff[2])^2    
# 
#  can be substituted by: 
# 
    left.eq <- t(par) %*% solve(S) %*% (par)
   
    right.eq <- (f.val / k.val)

    res <- mean((left.eq-right.eq)^2)
    
    if (is.na(res)){
        res<-Inf
    }
    
    
    return(res)
}


### Weibull model function to optimize

weibull <- function(params, t, X){
    
    alpha <- params[1]
    beta <- params[2]
    pred <- 1-exp(-alpha*t^(beta))
    res <- mean((X-pred)^2)
    
    if (is.na(res)){
        res<-Inf
    }
    
    return(res)
}


### Weibull model with fixed parameters of alpha and beta

weibull.fit <- function(t, alpha.calc,beta.calc){
    
    alpha <- alpha.calc
    beta <- beta.calc
    
    diss <- 1-exp(-alpha*t^(beta))
    
    return(diss)
}


### Two parameter functions

model.2params <- function(params, t, X, which.model){
  
  alpha <- params[1]
  beta <- params[2]
  
  which.model <- which.model
  
  if(which.model == "Weibull"){
    
    pred <- 100*(1-exp(-alpha*(t^(beta))))
    # y=100*(1-e^(-alpha*((t)^beta)))
    # Sathe PM, Tsong Y, Shah VP. In-vitro dissolution profile comparison: statistics and analysis, model dependent approach. Pharm Res. 1996;13:1799–803. doi: 10.1023/A:1016020822093.

  } else if (which.model == "Korsmeyer-Peppas"){
    
    pred <- alpha*t^beta
    pred <- pred * 100
    # y=k*t^n
    # Korsmeyer RW, Gurny R, Doelker E, Buri P, Peppas NA. Mechanisms of solute release from porous hydrophilic polymers. Int J Pharm. 1983;15:25–35. doi: 10.1016/0378-5173(83)90064-9
    
  } else if(which.model == "Peppas-Sahlin"){
    
    pred <- alpha*(t*60)^(0.5)+beta*(t*60)
    # y=k1*t^(0.5)+k2*t
    # Peppas NA, Sahlin JJ. A simple equation for the description of solute release III. Coupling of diffusion and relaxation. Int J Pharm. 1989;57:169–72. doi: 10.1016/0378-5173(89)90306-2.
    
  } else if (which.model=="Quadratic"){
    
    pred <- 100*(alpha*(t*60)^2+beta*(t*60))
    # y=100*(k1*t^2+k2*t)
    # Costa P, Sousa Lobo JM. Modeling and comparison of dissolution profiles. Eur J Pharm Sci. 2001;13:123–33. doi: 10.1016/S0928-0987(01)00095-1.
    
  } else if(which.model == "Logistic"){
    
    pred <- 100*((exp(alpha+beta*log10(t*60)))/(1+exp(alpha+beta*log10(t*60))))
    # y=100*((e^(alpha+beta*log(t)))/(1+e^(alpha+beta*log(t))))
    # Sathe PM, Tsong Y, Shah VP. In-vitro dissolution profile comparison: statistics and analysis, model dependent approach. Pharm Res. 1996;13:1799–803. doi: 10.1023/A:1016020822093.
    
  } else if(which.model == "Gompertz"){
    
    pred <- 100*exp(-alpha*exp(-beta*log10(t*60)))
    # y=100*e^(-alpha*e^(-beta*log(t)))
    # Tsong Y, Hammerstrom T, Chen JJ. Multipoint dissolution specification and acceptance sampling rule based on profile modeling and principal component analysis. J Biopharm Stat. 1997;7:423–39.
    
  # } else if (which.model=="Probit"){ 
  #   fi <- pnorm(time,mean=mean(time),sd=sd(time)) # should be normal distribution function of time, but it seems that smth is wrong here =========> switch-off
  #   pred <- 100*fi*(alpha+beta*log(t))
  #   # y = 100*fi*(alpha+beta*log(t))
  #   # Tsong Y, Hammerstrom T, Chen JJ. Multipoint dissolution specification and acceptance sampling rule based on profile modeling and principal component analysis. J Biopharm Stat. 1997;7:423–39.

      } else if (which.model == "Zero-order with Tlag"){
    
    pred <- alpha*(t*60 - beta)
    # y = alpha*(t - beta)
    # Borodkin S, Tucker FE. Linear drug release from laminated hydroxypropyl cellulose-polyvinyl acetate films. J Pharm Sci. 1975;64:1289–94.
  } else if(which.model == "Zero-order with F0"){
    
    pred <- alpha + beta*t*60
    # y = alpha + beta*t
    # Costa P, Sousa Lobo JM. Modeling and comparison of dissolution profiles. Eur J Pharm Sci. 2001;13:123–33.
  } else if(which.model == "First-order with Tlag"){
    
    pred <- 100*(1-exp(-alpha*(t*60-beta)))
    # y = 100*(1-exp(-alpha*(t-beta)))
    # Phaechamud T, Pitaksantayothin K, Kositwattanakoon P, Seehapong P, Jungvivatanavong S. Sustainable release of propranolol hydrochloride tablet using chitin as press-coating material. Silpakorn Univ Int J. 2002;2:147–59.
  } else if (which.model == "First-order with Fmax"){
    
    pred <- alpha*(1-exp(-beta*t*60))
    # y = alpha*(1-exp(-beta*t))
    # Tsong Y, Hammerstrom T, Chen JJ. Multipoint dissolution specification and acceptance sampling rule based on profile modeling and principal component analysis. J Biopharm Stat. 1997;7:423–39.
  } else if(which.model == "Higuchi with Tlag"){
    
    pred <- alpha*(t*60-beta)^0.5
    # y = alpha*(t-beta)^0.5
    # Tarvainen M, Peltonen S, Mikkonen H, Elovaara M, Tuunainen M, Paronen P et al. Aqueous starch acetate dispersion as a novel coating material for controlled release products. J Control Release. 2004;96:179–91.
  } else if(which.model == "Higuchi with F0"){
    
    pred <- alpha + beta*(t*60)^0.5
    # y = alpha + beta*t^0.5
    # Ford JL, Mitchell K, Rowe P, Armstrong DJ, Elliott PNC, Rostron C et al. Mathematical modelling of drug release from hydroxypropylmethylcellulose matrices: effect of temperature. Int J Pharm. 1991;71:95–104.
  } else if(which.model == "Hixon-Crowell with Tlag"){
    
    pred <- 100*(1-(1-alpha*(t*60-beta))^3)
    # y = 100*(1-(1-alpha*(t-beta))^3)
    # Mollo AR, Corrigan OI. An investigation of the mechanism of release of the amphoteric drug amoxycillin from poly( D , L -lactide-co-glycolide) matrices. Pharm Dev Technol. 2002;7:333–43.
  } else if(which.model == "Hopfenberg"){
    
    pred <- 100*(1-(1-alpha*t*60)^beta)
    # y = 100*(1-(1-alpha*t)^beta)
    # Enscore DJ, Hopfenberg HB, Stannett VT. Effect of particle size on the mechanism controlling n-hexane sorption in glassy polystyrene microspheres. Polymer. 1977;18:793–800.
  }
  
  # res <- sqrt(sum((pred-X)^2)/length(pred)) # RMSE
  res <- mean((X-pred)^2) # MSE
  
  
  if (is.na(res)){
    res<-Inf
  }
  
  closeAllConnections()
  return(res)
}


model.2params.fit <- function(t, alpha.calc, beta.calc, which.model){
  
  alpha <- alpha.calc
  beta <- beta.calc
  
  
  
  if(which.model == "Weibull"){
    
    diss <- 100*(1-exp(-alpha*(t^(beta))))
    # y=100*(1-e^(-alpha*((t)^beta)))
    # Sathe PM, Tsong Y, Shah VP. In-vitro dissolution profile comparison: statistics and analysis, model dependent approach. Pharm Res. 1996;13:1799–803. doi: 10.1023/A:1016020822093.
    
  } else if (which.model == "Korsmeyer-Peppas"){
    
    diss <- alpha*t^beta
    diss <- diss * 100
    # y=k*t^n
    # Korsmeyer RW, Gurny R, Doelker E, Buri P, Peppas NA. Mechanisms of solute release from porous hydrophilic polymers. Int J Pharm. 1983;15:25–35. doi: 10.1016/0378-5173(83)90064-9
    
  } else if(which.model == "Peppas-Sahlin"){
    
    diss <- (alpha*(t*60)^(0.5)+beta*(t*60))
    # y=k1*t^(0.5)+k2*t
    # Peppas NA, Sahlin JJ. A simple equation for the description of solute release III. Coupling of diffusion and relaxation. Int J Pharm. 1989;57:169–72. doi: 10.1016/0378-5173(89)90306-2.
    
  } else if (which.model=="Quadratic"){
    
    diss <- 100*(alpha*(t*60)^2+beta*(t*60))
    # y=100*(k1*t^2+k2*t)
    # Costa P, Sousa Lobo JM. Modeling and comparison of dissolution profiles. Eur J Pharm Sci. 2001;13:123–33. doi: 10.1016/S0928-0987(01)00095-1.
    
  } else if(which.model == "Logistic"){
    
    diss <- 100*((exp(alpha+beta*log10(t*60)))/(1+exp(alpha+beta*log10(t*60))))
    # y=100*((e^(alpha+beta*log(t)))/(1+e^(alpha+beta*log(t))))
    # Sathe PM, Tsong Y, Shah VP. In-vitro dissolution profile comparison: statistics and analysis, model dependent approach. Pharm Res. 1996;13:1799–803. doi: 10.1023/A:1016020822093.
    
  } else if(which.model == "Gompertz"){
    
    diss <- 100*exp(-alpha*exp(-beta*log10(t*60)))
    # y=100*e^(-alpha*e^(-beta*log(t)))
    # Tsong Y, Hammerstrom T, Chen JJ. Multipoint dissolution specification and acceptance sampling rule based on profile modeling and principal component analysis. J Biopharm Stat. 1997;7:423–39.
    
    # } else if (which.model=="Probit"){ 
    #   fi <- pnorm(time,mean=mean(time),sd=sd(time)) # should be normal distribution function of time, but it seems that smth is wrong here =========> switch-off
    #   diss <- 100*fi*(alpha+beta*log(t))
    #   # y = 100*fi*(alpha+beta*log(t))
    #   # Tsong Y, Hammerstrom T, Chen JJ. Multipoint dissolution specification and acceptance sampling rule based on profile modeling and principal component analysis. J Biopharm Stat. 1997;7:423–39.
  }  else if (which.model == "Zero-order with Tlag"){
    
    diss <- alpha*(t*60 - beta)
    # y = alpha*(t - beta)
    # Borodkin S, Tucker FE. Linear drug release from laminated hydroxypropyl cellulose-polyvinyl acetate films. J Pharm Sci. 1975;64:1289–94.
  } else if(which.model == "Zero-order with F0"){
    
    diss <- alpha + beta*t*60
    # y = alpha + beta*t
    # Costa P, Sousa Lobo JM. Modeling and comparison of dissolution profiles. Eur J Pharm Sci. 2001;13:123–33.
  } else if(which.model == "First-order with Tlag"){
    
    diss <- 100*(1-exp(-alpha*(t*60-beta)))
    # y = 100*(1-exp(-alpha*(t-beta)))
    # Phaechamud T, Pitaksantayothin K, Kositwattanakoon P, Seehapong P, Jungvivatanavong S. Sustainable release of propranolol hydrochloride tablet using chitin as press-coating material. Silpakorn Univ Int J. 2002;2:147–59.
  } else if (which.model == "First-order with Fmax"){
    
    diss <- alpha*(1-exp(-beta*t*60))
    # y = alpha*(1-exp(-beta*t))
    # Tsong Y, Hammerstrom T, Chen JJ. Multipoint dissolution specification and acceptance sampling rule based on profile modeling and principal component analysis. J Biopharm Stat. 1997;7:423–39.
  } else if(which.model == "Higuchi with Tlag"){
    
    diss <- alpha*(t*60-beta)^0.5
    # y = alpha*(t-beta)^0.5
    # Tarvainen M, Peltonen S, Mikkonen H, Elovaara M, Tuunainen M, Paronen P et al. Aqueous starch acetate dispersion as a novel coating material for controlled release products. J Control Release. 2004;96:179–91.
  } else if(which.model == "Higuchi with F0"){
    
    diss <- alpha + beta*(t*60)^0.5
    # y = alpha + beta*t^0.5
    # Ford JL, Mitchell K, Rowe P, Armstrong DJ, Elliott PNC, Rostron C et al. Mathematical modelling of drug release from hydroxypropylmethylcellulose matrices: effect of temperature. Int J Pharm. 1991;71:95–104.
  } else if(which.model == "Hixon-Crowell with Tlag"){
    
    diss <- 100*(1-(1-alpha*(t*60-beta))^3)
    # y = 100*(1-(1-alpha*(t-beta))^3)
    # Mollo AR, Corrigan OI. An investigation of the mechanism of release of the amphoteric drug amoxycillin from poly( D , L -lactide-co-glycolide) matrices. Pharm Dev Technol. 2002;7:333–43.
  } else if(which.model == "Hopfenberg"){
    
    diss <- 100*(1-(1-alpha*t*60)^beta)
    # y = 100*(1-(1-alpha*t)^beta)
    # Enscore DJ, Hopfenberg HB, Stannett VT. Effect of particle size on the mechanism controlling n-hexane sorption in glassy polystyrene microspheres. Polymer. 1977;18:793–800.
  }

  closeAllConnections()
  return(diss)
}


### Weibull model with fixed parameters of alpha and beta

weibull.fit <- function(t, alpha.calc,beta.calc){
  
  alpha <- alpha.calc
  beta <- beta.calc
  
  diss <- 1-exp(-alpha*t^(beta))
  
  return(diss)
}


### RMSE function # Function that returns Root Mean Squared Error
RMSE <- function(predicted, observed){
    
    error <- observed - predicted
    res <- sqrt(mean(error^2))
    
    return(res)
}


### Taken after package sparsediscrim, to calculate covariance matrix for multi-sample pooled covariance

cov_pool <- function (x, y){
    x <- as.matrix(x)
    y <- as.factor(y)
    n <- length(y)
    scatter_matrices <- tapply(seq_len(n), y, function(i) {
        (length(i) - 1) * cov(as.matrix(x[i, ]))
    })
    as.matrix(Reduce("+", scatter_matrices)/(n-nlevels(y))) ### MODIFIED PUT (n-nlevels(y)) to reflect the Sathe et al.
}
