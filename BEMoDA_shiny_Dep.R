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
# 

BEMoDA_shiny_Dep <- function(input_ref, input_test, input_std, optim.model.params, dissol.model, optim.ellipse.params){
#######################################################
### OPTIONS FOR TEST AND REFERENCE PRODUCTS        ####
#######################################################

# dissolution data for test & reference, format in columns t1 t2 t3 t4 ...
# rows contains ref_1 ref_2 ref_3 etc.
# TAB-delimited file with column names and rownames
  
print(input_ref)
print(input_test)
print(input_std)
print(optim.model.params)
print(dissol.model)
print(optim.ellipse.params)


mr <- input_ref
mt <- input_test
std.list <- input_std
std.no <- length(std.list)

plot.std <- list()
plot.ref <- list()
plot.test <- list()
rect.sr <- list()
res.table.std <- list()

###################################################################
###             OPTIONS FOR MODEL PARAMS OPTIMIZATION          ####
###################################################################

# Optimization method for mech. model fitting, same as for method parameter in optimx()

maxit <- optim.model.params$maxit

###################################################################
### OPTIONS FOR nloptr TO SEARCH FOR CR AND SR ELLIPSES        ####
###################################################################

params.no <- 2 # how many parameters in the equation
starting.params <- rnorm(params.no)/10 # generate starting params for optim functions
lower.boundary <- rep(-100*max(abs(starting.params)),times=params.no) # for alpha and beta
upper.boundary <- rep(100*max(abs(starting.params)),times=params.no) # for alpha and beta

sr.level <- optim.ellipse.params$sr_level
cr.level <- optim.ellipse.params$cr_level
ellipse.sr.npts <- optim.ellipse.params$sr_pts
ellipse.cr.npts <- optim.ellipse.params$cr_pts
maxit.ellipse <- optim.ellipse.params$iter
toler.ellipse <- optim.ellipse.params$toler


##############################
##### CHECK DATA FORMAT ######
#############################

nt <- nrow(mt)
nr <- nrow(mr)
ns <- nrow(std.list[[1]])

if(nt != nr){
    stop("Number of observations (dissolution profiles) in test and reference differ! Please check the data!")
}


p <- ncol(mt)

if(p != ncol(mr)){
    stop("Number of parameters (time points) in test and reference differ! Please check the data!")
}



# Prepare data frame for writing parameters, RMSE and R2
par.ref.df <- matrix(ncol = 5, nrow = nr)
par.test.df <- matrix(ncol = 5, nrow = nt)


for(i in 1:length(std.list)){

par.std.df.name <- paste("par.std.df",i,sep="")

assign(par.std.df.name, setNames(data.frame(matrix(ncol=3, nrow = ns)), c("Formulation","param_A","param_B")))

}


tmp.x <- c("Formulation","param_A","param_B","RMSE","R2")
colnames(par.ref.df) <- tmp.x
colnames(par.test.df) <- tmp.x


# Prepare matrix to save RMSE
diss.rmse.ref <- matrix(ncol = 2, nrow = nr)
diss.rmse.test <- matrix(ncol = 2, nrow = nt)

tmp.x <- c("Formulation","RMSE")
colnames(diss.rmse.ref) <- tmp.x
colnames(diss.rmse.test) <- tmp.x




##########################################################
###         MODEL FITTING FOR REF AND TEST            ####
##########################################################

for(i in 1:nt){

# Prepare new data frame
# Prepare data points - acquire them from the table header

dat.r <- extract.timepoints(mr)
dat.t <- extract.timepoints(mt)

dat.r <-rbind(t=dat.r,X=mr[i,])
dat.t <-rbind(t=dat.t,X=mt[i,])



# Transform time into hours
dat.r[1,] <- dat.r[1,]/60
dat.r[2,] <- dat.r[2,]

dat.t[1,] <- dat.t[1,]/60
dat.t[2,] <- dat.t[2,]


# Transpose data frame
dat.r <- as.data.frame(t(dat.r))
dat.t <- as.data.frame(t(dat.t))


# optimx optim method
if(optim.model.params$method == "optimx"){

optim.method <- optim.model.params$optimx_method
maxit <- optim.model.params$maxit

optimx.r <- optimx(starting.params, model.2params, t=dat.r[,1], X=dat.r[,2], which.model=dissol.model, method=optim.method, control=list(trace=FALSE,maxit=maxit))
optimx.t <- optimx(starting.params, model.2params, t=dat.t[,1], X=dat.t[,2], which.model=dissol.model, method=optim.method, control=list(trace=FALSE,maxit=maxit))

# Get the alpha and beta parameters from the model
r.alpha.calc <- optimx.r$p1
r.beta.calc <- optimx.r$p2

t.alpha.calc <- optimx.t$p1
t.beta.calc <- optimx.t$p2

}


# nloptr optim method

if(optim.model.params$method == "nloptr"){

maxit <- optim.model.params$maxit
optim_rel_tol <- optim.model.params$toler


nloptr.r <- nloptr(x0=starting.params, eval_f=model.2params, lb=lower.boundary, ub=upper.boundary,
                            t=dat.r[,1],X=dat.r[,2], which.model=dissol.model,
                            opts=list(algorithm="NLOPT_GN_CRS2_LM",xtol_rel=optim_rel_tol,maxeval=maxit,print_level=0,local_opts=list(algorithm="NLOPT_LD_MMA",xtol_rel=optim_rel_tol))
                            )

nloptr.t <- nloptr(x0=starting.params, eval_f=model.2params, lb=lower.boundary, ub=upper.boundary,
                            t=dat.t[,1],X=dat.t[,2], which.model=dissol.model,
                            opts=list(algorithm="NLOPT_GN_CRS2_LM",xtol_rel=optim_rel_tol,maxeval=maxit,print_level=0,local_opts=list(algorithm="NLOPT_LD_MMA",xtol_rel=optim_rel_tol))
                            )
                            
                            
  r.alpha.calc <- nloptr.r$solution[1]
  r.beta.calc <- nloptr.r$solution[2]
  
  t.alpha.calc <- nloptr.t$solution[1]
  t.beta.calc <- nloptr.t$solution[2]
                            
}

# GenSA optim method

if(optim.model.params$method == "genSA"){

  max_iter_gensa <- optim.model.params$maxit
  
  gensa.r <- GenSA(par=starting.params, fn=model.2params,
                            t=dat.r[,1],X=dat.r[,2], which.model=dissol.model, lower=lower.boundary,
                            upper=upper.boundary,
                            control=list(smooth=FALSE,maxit=max_iter_gensa,verbose=TRUE,nb.stop.improvement=max_iter_gensa)
                            )

  
  gensa.t <- GenSA(par=starting.params, fn=model.2params,
                            t=dat.t[,1],X=dat.t[,2], which.model=dissol.model, lower=lower.boundary,
                            upper=upper.boundary,
                            control=list(smooth=FALSE,maxit=max_iter_gensa,verbose=TRUE,nb.stop.improvement=max_iter_gensa)
                            )

  r.alpha.calc <- gensa.r$par[1] 
  r.beta.calc <- gensa.r$par[2]
  
  t.alpha.calc <- gensa.t$par[1] 
  t.beta.calc <- gensa.t$par[2]
  
}


# nls optim method
if(optim.model.params$method == "nls"){

 
  # Predict nonlinear least square model for different dissolution mech. models
  
  if(dissol.model == "Weibull"){
    
    fit.s <- nls(X ~ 100*(1-exp(-alpha*t^(beta))),start=list(alpha=1,beta=1),data=as.data.frame(dat.s),control = list(maxiter = maxit,warnOnly=T),trace=TRUE)
    
  } else if (dissol.model == "Korsmeyer-Peppas"){
    
    fit.s <- nls(X ~ 100*(alpha*(t^(beta))),start=list(alpha=1,beta=1),data=as.data.frame(dat.s),control = list(maxiter = maxit,warnOnly=T),trace=TRUE)
    
  } else if(dissol.model == "Peppas-Sahlin"){
    
    fit.s <- nls(X ~ (alpha*(t*60)^(0.5)+beta*(t*60)),start=list(alpha=1,beta=1),data=as.data.frame(dat.s),control = list(maxiter = maxit,warnOnly=T),trace=TRUE)
    
  } else if (dissol.model =="Quadratic"){
    
    fit.s <- nls(X ~ 100*(alpha*(t*60)^2+beta*(t*60)),start=list(alpha=1,beta=1),data=as.data.frame(dat.s),control = list(maxiter = maxit,warnOnly=T),trace=TRUE)
    
  } else if(dissol.model == "Logistic"){
    
    fit.s <- nls(X ~ 100*((exp(alpha+beta*log10(t*60)))/(1+exp(alpha+beta*log10(t*60)))),start=list(alpha=1,beta=1),data=as.data.frame(dat.s),control = list(maxiter = maxit,warnOnly=T),trace=TRUE)
    
  } else if(dissol.model == "Gompertz"){
    
    fit.s <- nls(X ~ 100*exp(-alpha*exp(-beta*log10(t*60))),start=list(alpha=1,beta=1),data=as.data.frame(dat.s),control = list(maxiter = maxit,warnOnly=T),trace=TRUE)
    
    # } else if (which.model=="Probit"){ 
    #   fi <- pnorm(time,mean=mean(time),sd=sd(time)) # should be normal distribution function of time, but it seems that smth is wrong here =========> switch-off
    #   pred <- 100*fi*(alpha+beta*log(t))
    #   # y = 100*fi*(alpha+beta*log(t))
    #   # Tsong Y, Hammerstrom T, Chen JJ. Multipoint dissolution specification and acceptance sampling rule based on profile modeling and principal component analysis. J Biopharm Stat. 1997;7:423–39.
    
  } else if(dissol.model == "Zero-order with Tlag"){
    
    fit.s <- nls(X ~ alpha*(t*60 - beta),start=list(alpha=1,beta=1),data=as.data.frame(dat.s),control = list(maxiter = maxit,warnOnly=T),trace=TRUE)
    
  } else if(dissol.model == "Zero-order with F0"){
    
    fit.s <- nls(X ~ alpha + beta*t*60,start=list(alpha=1,beta=1),data=as.data.frame(dat.s),control = list(maxiter = maxit,warnOnly=T),trace=TRUE)
    
  } else if(dissol.model == "First-order with Tlag"){
    
    fit.s <- nls(X ~ 100*(1-exp(-alpha*(t*60-beta))),start=list(alpha=1,beta=1),data=as.data.frame(dat.s),control = list(maxiter = maxit,warnOnly=T),trace=TRUE)
    
  } else if(dissol.model == "First-order with Fmax"){
    
    fit.s <- nls(X ~ alpha*(1-exp(-beta*t*60)),start=list(alpha=1,beta=1),data=as.data.frame(dat.s),control = list(maxiter = maxit,warnOnly=T),trace=TRUE)
    
  } else if(dissol.model == "Higuchi with Tlag"){
    
    fit.s <- nls(X ~ alpha*(t*60-beta)^0.5,start=list(alpha=1,beta=1),data=as.data.frame(dat.s),control = list(maxiter = maxit,warnOnly=T),trace=TRUE)
    
  } else if(dissol.model == "Higuchi with F0"){
    
    fit.s <- nls(X ~ alpha + beta*(t*60)^0.5,start=list(alpha=1,beta=1),data=as.data.frame(dat.s),control = list(maxiter = maxit,warnOnly=T),trace=TRUE)
    
  } else if(dissol.model == "Hixon-Crowell with Tlag"){
    
    fit.s <- nls(X ~ 100*(1-(1-alpha*(t*60-beta))^3),start=list(alpha=1,beta=1),data=as.data.frame(dat.s),control = list(maxiter = maxit,warnOnly=T),trace=TRUE)
    
  } else if(dissol.model == "Hopfenberg"){
    
    fit.s <- nls(X ~ 100*(1-(1-alpha*t*60)^beta),start=list(alpha=1,beta=1),data=as.data.frame(dat.s),control = list(maxiter = maxit,warnOnly=T),trace=TRUE)
    
  }
  


# Get the alpha and beta parameters from the model
r.alpha.calc <- coef(fit.r)[1]
r.beta.calc <- coef(fit.r)[2]

t.alpha.calc <- coef(fit.t)[1]
t.beta.calc <- coef(fit.t)[2]

}



# Prepare new data for time values
new.t.ref <- seq(0,max(dat.r[,1]),0.05)
new.t.test <- seq(0,max(dat.t[,1]),0.05)


# Calculate new Q values with new mech. models
new.X.ref <- model.2params.fit(new.t.ref, r.alpha.calc, r.beta.calc, dissol.model)
new.X.test <- model.2params.fit(new.t.test, t.alpha.calc, t.beta.calc, dissol.model)


# Calculate RMSE for optimized mech. models

# First calculate predicted dissolution based on mech. models
diss.ref <- model.2params.fit(dat.r[,1], r.alpha.calc, r.beta.calc, dissol.model)
diss.test <- model.2params.fit(dat.t[,1], t.alpha.calc, t.beta.calc, dissol.model)

# RMSE
rmse.ref <- RMSE(diss.ref, dat.r[,2])
rmse.test <- RMSE(diss.test, dat.t[,2])

# R squared
r2.ref <- summary(lm(diss.ref ~ dat.r[,2]))$r.squared
r2.test <- summary(lm(diss.test ~ dat.t[,2]))$r.squared

# Save parameters, RMSE and R2 in data frame par.ref.df and par.test.df
par.ref.df[i,1] <- i
par.ref.df[i,2] <- r.alpha.calc
par.ref.df[i,3] <- r.beta.calc
par.ref.df[i,4] <- rmse.ref
par.ref.df[i,5] <- r2.ref

par.test.df[i,1] <- i
par.test.df[i,2] <- t.alpha.calc
par.test.df[i,3] <- t.beta.calc
par.test.df[i,4] <- rmse.test
par.test.df[i,5] <- r2.test


# Save RMSE in matrix
diss.rmse.ref[i,1] <- i
diss.rmse.ref[i,2] <- rmse.ref

diss.rmse.test[i,1] <- i
diss.rmse.test[i,2] <- rmse.test
# 
# cat("\n")
# cat("=========================","\n")
# cat("ref ",i,"\n")
# cat("alpha = ", r.alpha.calc,"\n")
# cat("beta = ", r.beta.calc,"\n")
# cat("=========================","\n")
# 
# cat("\n")
# cat("=========================","\n")
# cat("test ",i,"\n")
# cat("alpha = ", t.alpha.calc,"\n")
# cat("beta = ", t.beta.calc,"\n")
# cat("=========================","\n")


# Plot the points and fitted model curve

plot.ref.tmp <- ggplot(data=dat.r, aes(x=t,y=X),color='blue') + geom_point() +
geom_line(color='red',data=data.frame(t=new.t.ref,X=new.X.ref),aes(x=t,y=X)) + 
  geom_text(hjust = 0, x = max(dat.r$t)*0.05, y = max(dat.r$X)*0.95, label=paste("RMSE = ", format(rmse.ref, digits = 3))) +
  geom_text(hjust = 0, x = max(dat.r$t)*0.05, y = max(dat.r$X)*0.90, label=paste("R^2 = ", format(r2.ref, digits = 3))) +
labs(x="Time [hours]",y="Q [%]") + ggtitle(paste("REFERENCE ", i , " " ,dissol.model," model")) + 
theme(plot.title = element_text(lineheight=.8, face="bold"))
# ggsave(paste("REF_",i,"_" ,dissol.model,"_plot.pdf",sep=""))

plot.test.tmp <- ggplot(data=dat.t, aes(x=t,y=X),color='blue') + geom_point() +
geom_line(color='red',data=data.frame(t=new.t.test,X=new.X.test),aes(x=t,y=X)) +
  geom_text(hjust = 0, x = max(dat.t$t)*0.05, y = max(dat.t$X)*0.95, label=paste("RMSE = ", format(rmse.test, digits = 3))) +
  geom_text(hjust = 0, x = max(dat.t$t)*0.05, y = max(dat.t$X)*0.90, label=paste("R^2 = ", format(r2.test, digits = 3))) +
labs(x="Time [hours]",y="Q [%]") + ggtitle(paste("TEST ", i , " " ,dissol.model," model")) + 
theme(plot.title = element_text(lineheight=.8, face="bold"))
# ggsave(paste("TEST_",i,"_",dissol.model,"_plot.pdf",sep=""))


plot.ref[[i]] <- plot.ref.tmp 
plot.test[[i]] <- plot.test.tmp 
    
}

res.table.test <- par.test.df
res.table.ref <- par.ref.df

par.test.df <- as.data.frame(par.test.df[,2:5])
par.ref.df <- as.data.frame(par.ref.df[,2:5])

cat("","\n")
cat("Model parameters for REFERENCE")
print(res.table.ref)
cat("\n")
cat("\n")
cat("Model parameters for TEST")
print(res.table.test)
cat("\n")
cat("\n")


################################################
###        MODEL FITTING FOR STD            ####
################################################

for(k in 1:length(std.list)){

plot.std.no <- list()
  
tmp.std <- std.list[[k]]   

tmp.ns <- nrow(tmp.std)

par.std.df.name <- paste("par.std.df",k,sep="")
assign(paste("tmp.par"),setNames(data.frame(matrix(nrow=tmp.ns,ncol=3)),c("Formulation","param_A","param_B")))

assign(paste("tmp.rmse.r2.table.std"),setNames(data.frame(matrix(nrow=tmp.ns,ncol=5)),c("Formulation","param_A","param_B","RMSE","R2")))


for(i in 1:ns){

plot.std.tmp <- list()  
  
# Prepare new data frame
# Prepare data points - acquire them from the table header

dat.s <- extract.timepoints(tmp.std)

dat.s <-rbind(t=dat.s,X=tmp.std[i,])

# Make fractions of Q and transform time into hours
dat.s[1,] <- dat.s[1,]/60
dat.s[2,] <- dat.s[2,]

# Transpose data frame
dat.s <- as.data.frame(t(dat.s))

# optimx optim method
if(optim.model.params$method == "optimx"){

optimx.s <- optimx(starting.params, model.2params, t=dat.s[,1],X=dat.s[,2], which.model=dissol.model, method=optim.method, control=list(trace=FALSE,maxit=maxit))

# Get the alpha and beta parameters from the model
s.alpha.calc <- optimx.s$p1
s.beta.calc <- optimx.s$p2

}


# nloptr optim method

if(optim.model.params$method == "nloptr"){

nloptr.s <- nloptr(x0=starting.params, eval_f=model.2params, lb=lower.boundary, ub=upper.boundary,
                            t=dat.s[,1],X=dat.s[,2], which.model=dissol.model,
                            opts=list(algorithm="NLOPT_GN_CRS2_LM",xtol_rel=optim_rel_tol,maxeval=maxit,print_level=0,local_opts=list(algorithm="NLOPT_LD_MMA",xtol_rel=optim_rel_tol))
                            )
            
  s.alpha.calc <- nloptr.s$solution[1]
  s.beta.calc <- nloptr.s$solution[2]
  
}

# GenSA optim method

if(optim.model.params$method == "genSA"){

  gensa.s <- GenSA(par=starting.params, fn=model.2params,
                            t=dat.s[,1], X=dat.s[,2], which.model = dissol.model, lower=lower.boundary,
                            upper=upper.boundary,
                            control=list(smooth=FALSE,maxit=maxit,verbose=TRUE,nb.stop.improvement=max_iter_gensa.weibull)
                            )

  s.alpha.calc <- gensa.s$par[1] 
  s.beta.calc <- gensa.s$par[2]
  
}


# nls optim method
if(optim.model.params$method == "nls"){

  
  # Predict nonlinear least square model for different dissolution mech. models
  
  if(dissol.model == "Weibull"){
    
    fit.s <- nls(X ~ 100*(1-exp(-alpha*t^(beta))),start=list(alpha=1,beta=1),data=as.data.frame(dat.s),control = list(maxiter = maxit,warnOnly=T),trace=TRUE)
    
  } else if (dissol.model == "Korsmeyer-Peppas"){
    
    fit.s <- nls(X ~ 100*(alpha*(t^(beta))),start=list(alpha=1,beta=1),data=as.data.frame(dat.s),control = list(maxiter = maxit,warnOnly=T),trace=TRUE)
    
  } else if(dissol.model == "Peppas-Sahlin"){
    
    fit.s <- nls(X ~ (alpha*(t*60)^(0.5)+beta*(t*60)),start=list(alpha=1,beta=1),data=as.data.frame(dat.s),control = list(maxiter = maxit,warnOnly=T),trace=TRUE)
    
  } else if (dissol.model =="Quadratic"){
    
    fit.s <- nls(X ~ 100*(alpha*(t*60)^2+beta*(t*60)),start=list(alpha=1,beta=1),data=as.data.frame(dat.s),control = list(maxiter = maxit,warnOnly=T),trace=TRUE)
    
  } else if(dissol.model == "Logistic"){
    
    fit.s <- nls(X ~ 100*((exp(alpha+beta*log10(t*60)))/(1+exp(alpha+beta*log10(t*60)))),start=list(alpha=1,beta=1),data=as.data.frame(dat.s),control = list(maxiter = maxit,warnOnly=T),trace=TRUE)
    
  } else if(dissol.model == "Gompertz"){
    
    fit.s <- nls(X ~ 100*exp(-alpha*exp(-beta*log10(t*60))),start=list(alpha=1,beta=1),data=as.data.frame(dat.s),control = list(maxiter = maxit,warnOnly=T),trace=TRUE)
    
    # } else if (which.model=="Probit"){ 
    #   fi <- pnorm(time,mean=mean(time),sd=sd(time)) # should be normal distribution function of time, but it seems that smth is wrong here =========> switch-off
    #   pred <- 100*fi*(alpha+beta*log(t))
    #   # y = 100*fi*(alpha+beta*log(t))
    #   # Tsong Y, Hammerstrom T, Chen JJ. Multipoint dissolution specification and acceptance sampling rule based on profile modeling and principal component analysis. J Biopharm Stat. 1997;7:423–39.
    
  } else if(dissol.model == "Zero-order with Tlag"){
    
    fit.s <- nls(X ~ alpha*(t*60 - beta),start=list(alpha=1,beta=1),data=as.data.frame(dat.s),control = list(maxiter = maxit,warnOnly=T),trace=TRUE)
    
  } else if(dissol.model == "Zero-order with F0"){
    
    fit.s <- nls(X ~ alpha + beta*t*60,start=list(alpha=1,beta=1),data=as.data.frame(dat.s),control = list(maxiter = maxit,warnOnly=T),trace=TRUE)
    
  } else if(dissol.model == "First-order with Tlag"){
    
    fit.s <- nls(X ~ 100*(1-exp(-alpha*(t*60-beta))),start=list(alpha=1,beta=1),data=as.data.frame(dat.s),control = list(maxiter = maxit,warnOnly=T),trace=TRUE)
    
  } else if(dissol.model == "First-order with Fmax"){
    
    fit.s <- nls(X ~ alpha*(1-exp(-beta*t*60)),start=list(alpha=1,beta=1),data=as.data.frame(dat.s),control = list(maxiter = maxit,warnOnly=T),trace=TRUE)
    
  } else if(dissol.model == "Higuchi with Tlag"){
    
    fit.s <- nls(X ~ alpha*(t*60-beta)^0.5,start=list(alpha=1,beta=1),data=as.data.frame(dat.s),control = list(maxiter = maxit,warnOnly=T),trace=TRUE)
    
  } else if(dissol.model == "Higuchi with F0"){
    
    fit.s <- nls(X ~ alpha + beta*(t*60)^0.5,start=list(alpha=1,beta=1),data=as.data.frame(dat.s),control = list(maxiter = maxit,warnOnly=T),trace=TRUE)
    
  } else if(dissol.model == "Hixon-Crowell with Tlag"){
    
    fit.s <- nls(X ~ 100*(1-(1-alpha*(t*60-beta))^3),start=list(alpha=1,beta=1),data=as.data.frame(dat.s),control = list(maxiter = maxit,warnOnly=T),trace=TRUE)
    
  } else if(dissol.model == "Hopfenberg"){
    
    fit.s <- nls(X ~ 100*(1-(1-alpha*t*60)^beta),start=list(alpha=1,beta=1),data=as.data.frame(dat.s),control = list(maxiter = maxit,warnOnly=T),trace=TRUE)
    
  }
  

# Get the alpha and beta parameters from the model
s.alpha.calc <- coef(fit.s)[1]
s.beta.calc <- coef(fit.s)[2]

}



# Prepare new data for time values
new.t.s <- seq(0,max(dat.s[,1]),0.05)

# Calculate new Q values with new model
new.X.s <- model.2params.fit(new.t.s, s.alpha.calc, s.beta.calc, dissol.model)

# Calculate RMSE for optimized equation

# First calculate predicted dissolution based on model
diss.std <- model.2params.fit(dat.s[,1], s.alpha.calc, s.beta.calc, dissol.model)

# Calculate % of dissolution
# diss.std <- 100*diss.std
# dat.s.diss <- 100*dat.s[,2]

# RMSE
rmse.std <- RMSE(diss.std, dat.s[,2])

#R squared
r2.std <- summary(lm(diss.std ~ dat.s[,2]))$r.squared


# Save prameters, RMSE and R2
tmp.par[i,1] <- i
tmp.par[i,2] <- s.alpha.calc
tmp.par[i,3] <- s.beta.calc
# tmp.par[i,4] <- rmse.std
# tmp.par[i,5] <- r2.std

tmp.rmse.r2.table.std[i,1] <- i
tmp.rmse.r2.table.std[i,2] <- s.alpha.calc
tmp.rmse.r2.table.std[i,3] <- s.beta.calc
tmp.rmse.r2.table.std[i,4] <- rmse.std
tmp.rmse.r2.table.std[i,5] <- r2.std



print("RMSE")
print(rmse.std)
# 
# # Save RMSE in matrix
# diss.rmse.std[i,1] <- i
# diss.rmse.std[i,2] <- rmse.std

cat("\n")
cat("=========================","\n")
cat("STD ",i,"\n")
cat("alpha = ", s.alpha.calc,"\n")
cat("beta = ", s.beta.calc,"\n")
cat("=========================","\n")



# Plot the points and fitted curve

# pdf(paste(i,"_mt_plot.pdf"))

plot.std.tmp <- ggplot(data=dat.s, aes(x=t,y=X),color='blue') + geom_point() +
geom_line(color='red',data=data.frame(t=new.t.s,X=new.X.s),aes(x=t,y=X)) +
  geom_text(hjust = 0, x = max(dat.s$t)*0.05, y = max(dat.s$X)*0.95, label=paste("RMSE = ", format(rmse.std, digits = 3))) +
  geom_text(hjust = 0, x = max(dat.s$t)*0.05, y = max(dat.s$X)*0.90, label=paste("R^2 = ", format(r2.std, digits = 3))) +
labs(x="Time [hours]",y="Q [%]") + ggtitle(paste("STANDARD ",i," batch no ", k, dissol.model," model")) + 
theme(plot.title = element_text(lineheight=.8, face="bold"))
# ggsave(paste("STD_",k,"_",i,"_",dissol.model,"_plot.pdf",sep=""))


plot.std.no[[i]] <- plot.std.tmp

}

plot.std[[k]] <- plot.std.no

res.table.std[[k]] <- tmp.rmse.r2.table.std[,2:5]


assign(par.std.df.name, tmp.par[,2:3])

cat("","\n")
cat(dissol.model," model parameters for STD_",k)
print(tmp.par)
cat("\n")

}


cat("","\n")
cat("res.table.std list")
print(res.table.std)
cat("\n")


# Calculate log-normal matrix for std batches
# OVERWRITTING STD.LIST !!!
if (optim.ellipse.params$log_trans == TRUE){
      
    for(i in 1:length(std.list)){
      
      par.std.df.name <- paste("par.std.df",i,sep="")
      
      assign(par.std.df.name, log(get(par.std.df.name)))
      
      std.list[[i]] <- get(paste("par.std.df",i,sep=""))
      
      }
} else if (optim.ellipse.params$log_trans == FALSE){
    
    for(i in 1:length(std.list)){
      
      par.std.df.name <- paste("par.std.df",i,sep="")
      
      assign(par.std.df.name, get(par.std.df.name))
      
      std.list[[i]] <- get(paste("par.std.df",i,sep=""))
    
  }
}

all.std.par <- bind_rows(std.list,.id="id")



cat("ALL.STD.LIST","\n")
print(all.std.par)
cat("\n")

### S.pooled - variance/covariance matrix of standard batches (std) 

S.sr <- cov_pool(all.std.par[,2:3],all.std.par[,1])

# Univariate SD.pooled
# In In Vitro-In Vivo Correlations, Eds.: David B. Young,John G. Devane,Jackie Butler
# In Vitro Dissolution Profile Comparison And IVIVR
# 
# # matrix of mean variances, intra-lot variances

sd.intra.lot.alpha <- var(tapply(all.std.par[,2], all.std.par[,1], mean))
sd.intra.lot.beta <- var(tapply(all.std.par[,3], all.std.par[,1], mean))

sd.log.std.alpha <- sapply(std.list, function(x) var(x[,1]))
sd.log.std.beta <- sapply(std.list, function(x) var(x[,2]))

sqr.sd.log.std.alpha <- sd.log.std.alpha^2
sqr.sd.log.std.beta <- sd.log.std.beta^2

SD.pooled.alpha <- sqrt((sum(sqr.sd.log.std.alpha)/length(std.list))+sd.intra.lot.alpha)
SD.pooled.beta <- sqrt((sum(sqr.sd.log.std.beta)/length(std.list))+sd.intra.lot.beta)


cat("SD.pooled.alpha","\n")
print(SD.pooled.alpha)
cat("SD.pooled.beta","\n")
print(SD.pooled.beta)

rect.sr1 <- data.frame(x1=SD.pooled.alpha, x2=-SD.pooled.alpha, y1=SD.pooled.beta, y2=-SD.pooled.beta)
rect.sr2 <- data.frame(x1=2*SD.pooled.alpha, x2=-2*SD.pooled.alpha, y1=2*SD.pooled.beta, y2=-2*SD.pooled.beta)
rect.sr3 <- data.frame(x1=3*SD.pooled.alpha, x2=-3*SD.pooled.alpha, y1=3*SD.pooled.beta, y2=-3*SD.pooled.beta)

rect.sr <- list(rect.sr1, rect.sr2, rect.sr3)

# About log-transformation
# After Zhang et al. An introduction to the approaches used by DDSolver
# (4) Because the MSD is calculated under the assumption that the model parameters are multivariate normally distributed, it should be determined whether the model parameters need to be log-transformed before the comparison is performed,
# which depends mainly on the statistical distribution properties of the compared parameters. Take the widely used Weibull model, for example: some studies have concluded similarity between two sets of parameters based on 
# natural-logarithm-transformed data (31), while some other studies have used untransformed data (32). To ensure rational application of model-dependent approaches, this issue obviously needs to be explicitly clarified.

par.test.df <- as.data.frame(par.test.df[,1:2])
par.ref.df <- as.data.frame(par.ref.df[,1:2])

if (optim.ellipse.params$log_trans == TRUE){
    
    log.mt <- log(par.test.df)
    log.mr <- log(par.ref.df)

  } else if (optim.ellipse.params$log_trans == FALSE ){
    
    log.mt <- par.test.df
    log.mr <- par.ref.df
  
}

# Check once again number of parameters and tablets
nt <- nrow(log.mt)
nr <- nrow(log.mr)
p <- ncol(log.mt)

if (optim.ellipse.params$log_trans == TRUE){
  
cat("log.mt: ","\n")
print(log.mt)
cat("\n")

cat("log.mr: ","\n")
print(log.mr)
cat("\n")

} else if (optim.ellipse.params$log_trans == FALSE){
  
  cat("mt: ","\n")
  print(log.mt)
  cat("\n")
  
  cat("mr: ","\n")
  print(log.mr)
  cat("\n")
  
} 

# Scaling factor
k.val <- ((nt+nr-p-1)/((nt+nr-2)*p))*((nt*nr)/(nt+nr)) 


# F-distribution value of critical region
F.cr <- qf(cr.level,p,(2*nt-p-1))

# F-distribution value of similarity region 
F.sr <- qf(sr.level,p,(2*nt-p-1))

# Covariance matrix for S_pooled 
S.cr <- (cov(log.mt)+cov(log.mr))/2

# Means of paramters for test and reference dissolution profiles
mean.log.mt <- colMeans(log.mt)
mean.log.mr <- colMeans(log.mr)
mean.log.diff <- mean.log.mt - mean.log.mr

# Difference of parameters
# log.diff <- log.mt - log.mr
log.diff <- log.mt - log.mr

# Mean of parameters difference
# mean.log.sr <- colMeans(rbind(colMeans(log.std1), colMeans(log.std2),colMeans(log.std3)))

###############################################################################################                            
#####     CALCULATE MAHALANOBIS DISTANCE AND HOTTELING T^2 FOR REFERENCE AND TEST          ####
###############################################################################################


# Mahalanobis distance, 'M' Distance
M.dist.cr <- sqrt(t(mean.log.diff) %*% solve(S.cr) %*% mean.log.diff)

# Hotelling's T^2, as Scaled 'M' Distance as in Sathe et. al - Hotelling's two-sample t-squared statistic
# https://en.wikipedia.org/wiki/Hotelling%27s_T-squared_distribution
# http://sites.stat.psu.edu/~ajw13/stat505/fa06/11_2sampHotel/01_2sampHotel.html
# 
M.dist.cr.scaled <- (nt*nr)/(nt+nr)*M.dist.cr^2



######################################################################################################################                            
#####     SEARCH FOR POINTS ON THE BORDER OF CRITICAL REGION (OPTIMIZE FUNCTION BY NLOPTR) AND SIMILARITY REGION  ####
######################################################################################################################


cat("S.sr: ","\n")
print(S.sr)
cat("\n")

cat("S.cr: ","\n")
print(S.cr)
cat("\n")

cat("M.dist.cr: ","\n")
print(M.dist.cr)
cat("\n")


cat("M.dist.cr.scaled: ","\n")
print(M.dist.cr.scaled)
cat("\n")



ellipse.sr <- matrix(nrow=ellipse.sr.npts,ncol=2,data=NA)
ellipse.cr <- matrix(nrow=ellipse.cr.npts,ncol=2,data=NA)


for (i in 1:ellipse.sr.npts){
fit <- nloptr(x0=starting.params, eval_f=sr.ellipse, lb=lower.boundary, ub=upper.boundary,
                            S=S.sr, k.val=k.val, p=p, nr=nr, level=sr.level, 
                            opts=list(algorithm="NLOPT_GN_CRS2_LM",xtol_rel=toler.ellipse, maxeval=maxit.ellipse,print_level=0,local_opts=list(algorithm="NLOPT_LD_MMA",xtol_rel=toler.ellipse))
                            )
ellipse.sr[i,] <- fit$solution
                            
                            
}
                            
for (i in 1:ellipse.cr.npts){
fit <- nloptr(x0=starting.params, eval_f=cr.ellipse, lb=lower.boundary, ub=upper.boundary,
                            S.cr=S.cr, mean.diff=mean.log.diff, k.val=k.val, p=p, nr=nr, level=cr.level, 
                            opts=list(algorithm="NLOPT_GN_CRS2_LM",xtol_rel=toler.ellipse,maxeval=maxit.ellipse,print_level=0,local_opts=list(algorithm="NLOPT_LD_MMA",xtol_rel=toler.ellipse))
                            )
ellipse.cr[i,] <- fit$solution
                            
                            
}


return.obj <- list(plot.ref, plot.test, plot.std, std.list, par.ref.df, par.test.df, ellipse.cr,   # 1 - 7
                   ellipse.sr, rect.sr, res.table.ref, res.table.test, res.table.std, S.cr, S.sr,  # 8 - 14
                   M.dist.cr, M.dist.cr.scaled, log.mt, log.mr, all.std.par, log.diff)             # 15 - 20

return(return.obj)

}
