# BEMoDA shiny v1.1.2 - BiowaivEr aid for Model Dependent-Independent Approach script for in-vitro dissolution profile comparison
# 
# Model Independent Approach script for in-vitro dissolution profile comparison as proposed by Tsong et al. in 1996
# (Tsong Y, Hammerstrom T, Sathe P, Shah VP. (1996) Statistical Assessment of Mean Differences between Two Dissolution Data Sets, Drug Info. J. 30:1105-1112).
# 
# The developed script was based on the code published at forum.bebac.at by Shuanghe (post: SAS/R code of MSD method for comparing dis­so­lution)
# 
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
# 
# dissolution data for test & reference, format in columns t1 t2 t3 t4 ...
# rows contains ref_1 ref_2 ref_3 etc.
# TAB-delimited file with column names and rownames

BEMoDA_shiny_InDep <- function(input_ref, input_test, dg, prob){

# dg is mean % of difference in dissolution between each time points. typically 10. 
# dg=15 was used to compare result with that of the article. 

dg <- dg


# prob is the vector of probabilities calculated via qf() function in R which is 
# Density, distribution function, quantile function and random
# generation for the F distribution with ‘df1’ and ‘df2’ degrees of
# freedom (and optional non-centrality parameter ‘ncp’).

prob <- prob


# nt and nr are number of unit, typically 12, Tsong's article only has 6 unit.
# p is number of time points, dg is mean % of difference in dissolution between each time points. typically 10. The assumption is that the SD is normally distributed across the variables. 
# dg=15 was used to compare result with that of the article.
# it's better to define those number firstly to be easier to modify the code

# read dissolution data
mt <- input_test # test
mr <- input_ref  # reference

nt <- nrow(mt)
nr <- nrow(mr)

tmp_p_mt <- ncol(mt)
tmp_p_mr <- ncol(mr)

if(nt != nr){
    stop("Number of observations (dissolution profiles) in test or reference differ! Please check the data!")
}


if(tmp_p_mt !=  tmp_p_mr){
    stop("Number of parameters (time points) in test and reference differ! Please check the data!")
}

p <- tmp_p_mt

# mean dissolution profile and difference
tbar <- colMeans(mt)
rbar <- colMeans(mr)
mdiff <- tbar-rbar

#matrix of covariance of test and reference and pooled covariance matrix together with it's inverse
st <- cov(mt)
sr <- cov(mr)
sp <- ((nt - 1) * st + (nr - 1) * sr) / (nt  + nr -2)
spinv <- solve(sp)

# Mahalanobis distance
dm <- sqrt(t(mdiff) %*% spinv %*% mdiff)

# k and fcrit
df <- nt + nr -p -1
k <- nr*nr/(nt+nr)*(nt+nr-p-1)/((nt+nr-2)*p)
fcrit <- qf(prob,p,df)

# create column vector with each element equal to predefined limit dg
mdg <- c(rep(dg, p))

# Mahalanobis distance for chosen criteria
dm_max <- sqrt(t(mdg) %*% spinv %*% mdg)

#Hotelling's T^2.
H2 <- nt*nr/(nt+nr)*dm^2

# check length
if(length(spinv) == 1){
   spinv <- rep(spinv, length(mdiff))
}

#vector on the boundary of confidence region
bound1 <- mdiff %*% (1+sqrt(fcrit/(k*t(mdiff) %*% spinv %*% mdiff)))
bound2 <- mdiff %*% (1-sqrt(fcrit/(k*t(mdiff) %*% spinv %*% mdiff)))

# 90% CI of Mahalanobis distance
ci1 <- sqrt(t(bound1) %*% spinv %*% bound1)
ci2 <- sqrt(t(bound2) %*% spinv %*% bound2)

# selecting lower and upper Mahalanobis distance range
dm_lower <- min(ci1, ci2)
dm_upper <- max(ci1, ci2)

  # F-stats
  f_stat_obs <- k * dm^2
  cat("Observed F-stat = ", as.numeric(f_stat_obs), "\n\n")
  cat("Observed F-crit = ", as.numeric(fcrit),"\n\n")
    
if(dm_upper < dm_max){
    
    conclusion <- "Similar"
    
        } else {
        
        conclusion <- "Not similar"
    }   

result <- data.frame(cbind(dm, dm_lower, dm_upper, dg, dm_max, H2, fcrit, f_stat_obs, conclusion))
names(result) <- c("Dm - 'M' Distance", "Dm_lower", "Dm_upper", "Dg (%)", "Dm_max", "Hotelling T Square - Scaled 'M' Distance",
                   "F-crit", "F-stat", "Conclusion")

cat("====================================","\n")
cat("Result matrix","\n")
cat("====================================","\n")

print(result)

cat("","\n")
cat("","\n")
cat("Conclusion. Profiles are: ","\n")
print(conclusion)

return(result)

}
