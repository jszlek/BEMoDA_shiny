# BEMoDA v1.0 - BiowaivEr aid for Model Dependent-Independent Approach script for in-vitro dissolution profile comparison
# 
# Model Independent Approach script for in-vitro dissolution profile comparison as proposed by Tsong et al. in 1996
# (Tsong Y, Hammerstrom T, Sathe P, Shah VP. (1996) Statistical Assessment of Mean Differences between Two Dissolution Data Sets, Drug Info. J. 30:1105-1112).
# 
# The developed script was based on the code published at forum.bebac.at by Shuanghe (post: SAS/R code of MSD method for comparing dis�so�lution)
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

BEMoDA_MANOVA <- function(filename_ref, filename_test, manova_test, manova_intercept, manova_tol){
  
  # filename of dissolution of test product
  filename_test <- filename_test
  
  # filename of dissolution of reference product
  filename_ref <- filename_ref
  
  
  # read dissolution data
  mt <- read.csv(file=filename_test,header=TRUE,row.names=1,sep="\t")
  mr <- read.csv(file=filename_ref,header=TRUE,row.names=1,sep="\t")
  
  my.data <- bind_rows(mr,mt,.id="Formulation")
  
  dep_part<- paste("~",colnames(my.data)[1],sep="")
  ind_part<- paste(colnames(my.data)[-1],collapse=",")
  dt_formula<- as.formula(paste("cbind(",ind_part,")", dep_part, sep=""))
  
  
  my.model <- manova(dt_formula, data = my.data)
  my.summary1 <- summary(my.model, test=manova_test, intercept= manova_intercept, tol = manova_tol)
  my.summary2 <- summary.aov(my.model)
  
  
  return(list(manova.summary=my.summary1, anova.summary=my.summary2))
}
