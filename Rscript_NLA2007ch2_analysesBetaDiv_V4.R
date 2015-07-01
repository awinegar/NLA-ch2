##################################################################################################
####NLA 2007 CHAPTER 2 PROJECT: DIATOM BETA-DIVERSIRY ACROSS THE US (historical/surface seds)    #
##################################################################################################

##In this script: Beta-diversity focused analyses and beta-diversity partitioning 
#Beta-diversity focused analyses 
#Previous script: Rscript_NLA2007ch2_analysesBetaDiv_V3.R
#R version: 3.1.2 (Pumpkin Helmet)

##Last update: July 1, 2015 
##Associated workspace: workspace_NLA2007ch2_analysesBetaDiv_V4.RData
##Associated markdown: 
##Associated .txt of R script: Rscript_NLA2007ch2_analysesBetaDiv_V4.txt
##Github: NLA-ch2 repository 

##################################################################################################

##################################################################################################
####CORRECTIONS TO MAKE                                                                          #
##################################################################################################

##As of July 1, 2015- once code/decisions finalized: 
#Remove 10 sites that are suspect in terms of dating and check for influence on results 
#(URT figures etc. can stay- just number adjustments). 

#Final grep cleaning of diatom species --> re-run for alpha/gamma/beta numbers and insert into
#manuscript. 

#Re-run of entire script with updated diatoms --> final check of all results (+replacing objects
#in R script). Use rm() to get rid of objects not needed. 

#Move all dataframes used as data input into script into repository, change file.choose() to 
#setwd() and push to github. 

#Take final code and adapt for rotifers (cleaned through network project) and soft phytoplankton 
#(cleaned through zUSA). 

##################################################################################################

##################################################################################################
####PACKAGES                                                                                     #
##################################################################################################
####Analyses/data manipulation
library(data.table) #For setnames() function
library(plyr) #For summary stats etc. 
library(vegan) #Ordination 
library(lme4) #Mixed-effect modelling 
library(sp) #For meauring distances between coordinates
library(rpart) #Regression trees
library(MASS) #Chi square 
library(reshape) #Manipulating long/wide format data 
library(reshape2) #Manipulating long/wide format data 
library(party) #Regression trees
library(rpart.plot) #Regression trees
library(partykit) #Regression trees
library(tree) #Regression trees
library(permute) #For rarefied richness
library(boot) #For rarefied richness
library(rich) #For rarefied richness
library(Iso) #For rarefied richness

##Visuals 
library(ggplot2) #Plotting 
library(maptools) #Maps
library(maps) #Maps 
library(rgeos) #Maps 
library(RColorBrewer) #Colour-schemes 
library(gridExtra) #Arranging plots 

##################################################################################################


##################################################################################################
####SOURCE FUNCTIONS                                                                             #
##################################################################################################

####BETA.DIV FUNCTION- for computing beta diversity (also computes LCBD and SCBD) (spatial)
##################################################################################################
beta.div <- function(Y, method="hellinger", sqrt.D=FALSE, samp=TRUE, nperm=999, save.D=FALSE, clock=FALSE)
  #
  # Compute estimates of total beta diversity as the total variance in Y, 
  # for 20 dissimilarity coefficients or analysis of raw data (not recommended). 
  # LCBD indices are tested by permutation within columns of Y.
  # This version includes direct calculation of the Jaccard, Sorensen and Ochiai 
  # coefficients for presence-absence data.
  #
  # Arguments --
  # 
  # Y : community composition data matrix.
  # method : name of one of the 20 dissimilarity coefficients, or "none" for
#          direct calculation on Y (also the case with method="euclidean").
# sqrt.D : If sqrt.D=TRUE, the distances in matrix D are square-rooted before 
#          computation of SStotal, BDtotal and LCBD. 
# samp : If samp=TRUE, the abundance-based distances (ab.jaccard, ab.sorensen,
#        ab.ochiai, ab.simpson) are computed for sample data. If samp=FALSE, 
#        they are computed for true population data.
# nperm : Number of permutations for test of LCBD.
# save.D : If save.D=TRUE, the distance matrix will appear in the output list.
# clock : If clock=TRUE, the computation time is printed in the R console.
#
# Reference --
#
# Legendre, P. and M. De Cáceres. 2013. Beta diversity as the variance of 
# community data: dissimilarity coefficients and partitioning. 
# Ecology Letters 16: 951-963. 
#
# License: GPL-2 
# Author:: Pierre Legendre, December 2012, April-May 2013
{
  ### Internal functions
  centre <- function(D,n)
    # Centre a square matrix D by matrix algebra
    # mat.cen = (I - 11'/n) D (I - 11'/n)
  {  One <- matrix(1,n,n)
  mat <- diag(n) - One/n
  mat.cen <- mat %*% D %*% mat
  }
  ###
  BD.group1 <- function(Y, method, save.D, per, n)
  {
    if(method=="profiles") Y = decostand(Y, "total")
    if(method=="hellinger") Y = decostand(Y, "hellinger")
    if(method=="chord") Y = decostand(Y, "norm")
    if(method=="chisquare") Y = decostand(Y, "chi.square")
    #
    s <- scale(Y, center=TRUE, scale=FALSE)^2   # eq. 1
    SStotal <- sum(s)          # eq. 2
    BDtotal <- SStotal/(n-1)   # eq. 3
    if(!per) { SCBD<-apply(s,2,sum)/SStotal }else{ SCBD<-NA }  # eqs. 4a and 4b
    LCBD <- apply(s, 1, sum)/SStotal  # eqs. 5a and 5b
    #
    D <- NA
    if(!per & save.D)   D <- dist(Y)
    #
    out <- list(SStotal_BDtotal=c(SStotal,BDtotal), SCBD=SCBD, LCBD=LCBD, 
                method=method, D=D)
  }
  ###
  BD.group2 <- function(Y, method, sqrt.D, n)
  {
    if(method == "divergence") {
      D = D11(Y)		
      
    } else if(any(method == 
                  c("jaccard","sorensen","ochiai"))) 
    {
      if(method=="jaccard") D = dist.binary(Y, method=1) # ade4 takes sqrt(D)
      if(method=="sorensen")  D = dist.binary(Y, method=5) #ade4 takes sqrt(D)
      if(method=="ochiai") D = dist.binary(Y, method=7) # ade4 takes sqrt(D)
      
    } else if(any(method == 
                  c("manhattan","canberra","whittaker","percentagedifference","wishart"))) 
    {
      if(method=="manhattan") D = vegdist(Y, "manhattan")
      if(method=="canberra")  D = vegdist(Y, "canberra")
      if(method=="whittaker") D = vegdist(decostand(Y,"total"),"manhattan")/2
      if(method=="percentagedifference") D = vegdist(Y, "bray")
      if(method=="wishart")   D = WishartD(Y)
    } else {
      if(method=="modmeanchardiff") D = D19(Y)
      if(method=="kulczynski")  D = vegdist(Y, "kulczynski")
      if(method=="ab.jaccard")  D = chao(Y, coeff="Jaccard", samp=samp)
      if(method=="ab.sorensen") D = chao(Y, coeff="Sorensen", samp=samp)
      if(method=="ab.ochiai")   D = chao(Y, coeff="Ochiai", samp=samp)
      if(method=="ab.simpson")  D = chao(Y, coeff="Simpson", samp=samp)
    }
    #
    if(sqrt.D) D = sqrt(D)
    SStotal <- sum(D^2)/n      # eq. 8
    BDtotal <- SStotal/(n-1)   # eq. 3
    delta1 <- centre(as.matrix(-0.5*D^2), n)   # eq. 9
    LCBD <- diag(delta1)/SStotal               # eq. 10b
    #
    out <- list(SStotal_BDtotal=c(SStotal,BDtotal), LCBD=LCBD, 
                method=method, D=D)
  }
  ###
  ###
  epsilon <- sqrt(.Machine$double.eps)
  method <- match.arg(method, c("euclidean", "manhattan", "modmeanchardiff", "profiles", "hellinger", "chord", "chisquare", "divergence", "canberra", "whittaker", "percentagedifference", "wishart", "kulczynski", "ab.jaccard", "ab.sorensen","ab.ochiai","ab.simpson","jaccard","sorensen","ochiai","none"))
  #
  if(any(method == c("profiles", "hellinger", "chord", "chisquare", "manhattan", "modmeanchardiff", "divergence", "canberra", "whittaker", "percentagedifference", "kulczynski"))) require(vegan)
  if(any(method == c("jaccard","sorensen","ochiai"))) require(ade4)
  #
  if(is.table(Y)) Y <- Y[1:nrow(Y),1:ncol(Y)]    # In case class(Y) is "table"
  n <- nrow(Y)
  if((n==2)&(dist(Y)[1]<epsilon)) stop("Y contains two identical rows, hence BDtotal = 0")
  #
  aa <- system.time({
    if(any(method == 
           c("euclidean", "profiles", "hellinger", "chord", "chisquare","none"))) {
      note <- "Info -- This coefficient is Euclidean"
      res <- BD.group1(Y, method, save.D, per=FALSE, n)
      #
      # Permutation test for LCBD indices, distances group 1
      if(nperm>0) {
        p <- ncol(Y)
        nGE.L = rep(1,n)
        for(iperm in 1:nperm) {
          Y.perm = apply(Y,2,sample)
          res.p <- BD.group1(Y.perm, method, save.D, per=TRUE, n)
          ge <- which(res.p$LCBD+epsilon >= res$LCBD)
          nGE.L[ge] <- nGE.L[ge] + 1
        }
        p.LCBD <- nGE.L/(nperm+1)
      } else { p.LCBD <- NA }
      #
      if(save.D) { D <- res$D } else { D <- NA }
      #
      out <- list(SStotal_BDtotal=res$SStotal_BDtotal, SCBD=res$SCBD, 
                  LCBD=res$LCBD, p.LCBD=p.LCBD, method=method, note=note, D=D)
      
    } else {
      #
      if(method == "divergence") {
        note = "Info -- This coefficient is Euclidean"
      } else if(any(method == c("jaccard","sorensen","ochiai"))) {
        note = c("Info -- This coefficient is Euclidean because dist.binary ",
                 "of ade4 computes it as sqrt(D). Use beta.div with option sqrt.D=FALSE")
      } else if(any(method == 
                    c("manhattan","canberra","whittaker","percentagedifference","wishart"))) {
        if(sqrt.D) {
          note = "Info -- This coefficient, in the form sqrt(D), is Euclidean"
        } else {
          note = c("Info -- For this coefficient, sqrt(D) would be Euclidean", 
                   "Use is.euclid(D) of ade4 to check Euclideanarity of this D matrix")
        }
      } else {
        note = c("Info -- This coefficient is not Euclidean", 
                 "Use is.euclid(D) of ade4 to check Euclideanarity of this D matrix")
      }
      #
      res <- BD.group2(Y, method, sqrt.D, n)
      #
      # Permutation test for LCBD indices, distances group 2
      if(nperm>0) {
        nGE.L = rep(1,n)
        for(iperm in 1:nperm) {
          Y.perm = apply(Y,2,sample)
          res.p <- BD.group2(Y.perm, method, sqrt.D, n)
          ge <- which(res.p$LCBD+epsilon >= res$LCBD)
          nGE.L[ge] <- nGE.L[ge] + 1
        }
        p.LCBD <- nGE.L/(nperm+1)
      } else { p.LCBD <- NA }
      #
      if(sqrt.D) note.sqrt.D<-"sqrt.D=TRUE"  else  note.sqrt.D<-"sqrt.D=FALSE"
      if(save.D) { D <- res$D } else { D <- NA }
      #
      out <- list(SStotal_BDtotal=res$SStotal_BDtotal, LCBD=res$LCBD,  
                  p.LCBD=p.LCBD, method=c(method,note.sqrt.D), note=note, D=D)
    }
    #
  })
  aa[3] <- sprintf("%2f",aa[3])
  if(clock) cat("Time for computation =",aa[3]," sec\n")
  #
  class(out) <- "beta.div"
  out
}

D11 <- function(Y, algo=1)
  #
  # Compute Clark's coefficient of divergence. 
  # Coefficient D11 in Legendre and Legendre (2012, eq. 7.51).
  #
  # License: GPL-2 
  # Author:: Pierre Legendre, April 2011
{
  Y <- as.matrix(Y)
  n <- nrow(Y)
  p <- ncol(Y)
  # Prepare to divide by pp = (p-d) = no. species present at both sites
  Y.ap <- 1 - decostand(Y, "pa")
  d <- Y.ap %*% t(Y.ap)
  pp <- p-d   # n. species present at the two compared sites
  #
  if(algo==1) {   # Faster algorithm
    D <- matrix(0, n, n)
    for(i in 2:n) {
      for(j in 1:(i-1)) {
        num <- (Y[i,]-Y[j,])
        den <- (Y[i,]+Y[j,])
        sel <- which(den > 0)
        D[i,j] = sqrt(sum((num[sel]/den[sel])^2)/pp[i,j])
      }
    }
    #
  } else {   # Slower algorithm 
    D <- matrix(0, n, n)
    for(i in 2:n) {
      for(j in 1:(i-1)) {
        temp = 0
        for(p2 in 1:p) {
          den = Y[i,p2] + Y[j,p2]
          if(den > 0) {
            temp = temp + ((Y[i,p2] - Y[j,p2])/den)^2
          }
        }
        D[i,j] = sqrt(temp/pp[i,j])
      }
    }
    #
  }	
  DD <- as.dist(D)
}

D19 <- function(Y)
  #
  # Compute the Modified mean character difference.
  # Coefficient D19 in Legendre and Legendre (2012, eq. 7.46).
  # Division is by pp = number of species present at the two compared sites
  #
  # License: GPL-2 
  # Author:: Pierre Legendre, April 2011
{
  Y <- as.matrix(Y)
  n <- nrow(Y)
  p <- ncol(Y)
  # Prepare to divide by pp = (p-d) = n. species present at both sites
  Y.ap <- 1 - decostand(Y, "pa")
  d <- Y.ap %*% t(Y.ap)
  pp <- p-d   # n. species present at the two compared sites
  #
  D <- vegdist(Y, "manhattan")
  DD <- as.dist(as.matrix(D)/pp)
}

WishartD <- function(Y)
  #
  # Compute dissimilarity - 1 - Wishart similarity ratio (Wishart 1969).
  #
  # License: GPL-2 
  # Author:: Pierre Legendre, August 2012
{
  CP = crossprod(t(Y))
  SS = apply(Y^2,1,sum)
  n = nrow(Y)
  mat.sq = matrix(0, n, n)
  for(i in 2:n) {
    for(j in 1:(n-1)) { mat.sq[i,j] = CP[i,j]/(SS[i] + SS[j] - CP[i,j]) }
  }
  mat = 1 - as.dist(mat.sq)
}

chao <- function(mat, coeff="Jaccard", samp=TRUE)
  #
  # Compute Chao et al. (2006) abundance-based indices.
  #
  # Arguments -
  # mat = data matrix, species abundances
  # coef = "Jaccard" : modified abundance-based Jaccard index
  #        "Sorensen": modified abundance-based Sørensen index
  #        "Ochiai"  : modified abundance-based Ochiai index
  #        "Simpson" : modified abundance-based Simpson index
  # samp=TRUE : Compute dissimilarities for sample data
  #     =FALSE: Compute dissimilarities for true population data
#
# Details -
# For coeff="Jaccard", the output values are identical to those
# produced by vegan's function vegdist(mat, "chao").
#
# Help received from A. Chao and T. C. Hsieh in July 2012 for the computation  
# of dissimilarities for true population data is gratefully acknowledged.
#
# Reference --
# Chao, A., R. L. Chazdon, R. K. Colwell and T. J. Shen. 2006. 
# Abundance-based similarity indices and their estimation when there 
# are unseen species in samples. Biometrics 62: 361–371.
#
# License: GPL-2 
# Author:: Pierre Legendre, July 2012
{
  require(vegan)
  nn = nrow(mat)
  res = matrix(0,nn,nn)
  if(samp) {   # First for sample data
    for(k in 2:nn) {
      for(j in 1:(k-1)) {
        #cat("k =",k,"  j =",j,"\n")
        v1 = mat[j,]   # Vector 1
        v2 = mat[k,]   # Vector 2
        v1.pa = decostand(v1,"pa")   # Vector 1 in presence-absence form
        v2.pa = decostand(v2,"pa")   # Vector 2 in presence-absence form
        N.j = sum(v1)   # Sum of abundances in vector 1
        N.k = sum(v2)   # Sum of abundances in vector 2
        shared.sp = v1.pa * v2.pa   # Vector of shared species ("pa")
        if(sum(shared.sp) == 0) { 
          res[k,j] = 1
        } else {
          C.j = sum(shared.sp * v1)   # Sum of shared sp. abundances in v1
          C.k = sum(shared.sp * v2)   # Sum of shared sp. abundances in v2
          # a1.j = sum(shared.sp * v1.pa)
          # a1.k = sum(shared.sp * v2.pa)
          a1.j = length(which((shared.sp * v2) == 1)) # Singletons in v2
          a1.k = length(which((shared.sp * v1) == 1)) # Singletons in v1
          a2.j = length(which((shared.sp * v2) == 2)) # Doubletons in v2
          if(a2.j == 0) a2.j <- 1
          a2.k = length(which((shared.sp * v1) == 2)) # Doubletons in v1
          if(a2.k == 0) a2.k <- 1
          # S.j = sum(v1[which(v2 == 1)]) # Sum abund. in v1 for singletons in v2
          # S.k = sum(v2[which(v1 == 1)]) # Sum abund. in v2 for singletons in v1
          sel2 = which(v2 == 1)
          sel1 = which(v1 == 1)
          if(length(sel2)>0) S.j = sum(v1[sel2]) else S.j = 0
          if(length(sel1)>0) S.k = sum(v2[sel1]) else S.k = 0
          
          U.j = (C.j/N.j) + ((N.k-1)/N.k) * (a1.j/(2*a2.j)) * (S.j/N.j) # Eq. 11
          if(U.j > 1) U.j <- 1
          U.k = (C.k/N.k) + ((N.j-1)/N.j) * (a1.k/(2*a2.k)) * (S.k/N.k) # Eq. 12
          if(U.k > 1) U.k <- 1
          
          if(coeff == "Jaccard") {                     # "Jaccard"
            res[k,j] = 1 - (U.j*U.k/(U.j + U.k - U.j*U.k))
          } else if(coeff == "Sorensen") {         # "Sorensen"
            res[k,j] = 1 - (2*U.j*U.k/(U.j + U.k))
          } else if(coeff == "Ochiai") {           # "Ochiai"
            res[k,j] = 1 - (sqrt(U.j*U.k))
          } else if(coeff == "Simpson") { 
            # Simpson (1943), or Lennon et al. (2001) in Chao et al. (2006)
            res[k,j] = 1 -
              (U.j*U.k/(U.j*U.k+min((U.j-U.j*U.k),(U.k-U.j*U.k))))
          } else { # 
            stop("Incorrect coefficient name")
          }
        }
      }
    }
    
  } else {   # Now for complete population data
    
    for(k in 2:nn) {
      for(j in 1:(k-1)) {
        v1 = mat[j,]   # Vector 1
        v2 = mat[k,]   # Vector 2
        v1.pa = decostand(v1,"pa")   # Vector 1 in presence-absence form
        v2.pa = decostand(v2,"pa")   # Vector 2 in presence-absence form
        shared.sp = v1.pa * v2.pa    # Vector of shared species ("pa")
        if(sum(shared.sp) == 0) { 
          res[k,j] = 1
        } else {
          N1 = sum(v1)   # Sum of abundances in vector 1
          N2 = sum(v2)   # Sum of abundances in vector 2
          U = sum(shared.sp * v1)/N1   # Sum of shared sp. abundances in v1
          V = sum(shared.sp * v2)/N2   # Sum of shared sp. abundances in v2
          
          if(coeff == "Jaccard") {                     # "Jaccard"
            res[k,j] = 1 - (U*V/(U + V - U*V))
          } else if(coeff == "Sorensen") {         # "Sorensen"
            res[k,j] = 1 - (2*U*V/(U + V))
          } else if(coeff == "Ochiai") {           # "Ochiai"
            res[k,j] = 1 - (sqrt(U*V))
          } else if(coeff == "Simpson") { # "Simpson"
            res[k,j] = 1 - (U*V/(U*V+min((U-U*V),(V-U*V)))) # Eq. ?
          } else { # 
            stop("Incorrect coefficient name")
          }
        }
      }
    }
  }
  res <- as.dist(res)
}

######## End of beta.div function
##################################################################################################

####BETA.DIV.COMP SOURCE FUNCTION- for partitioning beta diversity into components (spatial)
##################################################################################################
beta.div.comp <- function(mat, coef="J", quant=FALSE, save.abc=FALSE)
  #
  # Description --
  # 
  # Podani-family and Baselga-family decompositions of the Jaccard and Sørensen 
  # dissimilarity coefficients into replacement and richness difference 
  # components, for species presence-absence or abundance data, as described  
  # in Legendre (2014).
  #
  # Usage --
  #
  # beta.div.comp(mat, coef="J", quant=FALSE, save.abc=FALSE)
#
# Arguments --
#
# mat : Data in matrix or data.frame form.
# coef : Family of coefficients to be computed --
#        "S" or "Sorensen": Podani family, Sørensen-based indices
#        "J" or "Jaccard" : Podani family, Jaccard-based indices
#        "BS" : Baselga family, Sørensen-based indices
#        "BJ" : Baselga family, Sørensen-based indices
#        "N" : Podani & Schmera (2011) relativized nestedness index.
#        The quantitative form in Sørensen family is the percentage difference.
#        The quantitative form in the Jaccard family is the Ruzicka index.
#
# quant=TRUE : Compute the quantitative form of replacement, nestedness and D.
#      =FALSE: Compute the presence-absence form of the coefficients.
# save.abc=TRUE : Save the matrices of parameters a, b and c used in the
#      presence-absence calculations.
#
# Details --
#
#    For species presence-absence data, the distance coefficients are 
# Jaccard=(b+c)/(a+b+c) and Sørensen=(b+c)/(2*a+b+c) with usual abc notation.
#
#    For species abundance data, the distance coefficients are 
# the Ruzicka index = (B+C)/(A+B+C) and Odum's percentage difference 
# (incorrectly called Bray-Curtis) = (B+C)/(2A+B+C), where  
# A = sum of the intersections (or minima) of species abundances at two sites,
# B = sum at site 1 minus A, C = sum at site 2 minus A.
#
#    The binary (quant=FALSE) and quantitative (quant=TRUE) forms of the S and  
# J indices return the same values when computed for presence-absence data.
#
# Value --
#
# repl : Replacement matrix, class = 'dist'.
# rich : Richness/abundance difference or nestedness matrix, class = 'dist'.
#        With options "BJ", "BS" and "N", 'rich' contains nestedness indices.
#        With option "N", the 'repl' and 'rich' values do not add up to 'D'.
# D    : Dissimilarity matrix, class = 'dist'.
# part : Beta diversity partitioning -- 
#        1. Total beta div. = sum(D.ij)/(n*(n-1)) (Legendre & De Cáceres 2013)
#        2. Total replacement diversity 
#        3. Total richness difference diversity (or nestedness)
#        4. Total replacement div./Total beta div.
#        5. Total richness difference div. (or nestedness)/Total beta div.
# Note : Name of the dissimilarity coefficient.
#
# References --
#
# Baselga, A. (2010) Partitioning the turnover and nestedness components of beta 
# diversity. Global Ecology and Biogeography, 19, 134–143.
#
# Baselga, A. (2012) The relationship between species replacement, dissimilarity 
# derived from nestedness, and nestedness. Global Ecology and Biogeography, 21, 
# 1223–1232. 
#
# Baselga, A. (2013) Separating the two components of abundance-based 
# dissimilarity: balanced changes in abundance vs. abundance gradients. Methods 
# in Ecology and Evolution, 4, 552–557.
#
# Carvalho, J.C., Cardoso, P., Borges, P.A.V., Schmera, D. & Podani, J. (2013)
# Measuring fractions of beta diversity and their relationships to nestedness: 
# a theoretical and empirical comparison of novel approaches. Oikos, 122, 
# 825–834.
#
# Legendre, P. 2014. Interpreting the replacement and richness difference   
# components of beta diversity. Global Ecology and Biogeography 23: (in press).
#
# Podani, J., Ricotta, C. & Schmera, D. (2013) A general framework for analyzing 
# beta diversity, nestedness and related community-level phenomena based on 
# abundance data. Ecological Complexity, 15, 52-61.
#
# Podani, J. & Schmera, D. 2011. A new conceptual and methodological framework 
# for exploring and explaining pattern in presence-absence data. Oikos, 120, 
# 1625–1638.
#
# License: GPL-2 
# Author:: Pierre Legendre
{
  coef <- pmatch(coef, c("S", "J", "BS", "BJ", "N"))
  if(coef==5 & quant) stop("coef='N' and quant=TRUE: combination not programmed")
  mat <- as.matrix(mat)
  n <- nrow(mat)
  if(is.null(rownames(mat))) noms <- paste("Site",1:n,sep="")
  else noms <- rownames(mat)
  #
  if(!quant) {      # Binary data provided, or make the data binary
    if(coef==1) form="Podani family, Sorensen" 
    if(coef==2) form="Podani family, Jaccard"
    if(coef==3) form="Baselga family, Sorensen" 
    if(coef==4) form="Baselga family, Jaccard"
    if(coef==5) form="Podani & Schmera (2011) relativized nestedness"
    mat.b <- ifelse(mat>0, 1, 0)
    a <- mat.b %*% t(mat.b)
    b <- mat.b %*% (1 - t(mat.b))
    c <- (1 - mat.b) %*% t(mat.b)
    min.bc <- pmin(b,c)
    #
    if(coef==1 || coef==2) {
      repl <- 2*min.bc   # replacement, turnover, beta-3
      rich <- abs(b-c)   # nestedness, richness diff., beta-rich
      #
      # Add the denominators
      if(coef==1) {                # Sørensen-based components
        repl <- repl/(2*a+b+c)
        rich <- rich/(2*a+b+c)
        D <- (b+c)/(2*a+b+c)
      } else if(coef==2) {     # Jaccard-based components
        repl <- repl/(a+b+c)
        rich <- rich/(a+b+c)
        D <- (b+c)/(a+b+c)
      }
    } else if(coef==3) {     # Baselga 2010 components based on Sørensen
      D <- (b+c)/(2*a+b+c)             # Sørensen dissimilarity
      repl <- min.bc/(a+min.bc)        # replacement, turnover
      rich <- D-repl                   # richness difference
      
    } else if(coef==4) {      # Baselga 2012 components based on Jaccard
      D <- (b+c)/(a+b+c)               # Jaccard dissimilarity
      repl <- 2*min.bc/(a+2*min.bc)    # replacement, turnover
      rich <- D-repl                   # richness difference
    } else if(coef==5) {      # rich = Podani N = nestdness based on Jaccard
      repl <- 2*min.bc/(a+b+c)
      D <- (b+c)/(a+b+c)
      rich <- matrix(0,n,n)
      for(i in 2:n) {
        for(j in 1:(i-1)) {
          aa = a[i,j]; bb = b[i,j]; cc = c[i,j]
          if(a[i,j] == 0)  rich[i,j] <- 0  
          else  rich[i,j] <- (aa + abs(bb-cc))/(aa+bb+cc) 
        }
      }
    }
    
    rownames(repl) <- rownames(rich) <- rownames(D) <- noms
    D <- as.dist(D)
    repl <- as.dist(repl)
    rich <- as.dist(rich)
    total.div <- sum(D)/(n*(n-1))
    repl.div <- sum(repl)/(n*(n-1))
    rich.div <- sum(rich)/(n*(n-1))
    part <- c(total.div,repl.div,rich.div,repl.div/total.div,rich.div/total.div)
    #
    if(save.abc) {
      res <- list(repl=repl, rich=rich, D=D, part=part, Note=form, 
                  a=as.dist(a), b=as.dist(b), c=as.dist(c))
    } else { 
      res <- list(repl=repl, rich=rich, D=D, part=part, Note=form)
    }
    #
  } else {      # Quantitative data
    # Calculations based on individuals.within.species
    if(coef==1) form<-"Podani family, percentage difference" 
    if(coef==2) form<-"Podani family, Ruzicka"
    if(coef==3) form<-"Baselga family, percentage difference"
    if(coef==4) form<-"Baselga family, Ruzicka"
    # Baselga (2013) notation:
    # A = W = sum of minima in among-site comparisons
    # B = site.1 sum - W = K.1 - W
    # C = site.2 sum - W = K.2 - W
    K <- vector("numeric", n)   # site (row) sums
    W <- matrix(0,n,n)
    repl <- matrix(0,n,n)
    rich <- matrix(0,n,n)
    D <- matrix(0,n,n)
    rownames(repl) <- rownames(rich) <- rownames(D) <- noms
    K <- apply(mat,1,sum)         # Row sums
    for(i in 2:n) for(j in 1:(i-1)) W[i,j] <- sum(pmin(mat[i,], mat[j,]))
    #
    # Quantitative extensions of the S and J decompositions
    for(i in 2:n) {
      for(j in 1:(i-1)) {
        repl[i,j] <- 2*(min(K[i],K[j])-W[i,j]) # 2*min(B,C)
        rich[i,j] <- abs(K[i]-K[j])            # abs(B-C)
      }
    }
    #
    # Add the denominators
    if(coef==1) {         # Sørensen-based (% difference) components
      for(i in 2:n) {
        for(j in 1:(i-1)) {	                        # Baselga 2013 notation:
          repl[i,j] <- repl[i,j]/(K[i]+K[j])          # 2min(B,C)/(2A+B+C)
          rich[i,j] <- rich[i,j]/(K[i]+K[j])          # abs(B-C)/(2A+B+C)
          # cat(K[i], K[j], W[i,j],"\n")
          D[i,j] <- (K[i]+K[j]-2*W[i,j])/(K[i]+K[j])  # (B+C)/(2A+B+C)
        }
      }
    } else if(coef==2) {    # Jaccard-based (Ruzicka) components
      for(i in 2:n) {
        for(j in 1:(i-1)) {                         # Baselga 2013 notation:
          repl[i,j] <- repl[i,j]/(K[i]+K[j]-W[i,j])   # 2min(B,C)/(A+B+C)
          rich[i,j] <- rich[i,j]/(K[i]+K[j]-W[i,j])   # abs(B-C)/(A+B+C)
          # cat(K[i], K[j], W[i,j],"\n")
          D[i,j]<-(K[i]+K[j]-2*W[i,j])/(K[i]+K[j]-W[i,j]) # (B+C)/(A+B+C)
        }
      }
    }
    #
    # Baselga (2013): quantitative extensions of the Baselga (2010) indices
    if(coef==3) {   # Baselga (2013) indices decomposing percentage difference
      for(i in 2:n) {
        for(j in 1:(i-1)) {
          repl[i,j] <- (min(K[i],K[j])-W[i,j])/min(K[i],K[j])
          rich[i,j] <- abs(K[i]-K[j])*W[i,j]/((K[i]+K[j])*min(K[i],K[j]))
          # cat(K[i], K[j], W[i,j],"\n")
          D[i,j] <- (K[i]+K[j]-2*W[i,j])/(K[i]+K[j])
        }
      }
    }	
    if(coef==4) {   # Decomposing Ruzicka in the spirit of Baselga 2013
      for(i in 2:n) {
        for(j in 1:(i-1)) {
          repl[i,j] <- 
            2*(min(K[i],K[j])-W[i,j])/(2*min(K[i],K[j])-W[i,j])
          rich[i,j] <- abs(K[i]-K[j])*W[i,j]/
            ((K[i]+K[j]-W[i,j])*(2*min(K[i],K[j])-W[i,j]))
          # cat(K[i], K[j], W[i,j],"\n")
          D[i,j] <- (K[i]+K[j]-2*W[i,j])/(K[i]+K[j]-W[i,j])
        }
      }
    }	
    #
    repl <- as.dist(repl)
    rich <- as.dist(rich)
    D <- as.dist(D)
    repl.div <- sum(repl)/(n*(n-1))
    rich.div <- sum(rich)/(n*(n-1))
    total.div <- sum(D)/(n*(n-1))
    part <- c(total.div,repl.div,rich.div,repl.div/total.div,rich.div/total.div)
    #
    res <- list(repl=repl, rich=rich, D=D, part=part, Note=form)
  }
  res
}
########End of beta.div.comp function
##################################################################################################

####LCBD COMP SOURCE FUNCTION- partition LCBD into components (spatial)
##################################################################################################
##LCBD.comp() source function
LCBD.comp <- function(x, sqrt.x=TRUE)
{
  ### Internal function
  centre <- function(D,n)
    # Centre a square matrix D by matrix algebra
    # mat.cen = (I - 11'/n) D (I - 11'/n)
  {
    One <- matrix(1,n,n)
    mat <- diag(n) - One/n
    mat.cen <- mat %*% D %*% mat
  }
  ###
  n <- nrow(as.matrix(x))
  if(sqrt.x) {
    # x = sqrt(x)
    SStotal <- sum(x)/n # eq. 8
    BDtotal <- SStotal/(n-1) # eq. 3
    G <- centre(as.matrix(-0.5*x), n) # Gower-centred matrix
  } else {
    SStotal <- sum(x^2)/n # eq. 8
    BDtotal <- SStotal/(n-1) # eq. 3
    G <- centre(as.matrix(-0.5*x^2), n) # Gower-centred matrix
  }
  LCBD <- diag(G)/SStotal # Legendre & De Caceres (2013), eq. 10b
  out <- list(SStotal_BDtotal=c(SStotal,BDtotal), LCBD=LCBD, D=x)
} 

# Arguments --
#
# x : D or beta diversity component matrix, class=dist.
# sqrt.x : Take sqrt() of components before computing LCBD.comp. Use
# sqrt.x=TRUE for the replacement and richness/abundance difference indices
# computed by beta.div.comp(), as well as for the corresponding D matrices.

########End of LCBD.comp function 
##################################################################################################

####DECOMPOSE.D2 TEMPORAL BD SOURCE FUNCTION- for computing temporal beta diversity and gain and loss components
##################################################################################################
decompose.D2 <- function(Y1, Y2, den.type=2)
  # Compare two surveys:
  # Decompose the Ruzicka and percentage difference dissimilarities into A, B and C.
  #
  # Parameters --
  # 
  # Y1 : First survey data with sites in rows and species in columns. 
  # Y2 : Second survey data with sites in rows and species in columns. 
  # The sites and species must be the same in the two tables, and in the same order.
  # The files may contain species presence-absence or quantitative abundance data.
  # 
  # den.type -- Denominator type for the indices
#     1 : (A+B+C) as in the Ruzicka dissimilarity
#     2 : (2A+B+C) as in the Percentage difference dissimilarity (alias Bray-Curtis)
#
# Value (output of the function) --
#
# mat1 : A, B and C results, with sites in rows and A, B and C in columns.
# mat2 : A,B,C,D divided by a denominator [either (A+B+C) or (2A+B+C)]; D = (B+C).
#
# Details: Numerical results in output matrices --
# A -- aj is the part of the abundance of species j that is common to the two survey vectors: aj = min(y1j, y2j). A is the sum of the aj values for all species in the functional group under study.
# B -- bj is the part of the abundance of species j that is higher in survey 1 than in survey 2: bj = y1j – y2j. B is the sum of the bj values for all species in the functional group under study.
# C -- cj is the part of the abundance of species j that is higher in survey 2 than in survey 1: cj = y2j – y1j. C is the sum of the cj values for all species in the functional group under study.
#
# Example --
# test1 = matrix(runif(50,0,100),10,5)
# test2 = matrix(runif(50,0,100),10,5)
# (res = decompose.D2(test1, test2, den.type=1))# License: GPL-2
#
# License: GPL-2
# Author:: Pierre Legendre
{
  ### Internal function
  den <- function(A,B,C,den.type) if(den.type==1) den=(A+B+C) else den=(2*A+B+C) #specifying the demoninator to use
  ### End internal function
  #
  n = nrow(Y1)
  p = ncol(Y1)
  if(nrow(Y2)!=n) stop("The data tables do not have the same number of rows") #warning messages in case matrices don't match
  if(ncol(Y2)!=p) stop("The data tables do not have the same number of columns")
  #
  ABC = c("A","B","C")        # A = similarity #column headings for output matrix 1
  ABCD = c("A","B","C","D")   # D = dissimilarity #column headings for output matrix 2
  #
  mat1 = matrix(NA,n,3) #output matrix 1
  colnames(mat1) = ABC
  mat2 = matrix(NA,n,4) #output matrix 2 
  colnames(mat2) = ABCD
  if(!is.null(rownames(Y1))) { 
    rownames(mat1)=rownames(mat2)=rownames(Y1) 
  } else {
    rownames(mat1)=rownames(mat2)=paste("Site",1:n,sep=".") }
  #
  for(i in 1:n) {
    YY = rbind(Y1[i,], Y2[i,])
    A = sum(apply(YY,2,min))
    tmp = YY[1,] - YY[2,]
    B = sum(tmp[tmp>0])
    C = -sum(tmp[tmp<0])
    D = B+C
    mat1[i,] = c(A,B,C)
    mat2[i,] = c(A,B,C,D)/den(A,B,C,den.type) #uses the internal function from above
  }
  #
  list(mat1=mat1, mat2=mat2)
}
########End of decompose.D2 function 
##################################################################################################

####PAIRED.DIFF2_2.R TEMPORAL EXCEPTIONAL SITES FUNCTION (with the BDC matrix)- calculating significant temporal change
##################################################################################################
paired.diff2 <- function(mat1,mat2,method="hellinger", pa.tr=FALSE, nperm=99, permute.sp=1, BCD=TRUE, replace=FALSE, clock=FALSE)
  #
  # Compute and test differences between pairs of data vectors of observations at T1 and T2.
  #
  # Test hypothesis (H0) that an object is not exceptionally different between T1 and T2.
  # Example in palaeoecology: ancient and modern diatom communities in sediment cores.
  # Example in sequence data: "regions" before and after treatment.
  # 
  # Arguments --
  # mat1, mat2: two matrices (or data.frames) with the same number of rows and columns. 
  # The rows must correspond to the same objects (e.g. sites) and the colums to the same 
  # variables (e.g. species).
# method={"hellinger", "chord", "ruzicka", "%difference", "euclidean"}. 
#   Methods {"hellinger", "chord"} are obtained by transformation of the  
#     species data followed by calculation of the Euclidean distance. These distances 
#     have the Euclidean property. 
#     If pa.tr=TRUE, sqrt(2)*sqrt(1-Ochiai) is computed.
#   Methods {"ruzicka", "%difference"} are obtained by computing a dissimilarity function. 
#     It is recommended to take the square root of these dissimilarities before computing 
#     ordinations by principal coordinate analysis. However, that precaution is not 
#     important here; the results of the permutation tests will be the same for these
#     dissimilarities square-rooted or not.
#     If pa.tr=TRUE, either the Jaccard or the Sørensen coefficient is computed.
# pa.tr=FALSE: do NOT transform the data to presence-absence.
#      =TRUE : transform the data to binary (i.e. presence-absence) form.
# nperm = number of permutations for the permutation test.
##
# This version of the function contains three permutation methods --
# permute.sp=1 : permute data separately in each column, both matrices in the same way.
#           =2 : permute data separately in each column. Do not force the permutations to 
#  			 start at the same point in the two matrices.
#           =3 : permute entire rows in each matrix separately (suggestion D. Borcard).
##
# BCD=TRUE  : Compute and save the B and C components of the %difference and Ruzicka D.
#             For the %difference, they are expressed as B/(2A+B+C) and C/(2A+B+C).
#             For the Ruzicka D, they are expressed as B/(A+B+C) and C/(A+B+C).
# BCD=FALSE : Do not compute the components. BCD=FALSE for D other than %diff and Ruzicka.
# replace=FALSE : sampling without replacement for regular permutation test.
#        =TRUE  : sampling with replacement. The testing method is then bootstrapping.
# clock=FALSE : Do not print the computation time.
#      =TRUE  : Print the time (in sec) used for computation.
#
# Details --
# H0: in each matrix (e.g. each time), the sites do not differ in species composition. 
#     They only differ by random sampling of each species' statistical population.
# H1: Some sites are exceptionally different between T1 and T2.
# 
# The randomization procedures are the following:
# 1. In each matrix, the original values (e.g. species abundances) are permuted at random, 
# independently in each column. Permutation of the two matrices is started with the same 
# random seed, so that the values in each column (e.g. species) are permuted in the same 
# way in mat1.perm and mat2.perm. 
# 2. The transformation, if any, is recomputed on the permuted data matrices. This is 
# necessary to make sure that the permuted data are transformed in the same way as the 
# initial data, with row sums or row lengths of 1. In this way, the D of the permuted data 
# will be comparable to the reference D.
# 3. The distances between T1 and T2 are recomputed, for each site separately.
#
# For presence-absence data, this function computes the binary forms of the quantitative 
# coefficients listed under the 'method' parameter. The "hellinger" and "chord" 
# transformations produce the Ochiai distance, or more precisely: 
# D.Hellinger = D.chord = sqrt(2) * sqrt(1 - S.Ochiai) 
# where "S.Ochiai" designates the Ochiai similarity coefficient.  
# The "%difference" dissimilarity produces (1 – S.Sørensen) 
# whereas the "ruzicka" dissimilarity produces (1 – S.Jaccard).
#
# Community composition data could be log-transformed prior to analysis. Only the 
# Euclidean distance option should be used with log-transformed data. It is meaningless to 
# subject log-transformed data to the {"hellinger", "chord"} transformations 
# available in this function. - One can use either the log(y+1 transformation (log1p() 
# function of {base}), or Anderson et al. (2006) special log transformation available in 
# {vegan}: decostand(mat, "log", logbase=10).
#
# Value --
# A list containing the vector of distances between T1 and T2 for each object and a 
# corresponding vector of p-values. The significant p-values (e.g. p.dist ≤ 0.05) indicate 
# exceptional objects for the difference of their species composition. The p-values should be corrected for multiple testing using function p.adjust() of {stats}. A good general choice is method="holm", which is the default option of the function.
# An output table containing B, C and D.
#
# Author:: Pierre Legendre
# License: GPL-2 
{
  ### Internal functions
  dissim <- function(mat1, mat2, n, method, tr=TRUE, BCD, ref)
  {
    vecD = vector(mode="numeric",length=n)
    if(BCD) { 
      vecB = vector(mode="numeric",length=n)
      vecC = vector(mode="numeric",length=n)
      vecD = vector(mode="numeric",length=n)
    } else { vecB=NA; vecC=NA; vecD=NA }
    #
    # Compute the dissimilarity between T1 and T2 for each object (site)
    # 1. If method = {"hellinger", "chord"}, tr is TRUE
    if(tr) for(i in 1:n) vecD[i] = dist(rbind(mat1[i,], mat2[i,]))
    #
    # 2. Compute the Euclidean distance
    if(method == "euclidean")  
      for(i in 1:n) vecD[i] = dist(rbind(mat1[i,], mat2[i,])) 
      # 3. Compute the Ruzicka or %difference dissimilarity 
      if(method == "ruzicka") dissimil=1       # Quantitative form of Jaccard
      if(method == "%difference") dissimil=2   # Quantitative form of Sørensen
      if(any(method == c("ruzicka", "%difference"))) { 
        for(i in 1:n) {
          tmp = RuzickaD(mat1[i,], mat2[i,], method=method, BCD=BCD, ref=ref) 
          if(BCD) {
            vecB[i] <- tmp$B
            vecC[i] <- tmp$C }
          vecD[i] <- tmp$D
        }
      }
      # Alternative method (not used here) to compute the %difference dissimilarity:
      #	for(i in 1:n) vecD[i] = vegdist(rbind(mat1[i,], mat2[i,]), "bray")         #Slower
      list(vecB=vecB, vecC=vecC, vecD=vecD)
  }
  ###
  transform <- function(mat, method)
  {
    if(method=="hellinger") mat = decostand(mat, "hellinger")
    if(method=="chord")     mat = decostand(mat, "norm")
    mat
  }
  ### End internal functions
  ###
  A <- system.time({
    
    epsilon <- sqrt(.Machine$double.eps)
    method <- match.arg(method, c("euclidean", "hellinger", "chord", "ruzicka", "%difference")) 
    n = nrow(mat1)
    p = ncol(mat1)
    if((nrow(mat2)!=n) | (ncol(mat2)!=p)) stop("The matrices are not of the same size.")
    #
    if(pa.tr) {
      mat1 <- ifelse(mat1>0, 1, 0)
      mat2 <- ifelse(mat2>0, 1, 0) }
    if(any(method == c("hellinger", "chord"))) {
      tr <- TRUE
      require(vegan)
    } else { tr <- FALSE }
    if( (any(method == c("ruzicka", "%difference"))) & BCD) { 
      BCD.mat <- matrix(0,n,3)
      if(method=="ruzicka")    colnames(BCD.mat) <- 
          c("B/(A+B+C)","C/(A+B+C)","D=(B+C)/(A+B+C)")
      if(method=="%difference") colnames(BCD.mat) <- 
          c("B/(2A+B+C)","C/(2A+B+C)","D=(B+C)/(2A+B+C)")
      rownames(BCD.mat) <- paste("Obj",1:n,sep=".")
    } else {
      BCD <- FALSE 
      BCD.mat <- NA }
    ###
    # 1. Compute the reference D for each object from corresponding vectors in the 2 matrices.
    if(tr) { 
      tmp <-dissim(transform(mat1,method), transform(mat2,method),n,method,tr,BCD,ref=FALSE)
    } else { tmp <- dissim(mat1, mat2, n, method, tr, BCD, ref=TRUE) }
    vecD.ref <- tmp$vecD
    if(BCD) { BCD.mat[,1]<-tmp$vecB ; BCD.mat[,2]<-tmp$vecC ; BCD.mat[,3]<-tmp$vecD }
    ###
    if(permute.sp!=3) {   # Permute the data separately in each column.
      # 2. Permutation methods 1 and 2 --
      # Permute *the raw data* by columns. Permute the two matrices in the same way, saving the seed before the two sets of permutations through sample(). 
      # Permutation test for each distance in vector D
      # seed: seed for random number generator, used by the permutation function 
      #       sample(). It is reset to that same value before permuting the values in the  
      #       columns of the second matrix. 
      if(nperm>0) {
        nGE.D = rep(1,n)
        for(iperm in 1:nperm) {
          BCD <- FALSE
          if(permute.sp==1) {    # Permutation methods 1
            seed <- ceiling(runif(1,max=100000))
            # cat("seed =",seed,'\n')
            set.seed(seed)
            mat1.perm <- apply(mat1,2,sample)
            set.seed(seed)
            mat2.perm <- apply(mat2,2,sample)
          } else {  # Permutation methods 2 - Do not force the permutations 
            # to start at the same point in the two matrices.
            mat1.perm <- apply(mat1,2,sample)
            mat2.perm <- apply(mat2,2,sample)
          }
          # 3. Recompute transformations of the matrices and the D values of the paired vectors.
          if(tr) { tmp <- dissim(transform(mat1.perm,method), 
                                 transform(mat2.perm,method), n, method, tr, BCD, ref=FALSE)
          } else { tmp <- dissim(mat1.perm, mat2.perm, n, method, tr, BCD, ref=FALSE) }
          vecD.perm <- tmp$vecD
          ge <- which(vecD.perm+epsilon >= vecD.ref)
          nGE.D[ge] <- nGE.D[ge] + 1
        }
        # 4. Compute the p-value associated with each distance.
        p.dist <- nGE.D/(nperm+1)
      } else { p.dist <- NA }   # if nperm=0
      
    } else if(permute.sp==3) {   
      # 2.bis  Permutation method 3 -- 
      # Permute entire rows in each matrix separately.
      if(nperm>0) {
        seed <- ceiling(runif(1,max=100000))
        set.seed(seed)
        nGE.D = rep(1,n)
        for(iperm in 1:nperm) {
          BCD <- FALSE
          mat1.perm <- mat1[sample(n),]
          mat2.perm <- mat2[sample(n),]
          #
          # 3.bis Recompute the D values of the paired vectors.
          if(tr) { tmp <- dissim(transform(mat1.perm,method), 
                                 transform(mat2.perm,method), n, method, tr, BCD, ref=FALSE)
          } else { tmp <- dissim(mat1.perm, mat2.perm, n, method, tr, BCD, ref=FALSE) }
          vecD.perm <- tmp$vecD
          ge <- which(vecD.perm+epsilon >= vecD.ref)
          nGE.D[ge] <- nGE.D[ge] + 1
        }
        # 4.bis Compute the p-value associated with each distance.
        p.dist <- nGE.D/(nperm+1)
      } else { p.dist <- NA }   # if nperm=0
    }
    p.adj <- p.adjust(p.dist,"holm")
  })
  A[3] <- sprintf("%2f",A[3])
  if(clock) cat("Computation time =",A[3]," sec",'\n')
  #
  list(vecD.ref=vecD.ref, p.dist=p.dist, p.adj=p.adj, BCD.mat=BCD.mat)
}

RuzickaD <- function(vec1, vec2, method="ruzicka", BCD=FALSE, ref=TRUE)
  #
  # Compute the Ruzicka dissimilarity (quantitative form of the Jaccard dissimilarity)
  # or the percentage difference (quantitative form of the Sørensen dissimilarity).
  # A single dissimilarity is computed because there are only two vectors in this function.
  #
  # Arguments --
  # vec1, vec2 : data vectors (species abundance or presence-absence data)
  # method == c("ruzicka", "%difference")
  # BCD=TRUE  : Compute and save the B and C components of the %difference and Ruzicka D.
  #             For the %difference, they are B/(2A+B+C), C/(2A+B+C), D/(2A+B+C).
  #             For the Ruzicka D, they are B/(A+B+C), C/(A+B+C), D/(A+B+C).
# BCD=FALSE : Do not compute the components. BCD=FALSE for D other than %diff and Ruzicka.
# ref=TRUE  : Compute the reference values of D, B and C
#    =FALSE : Under permutation, compute only the value of D. Use separate code (shorter).
#
# License: GPL-2 
# Author:: Pierre Legendre, April 2015
{
  # An algorithm applicable to matrices Y containing two data vectors only
  #
  A <- sum(pmin(vec1, vec2))          # A = sum of minima from comparison of the 2 vectors
  sum.Y <- sum(vec1, vec2)            # Sum of all values in the two vectors, (2A+B+C)
  #
  if(ref) {    # Compute the reference values of statistics D, B and C
    tmp = vec1 - vec2
    B = sum(tmp[tmp>0])                 # Sum of the species losses between T1 and T2
    C = -sum(tmp[tmp<0])                # Sum of the species gains between T1 and T2
    D = B+C                             # Dissimilarity
    
    # Under permutation, compute only the value of D. - Shorter computation time.
  } else { 
    D <- sum.Y-2*A                      # (B+C)
  }
  # Compute the denominator (den) of the Ruzicka or %difference index
  if(method == "ruzicka") { den <-(sum.Y-A)  # den = (A+B+C)
  } else { den <- sum.Y }                # den = (2A+B+C)
  if(!BCD) { B <- NA ; C <- NA }
  list(B=B/den, C=C/den, D=D/den)
}

# Examples -- 
# data(mite)
# res1 = paired.diff(mite[1:10,],mite[61:70,],method="hellinger",nperm=999,permute.sp=1)
# Computation time = 1.971000  sec 
# res2 = paired.diff(mite[1:10,],mite[61:70,],method="hellinger",nperm=999,permute.sp=3)

########End of paired.diff2 function 
##################################################################################################


##################################################################################################
####DATA SET-UP                                                                                  #
##################################################################################################

####Master abiotic/site data####
##General data
nla.data2<- read.csv(file.choose()) #nla2007_lakes_topbotsamples_UpdateNov2014_V2.csv
#McGill/PhD Chapters and projects/EPA National lakes assessment/Chapter 2/ NLA data
#Only cores greater than 30cm (249, 186 of those High confidence)
#10 cores that could potentially have been excluded due to radiometric dating results are still included here. 

##Subsetting of data required for matching/merging dataframes

#In nla.data2 --> retain only high confidence cores
nla.data2.highC<- subset(nla.data2, HIGH_C_CORE == "YES", drop=T) #372 rows (186 top, 186 bottom)

#In nla.data2.highC --> subset into surface and historical sediments 
nla.highC.surf<- subset(nla.data2.highC, SED_TYPE == "Surface") #Surface/top sediments #186
nla.highC.hist<- subset(nla.data2.highC, SED_TYPE == "Historical") #Historical/bottom sediments #186

####Add water quality data#### 
##Extract water quality data columns: 
extradat2<- read.csv(file.choose()) #Beaulieu&Taranu2014_NLAdataset.csv
#File path: C:\Users\Winegardner\Documents\MCGILL\PhD chapters and projects\EPA National lakes Assessment\Chapter 2
#All visit 1

#Keep only SECMEAN, CHLA, PTL, NTL, MEAN_T, COND, PH_FIELD
waterqual.dat<- as.data.frame(subset(extradat2, select=c("SITE_ID", "SECMEAN", "CHLA", "PTL", "NTL", "MEAN_T", "COND", "PH_FIELD")))

#Remove sites from waterqual.dat that aren't in nla.highC.surf (using surface as template)
waterqual.dat.selected<- (waterqual.dat$SITE_ID %in% nla.highC.surf$SITE_ID)
waterqual.dat.red<- waterqual.dat[waterqual.dat.selected,] #ok, same sites(CHECK), but in different order than in nla.highC.surf
#Need to order so that will be able to copy env distance into a column. 

#Re-order so can check that same- yes. 
waterqual.dat.red<- waterqual.dat.red[ order(match(waterqual.dat.red$SITE_ID, nla.highC.surf$SITE_ID)), ]#ok

#Add WSA_ECOREGION
waterqual.dat.red<- as.data.frame(cbind(waterqual.dat.red, nla.highC.surf$WSA_ECOREGION))
colnames(waterqual.dat.red)[9]<- 'WSA_ECOREGION'

####Spatial distances between sites#### 
##Use latitude and longitude from nla.highC.surf $LONG_DD, $LAT_DD
latlong.coords<- as.data.frame(subset(nla.highC.surf, select=c("SITE_ID", "LONG_DD", "LAT_DD")))
rownames(latlong.coords)<- as.character(latlong.coords[,1])
colnames(latlong.coords)[2]<- 'long' #changing name to long so works with function
colnames(latlong.coords) [3]<- 'lat' #changing name to lat so works with function

#Generate a matrix containing distance between all pairs of coordinates in a matrix (in Kms)
km.matrix<- apply(latlong.coords[,2:3], 1, function(eachPoint) spDistsN1(as.matrix(latlong.coords[,2:3]), eachPoint, longlat=TRUE))
as.data.frame(km.matrix)
#Write to .csv so that can create master matrix of distances (for: nla2007_spatialdistances_April2015.csv)
#write.csv(km.matrix, "km.matrix.csv")

####Environmental distances between sites#### 
##PCA of environmental variables:
#NAs in waterqual.dat.red- replace with median values.
#The NAs are only in SECMEAN column. --> so replace those NAs with the median value from the SECMEAN column. 
waterqual.dat.red$SECMEAN[is.na(waterqual.dat.red$SECMEAN)] <- median(waterqual.dat.red$SECMEAN, na.rm=TRUE)

#PCA
rownames(waterqual.dat.red)<- as.character(waterqual.dat.red[,1])
waterqual.pca<- rda(decostand(waterqual.dat.red[,2:8], method="standardize"))
summary(waterqual.pca)
#Take a look at the biplot
plot(waterqual.pca, dis=c("sp", "sites"))
plot(waterqual.pca, type="n")
points(waterqual.pca, dis="sites")
text(waterqual.pca, dis="sp")

#Extract PC1 and PC2
waterqual.scores<- as.data.frame(scores(waterqual.pca, dis="sites", choices=c(1:2))) #Site scores 
#2 column dataframe, PC1 and PC2, where rownames show up as SITE_ID
waterqual.sp.scores<- as.data.frame(scores(waterqual.pca, dis="sp", choices=c(1:2))) #Species scores

#Bind SITE_ID to waterqual.scores
waterqual.scores<- as.data.frame(cbind(waterqual.dat.red$SITE_ID, waterqual.scores))
colnames(waterqual.scores)[1]<- 'SITE_ID'

#Site by site matrix for environmental distance for the waterqual.scores dataframe. 
waterqual.dist<- vegdist(waterqual.scores[,2:3], method="bray") #Base dissimilarities on PC1 and PC2 columns
#best choice of dissimilarity metric? Bray for now. NEED TO DECIDE ON THIS. ISSUE WITH NEGATIVE ENTRIES?

#write.csv(as.matrix(waterqual.dist), "waterqual.dist.csv")
#I used this to fill in nla2007_spatialdistances_April2015.csv 

dist.data<- read.csv(file.choose()) #nla2007_spatialdistances_April2015.csv 
#File path: C:\Users\Winegardner\Documents\MCGILL\PhD chapters and projects\EPA National lakes Assessment\Chapter 2\NLA data\Env and spatial distance

hist(dist.data$km_between) #Spatial distances between sites 
hist(dist.data$Env_distance_bray) #Environmental distances between sites 

##Make a clean PCA plot for manuscript Supplementary data 2
watqual.pca.plot<-ggplot()
watqual.pca.plot<- watqual.pca.plot + geom_vline(x=0,colour="grey50") 
watqual.pca.plot<- watqual.pca.plot+ geom_hline(y=0,colour="grey50") 
watqual.pca.plot<- watqual.pca.plot + labs(x= "PC1 (54% var explained)", y= "PC2 (18% var exp)") + theme_bw()
watqual.pca.plot<- watqual.pca.plot + geom_point(data = waterqual.scores, aes(x = PC1, y = PC2), size=4) 
watqual.pca.plot<- watqual.pca.plot + geom_text(data = waterqual.sp.scores, aes(x = PC1, y = PC2, label=rownames(waterqual.sp.scores)), size=6) 
watqual.pca.plot<- watqual.pca.plot + theme(axis.text.x = element_text(colour="black", size=16))
watqual.pca.plot<- watqual.pca.plot + theme(axis.text.y = element_text(colour="black", size=16))
watqual.pca.plot<- watqual.pca.plot + theme(axis.title.x = element_text(size = rel(2), angle=00))
watqual.pca.plot<- watqual.pca.plot + theme(axis.title.y = element_text(size = rel(2), angle=90))

####Presence-absence diatoms#### 
##Diatom data
nla.diat<- read.csv(file.choose()) #combo.red30highC.diatoms_matchedpa_FIXED.csv
#File path: McGill/PhD Chapters and projects/EPA National lakes assessment/Chapter 2/ NLA data/ Files for analysesNov2014Rscript

#In nla.data2.highC --> subset into surface and historical sediments 
nla.highC.surf<- subset(nla.data2.highC, SED_TYPE == "Surface") #Surface/top sediments #186
nla.highC.hist<- subset(nla.data2.highC, SED_TYPE == "Historical") #Historical/bottom sediments #186

#In nla.diat --> subset into surface and historical sediments 
nla.diat.surf<- subset(nla.diat, SED_TYPE == "Surface") #Surface/top sediments #186 
nla.diat.hist<- subset(nla.diat, SED_TYPE == "Historical") #Historical/bottom sediments #186

##Surface diatoms 
#In nla.diat.surf --> retain only sites (rows) found in nla.highC.surf
#Think that they might actually match already- check. 
#Sites in nla.diat.surf that are found in nla.highC.surf
surf.check<- (nla.diat.surf$SITE_ID %in% nla.highC.surf$SITE_ID) #All match

#Add MY_SAMPLE_ID column to the diatom file
#Add WSA_ECOREGION column to the diatom file
#Add LONG_DD
#Add LAT_DD
nla.diat.surf<- as.data.frame(cbind(nla.diat.surf, nla.highC.surf[,2], nla.highC.surf[,11], nla.highC.surf[,4], nla.highC.surf[,5]))
summary(nla.diat.surf) #last 4 columns are now "MY_SAMPLE ID", "WSA_ECOREGION", "LONG_DD", "LAT_DD" (but labelled as column numbers)
#Rename those last two columns as "MY_SAMPLE_ID", "WSA_ECOREGION", "LONG_DD", "LAT_DD"
colnames(nla.diat.surf)[1248]<- 'MY_SAMPLE_ID'
colnames(nla.diat.surf)[1249]<- 'WSA_ECOREGION'
colnames(nla.diat.surf)[1250]<- 'LONG_DD'
colnames(nla.diat.surf)[1251]<- 'LAT_DD'

#Ok- breakdown for new nla.diat.surf: 
#[,1]: SITE_ID
#[,2]: SED_TYPE
#[,3:1247]: diatom species columns
#[,1248]: MY_SAMPLE_ID (unique identifier that can match to abiotic file as well)
#[,1249]: WSA_ECOREGION 
#[,1250]: LONG_DD
#[,1251]: LAT_DD
#So will be able to subsample the diatom files based on ecoregion later on. 

##Historical diatoms 
#In nla.diat.hist --> retain only sites (rows) found in nla.highC.hist
#Think that they might actually match already- check. 
#Sites in nla.diat.surf that are found in nla.highC.surf
hist.check<- (nla.diat.hist$SITE_ID %in% nla.highC.hist$SITE_ID) #All match 

#Add MY_SAMPLE_ID column to the diatom file
#Add WSA_ECOREGION column to the diatom file
#Add LONG_DD
#Add LAT_DD
nla.diat.hist<- as.data.frame(cbind(nla.diat.hist, nla.highC.hist[,2], nla.highC.hist[,11], nla.highC.hist[,4], nla.highC.hist[,5]))
summary(nla.diat.hist) #last 4 columns are now "MY_SAMPLE ID", "WSA_ECOREGION", "LONG_DD", "LAT_DD" (but labelled as column numbers)
#Rename those last two columns as "MY_SAMPLE_ID", "WSA_ECOREGION", "LONG_DD", "LAT_DD"
colnames(nla.diat.hist)[1248]<- 'MY_SAMPLE_ID'
colnames(nla.diat.hist)[1249]<- 'WSA_ECOREGION'
colnames(nla.diat.hist)[1250]<- 'LONG_DD'
colnames(nla.diat.hist)[1251]<- 'LAT_DD'

#Ok- breakdown for new nla.diat.hist: 
#[,1]: SITE_ID
#[,2]: SED_TYPE
#[,3:1247]: diatom species columns
#[,1248]: MY_SAMPLE_ID (unique identifier that can match to abiotic file as well)
#[,1249]: WSA_ECOREGION 
#[,1250]: LONG_DD
#[,1251]: LAT_DD
#So will be able to subsample the diatom files based on ecoregion later on. 

####Quantitative diatoms#### 
##Surface diatoms 
surf.diat.abund<- read.csv(file.choose()) #topv1.red30highC.diatoms_abund
#File path: C:\Users\Winegardner\Documents\MCGILL\PhD chapters and projects\EPA National lakes Assessment\Chapter 2\NLA data\Files for analysesNov2014R script

##Historical diatoms 
hist.diat.abund<- read.csv(file.choose()) #botv1.red30highC.diatoms_abund
#File path: C:\Users\Winegardner\Documents\MCGILL\PhD chapters and projects\EPA National lakes Assessment\Chapter 2\NLA data\Files for analysesNov2014R script

####Radiometric dating data#### 
##Data as of April 18, 2015 (still waiting on a few samples)
dates.dat<- read.csv(file.choose()) #nla_dating_binary.csv
#File path: C:\Users\Winegardner\Documents\MCGILL\PhD chapters and projects\EPA National lakes Assessment\Chapter 2\Radiometric dating
#1 means sample is thought to be pre-1850, 0 means sample NOT pre-1850

####Final data available####
nla.highC.surf #Site info/land use data- surface sediments (high C sites only)
nla.highC.hist #Site info/land use data- historical sediments (high C sites only)
waterqual.dat.red #Water qulaity data (contemporary) for the 186 sites 
km.matrix #Spatial distances bwteen sites 
waterqual.dist #Water-quality (env) distances between sites 
dist.data #Dataframe with both spatial and env distances between sites
nla.diat.surf #Presence-absence surface diatoms
nla.diat.hist #Presence-absence historical diatoms 
surf.diat.abund #Quantitative surface diatoms
hist.diat.abund #Quantitative historical diatoms
dates.dat #Dating data with  binary assignments for pre-1850 based on radiometric dating 

##################################################################################################

##################################################################################################
####RADIOMETRIC DATING LOGISTIC REGRESSION AND CHI-SQAURE TESTS                                  #
##################################################################################################

####Logistic regression on binary results for 210Pb ratios####
#Use dates.dat 
#Looking for logistic relationship between core height and whether a core was pre-1850 as 
#evaluated by radiometric dating. 
dates.logr <- glm(Pre_1850_binary ~ Core_height, data=dates.dat, family=binomial(link="logit"))
dates.logr
summary(dates.logr)
logr.plot<- ggplot(dates.dat, aes(x=Core_height, y=Pre_1850_binary)) + geom_point() + stat_smooth(method="glm", family="binomial", se=FALSE, size=1)

##Chi-square test to test the hypothesis that age assignment is independent of core length
ctbl<- table(dates.dat$Core_height_description, dates.dat$Pre_1850) #contingency table #May be an issue with NAs? (fix dating.binary as.is=T)

chisq.test(ctbl) #p-value = 0.9586
#cannot reject null that core length is INDEPENDENT of pre-1850 or not
#so they could be dependent on each other. 

##Chi-sqaure of pre 1850 using radiometric dating and pre-1850 as predicted by Brothers eqn ages
ctbl2<- table(dates.dat$Pre_1850_binary_Brothers, dates.dat$Pre_1850_binary) #contingency table

chisq.test(ctbl2) #p-value = ~0.6
#cannot reject null that radiometric dates are INDEPENDENT of Brothers eqn dates.
#so could be related. 

##################################################################################################

##################################################################################################
####SPATIAL BETA-DIVERSITY (LEGENDRE FUNCTIONS)                                                  #
##################################################################################################

####Use beta.div to calculate total spatial BD on appropriately transformed data#### 
##NB: beta.div uses non-transformed quantitative data 

####Using quantitative data for both Historical and Surface sediments
#Sorensen index. 

##Prepping surface sediments - Quantitative, Sorensen
surf.diat.abund

#Non-transformed 
#Bind the non-transformed species columns to site, ecoregion and lat/long columns. 
surf.diat.abundnon<- as.data.frame(cbind(nla.diat.surf[,1], nla.diat.surf[,1249], nla.diat.surf[,1250], nla.diat.surf[,1251], surf.diat.abund[,4:953]))
colnames(surf.diat.abundnon)[1]<- 'SITE_ID'
colnames(surf.diat.abundnon)[2]<- 'WSA_ECOREGION'
colnames(surf.diat.abundnon)[3]<- 'LONG_DD'
colnames(surf.diat.abundnon)[4]<- 'LAT_DD'

#Hellinger-transformed 
surf.diat.abundh<- decostand(surf.diat.abund[,4:953], method="hell") #Hellinger on species columns. 
#Bind the Hellinger transformed species columns to site, ecoregion and lat/long columns. 
surf.diat.abundh<- as.data.frame(cbind(nla.diat.surf[,1], nla.diat.surf[,1249], nla.diat.surf[,1250], nla.diat.surf[,1251], surf.diat.abundh))
colnames(surf.diat.abundh)[1]<- 'SITE_ID'
colnames(surf.diat.abundh)[2]<- 'WSA_ECOREGION'
colnames(surf.diat.abundh)[3]<- 'LONG_DD'
colnames(surf.diat.abundh)[4]<- 'LAT_DD'

##Prepping historical sediments - Quantitative, Sorensen 
hist.diat.abund

#Non-transformed
#Bind the non-transformed species columns to site, ecoregion and lat/long columns. 
hist.diat.abundnon<- as.data.frame(cbind(nla.diat.hist[,1], nla.diat.hist[,1249], nla.diat.hist[,1250], nla.diat.hist[,1251], hist.diat.abund[,3:1033]))
colnames(hist.diat.abundnon)[1]<- 'SITE_ID'
colnames(hist.diat.abundnon)[2]<- 'WSA_ECOREGION'
colnames(hist.diat.abundnon)[3]<- 'LONG_DD'
colnames(hist.diat.abundnon)[4]<- 'LAT_DD'

#Hellinger-transformed 
hist.diat.abundh<- decostand(hist.diat.abund[,3:1033], method="hell") #Hellinger on species columns. 
#Bind the Hellinger transformed species columns to site, ecoregion and lat/long columns. 
hist.diat.abundh<- as.data.frame(cbind(nla.diat.hist[,1], nla.diat.hist[,1249], nla.diat.hist[,1250], nla.diat.hist[,1251], hist.diat.abundh))
colnames(hist.diat.abundh)[1]<- 'SITE_ID'
colnames(hist.diat.abundh)[2]<- 'WSA_ECOREGION'
colnames(hist.diat.abundh)[3]<- 'LONG_DD'
colnames(hist.diat.abundh)[4]<- 'LAT_DD'

##Surface - beta.div
surf.beta.div <- beta.div(surf.diat.abundnon[,5:954], method="hellinger", sqrt.D=FALSE, samp=TRUE, nperm=999, save.D=FALSE, clock=FALSE)
#Hellinger transformation for calculating total BD. 
summary(surf.beta.div)
surf.beta.div$SStotal_BDtotal
surf.beta.div$SCBD
surf.beta.div$LCBD
surf.beta.div$p.LCBD

#Without Hellinger transformation
surf.beta.div2 <- beta.div(surf.diat.abundnon[,5:954], sqrt.D=FALSE, samp=TRUE, nperm=999, save.D=FALSE, clock=FALSE)
summary(surf.beta.div2)
surf.beta.div2$SStotal_BDtotal
surf.beta.div2$SCBD
surf.beta.div2$LCBD
surf.beta.div2$p.LCBD

#Look at SCBD specifically
surf.SCBD.summary<- as.data.frame(surf.beta.div2$SCBD)
#write.csv(surf.SCBD.summary, "surf.SCBD.summary.csv") #so can order by magnitude and look at highest contributing species.

#Use percentage difference 
surf.beta.div3 <- beta.div(surf.diat.abundnon[,5:954], method="percentagedifference", sqrt.D=FALSE, samp=TRUE, nperm=999, save.D=FALSE, clock=FALSE)
summary(surf.beta.div3)
surf.beta.div3$SStotal_BDtotal
surf.beta.div3$SCBD
surf.beta.div3$LCBD
surf.beta.div3$p.LCBD

#Make a dataframe with LCBD information for surface data using percentage difference
surf.LCBD.values<- as.data.frame(surf.beta.div3$LCBD)
surf.LCBD.p<- as.data.frame(surf.beta.div3$p.LCBD)
surf.LCBD.summary<- as.data.frame(cbind(surf.diat.abundnon$SITE_ID, surf.diat.abundnon$WSA_ECOREGION, surf.diat.abundnon$LONG_DD, surf.diat.abundnon$LAT_DD, surf.LCBD.values, surf.LCBD.p))
colnames(surf.LCBD.summary)[1]<- 'SITE_ID'
colnames(surf.LCBD.summary)[2]<- 'WSA_ECOREGION'
colnames(surf.LCBD.summary)[3]<- 'LONG_DD' 
colnames(surf.LCBD.summary)[4]<- 'LAT_DD'
colnames(surf.LCBD.summary)[5]<- 'LCBD'
colnames(surf.LCBD.summary)[6]<- 'LCBD.P'

#Look at sites with significant (p<0.05)
surf.LCBD.summary.sig<- as.data.frame(subset(surf.LCBD.summary, LCBD.P<0.05))

##Historical - beta.div
hist.beta.div <- beta.div(hist.diat.abundnon[,5:1035], method="hellinger", sqrt.D=FALSE, samp=TRUE, nperm=999, save.D=FALSE, clock=FALSE)
summary(hist.beta.div)
hist.beta.div$SStotal_BDtotal
hist.beta.div$SCBD
hist.beta.div$LCBD
hist.beta.div$p.LCBD

#Without Hellinger transformation
hist.beta.div2 <- beta.div(hist.diat.abundnon[,5:1035], sqrt.D=FALSE, samp=TRUE, nperm=999, save.D=FALSE, clock=FALSE)
summary(hist.beta.div2)
hist.beta.div2$SStotal_BDtotal
hist.beta.div2$SCBD
hist.beta.div2$LCBD
hist.beta.div2$p.LCBD

#Look at SCBD specifically
hist.SCBD.summary<- as.data.frame(hist.beta.div2$SCBD)
#write.csv(hist.SCBD.summary, "hist.SCBD.summary.csv") #so can order by magnitude and look at highest contributing species.

#Use percentage difference instead of Hellinger
hist.beta.div3 <- beta.div(hist.diat.abundnon[,5:1035], method="percentagedifference", sqrt.D=FALSE, samp=TRUE, nperm=999, save.D=FALSE, clock=FALSE)
summary(hist.beta.div3)
hist.beta.div3$SStotal_BDtotal
hist.beta.div3$SCBD
hist.beta.div3$LCBD
hist.beta.div3$p.LCBD

#Make a dataframe with LCBD information for historical data using percentage difference
hist.LCBD.values<- as.data.frame(hist.beta.div3$LCBD)
hist.LCBD.p<- as.data.frame(hist.beta.div3$p.LCBD)
hist.LCBD.summary<- as.data.frame(cbind(hist.diat.abundnon$SITE_ID, hist.diat.abundnon$WSA_ECOREGION, hist.diat.abundnon$LONG_DD, hist.diat.abundnon$LAT_DD, hist.LCBD.values, hist.LCBD.p))
colnames(hist.LCBD.summary)[1]<- 'SITE_ID'
colnames(hist.LCBD.summary)[2]<- 'WSA_ECOREGION'
colnames(hist.LCBD.summary)[3]<- 'LONG_DD' 
colnames(hist.LCBD.summary)[4]<- 'LAT_DD'
colnames(hist.LCBD.summary)[5]<- 'LCBD'
colnames(hist.LCBD.summary)[6]<- 'LCBD.P'

#Look at sites with significant (p<0.05)
hist.LCBD.summary.sig<- as.data.frame(subset(hist.LCBD.summary, LCBD.P<0.05))

##################################################################################################

##################################################################################################
####FIGURES RELATED TO SPATIAL BETA-DIVERSITY                                                    #
##################################################################################################

####Use results from beta.div to create figures for manuscript####
##Map LCBD values across the landscape, coding by signifance (p-value)
#Base map
all.states <- map_data("state")
us.map<- ggplot()
us.map<- us.map + geom_polygon(data=all.states, aes(x=long, y=lat, group=group), colour="black", fill="grey50")

##FIGURE 2, map of sites by ecoregion along with lake water quality variables
#(A) map
us.map
us.map.ecoregion<- us.map + geom_point(data=nla.highC.surf, aes(x=LONG_DD, y=LAT_DD, colour=WSA_ECOREGION), size=4)
us.map.ecoregion<- us.map.ecoregion + theme_bw() + labs(x= "Longitude (DD)", y="Latitude (DD)")
us.map.ecoregion<- us.map.ecoregion + theme(axis.text.x = element_text(colour="black", size=16))
us.map.ecoregion<- us.map.ecoregion + theme(axis.text.y = element_text(colour="black", size=16))
us.map.ecoregion<- us.map.ecoregion + theme(axis.title.x = element_text(size = rel(2), angle=00))
us.map.ecoregion<- us.map.ecoregion + theme(axis.title.y = element_text(size = rel(2), angle=90))
us.map.ecoregion<- us.map.ecoregion + scale_colour_brewer(type="qual", palette="Dark2", "Ecoregion (WSA)", breaks=c("CPL", "NAP", "NPL", "SPL", "TPL", "UMW", "WMT", "XER"), labels=c("Coastal Plains", "N. Appalachians", "N. Plains", "S. Plains", "Temperate Plains", "Upper MidWest", "W. Mountains", "Xeric"))  
#us.map.ecoregion<- us.map.ecoregion + annotate("text", x=-120, y=50, label = "n = 35") #WMT
#us.map.ecoregion<- us.map.ecoregion + annotate("text", x=-118, y=42, label = "n = 1") #Xeric
#us.map.ecoregion<- us.map.ecoregion + annotate("text", x=-100, y=43, label = "n = 2") #S Plains
#us.map.ecoregion<- us.map.ecoregion + annotate("text", x=-93, y=43, label = "n = 6") #Temperate Plains
#us.map.ecoregion<- us.map.ecoregion + annotate("text", x=-94, y=48, label = "n = 69") #UMW
#us.map.ecoregion<- us.map.ecoregion + annotate("text", x=-75, y=46, label = "n = 57") #N App
#us.map.ecoregion<- us.map.ecoregion + annotate("text", x=-78, y=30, label = "n = 10") #Coastal plains
us.map.ecoregion<- us.map.ecoregion + theme(legend.title=element_text(size=16))
us.map.ecoregion<- us.map.ecoregion + theme(legend.text=element_text(size=16))
us.map.ecoregion<- us.map.ecoregion + annotate("text", x=-120, y=52, label="(a)", size=10)
#export 

#(B) Boxplots of lake area 
lakearea.bp<- ggplot(nla.highC.surf, aes(x=WSA_ECOREGION, y=LAKEAREA_KM2, fill=WSA_ECOREGION)) + geom_boxplot()
lakearea.bp<- lakearea.bp + theme_bw() + labs(x= "Ecoregion", y="Lake Area (Km2)")
lakearea.bp<- lakearea.bp + theme(axis.text.x = element_text(colour="black", size=16))
lakearea.bp<- lakearea.bp + theme(axis.text.y = element_text(colour="black", size=16))
lakearea.bp<- lakearea.bp + theme(axis.title.x = element_text(size = rel(2), angle=00))
lakearea.bp<- lakearea.bp + theme(axis.title.y = element_text(size = rel(2), angle=90))
lakearea.bp<- lakearea.bp + scale_fill_brewer(type="qual", palette="Dark2", "Ecoregion (WSA)", breaks=c("CPL", "NAP", "NPL", "SPL", "TPL", "UMW", "WMT", "XER"), labels=c("Coastal Plains", "N. Appalachians", "N. Plains", "S. Plains", "Temperate Plains", "Upper MidWest", "W. Mountains", "Xeric"))  
lakearea.bp<- lakearea.bp + theme(legend.position="none")
lakearea.bp<- lakearea.bp + annotate("text", x=1, y=120, label="(b)", size=10)

#(C) Boxplots of Z_max
zmax.bp<- ggplot(nla.highC.surf, aes(x=WSA_ECOREGION, y=Z_MAX, fill=WSA_ECOREGION)) + geom_boxplot()
zmax.bp<- zmax.bp + theme_bw() + labs(x= "Ecoregion", y="Observed maximum depth (m)")
zmax.bp<- zmax.bp + theme(axis.text.x = element_text(colour="black", size=16))
zmax.bp<- zmax.bp + theme(axis.text.y = element_text(colour="black", size=16))
zmax.bp<- zmax.bp + theme(axis.title.x = element_text(size = rel(2), angle=00))
zmax.bp<- zmax.bp + theme(axis.title.y = element_text(size = rel(2), angle=90))
zmax.bp<- zmax.bp + scale_fill_brewer(type="qual", palette="Dark2", "Ecoregion (WSA)", breaks=c("CPL", "NAP", "NPL", "SPL", "TPL", "UMW", "WMT", "XER"), labels=c("Coastal Plains", "N. Appalachians", "N. Plains", "S. Plains", "Temperate Plains", "Upper MidWest", "W. Mountains", "Xeric"))  
zmax.bp<- zmax.bp + theme(legend.position="none")
zmax.bp<- zmax.bp + annotate("text", x=1, y=60, label="(c)", size=10)

#(D) Boxplots of pH
ph.bp<- ggplot(waterqual.dat.red, aes(x=WSA_ECOREGION, y=PH_FIELD, fill=WSA_ECOREGION)) + geom_boxplot()
ph.bp<- ph.bp + theme_bw() + labs(x= "Ecoregion", y="pH")
ph.bp<- ph.bp + theme(axis.text.x = element_text(colour="black", size=16))
ph.bp<- ph.bp + theme(axis.text.y = element_text(colour="black", size=16))
ph.bp<- ph.bp + theme(axis.title.x = element_text(size = rel(2), angle=00))
ph.bp<- ph.bp + theme(axis.title.y = element_text(size = rel(2), angle=90))
ph.bp<- ph.bp + scale_fill_brewer(type="qual", palette="Dark2", "Ecoregion (WSA)", breaks=c("CPL", "NAP", "NPL", "SPL", "TPL", "UMW", "WMT", "XER"), labels=c("Coastal Plains", "N. Appalachians", "N. Plains", "S. Plains", "Temperate Plains", "Upper MidWest", "W. Mountains", "Xeric"))  
ph.bp<- ph.bp + theme(legend.position="none")
ph.bp<- ph.bp + annotate("text", x=1, y=10, label="(d)", size=10)

#(E) Boxplots of TP
tp.bp<- ggplot(waterqual.dat.red, aes(x=WSA_ECOREGION, y=PTL, fill=WSA_ECOREGION)) + geom_boxplot()
tp.bp<- tp.bp + theme_bw() + labs(x= "Ecoregion", y="Total phosphorus (ug/L)")
tp.bp<- tp.bp + theme(axis.text.x = element_text(colour="black", size=16))
tp.bp<- tp.bp + theme(axis.text.y = element_text(colour="black", size=16))
tp.bp<- tp.bp + theme(axis.title.x = element_text(size = rel(2), angle=00))
tp.bp<- tp.bp + theme(axis.title.y = element_text(size = rel(2), angle=90))
tp.bp<- tp.bp + scale_fill_brewer(type="qual", palette="Dark2", "Ecoregion (WSA)", breaks=c("CPL", "NAP", "NPL", "SPL", "TPL", "UMW", "WMT", "XER"), labels=c("Coastal Plains", "N. Appalachians", "N. Plains", "S. Plains", "Temperate Plains", "Upper MidWest", "W. Mountains", "Xeric"))  
tp.bp<- tp.bp + theme(legend.position="none")
tp.bp<- tp.bp + annotate("text", x=1, y=1500, label="(e)", size=10)

#Combo of B-E
be.panel<- grid.arrange(lakearea.bp, zmax.bp, ph.bp, tp.bp, nrow=2) #export 

#Combine map
fig2<- grid.arrange(us.map.ecoregion, lakearea.bp, zmax.bp, ph.bp, tp.bp, nrow=3)

##FIGURE 3, maps of LCBD and temporal exceptional sites 
#(A) Historical LCBD
#Map LCBD values across the landscape, coding by signifance (p-value)
us.map.histLCBD<- us.map + geom_point(data=hist.LCBD.summary, aes(x=LONG_DD, y=LAT_DD, size=LCBD, colour=LCBD.P<0.05))
us.map.histLCBD<- us.map.histLCBD + scale_size(name="LCBD") 
#us.map.histLCBD<- us.map.histLCBD + ggtitle("Historical LCBD")
#Map of LCBD values across landscape highlighting sites that make a SIGNIFICANT (p<0.05) contribution to BD. 
us.map.histLCBD<- us.map.histLCBD + labs(x="Longitude", y="Latitude") 
us.map.histLCBD<- us.map.histLCBD + theme_bw()
us.map.histLCBD<- us.map.histLCBD + theme(axis.text.x = element_text(colour="black", size=16))
us.map.histLCBD<- us.map.histLCBD + theme(axis.text.y = element_text(colour="black", size=16))
us.map.histLCBD<- us.map.histLCBD + theme(axis.title.x = element_text(size = rel(2), angle=00))
us.map.histLCBD<- us.map.histLCBD + theme(axis.title.y = element_text(size = rel(2), angle=90))
us.map.histLCBD<- us.map.histLCBD + theme(legend.title=element_text(size=16))
us.map.histLCBD<- us.map.histLCBD + theme(legend.text=element_text(size=16))
us.map.histLCBD<- us.map.histLCBD + annotate("text", x=-118, y=51, label="(a)", size=10) #for using in a panel with other LCBD map

#(B) Surface (2007) LCBD
us.map.surfLCBD<- us.map + geom_point(data=surf.LCBD.summary, aes(x=LONG_DD, y=LAT_DD, size=LCBD, colour=LCBD.P<0.05))
us.map.surfLCBD<- us.map.surfLCBD + scale_size(name="LCBD") 
#us.map.surfLCBD<- us.map.surfLCBD + ggtitle("2007 LCBD")
#Map of LCBD values across landscape highlighting sites that make a SIGNIFICANT (p<0.05) contribution to BD. 
us.map.surfLCBD<- us.map.surfLCBD + labs(x="Longitude", y="Latitude") 
us.map.surfLCBD<- us.map.surfLCBD + theme_bw()
us.map.surfLCBD<- us.map.surfLCBD + theme(axis.text.x = element_text(colour="black", size=16))
us.map.surfLCBD<- us.map.surfLCBD + theme(axis.text.y = element_text(colour="black", size=16))
us.map.surfLCBD<- us.map.surfLCBD + theme(axis.title.x = element_text(size = rel(2), angle=00))
us.map.surfLCBD<- us.map.surfLCBD + theme(axis.title.y = element_text(size = rel(2), angle=90))
us.map.surfLCBD<- us.map.surfLCBD + theme(legend.title=element_text(size=16))
us.map.surfLCBD<- us.map.surfLCBD + theme(legend.text=element_text(size=16))
us.map.surfLCBD<- us.map.surfLCBD + annotate("text", x=-118, y=51, label="(b)", size=9) #for using in a panel with other LCBD map

#(C) Temporal exceptional sites #temporalBD.BCD has to be loaded first (see below). 
us.map.tempsig<- us.map + geom_point(data=temporalBD.BCD, aes(x=LONG_DD, y=LAT_DD, size=Total_beta, colour=p.adj<0.05))
us.map.tempsig<- us.map.tempsig + scale_size(name="Temporal beta") 
us.map.tempsig<- us.map.tempsig + labs(x="Longitude", y="Latitude") 
us.map.tempsig<- us.map.tempsig + theme_bw()
us.map.tempsig<- us.map.tempsig + theme(axis.text.x = element_text(colour="black", size=16))
us.map.tempsig<- us.map.tempsig + theme(axis.text.y = element_text(colour="black", size=16))
us.map.tempsig<- us.map.tempsig + theme(axis.title.x = element_text(size = rel(2), angle=00))
us.map.tempsig<- us.map.tempsig + theme(axis.title.y = element_text(size = rel(2), angle=90))
us.map.tempsig<- us.map.tempsig + theme(legend.title=element_text(size=16))
us.map.tempsig<- us.map.tempsig + theme(legend.text=element_text(size=16))
us.map.tempsig<- us.map.tempsig + annotate("text", x=-118, y=51, label="(c)", size=9) #for using in a panel with other LCBD maps

##Combo panel for figure 3
LCBD.combo<- grid.arrange(us.map.histLCBD, us.map.surfLCBD, us.map.tempsig, nrow=3)

##################################################################################################

##################################################################################################
####SPATIAL BETA-DIVERSITY PARTITIONING (LEGENDRE FUNCTIONS)                                     #
##################################################################################################

####Use non-transformed data to partition spatial beta-diversity into replacement, AbDiff components#### 
##Surface diatoms
surf.abund.bdc<- beta.div.comp(surf.diat.abundnon[,5:954], coef="S", quant=TRUE, save.abc=FALSE)
summary(surf.abund.bdc$repl) #Replacement
summary(surf.abund.bdc$rich) #Abundance difference 
surf.abund.bdc$part
#1st value = Total beta div  
#2nd value = Total replacement diversity
#3rd value = Total abundance difference diversity 
#4th value = Total replcement/Total beta diversity
#5th value = Total abundance difference/Total beta diversity 

#write surf.abund.bdc to a .csv matrix so that can put pairwise comparisons into nla2007_spatialdistances_April2015.csv
surf.D<- as.matrix(surf.abund.bdc$D)
#write.csv(surf.D, "surface_spatialbeta.csv")

##For ecoregions (had to bind new columns)
#CPL
surf.diat.abundnonCPL<- as.data.frame(subset(surf.diat.abundnon, WSA_ECOREGION == "CPL", drop=T))
surf.abund.bdcCPL<- beta.div.comp(surf.diat.abundnonCPL[,5:954], coef="S", quant=TRUE, save.abc=FALSE)
surf.abund.bdcCPL$part

#NAP
surf.diat.abundnonNAP<- as.data.frame(subset(surf.diat.abundnon, WSA_ECOREGION == "NAP", drop=T))
surf.abund.bdcNAP<- beta.div.comp(surf.diat.abundnonNAP[,5:954], coef="S", quant=TRUE, save.abc=FALSE)
surf.abund.bdcNAP$part

#SPL
surf.diat.abundnonSPL<- as.data.frame(subset(surf.diat.abundnon, WSA_ECOREGION == "SPL", drop=T))
surf.abund.bdcSPL<- beta.div.comp(surf.diat.abundnonSPL[,5:954], coef="S", quant=TRUE, save.abc=FALSE)
surf.abund.bdcSPL$part

#TPL
surf.diat.abundnonTPL<- as.data.frame(subset(surf.diat.abundnon, WSA_ECOREGION == "TPL", drop=T))
surf.abund.bdcTPL<- beta.div.comp(surf.diat.abundnonTPL[,5:954], coef="S", quant=TRUE, save.abc=FALSE)
surf.abund.bdcTPL$part

#UMW
surf.diat.abundnonUMW<- as.data.frame(subset(surf.diat.abundnon, WSA_ECOREGION == "UMW", drop=T))
surf.abund.bdcUMW<- beta.div.comp(surf.diat.abundnonUMW[,5:954], coef="S", quant=TRUE, save.abc=FALSE)
surf.abund.bdcUMW$part

#WMT
surf.diat.abundnonWMT<- as.data.frame(subset(surf.diat.abundnon, WSA_ECOREGION == "WMT", drop=T))
surf.abund.bdcWMT<- beta.div.comp(surf.diat.abundnonWMT[,5:954], coef="S", quant=TRUE, save.abc=FALSE)
surf.abund.bdcWMT$part

##Historical diatoms 
hist.abund.bdc<- beta.div.comp(hist.diat.abundnon[,5:1035], coef="S", quant=TRUE, save.abc=FALSE)
summary(hist.abund.bdc$repl) #Replacement
summary(hist.abund.bdc$rich) #Abundance difference 
hist.abund.bdc$part

#write hist.abund.bdc to a .csv matrix so that can put pairwise comparisons into nla2007_spatialdistances_April2015.csv
hist.D<- as.matrix(hist.abund.bdc$D)
#write.csv(hist.D, "historical_spatialbeta.csv")

##For ecoregions (had to bind new columns)
#CPL
hist.diat.abundnonCPL<- as.data.frame(subset(hist.diat.abundnon, WSA_ECOREGION == "CPL", drop=T))
hist.abund.bdcCPL<- beta.div.comp(hist.diat.abundnonCPL[,5:1035], coef="S", quant=TRUE, save.abc=FALSE)
hist.abund.bdcCPL$part

#NAP
hist.diat.abundnonNAP<- as.data.frame(subset(hist.diat.abundnon, WSA_ECOREGION == "NAP", drop=T))
hist.abund.bdcNAP<- beta.div.comp(hist.diat.abundnonNAP[,5:1035], coef="S", quant=TRUE, save.abc=FALSE)
hist.abund.bdcNAP$part

#SPL
hist.diat.abundnonSPL<- as.data.frame(subset(hist.diat.abundnon, WSA_ECOREGION == "SPL", drop=T))
hist.abund.bdcSPL<- beta.div.comp(hist.diat.abundnonSPL[,5:1035], coef="S", quant=TRUE, save.abc=FALSE)
hist.abund.bdcSPL$part

#TPL
hist.diat.abundnonTPL<- as.data.frame(subset(hist.diat.abundnon, WSA_ECOREGION == "TPL", drop=T))
hist.abund.bdcTPL<- beta.div.comp(hist.diat.abundnonTPL[,5:1035], coef="S", quant=TRUE, save.abc=FALSE)
hist.abund.bdcTPL$part

#UMW
hist.diat.abundnonUMW<- as.data.frame(subset(hist.diat.abundnon, WSA_ECOREGION == "UMW", drop=T))
hist.abund.bdcUMW<- beta.div.comp(hist.diat.abundnonUMW[,5:1035], coef="S", quant=TRUE, save.abc=FALSE)
hist.abund.bdcUMW$part

#WMT
hist.diat.abundnonWMT<- as.data.frame(subset(hist.diat.abundnon, WSA_ECOREGION == "WMT", drop=T))
hist.abund.bdcWMT<- beta.div.comp(hist.diat.abundnonWMT[,5:1035], coef="S", quant=TRUE, save.abc=FALSE)
hist.abund.bdcWMT$part


##Scatterplots of the explanatory components

all.bdc<- read.csv(file.choose()) #total_bd_partitioning.csv ##Should move up to DATA SET-UP!!
#File path: C:\Users\Winegardner\Documents\MCGILL\PhD chapters and projects\EPA National lakes Assessment\Chapter 2
#Note: Made this .csv with results from partitioning, but "Total BD" above is based on total BD
#from beta.div.comp() and differs from beta.div() for some reason

##################################################################################################

##################################################################################################
####SPATIAL BETA-DIVERSITY LCBD PARTITIONING                                                     #
##################################################################################################

####Use LCDB.comp to partition spatial LCBD (total LCBD should already have been calculated####
#using beta.div)
##Surface diatoms - Repl
#surf.abund.bdc$repl
surf.repl.LCBD<-LCBD.comp(surf.abund.bdc$repl, sqrt.x=TRUE)
surf.repl.LCBD
surf.repl.LCBD$SStotal_BDtotal #Total SS and total BD 
surf.repl.LCBD$LCBD #Vector of local contributions to beta diversity for sites
surf.repl.LCBD$D

#Map LCBD vector to site names. 
surf.repl.LCBDmat<- as.data.frame(surf.repl.LCBD$LCBD)
colnames(surf.repl.LCBDmat)[1]<- 'Surface_Repl_LCBD'

##Surface diatoms - AbDiff
#surf.abund.bdc$rich
surf.abdiff.LCBD<-LCBD.comp(surf.abund.bdc$rich, sqrt.x=TRUE)
surf.abdiff.LCBD
surf.abdiff.LCBD$SStotal_BDtotal #Total SS and total BD 
surf.abdiff.LCBD$LCBD #Vector of local contributions to beta diversity for sites
surf.abdiff.LCBD$D

#Map LCBD vector to site names. 
surf.abdiff.LCBDmat<- as.data.frame(surf.abdiff.LCBD$LCBD)
colnames(surf.abdiff.LCBDmat)[1]<- 'Surface_AbDiff_LCBD'

##Historical diatoms - Repl
#hist.abund.bdc$repl
hist.repl.LCBD<-LCBD.comp(hist.abund.bdc$repl, sqrt.x=TRUE)
hist.repl.LCBD
hist.repl.LCBD$SStotal_BDtotal #Total SS and total BD 
hist.repl.LCBD$LCBD #Vector of local contributions to beta diversity for sites
hist.repl.LCBD$D

#Map LCBD vector to site names. 
hist.repl.LCBDmat<- as.data.frame(hist.repl.LCBD$LCBD)
colnames(hist.repl.LCBDmat)[1]<- 'Historical_Repl_LCBD'

##Historical diatoms - AbDiff
#hist.abund.bdc$rich 
hist.abdiff.LCBD<-LCBD.comp(hist.abund.bdc$rich, sqrt.x=TRUE)
hist.abdiff.LCBD
hist.abdiff.LCBD$SStotal_BDtotal #Total SS and total BD 
hist.abdiff.LCBD$LCBD #Vector of local contributions to beta diversity for sites
hist.abdiff.LCBD$D

#Map LCBD vector to site names. 
hist.abdiff.LCBDmat<- as.data.frame(hist.abdiff.LCBD$LCBD)
colnames(hist.abdiff.LCBDmat)[1]<- 'Historical_AbDiff_LCBD'

##Combine into one dataframe with SITE_ID, WSA_ECOREGION and LAT/LONG
LCBDpart.combos<- as.data.frame(cbind(surf.diat.abundh$SITE_ID, surf.diat.abundh$WSA_ECOREGION, surf.diat.abundh$LONG_DD, surf.diat.abundh$LAT_DD, nla.highC.surf$NLA_WGT, surf.repl.LCBDmat, surf.abdiff.LCBDmat, hist.repl.LCBDmat, hist.abdiff.LCBDmat))
colnames(LCBDpart.combos)[1]<- 'SITE_ID'
colnames(LCBDpart.combos)[2]<- 'WSA_ECOREGION'
colnames(LCBDpart.combos)[3]<- 'LONG_DD'
colnames(LCBDpart.combos)[4]<- 'LAT_DD'
colnames(LCBDpart.combos)[5]<- 'NLA_WGT'

##################################################################################################

##################################################################################################
####TEMPORAL BETA-DIVERSITY (DECOMPOSE.DR2 FUNCTION)                                             #
##################################################################################################

####Using new decompose.D2 function which allows for directional (time point 1 to time point 2 comparisons)####
#Create quantitative data matrices with same species
nla.diat.new<- as.data.frame(cbind(nla.diat, nla.data2.highC$MY_SAMPLE_ID, nla.data2.highC$WSA_ECOREGION, nla.data2.highC$LONG_DD, nla.data2.highC$LAT_DD))
summary(nla.diat.new)
colnames(nla.diat.new)[1248]<- 'MY_SAMPLE_ID'
colnames(nla.diat.new)[1249]<- 'WSA_ECOREGION'
colnames(nla.diat.new)[1250]<- 'LONG_DD'
colnames(nla.diat.new)[1251]<- 'LAT_DD'

#Now subset that matrix again by SED_TYPE so that ensure row names same betwen surface and historical. 
nla.diat.new.surf<- subset(nla.diat.new, SED_TYPE == "Surface", drop=T)
nla.diat.new.hist<- subset(nla.diat.new, SED_TYPE == "Historical", drop=T)

#Check that site names match between the subsets
#Sites in nla.diat.new.hist that are found in nla.diat.new.surf
site.check<- (nla.diat.new.hist$SITE_ID %in% nla.diat.new.surf$SITE_ID) #All match

#Set row names as SITE_ID
rownames(nla.diat.new.surf)<- as.character(nla.diat.new.surf[,1])
rownames(nla.diat.new.hist)<- as.character(nla.diat.new.hist[,1])

##Temporal BD between historical and surface using decompose.D2 (PRESENCE-ABSENCE DATA)
temp.D2.pa<- decompose.D2(nla.diat.new.hist[,3:1247], nla.diat.new.surf[,3:1247], den.type=1)
#using den.type=1, where denomination = (A+B+C) because presence-absence data

#Looking at output
temp.D2.pa
temp.D2pa.mat1<- as.data.frame(temp.D2.pa$mat1)
temp.D2pa.mat2<- as.data.frame(temp.D2.pa$mat2)
summary(temp.D2pa.mat2)

##Temporal BD between historical and surface using decompose.D2 (QUANTITATIVE DATA)
#Need non-transformed quantitative data with same species. 

#Matrices to stack: 
#surf.diat.abundnon #surface diatoms, quantitative, non-transformed, species col 5:954
#hist.diat.abundnon #historical diatoms, quantitative, non-transformed, species col 5:1035

#Add a column with Surface or Historical to each matrix
#Surface
surf.diat.abundnon.label<- surf.diat.abundnon #create new dataframe that can use for manipulation
surf.diat.abundnon.label$TYPE<- rep("Surface", nrow(surf.diat.abundnon.label))

#Historical
hist.diat.abundnon.label<- hist.diat.abundnon
hist.diat.abundnon.label$TYPE<- rep("Historical", nrow(hist.diat.abundnon.label))

#Melt each so in long format: 
surf.long<- melt(surf.diat.abundnon.label, id.vars=c("SITE_ID", "WSA_ECOREGION", "LONG_DD", "LAT_DD", "TYPE"))
hist.long<- melt(hist.diat.abundnon.label, id.vars=c("SITE_ID", "WSA_ECOREGION", "LONG_DD", "LAT_DD", "TYPE"))

#Stack the long format dataframes on top of eachother: 
combo.long<- rbind(surf.long, hist.long)

#Cast into wideformat
combo.wide<- dcast(combo.long, SITE_ID + WSA_ECOREGION + LONG_DD + LAT_DD + TYPE ~ variable)
#Now for each site there is a surface and historical abundance measurement, common set of species 
#NAs for species not found in one of the matrices. 
#Replace NAs with zeros. 
combo.wide[is.na(combo.wide)]<-0

#Resubset into separate matrices
surf.abund.fortempBD<- as.data.frame(subset(combo.wide, TYPE == "Surface", drop=T)) #6:1251 are columns
#Make rownames the SITE_ID
rownames(surf.abund.fortempBD)<- as.character(surf.abund.fortempBD[,1])

hist.abund.fortempBD<- as.data.frame(subset(combo.wide, TYPE == "Historical", drop=T)) #6:1251 are columns
rownames(hist.abund.fortempBD)<- as.character(hist.abund.fortempBD[,1])

#decompose.D2
temp.D2.quant<- decompose.D2(hist.abund.fortempBD[,6:1251], surf.abund.fortempBD[,6:1251], den.type=2)

#Looking at output
temp.D2.quant
temp.D2quant.mat1<- as.data.frame(temp.D2.quant$mat1)
temp.D2quant.mat2<- as.data.frame(temp.D2.quant$mat2)

#Bind this quantitative output back with ecoregion data
temp.D2quant.mat3<- as.data.frame(cbind(surf.abund.fortempBD[,1:4], temp.D2quant.mat2)) #Site data etc with the output matrix 2, don't include "Surface" label                               
colnames(temp.D2quant.mat3)[5]<- 'A' #Similarity 
colnames(temp.D2quant.mat3)[6]<- 'B' #Species loss from T1
colnames(temp.D2quant.mat3)[7]<- 'C' #Species gain in T2
colnames(temp.D2quant.mat3)[8]<- 'D' #Dissimilarity 

summary(temp.D2quant.mat3)

temp.D2quant.summary<- as.data.frame(ddply(temp.D2quant.mat3, "WSA_ECOREGION", summarize, 
                                           mean_A = mean(A, na.rm=T), 
                                           sd_A = sd(A, na.rm=T),
                                           mean_B = mean(B, na.rm=T),
                                           sd_B = sd(B, na.rm=T),
                                           mean_C = mean(C, na.rm=T),
                                           sd_C = sd(C, na.rm=T),
                                           mean_D = mean(D, na.rm=T), 
                                           sd_D = sd(D, na.rm=T)))

##################################################################################################

##################################################################################################
####FINDING SIGNIFICANT SITES IN TEMPORAL BETA-DIVERSITY (paired.diff2 FUNCTION)                 #
##################################################################################################

####Using paired.diff2_2.R to look for exceptional sites in temporal beta diversity####
##Using percentage difference on abundance data, testing the different permutation options.
##H0 is that there are no real difference between T1 and T2 for a site. (sig values indicate exceptional sites)
##This revised function produces matrix of B and C values as well. 

##Abundance data- using permutation method 1. 
#Permutations n=9999, % difference 
except.sites2.BCD<- paired.diff2(hist.abund.fortempBD[,6:1251], surf.abund.fortempBD[,6:1251], method="%difference", pa.tr=FALSE, nperm=9999, permute.sp=1, BCD=TRUE, replace=FALSE, clock=FALSE)

except.sites2.BCD
except.sites2.BCD$vecD.ref
except.sites2.BCD$p.dist
except.sites2.BCD$p.adj
except.sites2.BCD$BCD.mat #need to subset this matrix to highlight the "objects" with significant p-values. 

summary(except.sites2.BCD$p.dist < 0.05) #44 sig
summary(except.sites2.BCD$p.adj < 0.05) #adjusted for multiple testing, 14 sig, 40 when ran again? #CHECK THIS

##Put the BCD components of except sites into a matrix/dataframe with other descriptors.
#Basic matrix 
temporalBD.BCD<- as.data.frame(cbind(except.sites2.BCD$BCD.mat, except.sites2.BCD$p.dist, except.sites2.BCD$p.adj))
colnames(temporalBD.BCD)[1]<- 'Species_loss_component' #B/(2A+B+C)
colnames(temporalBD.BCD)[2]<- 'Species_gain_component' #C/(2A+B+C)
colnames(temporalBD.BCD)[3]<- 'Total_beta' #D=(B+C)/(2A+B+C)
colnames(temporalBD.BCD)[4]<- 'p.dist' #non-adjusted p-value, p<0.05, sig D between time 1 and 2
colnames(temporalBD.BCD)[5]<- 'p.adj' #p-value adjusted for multiple testing (holm method)

#Add columns for binary designations of the signifigance of D. 
temporalBD.BCD$p.dist.binary[temporalBD.BCD$p.dist < 0.05]<- 1 #add a column for binary designation of p.dist, 1= p<0.05 (sig), 0=p>0.05 (not sig.) 
temporalBD.BCD$p.dist.binary[temporalBD.BCD$p.dist > 0.05]<- 0

temporalBD.BCD$p.adj.binary[temporalBD.BCD$p.adj < 0.05]<- 1 #add a column for binary designation of p.adj, 1= p<0.05 (sig), 0=p>0.05 (not sig.)  
temporalBD.BCD$p.adj.binary[temporalBD.BCD$p.adj > 0.05]<- 0

#Add columns for site ID, ecoregion, lake origin etc, and land use variables from 2007 representing percent development in a basin and percent agriculture in a basin
temporalBD.BCD<- as.data.frame(cbind(surf.diat.abundnon.label$SITE_ID, surf.diat.abundnon.label$WSA_ECOREGION, surf.diat.abundnon.label$LONG_DD, surf.diat.abundnon.label$LAT_DD, nla.highC.surf$ORIGIN, nla.highC.surf$LAKE_TYPE, nla.highC.surf$PCT_DEVELOPED_BSN, nla.highC.surf$PCT_AGRIC_BSN, temporalBD.BCD))
colnames(temporalBD.BCD)[1]<- 'SITE_ID'
colnames(temporalBD.BCD)[2]<- 'WSA_ECOREGION'
colnames(temporalBD.BCD)[3]<- 'LONG_DD'
colnames(temporalBD.BCD)[4]<- 'LAT_DD'
colnames(temporalBD.BCD)[5]<- 'ORIGIN' 
colnames(temporalBD.BCD)[6]<- 'LAKE_TYPE' 
colnames(temporalBD.BCD)[7]<- 'PCT_DEVELOPED_BSN'
colnames(temporalBD.BCD)[8]<- 'PCT_AGRIC_BSN'



##Code for when trying to do double access, don't use###
#write.csv(temporalBD.BCD, "temporalBD.BCD.csv") #modified and re-read in. 

fig5.dat<- read.csv(file.choose()) #temporalBD.BCD_forfig4.csv
#File path: C:\Users\Winegardner\Documents\MCGILL\PhD chapters and projects\EPA National lakes Assessment\Chapter 2\Manuscript figures
#Stacked matrix, where species loss and gainn component combined into one column called "Species_component"
#where species gain are positive #s, species loss are negative #s. 

##################################################################################################

##################################################################################################
####FIGURES RELATED TO TEMPORAL BETA-DIVERSITY COMPONENTS                                        #
##################################################################################################

####Use temporal beta-diversity figures for manuscript, figure 5#### 
##FIGURE 5 - Temporal beta-diversity and components, logistic regression of B and C components as predictors, and sig. D value as outcome.

#(5A) Species loss component 
#Species loss component (B)-Adjusted P-value
temporalB2.logr <- glm(p.adj.binary ~ Species_loss_component, data=temporalBD.BCD, family=binomial(link="logit"))
summary(temporalB2.logr)
temporalB2logr.plot<- ggplot(temporalBD.BCD, aes(x=Species_loss_component, y=p.adj.binary)) + geom_point() + stat_smooth(method="glm", family="binomial", se=FALSE, size=1)
temporalB2logr.plot<- temporalB2logr.plot + labs(x="Species loss component", y="Significance of BD") 
temporalB2logr.plot<- temporalB2logr.plot + theme_bw()
temporalB2logr.plot<- temporalB2logr.plot + theme(axis.text.x = element_text(colour="black", size=16))
temporalB2logr.plot<- temporalB2logr.plot + theme(axis.text.y = element_text(colour="black", size=16))
temporalB2logr.plot<- temporalB2logr.plot + theme(axis.title.x = element_text(size = rel(2), angle=00))
temporalB2logr.plot<- temporalB2logr.plot + theme(axis.title.y = element_text(size = rel(2), angle=90))
temporalB2logr.plot<- temporalB2logr.plot + annotate("text", x=0.1, y=1.05, label="(a)", size=10)

#(5B)
#Species gain component (C)- Adjusted P-value 
temporalC2.logr <- glm(p.adj.binary ~ Species_gain_component, data=temporalBD.BCD, family=binomial(link="logit"))
summary(temporalC2.logr)
temporalC2logr.plot<- ggplot(temporalBD.BCD, aes(x=Species_gain_component, y=p.adj.binary)) + geom_point() + stat_smooth(method="glm", family="binomial", se=FALSE, size=1)
temporalC2logr.plot<- temporalC2logr.plot + labs(x="Species gain component", y="Significance of BD") 
temporalC2logr.plot<- temporalC2logr.plot + theme_bw()
temporalC2logr.plot<- temporalC2logr.plot + theme(axis.text.x = element_text(colour="black", size=16))
temporalC2logr.plot<- temporalC2logr.plot + theme(axis.text.y = element_text(colour="black", size=16))
temporalC2logr.plot<- temporalC2logr.plot + theme(axis.title.x = element_text(size = rel(2), angle=00))
temporalC2logr.plot<- temporalC2logr.plot + theme(axis.title.y = element_text(size = rel(2), angle=90))
temporalC2logr.plot<- temporalC2logr.plot + annotate("text", x=0.11, y=1.05, label="(b)", size=10)

#Combo
temporallogr.combo<- grid.arrange(temporalB2logr.plot, temporalC2logr.plot, nrow=1)

#(5C)
temporalBD.BCD.gainmean<- as.data.frame(ddply(temporalBD.BCD, "WSA_ECOREGION", summarize, 
                                              mean_Totalbeta = mean(Total_beta, na.rm=T),
                                              mean_Speciesgain = mean(Species_gain_component, na.rm=T)))

fig5c.plot<- ggplot(temporalBD.BCD.gainmean, aes(x=mean_Totalbeta, y=mean_Speciesgain, colour=WSA_ECOREGION)) + geom_point(size=4)
fig5c.plot<- fig5c.plot + labs(x="Mean temporal BD", y="Mean species gain") 
fig5c.plot<- fig5c.plot + theme_bw()
fig5c.plot<- fig5c.plot + theme(axis.text.x = element_text(colour="black", size=16))
fig5c.plot<- fig5c.plot + theme(axis.text.y = element_text(colour="black", size=16))
fig5c.plot<- fig5c.plot + theme(axis.title.x = element_text(size = rel(2), angle=00))
fig5c.plot<- fig5c.plot + theme(axis.title.y = element_text(size = rel(2), angle=90))
fig5c.plot<- fig5c.plot + annotate("text", x=0.27, y=0.35, label="(c)", size=10)
fig5c.plot<- fig5c.plot + theme(legend.title=element_text(size=16))
fig5c.plot<- fig5c.plot + theme(legend.text=element_text(size=16))
fig5c.plot<- fig5c.plot + scale_colour_brewer(type="qual", palette="Dark2", "Ecoregion (WSA)", breaks=c("CPL", "NAP", "NPL", "SPL", "TPL", "UMW", "WMT", "XER"), labels=c("Coastal Plains", "N. Appalachians", "N. Plains", "S. Plains", "Temperate Plains", "Upper MidWest", "W. Mountains", "Xeric"))  

##################################################################################################

##################################################################################################
####RELATIONSHIPS WITH SPATIAL AND ENVIRONMENTAL DISTANCE AND BETA-DIVERSITY                     #
##################################################################################################

##NEED TO CHECK THIS SECTION. 

####Relationship between spatial distance and environmental heterogeneity- across whole landscape####
km.env<- lm(Env_distance_bray~km_between, data=dist.data) #Need to actually go through the linear model. 
plot(km.env)

#Subset for comparisons within the same ecoregion so have an estimate of environmental het comparisons within each ecoregion. 
cpl.comparisons<- subset(dist.data, Ecoregion_sit1 == "CPL" & Ecoregion_sit2 == "CPL") #45 comparisons

nap.comparisons<- subset(dist.data, Ecoregion_sit1 == "NAP" & Ecoregion_sit2 == "NAP") #1827 comparisons

spl.comparisons<- subset(dist.data, Ecoregion_sit1 == "SPL" & Ecoregion_sit2 == "SPL") #10 comparisons 

tpl.comparisons<- subset(dist.data, Ecoregion_sit1 == "TPL" & Ecoregion_sit2 == "TPL") #15 comparisons 

umw.comparisons<- subset(dist.data, Ecoregion_sit1 == "UMW" & Ecoregion_sit2 == "UMW") #2346 comparisons 

wmt.comparisons<- subset(dist.data, Ecoregion_sit1 == "WMT" & Ecoregion_sit2 == "WMT") #561 comparisons

xer.comparisons<- subset(dist.data, Ecoregion_sit1 == "XER" & Ecoregion_sit2 == "XER") #no comparisons 

#Relationship between spatial distance and env het for each ecoregion
km.env.cpl<- lm(Env_distance_bray~km_between, data=cpl.comparisons)
summary(km.env.cpl)

km.env.nap<- lm(Env_distance_bray~km_between, data=nap.comparisons)
summary(km.env.nap)

km.env.spl<- lm(Env_distance_bray~km_between, data=spl.comparisons)
summary(km.env.spl)

km.env.tpl<- lm(Env_distance_bray~km_between, data=tpl.comparisons)
summary(km.env.tpl)

km.env.umw<- lm(Env_distance_bray~km_between, data=umw.comparisons)
summary(km.env.umw)

km.env.wmt<- lm(Env_distance_bray~km_between, data=wmt.comparisons)
summary(km.env.wmt)

#Similar to other results in the literature where sites can be very different from other
#sites that are both close and far. 

#Mean env het and km distance for each ecoregion

#CPL
summary(cpl.comparisons$km_between) #average distance between lakes in CPL = 893km
summary(cpl.comparisons$Env_distance_bray) #mean env distance = 0.35

#NAP
summary(nap.comparisons$km_between) #average distance between lakes in NAP = 455km
summary(nap.comparisons$Env_distance_bray) #mean env distance = 1.88

#SPL
summary(spl.comparisons$km_between) #average distance between lakes in SPL = 106.4km
summary(spl.comparisons$Env_distance_bray) #mean env distance = 0.68

#TPL
summary(tpl.comparisons$km_between) #average distance between lakes in TPL = 488.8km
summary(tpl.comparisons$Env_distance_bray) #mean env distance = 0.56

#UMW
summary(umw.comparisons$km_between) #average distance between lakes in UMW = 487.7km
summary(umw.comparisons$Env_distance_bray) #mean env distance = 2.9

#WMT
summary(wmt.comparisons$km_between) #average distance between lakes in WMT = 615.8km 
summary(wmt.comparisons$Env_distance_bray) #mean env distance = 0.33 

##Relationship between spatial Beta-diversity D-matrices and spatial and env distance 
km.surfacebeta<- lm(Surface_spatialBeta~km_between, data=dist.data)
summary(km.surfacebeta)

km.historicalbeta<- lm(Historical_spatialBeta~km_between, data=dist.data)
summary(km.historicalbeta)

env.surfacebeta<- lm(Surface_spatialBeta~Env_distance_bray, data=dist.data)
summary(env.surfacebeta)

env.historicalbeta<- lm(Historical_spatialBeta~Env_distance_bray, data=dist.data)
summary(env.historicalbeta)

surfbeta.histbeta<- lm(Historical_spatialBeta~Surface_spatialBeta, data=dist.data)
summary(surfbeta.histbeta)

surfbeta.histbeta.plot<- ggplot(dist.data, aes(x=Surface_spatialBeta, y=Historical_spatialBeta)) + geom_point() + geom_smooth(method="lm", se=FALSE, size=1)
surfbeta.histbeta.plot<- surfbeta.histbeta.plot + labs(x="2007 spatial beta diversity", y="Historical spatial beta diversity") 
surfbeta.histbeta.plot<- surfbeta.histbeta.plot + theme_bw()
surfbeta.histbeta.plot<- surfbeta.histbeta.plot + theme(axis.text.x = element_text(colour="black", size=16))
surfbeta.histbeta.plot<- surfbeta.histbeta.plot + theme(axis.text.y = element_text(colour="black", size=16))
surfbeta.histbeta.plot<- surfbeta.histbeta.plot + theme(axis.title.x = element_text(size = rel(2), angle=00))
surfbeta.histbeta.plot<- surfbeta.histbeta.plot + theme(axis.title.y = element_text(size = rel(2), angle=90))

##################################################################################################

##################################################################################################
####BETA-DIVERISTY AND URTs of LANDUSE/WATER QUALITY                                             #
##################################################################################################

####Making URTs with LCBD and temporal BD for water quality and land cover variables####
##Surface LCBD and water quality
surf.LCBD.env<- as.data.frame(cbind(surf.LCBD.summary, waterqual.dat.red[,2:8], nla.highC.surf$ORIGIN))
colnames(surf.LCBD.env)[14]<- 'ORIGIN'  #though only a few man-made, majority are natural. 

surf.fit2<- rpart(LCBD~ORIGIN + WSA_ECOREGION + LONG_DD + LAT_DD + SECMEAN + CHLA + PTL + NTL + MEAN_T + COND + PH_FIELD, method="anova", data=surf.LCBD.env)
printcp(surf.fit2) #Variables actually used in construction: Chla, Cond, Lat_DD, Long_DD, Mean temp, TN, TP, Secchi
plot(surf.fit2)
summary(surf.fit2)
surf.fit2prune<-  prune(surf.fit2, cp=surf.fit2$cptable[which.min(surf.fit2$cptable[,"xerror"]),"CP"])
summary(surf.fit2prune)

par(mfrow=c(1,2)) 
rsq.rpart(surf.fit2prune)  # visualize cross-validation result
#R2 is 1-Rel Error
# R2      
1-surf.fit2$cptable[7,3]   #column 3, row 7 in cptable --> this gives you the adjusted R2 of the model. 

##Surface LCBD and land cover
surf.LCBD.land<- as.data.frame(cbind(surf.LCBD.summary, nla.highC.surf[,25:42]))

surf.fit3<-  rpart(LCBD~PCT_DEVELOPED_BSN + PCT_BARREN_BSN + PCT_FOREST_BSN + PCT_AGRIC_BSN + PCT_WETLAND_BSN, method="anova", data=surf.LCBD.land)
printcp(surf.fit3) #Variables actually used in construction: AGRIC, DEVELOPED, FOREST, WETLAND
plot(surf.fit3)
summary(surf.fit3)
surf.fit3prune<-  prune(surf.fit3, cp=surf.fit3$cptable[which.min(surf.fit3$cptable[,"xerror"]),"CP"])
summary(surf.fit3prune)

par(mfrow=c(1,2)) 
rsq.rpart(surf.fit2prune)  # visualize cross-validation result
#R2 is 1-Rel Error
# R2      
1-surf.fit3$cptable[9,3]   #column 3, row 9 in cptable --> this gives you the adjusted R2 of the model. 

##Temporal BD and land cover
temporalBD.landuse<- as.data.frame(cbind(surf.diat.abundnon.label$SITE_ID, surf.diat.abundnon.label$WSA_ECOREGION, surf.diat.abundnon.label$LONG_DD, surf.diat.abundnon.label$LAT_DD, nla.highC.surf$ORIGIN, nla.highC.surf$LAKE_TYPE, nla.highC.surf$PCT_DEVELOPED_BSN, nla.highC.surf$PCT_FOREST_BSN, nla.highC.surf$PCT_SHRUBLAND_BSN, nla.highC.surf$PCT_GRASS_BSN, nla.highC.surf$PCT_AGRIC_BSN, nla.highC.surf$PCT_WETLAND_BSN, temporalBD.BCD[,9:15]))
colnames(temporalBD.landuse)[1]<- 'SITE_ID'
colnames(temporalBD.landuse)[2]<- 'WSA_ECOREGION'
colnames(temporalBD.landuse)[3]<- 'LONG_DD'
colnames(temporalBD.landuse)[4]<- 'LAT_DD'
colnames(temporalBD.landuse)[5]<- 'ORIGIN' 
colnames(temporalBD.landuse)[6]<- 'LAKE_TYPE' 
colnames(temporalBD.landuse)[7]<- 'PCT_DEVELOPED_BSN'
colnames(temporalBD.landuse)[8]<- 'PCT_FOREST_BSN'
colnames(temporalBD.landuse)[9]<- 'PCT_SHRUBLAND_BSN'
colnames(temporalBD.landuse)[10]<- 'PCT_GRASS_BSN'
colnames(temporalBD.landuse)[11]<- 'PCT_AGRIC_BSN'
colnames(temporalBD.landuse)[12]<- 'PCT_WETLAND_BSN'
#All same units (percent)- don't need to standardize. 

temporalland.fit1<- rpart(Total_beta~PCT_DEVELOPED_BSN + PCT_FOREST_BSN + PCT_SHRUBLAND_BSN + PCT_GRASS_BSN + PCT_AGRIC_BSN + PCT_WETLAND_BSN, method="anova", data=temporalBD.landuse) 
printcp(temporalland.fit1) #Variables actually used in construction: AGRIC, DEVELOPED, GRASS, FOREST, SHRUBLAND, WETLAND
summary(temporalland.fit1)
temporal.land1prune<-  prune(temporalland.fit1, cp=temporalland.fit1$cptable[which.min(temporalland.fit1$cptable[,"xerror"]),"CP"])
summary(temporal.land1prune)
#R2 is 1-Rel Error
# R2      
1-temporalland.fit1$cptable[10,3]  

##Temporal species gain component and land cover 
temporalland.fit2<- rpart(Species_gain_component~PCT_DEVELOPED_BSN + PCT_FOREST_BSN + PCT_SHRUBLAND_BSN + PCT_GRASS_BSN + PCT_AGRIC_BSN + PCT_WETLAND_BSN, method="anova", data=temporalBD.landuse) 
printcp(temporalland.fit2) #Variables actually used in construction: AGRIC, DEVELOPED, FOREST, SHRUBLAND, WETLAND
summary(temporalland.fit2)
temporal.land2prune<-  prune(temporalland.fit2, cp=temporalland.fit2$cptable[which.min(temporalland.fit2$cptable[,"xerror"]),"CP"])
summary(temporal.land2prune)
#R2 is 1-Rel Error
# R2      
1-temporalland.fit2$cptable[7,3]  

##################################################################################################

##################################################################################################
####URT FIGURES                                                                                  #
##################################################################################################

##NEED TO CHECK THIS SECTION- FOR CONSISTENCY WITH OTHER RESULTS. 

####Using URT results to make rough figures for manuscript####
##FIGURE 4, LCBD trees
#(4A)- Water quality
plot(as.party(surf.fit2prune), tp_args = list(id = FALSE))

#(4B)- Land cover
plot(as.party(surf.fit3prune), tp_args = list(id = FALSE)) 

##FIGURE 6, temporal BD trees
#(6A)- Total BD and land cover
plot(as.party(temporal.land1prune), tp_args = list(id = FALSE)) #Pruned results only in forest

#(6B)- Species gain component and land cover 
plot(as.party(temporal.land2prune), tp_args = list(id = FALSE)) #Pruned results: forest, shurbland, wetland

##################################################################################################

##################################################################################################
####PLOTS OF CERTAIN DIATOM SPECIES ABUNDANCE WITH HIGH SCBD FOR S2                              #
##################################################################################################

####Diatom plots for SCBD for Supplementary Information 2####
##Plot abundance between historical and 2007 sediments of diatom taxa with highest SCBD in both sediment sets. 

A.formosa<- as.data.frame(cbind(hist.diat.abund$SITE_ID, hist.diat.abund$Asterionella.formosa.Hassal, surf.diat.abund$Asterionella.formosa.Hassal))
colnames(A.formosa) [1]<- 'SITE_ID'
colnames(A.formosa) [2]<- 'Hist_A.formosa'
colnames(A.formosa) [3]<- 'Surf_A.formosa'

A.formosa.long<- melt(A.formosa, id.vars=c("SITE_ID"))

A.ambigua<- as.data.frame(cbind(hist.diat.abund$SITE_ID, hist.diat.abund$Aulacoseira.ambigua..Grunow..Simonsen, surf.diat.abund$Aulacoseira.ambigua..Grunow..Simonsen))
colnames(A.ambigua) [1]<- 'SITE_ID'
colnames(A.ambigua) [2]<- 'Hist_A.ambigua'
colnames(A.ambigua) [3]<- 'Surf_A.ambigua'

A.ambigua.long<- melt(A.ambigua, id.vars=c("SITE_ID"))

D.pseudostelligera<- as.data.frame(cbind(hist.diat.abund$SITE_ID, hist.diat.abund$Discostella.pseudostelligera..Hustedt..Houk.et.Klee, surf.diat.abund$Discostella.pseudostelligera..Hustedt..Houk.et.Klee))
colnames(D.pseudostelligera) [1]<- 'SITE_ID'
colnames(D.pseudostelligera) [2]<- 'Hist_D.pseudo'
colnames(D.pseudostelligera) [3]<- 'Surf_D.pseudo'

D.pseudostelligera.long<- melt(D.pseudostelligera, id.vars=c("SITE_ID"))

T.flocculosa<- as.data.frame(cbind(hist.diat.abund$SITE_ID, hist.diat.abund$Tabellaria.flocculosa..Roth..K�tzing, surf.diat.abund$Tabellaria.flocculosa..Roth..K�tzing))
colnames(T.flocculosa) [1]<- 'SITE_ID'
colnames(T.flocculosa) [2]<- 'Hist_T.flocculosa'
colnames(T.flocculosa) [3]<- 'Surf_T.flocculosa'

T.flocculosa.long<- melt(T.flocculosa, id.vars=c("SITE_ID"))

S.construens<- as.data.frame(cbind(hist.diat.abund$SITE_ID, hist.diat.abund$Staurosira.construens.Ehrenberg, surf.diat.abund$Staurosira.construens.Ehrenberg))
colnames(S.construens) [1]<- 'SITE_ID'
colnames(S.construens) [2]<- 'Hist_S.construens'
colnames(S.construens) [3]<- 'Surf_S.construens'

S.construens.long<- melt(S.construens, id.vars=c("SITE_ID"))

##Plots for Supplementary info 1
A.formosa.plot<- ggplot(A.formosa.long, aes(x=variable, y=value)) + geom_boxplot()
A.formosa.plot<- A.formosa.plot + labs(x="Sediment type", y="Abundance (count)") 
A.formosa.plot<- A.formosa.plot + theme_bw()
A.formosa.plot<- A.formosa.plot + theme(axis.text.x = element_text(colour="black", size=16))
A.formosa.plot<- A.formosa.plot + theme(axis.text.y = element_text(colour="black", size=16))
A.formosa.plot<- A.formosa.plot + theme(axis.title.x = element_text(size = rel(2), angle=00))
A.formosa.plot<- A.formosa.plot + theme(axis.title.y = element_text(size = rel(2), angle=90))
#A.formosa.plot<- A.formosa.plot + annotate("text", x=0.5, y=450, label="(a)", size=10)

A.ambigua.plot<- ggplot(A.ambigua.long, aes(x=variable, y=value)) + geom_boxplot()
A.ambigua.plot<- A.ambigua.plot + labs(x="Sediment type", y="Abundance (count)") 
A.ambigua.plot<- A.ambigua.plot + theme_bw()
A.ambigua.plot<- A.ambigua.plot + theme(axis.text.x = element_text(colour="black", size=16))
A.ambigua.plot<- A.ambigua.plot + theme(axis.text.y = element_text(colour="black", size=16))
A.ambigua.plot<- A.ambigua.plot + theme(axis.title.x = element_text(size = rel(2), angle=00))
A.ambigua.plot<- A.ambigua.plot + theme(axis.title.y = element_text(size = rel(2), angle=90))
#A.ambigua.plot<- A.ambigua.plot + annotate("text", x=0.5, y=450, label="(b)", size=10)

D.pseudostelligera.plot<- ggplot(D.pseudostelligera.long, aes(x=variable, y=value)) + geom_boxplot()
D.pseudostelligera.plot<- D.pseudostelligera.plot + labs(x="Sediment type", y="Abundance (count)") 
D.pseudostelligera.plot<- D.pseudostelligera.plot + theme_bw()
D.pseudostelligera.plot<- D.pseudostelligera.plot + theme(axis.text.x = element_text(colour="black", size=16))
D.pseudostelligera.plot<- D.pseudostelligera.plot + theme(axis.text.y = element_text(colour="black", size=16))
D.pseudostelligera.plot<- D.pseudostelligera.plot + theme(axis.title.x = element_text(size = rel(2), angle=00))
D.pseudostelligera.plot<- D.pseudostelligera.plot + theme(axis.title.y = element_text(size = rel(2), angle=90))
#D.pseudostelligera.plot<- D.pseudostelligera.plot + annotate("text", x=0.5, y=450, label="(c)", size=10)

T.flocculosa.plot<- ggplot(T.flocculosa.long, aes(x=variable, y=value)) + geom_boxplot()
T.flocculosa.plot<- T.flocculosa.plot + labs(x="Sediment type", y="Abundance (count)") 
T.flocculosa.plot<- T.flocculosa.plot + theme_bw()
T.flocculosa.plot<- T.flocculosa.plot + theme(axis.text.x = element_text(colour="black", size=16))
T.flocculosa.plot<- T.flocculosa.plot + theme(axis.text.y = element_text(colour="black", size=16))
T.flocculosa.plot<- T.flocculosa.plot + theme(axis.title.x = element_text(size = rel(2), angle=00))
T.flocculosa.plot<- T.flocculosa.plot + theme(axis.title.y = element_text(size = rel(2), angle=90))
#T.flocculosa.plot<- T.flocculosa.plot + annotate("text", x=0.5, y=450, label="(d)", size=10)

S.construens.plot<- ggplot(S.construens.long, aes(x=variable, y=value)) + geom_boxplot()
S.construens.plot<- S.construens.plot + labs(x="Sediment type", y="Abundance (count)") 
S.construens.plot<- S.construens.plot + theme_bw()
S.construens.plot<- S.construens.plot + theme(axis.text.x = element_text(colour="black", size=16))
S.construens.plot<- S.construens.plot + theme(axis.text.y = element_text(colour="black", size=16))
S.construens.plot<- S.construens.plot + theme(axis.title.x = element_text(size = rel(2), angle=00))
S.construens.plot<- S.construens.plot + theme(axis.title.y = element_text(size = rel(2), angle=90))
#S.construens.plot<- S.construens.plot + annotate("text", x=0.5, y=450, label="(e)", size=10)

##################################################################################################

##################################################################################################
####LAND COVER PCA FOR S1                                                                        #
##################################################################################################

####Use Land cover data to make other PCA for Supplementary Info 1####
landcover.pca<- rda(decostand(temporalBD.landuse[,7:12], method="standardize"))
summary(landcover.pca)
#Extract PC1 and PC2
landcover.scores<- as.data.frame(scores(landcover.pca, dis="sites", choices=c(1:2))) #Site scores 
landcover.sp.scores<- as.data.frame(scores(landcover.pca, dis="sp", choices=c(1:2))) #Species scores

landcover.pca.plot<-ggplot()
landcover.pca.plot<- landcover.pca.plot + geom_vline(x=0,colour="grey50") 
landcover.pca.plot<- landcover.pca.plot+ geom_hline(y=0,colour="grey50") 
landcover.pca.plot<- landcover.pca.plot + labs(x= "PC1 (34% var explained)", y= "PC2 (23% var exp)") + theme_bw()
landcover.pca.plot<- landcover.pca.plot + geom_point(data = landcover.scores, aes(x = PC1, y = PC2), size=4) 
landcover.pca.plot<- landcover.pca.plot + geom_text(data = landcover.sp.scores, aes(x = PC1, y = PC2, label=rownames(landcover.sp.scores)), size=6) 
landcover.pca.plot<- landcover.pca.plot + theme(axis.text.x = element_text(colour="black", size=16))
landcover.pca.plot<- landcover.pca.plot + theme(axis.text.y = element_text(colour="black", size=16))
landcover.pca.plot<- landcover.pca.plot + theme(axis.title.x = element_text(size = rel(2), angle=00))
landcover.pca.plot<- landcover.pca.plot + theme(axis.title.y = element_text(size = rel(2), angle=90))

##################################################################################################

##################################################################################################
####ALPHA AND GAMMA DIVERSITY                                                                    #
##################################################################################################

####Ongoing calculations of alpha and gamma diversity####
#Use 
#surf.diat.abund   #quantitative data for surface diatoms (~953 species)
#hist.diat.abund  #quantitative data for historical diatoms (~1033 species)

#Add an ecoregion column for subsetting using:
#nla.highC.surf
#nla.highC.hist 

########################This code to clean up the diatoms is not working- fix and then re-do diversity analyses
##NOTE that could use code from Julien to parse the names and make a better cleaned up list for these species
#if put into long form for each- should do? But check beta-diversity estimates too. 

##Historical sediment diatoms
#hist.diat.abund    #in raw form: 1031 species columns + 2 columns 
#Rename to an object so not overwriting any other objects in the workspace
hd.abund<- hist.diat.abund
hd.abund<- hd.abund[,-2] #Remove column 2 ("variable" column, so only SITE_ID and species columns left)

#Melt to make long
hd.abund.long<- melt(hd.abund, id.vars = c("SITE_ID")) 
colnames(hd.abund.long) [2] <- 'TAXANAME'
colnames(hd.abund.long) [3] <- 'ABUND'

#Now clean up species in the TAXANAME column. 
# Create diatom dataset
hist.taxanames <- hd.abund.long$TAXANAME
hist.taxanames.str<-toString(hist.taxanames) #not sure about this. 

## Standardize all taxanames
# Parsing functions
parse_tax_cf <- function(hist.taxanames.str) {
  # Removes "cf." from taxaname
  hist.taxanames.str <- gsub("cf..", " ", hist.taxanames.str) #Remove cf.. from names 
  hist.taxanames.str <- gsub("cff..", " ", hist.taxanames.str) #Remove cff.. from names (don't think there are any)
  hist.taxanames.str <- gsub("var..", " ", hist.taxanames.str) #Remove var.. from names 
  hist.taxanames.str
}
parse_tax_sp <- function(hist.taxanames.str) {
  # Combines sp, spp, sp1, etc. variants
  hist.taxanames.str <- gsub("sp..[0-9].", " sp.", hist.taxanames.str) #Replace sp..# with sp.
  hist.taxanames.str <- gsub("spp.", " sp.", hist.taxanames.str) #Replace spp. with sp.
  hist.taxanames.str <- gsub("sp.[0-9].", " sp.", hist.taxanames.str) #Replace sp.# with sp. 
  hist.taxanames.str <- gsub("spp..", " sp. ", hist.taxanames.str) #Replace spp.. with sp. 
  hist.taxanames.str <- gsub("[?]", "", hist.taxanames.str)
  hist.taxanames.str
}
# Parse taxaname function
parse_tax <- function(hist.taxanames.str) {
  # Combines parsing sub-functions
  hist.taxanames.str <- parse_tax_cf(hist.taxanames.str)
  hist.taxanames.str <- parse_tax_sp(hist.taxanames.str)
  hist.taxanames.str
}
# Parse all taxanames
hist.taxanames_new <- unlist(sapply(hist.taxanames.str, FUN = parse_tax))

# Split and remake taxanames (removing non species epithets)
split_join_tax <- function(hist.taxanames.str) {
  split_string <- strsplit(as.character(hist.taxanames.str), " ")[[1]]
  hist.taxanames.str <- paste(split_string[1], split_string[2], sep=".")
  hist.taxanames.str
}
# Split and join taxanames
hist.taxanames_new <- unlist(sapply(taxanames_new, FUN = split_join_tax))

# Replace diatom taxanames with new taxanames
hd.abund.long$TAXANAME <- hist.taxanames_new

#Recast back into wide format

##Surface sediment diatoms 
#Would use parsing subfunctions as above but surf.taxanames 

sd.abund<- surf.diat.abund
sd.abund<- surf.abund[,-2] 

#Melt to make long
sd.abund.long<- melt(sd.abund, id.vars = c("SITE_ID")) 
colnames(sd.abund.long) [2] <- 'TAXANAME'
colnames(sd.abund.long) [3] <- 'ABUND'

surf.taxanames <- sd.abund.long$TAXANAME
surf.taxanames.str<-toString(hist.taxanames) 

################################################

##Diversity analyses steps ##RE_DO ONCE PARSING CODE FIGURED OUT ABOVE. 

#Take wide-format surface and historical diatom dataframes with quantitative data and bind in ecoregion column
#surf.diat.abund   #quantitative data for surface diatoms (~953 species) #Use for now
#hist.diat.abund  #quantitative data for historical diatoms (~1033 species) #Use for now 

#Add an ecoregion column for subsetting using:
#nla.highC.surf
#nla.highC.hist 

surf.eco<- as.data.frame(cbind(surf.diat.abund, nla.highC.surf$WSA_ECOREGION))
colnames(surf.eco) [954] <- 'WSA_ECOREGION'

hist.eco<- as.data.frame(cbind(hist.diat.abund, nla.highC.hist$WSA_ECOREGION))
colnames(hist.eco) [1034] <- 'WSA_ECOREGION'

#Subset historical sediment dataframe into ecoregions
hist.eco.CPL<- as.data.frame(subset(hist.eco, WSA_ECOREGION == "CPL", drop=T))
hist.eco.NAP<- as.data.frame(subset(hist.eco, WSA_ECOREGION == "NAP", drop=T))
hist.eco.SPL<- as.data.frame(subset(hist.eco, WSA_ECOREGION == "SPL", drop=T))
hist.eco.TPL<- as.data.frame(subset(hist.eco, WSA_ECOREGION == "TPL", drop=T))
hist.eco.UMW<- as.data.frame(subset(hist.eco, WSA_ECOREGION == "UMW", drop=T))
hist.eco.WMT<- as.data.frame(subset(hist.eco, WSA_ECOREGION == "WMT", drop=T))

#Subset surface sediment dataframe into ecoregions
surf.eco.CPL<- as.data.frame(subset(surf.eco, WSA_ECOREGION == "CPL", drop=T))
surf.eco.NAP<- as.data.frame(subset(surf.eco, WSA_ECOREGION == "NAP", drop=T))
surf.eco.SPL<- as.data.frame(subset(surf.eco, WSA_ECOREGION == "SPL", drop=T))
surf.eco.TPL<- as.data.frame(subset(surf.eco, WSA_ECOREGION == "TPL", drop=T))
surf.eco.UMW<- as.data.frame(subset(surf.eco, WSA_ECOREGION == "UMW", drop=T))
surf.eco.WMT<- as.data.frame(subset(surf.eco, WSA_ECOREGION == "WMT", drop=T))

#Rarefied species richness for historical sediment subsets
#All
hist.Srar <- rarefy(hist.eco[,3:1033], min(rowSums(hist.eco[,3:1033]))) 
hist.Srar #String
hist.Srar<- as.data.frame(hist.Srar)
summary(hist.Srar) #can get mean 

#CPL
hist.CPL.Srar <- rarefy(hist.eco.CPL[,3:1033], min(rowSums(hist.eco.CPL[,3:1033]))) 
hist.CPL.Srar #String
hist.CPL.Srar<- as.data.frame(hist.CPL.Srar)
summary(hist.CPL.Srar) #can get mean 

#NAP
hist.NAP.Srar <- rarefy(hist.eco.NAP[,3:1033], min(rowSums(hist.eco.NAP[,3:1033])))
hist.NAP.Srar #String
hist.NAP.Srar<- as.data.frame(hist.NAP.Srar)
summary(hist.NAP.Srar) #can get mean 

#SPL
hist.SPL.Srar <- rarefy(hist.eco.SPL[,3:1033], min(rowSums(hist.eco.SPL[,3:1033])))
hist.SPL.Srar #String
hist.SPL.Srar<- as.data.frame(hist.SPL.Srar)
summary(hist.SPL.Srar) #can get mean 

#TPL
hist.TPL.Srar <- rarefy(hist.eco.TPL[,3:1033], min(rowSums(hist.eco.TPL[,3:1033])))
hist.TPL.Srar #String
hist.TPL.Srar<- as.data.frame(hist.TPL.Srar)
summary(hist.TPL.Srar) #can get mean 

#UMW
hist.UMW.Srar <- rarefy(hist.eco.UMW[,3:1033], min(rowSums(hist.eco.UMW[,3:1033])))
hist.UMW.Srar #String
hist.UMW.Srar<- as.data.frame(hist.UMW.Srar)
summary(hist.UMW.Srar) #can get mean 

#WMT
hist.WMT.Srar <- rarefy(hist.eco.WMT[,3:1033], min(rowSums(hist.eco.WMT[,3:1033])))
hist.WMT.Srar #String
hist.WMT.Srar<- as.data.frame(hist.WMT.Srar)
summary(hist.WMT.Srar) #can get mean 

#Rarefied species richness for surface sediment subsets 
#All
surf.Srar <- rarefy(surf.eco[,4:953], min(rowSums(surf.eco[,4:953]))) #extra column in front of Site_ID
surf.Srar #String
surf.Srar<- as.data.frame(surf.Srar)
summary(surf.Srar) #can get mean 

#CPL
surf.CPL.Srar <- rarefy(surf.eco.CPL[,4:953], min(rowSums(surf.eco.CPL[,4:953]))) #extra column in front of Site_ID
surf.CPL.Srar #String
surf.CPL.Srar<- as.data.frame(surf.CPL.Srar)
summary(surf.CPL.Srar) #can get mean 

#NAP
surf.NAP.Srar <- rarefy(surf.eco.NAP[,4:953], min(rowSums(surf.eco.NAP[,4:953]))) #extra column in front of Site_ID
surf.NAP.Srar #String
surf.NAP.Srar<- as.data.frame(surf.NAP.Srar)
summary(surf.NAP.Srar) #can get mean 

#SPL
surf.SPL.Srar <- rarefy(surf.eco.SPL[,4:953], min(rowSums(surf.eco.SPL[,4:953]))) #extra column in front of Site_ID
surf.SPL.Srar #String
surf.SPL.Srar<- as.data.frame(surf.SPL.Srar)
summary(surf.SPL.Srar) #can get mean 

#TPL
surf.TPL.Srar <- rarefy(surf.eco.TPL[,4:953], min(rowSums(surf.eco.TPL[,4:953]))) #extra column in front of Site_ID
surf.TPL.Srar #String
surf.TPL.Srar<- as.data.frame(surf.TPL.Srar)
summary(surf.TPL.Srar) #can get mean 

#UMW
surf.UMW.Srar <- rarefy(surf.eco.UMW[,4:953], min(rowSums(surf.eco.UMW[,4:953]))) #extra column in front of Site_ID
surf.UMW.Srar #String
surf.UMW.Srar<- as.data.frame(surf.UMW.Srar)
summary(surf.UMW.Srar) #can get mean 

#WMT
surf.WMT.Srar <- rarefy(surf.eco.WMT[,4:953], min(rowSums(surf.eco.WMT[,4:953]))) #extra column in front of Site_ID
surf.WMT.Srar #String
surf.WMT.Srar<- as.data.frame(surf.WMT.Srar)
summary(surf.WMT.Srar) #can get mean 


#Shannon alpha diversity for historical sediment subsets
#All
hist.Shannon<- diversity(hist.eco[,3:1033], index = "shannon")
summary(hist.Shannon)

#CPL
hist.CPL.Shannon<- diversity(hist.eco.CPL[,3:1033], index = "shannon")
summary(hist.CPL.Shannon)

#NAP
hist.NAP.Shannon<- diversity(hist.eco.NAP[,3:1033], index = "shannon")
summary(hist.NAP.Shannon)

#SPL
hist.SPL.Shannon<- diversity(hist.eco.SPL[,3:1033], index = "shannon")
summary(hist.SPL.Shannon)

#TPL
hist.TPL.Shannon<- diversity(hist.eco.TPL[,3:1033], index = "shannon")
summary(hist.TPL.Shannon)

#UMW
hist.UMW.Shannon<- diversity(hist.eco.UMW[,3:1033], index = "shannon")
summary(hist.UMW.Shannon)

#WMT
hist.WMT.Shannon<- diversity(hist.eco.WMT[,3:1033], index = "shannon")
summary(hist.WMT.Shannon)

#Shannon alpha diversity for surface sediment subsets
#All
surf.Shannon<- diversity(surf.eco[,4:953], index = "shannon")
summary(surf.Shannon)

#CPL
surf.CPL.Shannon<- diversity(surf.eco.CPL[,4:953], index = "shannon")
summary(surf.CPL.Shannon)

#NAP
surf.NAP.Shannon<- diversity(surf.eco.NAP[,4:953], index = "shannon")
summary(surf.NAP.Shannon)

#SPL
surf.SPL.Shannon<- diversity(surf.eco.SPL[,4:953], index = "shannon")
summary(surf.SPL.Shannon)

#TPL
surf.TPL.Shannon<- diversity(surf.eco.TPL[,4:953], index = "shannon")
summary(surf.TPL.Shannon)

#UMW
surf.UMW.Shannon<- diversity(surf.eco.UMW[,4:953], index = "shannon")
summary(surf.UMW.Shannon)

#WMT
surf.WMT.Shannon<- diversity(surf.eco.WMT[,4:953], index = "shannon")
summary(surf.WMT.Shannon)


#Simpson alpha for historical sediment subsets 
#All
hist.Simpson<- diversity(hist.eco[,3:1033], index = "simpson")
summary(hist.Simpson)

#CPL
hist.CPL.Simpson<- diversity(hist.eco.CPL[,3:1033], index = "simpson")
summary(hist.CPL.Simpson)

#NAP
hist.NAP.Simpson<- diversity(hist.eco.NAP[,3:1033], index = "simpson")
summary(hist.NAP.Simpson)

#SPL
hist.SPL.Simpson<- diversity(hist.eco.SPL[,3:1033], index = "simpson")
summary(hist.SPL.Simpson)

#TPL
hist.TPL.Simpson<- diversity(hist.eco.TPL[,3:1033], index = "simpson")
summary(hist.TPL.Simpson)

#UMW
hist.UMW.Simpson<- diversity(hist.eco.UMW[,3:1033], index = "simpson")
summary(hist.UMW.Simpson)

#WMT
hist.WMT.Simpson<- diversity(hist.eco.WMT[,3:1033], index = "simpson")
summary(hist.WMT.Simpson)

#Simpson slpha for surface sediment subsets
#All
surf.Simpson<- diversity(surf.eco[,4:953], index = "simpson")
summary(surf.Simpson)

#CPL
surf.CPL.Simpson<- diversity(surf.eco.CPL[,4:953], index = "simpson")
summary(surf.CPL.Simpson)

#NAP
surf.NAP.Simpson<- diversity(surf.eco.NAP[,4:953], index = "simpson")
summary(surf.NAP.Simpson)

#SPL
surf.SPL.Simpson<- diversity(surf.eco.SPL[,4:953], index = "simpson")
summary(surf.SPL.Simpson)

#TPL
surf.TPL.Simpson<- diversity(surf.eco.TPL[,4:953], index = "simpson")
summary(surf.TPL.Simpson)

#UMW
surf.UMW.Simpson<- diversity(surf.eco.UMW[,4:953], index = "simpson")
summary(surf.UMW.Simpson)

#WMT
surf.WMT.Simpson<- diversity(surf.eco.WMT[,4:953], index = "simpson")
summary(surf.WMT.Simpson)


#Gamma diversity for historical sediment subsets (Alpha-Shannon * BD) #so do not need to compute. 
#Gamma diversity for surface sediment subsets 

##################################################################################################

##################################################################################################
####URTs WITH ALPHA DIVERSITY                                                                    #
##################################################################################################

####Use alpha diversity from surface sediments as a responese variable and do URTs with land cover and water quality####
surf.Shannon

surf.Shannon.dframe<- as.data.frame(surf.Shannon)

#Water quality
surf.Shannon.env<- as.data.frame(cbind(surf.Shannon.dframe, waterqual.dat.red[,2:8])) #Didn't include latitude, ecoregion etc. 
surf.shannon.env2<- rpart(surf.Shannon~SECMEAN + CHLA + PTL + NTL + MEAN_T + COND + PH_FIELD, method="anova", data=surf.Shannon.env) 
printcp(surf.shannon.env2) #Variables actually used in construction: Cond, Mean_T, NTL, pH, PTL, Secchi
plot(surf.shannon.env2)
summary(surf.shannon.env2)
surf.shannonenvprune<-  prune(surf.shannon.env2, cp=surf.shannon.env2$cptable[which.min(surf.shannon.env2$cptable[,"xerror"]),"CP"])
summary(surf.shannonenvprune)

par(mfrow=c(1,2)) 
rsq.rpart(surf.shannonenvprune)  # visualize cross-validation result
#R2 is 1-Rel Error
# R2      
1-surf.shannonenvprune$cptable[2,3]   #column 3, row 2 in cptable --> this gives you the adjusted R2 of the model. 
#0.15
plot(as.party(surf.shannonenvprune), tp_args = list(id = FALSE))

#Land cover
surf.Shannon.land<- as.data.frame(cbind(surf.Shannon.dframe, nla.highC.surf[,25:42]))
surf.shannon.land2<-  rpart(surf.Shannon~PCT_DEVELOPED_BSN + PCT_BARREN_BSN + PCT_FOREST_BSN + PCT_AGRIC_BSN + PCT_WETLAND_BSN, method="anova", data=surf.Shannon.land)
printcp(surf.shannon.land2) #Variables actually used in construction: AGRIC, BARREN, DEVELOPED, FOREST, WETLAND
plot(surf.shannon.land2)
summary(surf.shannon.land2)
surf.shannonlandprune<-  prune(surf.shannon.land2, cp=surf.shannon.land2$cptable[which.min(surf.shannon.land2$cptable[,"xerror"]),"CP"])
summary(surf.shannonlandprune)

par(mfrow=c(1,2)) 
rsq.rpart(surf.shannonlandprune)  # visualize cross-validation result
#R2 is 1-Rel Error
# R2      
1-surf.shannonlandprune$cptable[8,3]   #column 3, row 8 in cptable --> this gives you the adjusted R2 of the model. 
#0.33
plot(as.party(surf.shannonlandprune), tp_args = list(id = FALSE))
##################################################################################################