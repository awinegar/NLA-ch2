##################################################################################################
####NLA 2007 CHAPTER 2 PROJECT: DIATOM BETA-DIVERSITY ACROSS THE US (historical/surface seds)    #
##################################################################################################

##In this script: Beta-diversity focused analyses and beta-diversity partitioning 
#Beta-diversity focused analyses 
#Previous script: Rscript_NLA2007ch2_analysesBetaDiv_V7.R (and V8)
#R version: 3.1.2 (Pumpkin Helmet)

##Last update: April 24, 2017
##Associated workspace: workspace_NLA2007ch2_analysesBetaDiv_V9.RData 
##Associated markdown: 
##Github: 

##NOTE: THIS SCRIPT IS FOR RE-DO OF ANALYSES WITHOUT FLORIDA CORES FOR APRIL 2017 RESUBMISSION
#Load V7 workspace and then resave as V9. 

##################################################################################################
setwd("C:/Users/Winegardner/Documents/MCGILL/PhD chapters and projects/EPA National lakes Assessment/Chapter 2/Final_analyses_data")

#Make sure have loaded V9 workspace. 

##Packages
library(data.table) 
library(plyr) 
library(vegan) 
library(lme4)  
library(sp) 
library(rpart) 
library(MASS) 
library(reshape) 
library(reshape2) 
library(party) 
library(rpart.plot) 
library(partykit)
library(tree) 
library(permute) 
library(boot) 
library(rich) 
library(Iso) 
library(corrplot) 

##Visuals 
library(ggplot2) 
library(maptools) 
library(maps)  
library(rgeos)  
library(RColorBrewer) 
library(gridExtra) 

##Baselga analyses
library(betapart)


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

####BETA.DIV.COMP2 SOURCE FUNCTION- updated by Pierre in May 2017 and will replace BETA.DIV.COMP. Simpler code, same results
##################################################################################################
beta.div.comp2 <- function(mat, coef="J", quant=FALSE, save.abc=FALSE)
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
#        "BJ" : Baselga family, Jaccard-based indices
#        "N" : Podani & Schmera (2011) relativized nestedness index.
#        The quantitative form of Sørensen index is the percentage difference.
#        The quantitative form of the Jaccard index is the Ruzicka index.
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
# B = sum of abundances at site 1 minus A, C = sum of abundances at site 2 minus A.
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
# components of beta diversity. Global Ecology and Biogeography, 23, 1324-1334.
#
# Legendre, P. and M. De Cáceres. 2013. Beta diversity as the variance of community data: 
# dissimilarity coefficients and partitioning. Ecology Letters 16: 951-963.
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
                  a=as.dist(a), b=as.dist(t(b)), c=as.dist(t(c)))
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
    # A = sum of minima in among-site comparisons (= W in L&L 2012, p. 285 and 311)
    # B = site.1 sum - A
    # C = site.2 sum - A
    repl <- matrix(0,n,n)
    rich <- matrix(0,n,n)
    D <- matrix(0,n,n)
    rownames(repl) <- rownames(rich) <- rownames(D) <- noms
    
    for(i in 2:n) {
      for(j in 1:(i-1)) {
        tmp = mat[i,] - mat[j,]
        A = sum(pmin(mat[i,], mat[j,]))
        B = sum(tmp[tmp>0])       # Sum of the species losses between T1 and T2
        C = -sum(tmp[tmp<0])      # Sum of the species gains between T1 and T2
        if(coef==1|| coef==3) { 
          den <- (2*A+B+C)      # Sørensen-based (percentage difference) components
        } else if(coef==2|| coef==4) {
          den <- (A+B+C)        # Jaccard-based (Ruzicka) components
        }
        #
        if(coef==1 || coef==2) {  # Podani indices: percentage difference, Ruzicka
          repl[i,j] <- 2*(min(B,C))/den		# 2*min(B,C)/den
          rich[i,j] <- abs(B-C)/den           # abs(B-C)/den
          D[i,j] <- (B+C)/den     			# (B+C)/den
        }
        #
        # Baselga (2013): quantitative extensions of the Baselga (2010) indices
        if(coef==3) {  # Baselga indices: percentage difference
          repl[i,j] <- (min(B,C))/(A+min(B,C))
          rich[i,j] <- abs(B-C)*A/(den*(A+min(B,C)))
          D[i,j] <- (B+C)/den     			# (B+C)/den
        }
        if(coef==4) {  # Baselga indices: Ruzicka
          repl[i,j] <- 2*(min(B,C))/(A+2*min(B,C))
          rich[i,j] <- abs(B-C)*A/(den*(A+2*min(B,C)))
          D[i,j] <- (B+C)/den     			# (B+C)/den
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
########End of beta.div.comp2 function
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

####TBI SOURCE FUNCTION- for computing temporal beta diversity and gain and loss components- replaces decompose.D2() and paired.diff()
##################################################################################################
TBI <- function(mat1,mat2,method="%difference", pa.tr=FALSE, nperm=99, permute.sp=1, BCD=TRUE, replace=FALSE, clock=FALSE)
{
  ### Internal functions
  dissim <- function(mat1, mat2, n, method, tr=TRUE, BCD, ref)
    # tr=TRUE  : The species data have been transformed by decostand in transform()
    # BCD=TRUE : Method is {"ruzicka", "%difference"} and output table BCD was requested
    # ref=TRUE : The function is called to compute the reference values of the TBI dissimil.
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
    #
    epsilon <- sqrt(.Machine$double.eps)
    method <- match.arg(method, c("%difference", "ruzicka", "hellinger", "chord", "jaccard", "sorensen", "ochiai", "euclidean")) 
    n = nrow(mat1)
    p = ncol(mat1)
    if((nrow(mat2)!=n) | (ncol(mat2)!=p)) stop("The matrices are not of the same size.")
    #
    if(method=="jaccard") { pa.tr=TRUE; method="ruzicka" }
    if(method=="sorensen") { pa.tr=TRUE; method="%difference" }
    if(method=="ochiai") { pa.tr=TRUE; method="hellinger" }
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
      if(method=="%difference") colnames(BCD.mat) <- 
          c("B/(2A+B+C)","C/(2A+B+C)","D=(B+C)/(2A+B+C)")
      if(method=="ruzicka")    colnames(BCD.mat) <- 
          c("B/(A+B+C)","C/(A+B+C)","D=(B+C)/(A+B+C)")
      rownames(BCD.mat) <- paste("Obj",1:n,sep=".")
    } else {
      BCD <- FALSE 
      BCD.mat <- NA }
    ###
    # 1. Compute the reference D for each object from corresponding vectors in the 2 matrices.
    if(tr) { 
      tmp <- dissim(transform(mat1,method),transform(mat2,method),n,method,tr,BCD,ref=FALSE)
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
          if(permute.sp==1) {    # Permutation method 1
            seed <- ceiling(runif(1,max=100000))
            # cat("seed =",seed,'\n')
            set.seed(seed)
            mat1.perm <- apply(mat1,2,sample)
            set.seed(seed)
            mat2.perm <- apply(mat2,2,sample)
          } else {  # Permutation method 2 - Do not force the permutations 
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
  list(TBI=vecD.ref, p.TBI=p.dist, p.adj=p.adj, BCD.mat=BCD.mat)
}

RuzickaD <- function(vec1, vec2, method="ruzicka", BCD=FALSE, ref=TRUE)
  #
  # Compute the Ruzicka dissimilarity (quantitative form of the Jaccard dissimilarity)
  # or the percentage difference (quantitative form of the Sørensen dissimilarity).
  # A single dissimilarity is computed because there are only two data vectors.
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
# res1 = paired.diff(mite[1:10,],mite[61:70,],method="%diff",nperm=999,permute.sp=1)
# 
# Example using permute.sp=3. This method is not recommended (low power).
# res2 = paired.diff(mite[1:10,],mite[61:70,],method="hellinger",nperm=999,permute.sp=3)
##################################################################################################


####DATA SET-UP

##Use V7 objects and remove additional sites (Florida cores)
#Remove 9 cores from Florida 
remove.florida<- c("NLA06608-0161", "NLA06608-0325", "NLA06608-0709", "NLA06608-1413", "NLA06608-1569", "NLA06608-1733", "NLA06608-2117", "NLA06608-FL:18261987", "NLA06608-FL:16674741", "NLA06608-2074", "NLA06608-0069", "NLA06608-0453", "NLA06608-0609", "NLA06608-0865", "NLA06608-1185", "NLA06608-1349", "NLA06608-1825", "NLA06608-1861", "NLA06608-2565", "NLA06608-2593", "NLA06608-2629", "NLA06608-2657", "NLA06608-3169", "NLA06608-FL:107895579", "NLA06608-FL:99324403", "NLA06608-FL:99344895")

#ABIOTIC DATA 
nla.surf #n = 176
nla.hist #n = 176

rownames(nla.surf)<- as.character(nla.surf[,1])
rownames(nla.hist)<- as.character(nla.hist[,1])

nla.surf<- nla.surf[ !(rownames(nla.surf) %in% remove.florida),] #now n= 169 sites
nla.hist<- nla.hist[ !(rownames(nla.hist) %in% remove.florida),] #now n= 169 sites

#WATER QUALITY DATA
waterqual.dat.red

rownames(waterqual.dat.red)<- as.character(waterqual.dat.red[,1])

waterqual.dat.red<- waterqual.dat.red[ !(rownames(waterqual.dat.red) %in% remove.florida),] #now n= 1769 sites

#GENUS LEVEL QUANTITATIVE DIATOMS
#Surface
surf.genus
rownames(surf.genus)<- as.character(surf.genus[,1])
surf.genus<- surf.genus[ !(rownames(surf.genus) %in% remove.florida),] 
#n = 169
#97 genera

#Historical 
hist.genus
rownames(hist.genus)<- as.character(hist.genus[,1])
hist.genus<- hist.genus[ !(rownames(hist.genus) %in% remove.florida),]
#n = 170
#106 genera 

#SPECIES LEVEL QUANTITATIVE DIATOMS
#Surface
surf.sp.red
rownames(surf.sp.red)<- as.character(surf.sp.red[,2])
surf.sp.red<- surf.sp.red[ !(rownames(surf.sp.red) %in% remove.florida),] 
#52 sites
#927 species

#Historical
hist.sp.red
rownames(hist.sp.red)<- as.character(hist.sp.red[,1])
hist.sp.red<- hist.sp.red[ !(rownames(hist.sp.red) %in% remove.florida),]
#52 sites
#1002 species

#SPECIES REDUCED SITES ##NEEDS FIXING 
waterqual.sp.selected<- (waterqual.dat$SITE_ID %in% nla.surf.charles$SITE_ID)
waterqual.sp.red<- waterqual.dat[waterqual.sp.selected,] #n = 59

#Re-order so can check that same- yes. 
waterqual.sp.red<- waterqual.sp.red[ order(match(waterqual.sp.red$SITE_ID, nla.surf.charles$SITE_ID)), ]#ok

#Add WSA_ECOREGION
waterqual.sp.red<- as.data.frame(cbind(waterqual.sp.red, nla.surf.charles$WSA_ECOREGION))
colnames(waterqual.sp.red)[11]<- 'WSA_ECOREGION'

waterqual.sp.red<- waterqual.sp.red[ !(rownames(waterqual.sp.red) %in% remove.florida),]


##DATA
nla.surf #Site info/land use data- surface sediments (high C sites only)   
nla.hist #Site info/land use data- historical sediments (high C sites only) 
waterqual.dat.red #Water qulaity data (contemporary) for the 186 sites 
surf.genus #Quantitative surface diatoms - GENUS LEVEL
hist.genus #Quantitative historical diatoms - GENUS LEVEL
surf.sp.red #Quantitative surface diatoms - SPECIES LEVEL
hist.sp.red #Quantitative historical diatoms - SPECIES LEVEL


####RE-DO OF ANALYSES 

##SPATIAL BETA DIVERSITY 

#GENUS LEVEL- SURFACE
#Bind the non-transformed species columns to site, ecoregion and lat/long columns. 
surf.genus.abundnon<- as.data.frame(cbind(nla.surf[,1], nla.surf[,11], nla.surf[,4], nla.surf[,5], surf.genus[,3:99]))
colnames(surf.genus.abundnon)[1]<- 'SITE_ID'
colnames(surf.genus.abundnon)[2]<- 'WSA_ECOREGION'
colnames(surf.genus.abundnon)[3]<- 'LONG_DD'
colnames(surf.genus.abundnon)[4]<- 'LAT_DD'

#Use percentage difference 
surf.beta.div <- beta.div(surf.genus.abundnon[,5:99], method="percentagedifference", sqrt.D=FALSE, samp=TRUE, nperm=999, save.D=FALSE, clock=FALSE)
summary(surf.beta.div)
surf.beta.div$SStotal_BDtotal #SS total = 44.3, BD total = 0.264
surf.beta.div$SCBD #Null 
surf.beta.div$LCBD
surf.beta.div$p.LCBD

#Make a dataframe with LCBD information for surface data using percentage difference
surf.LCBD.values<- as.data.frame(surf.beta.div$LCBD)
surf.LCBD.p<- as.data.frame(surf.beta.div$p.LCBD)
surf.LCBD.summary<- as.data.frame(cbind(surf.genus.abundnon$SITE_ID, surf.genus.abundnon$WSA_ECOREGION, surf.genus.abundnon$LONG_DD, surf.genus.abundnon$LAT_DD, surf.LCBD.values, surf.LCBD.p))
colnames(surf.LCBD.summary)[1]<- 'SITE_ID'
colnames(surf.LCBD.summary)[2]<- 'WSA_ECOREGION'
colnames(surf.LCBD.summary)[3]<- 'LONG_DD' 
colnames(surf.LCBD.summary)[4]<- 'LAT_DD'
colnames(surf.LCBD.summary)[5]<- 'LCBD'
colnames(surf.LCBD.summary)[6]<- 'LCBD.P'

#Look at sites with significant (p<0.05)
surf.LCBD.summary.sig<- as.data.frame(subset(surf.LCBD.summary, LCBD.P<0.05)) #17 sites

#GENUS LEVEL- HISTORICAL 
#Bind the non-transformed species columns to site, ecoregion and lat/long columns. 
hist.genus.abundnon<- as.data.frame(cbind(nla.hist[,1], nla.hist[,11], nla.hist[,4], nla.hist[,5], hist.genus[,3:108]))
colnames(hist.genus.abundnon)[1]<- 'SITE_ID'
colnames(hist.genus.abundnon)[2]<- 'WSA_ECOREGION'
colnames(hist.genus.abundnon)[3]<- 'LONG_DD'
colnames(hist.genus.abundnon)[4]<- 'LAT_DD'

#Use percentage difference 
hist.beta.div <- beta.div(hist.genus.abundnon[,5:110], method="percentagedifference", sqrt.D=FALSE, samp=TRUE, nperm=999, save.D=FALSE, clock=FALSE)
summary(hist.beta.div)
hist.beta.div$SStotal_BDtotal #SS total = 44.3, BD total = 0.26
hist.beta.div$SCBD #Null 
hist.beta.div$LCBD
hist.beta.div$p.LCBD

#Make a dataframe with LCBD information for surface data using percentage difference
hist.LCBD.values<- as.data.frame(hist.beta.div$LCBD)
hist.LCBD.p<- as.data.frame(hist.beta.div$p.LCBD)
hist.LCBD.summary<- as.data.frame(cbind(hist.genus.abundnon$SITE_ID, hist.genus.abundnon$WSA_ECOREGION, hist.genus.abundnon$LONG_DD, hist.genus.abundnon$LAT_DD, hist.LCBD.values, hist.LCBD.p))
colnames(hist.LCBD.summary)[1]<- 'SITE_ID'
colnames(hist.LCBD.summary)[2]<- 'WSA_ECOREGION'
colnames(hist.LCBD.summary)[3]<- 'LONG_DD' 
colnames(hist.LCBD.summary)[4]<- 'LAT_DD'
colnames(hist.LCBD.summary)[5]<- 'LCBD'
colnames(hist.LCBD.summary)[6]<- 'LCBD.P'

#Look at sites with significant (p<0.05)
hist.LCBD.summary.sig<- as.data.frame(subset(hist.LCBD.summary, LCBD.P<0.05)) #20 sites


##SPATIAL BETA DIVERSITY PARTITIONING

#GENUS LEVEL - SURFACE
####Use non-transformed data to partition spatial beta-diversity into replacement, AbDiff components#### 
##Surface diatoms
surf.abund.bdc<- beta.div.comp(surf.genus.abundnon[,5:101], coef="S", quant=TRUE, save.abc=FALSE)
summary(surf.abund.bdc$repl) #Replacement
summary(surf.abund.bdc$rich) #Abundance difference 
surf.abund.bdc$part
#1st value = Total beta div  
#2nd value = Total replacement diversity
#3rd value = Total abundance difference diversity 
#4th value = Total replcement/Total beta diversity
#5th value = Total abundance difference/Total beta diversity

#Checking using new version of beta.div.comp (condensed)
surf.abund.bdc2<- beta.div.comp2(surf.genus.abundnon[,5:101], coef="S", quant=TRUE, save.abc=FALSE)
summary(surf.abund.bdc2$repl) #Replacement
summary(surf.abund.bdc2$rich) #Abundance difference 
surf.abund.bdcBS$part
#Note: results are the same! Good. 

#Run with Baselga S coefficient
surf.abund.bdcBS<- beta.div.comp2(surf.genus.abundnon[,5:101], coef="BS", quant=TRUE, save.abc=FALSE)
summary(surf.abund.bdcBS$repl) #Replacement
summary(surf.abund.bdcBS$rich) #Abundance difference 
surf.abund.bdcBS$part

##For ecoregions (had to bind new columns)
#CPL
surf.genus.abundnonCPL<- as.data.frame(subset(surf.genus.abundnon, WSA_ECOREGION == "CPL", drop=T))
surf.genus.bdcCPL<- beta.div.comp(surf.genus.abundnonCPL[,5:101], coef="S", quant=TRUE, save.abc=FALSE)
surf.genus.bdcCPL$part

#NAP
surf.genus.abundnonNAP<- as.data.frame(subset(surf.genus.abundnon, WSA_ECOREGION == "NAP", drop=T))
surf.genus.bdcNAP<- beta.div.comp(surf.genus.abundnonNAP[,5:101], coef="S", quant=TRUE, save.abc=FALSE)
surf.genus.bdcNAP$part

#SPL
surf.genus.abundnonSPL<- as.data.frame(subset(surf.genus.abundnon, WSA_ECOREGION == "SPL", drop=T))
surf.genus.bdcSPL<- beta.div.comp(surf.genus.abundnonSPL[,5:101], coef="S", quant=TRUE, save.abc=FALSE)
surf.genus.bdcSPL$part

#TPL
surf.genus.abundnonTPL<- as.data.frame(subset(surf.genus.abundnon, WSA_ECOREGION == "TPL", drop=T))
surf.genus.bdcTPL<- beta.div.comp(surf.genus.abundnonTPL[,5:101], coef="S", quant=TRUE, save.abc=FALSE)
surf.genus.bdcTPL$part

#UMW
surf.genus.abundnonUMW<- as.data.frame(subset(surf.genus.abundnon, WSA_ECOREGION == "UMW", drop=T))
surf.genus.bdcUMW<- beta.div.comp(surf.genus.abundnonUMW[,5:101], coef="S", quant=TRUE, save.abc=FALSE)
surf.genus.bdcUMW$part

#WMT
surf.genus.abundnonWMT<- as.data.frame(subset(surf.genus.abundnon, WSA_ECOREGION == "WMT", drop=T))
surf.genus.bdcWMT<- beta.div.comp(surf.genus.abundnonWMT[,5:101], coef="S", quant=TRUE, save.abc=FALSE)
surf.genus.bdcWMT$part


#GENUS LEVEL - HISTORICAL 
hist.genus.bdc<- beta.div.comp(hist.genus.abundnon[,5:110], coef="S", quant=TRUE, save.abc=FALSE)
summary(hist.genus.bdc$repl) #Replacement
summary(hist.genus.bdc$rich) #Abundance difference 
hist.genus.bdc$part

#Checking with beta.div.comp2
hist.genus.bdc2<- beta.div.comp2(hist.genus.abundnon[,5:110], coef="S", quant=TRUE, save.abc=FALSE)
summary(hist.genus.bdc2$repl) #Replacement
summary(hist.genus.bdc2$rich) #Abundance difference 
hist.genus.bdc2$part

#Running with Baselga S
hist.genus.bdcBS<- beta.div.comp2(hist.genus.abundnon[,5:110], coef="BS", quant=TRUE, save.abc=FALSE)
summary(hist.genus.bdcBS$repl) #Replacement
summary(hist.genus.bdcBS$rich) #Abundance difference 
hist.genus.bdcBS$part



##For ecoregions (had to bind new columns)
#CPL
hist.genus.abundnonCPL<- as.data.frame(subset(hist.genus.abundnon, WSA_ECOREGION == "CPL", drop=T))
hist.genus.bdcCPL<- beta.div.comp(hist.genus.abundnonCPL[,5:110], coef="S", quant=TRUE, save.abc=FALSE)
hist.genus.bdcCPL$part

#NAP
hist.genus.abundnonNAP<- as.data.frame(subset(hist.genus.abundnon, WSA_ECOREGION == "NAP", drop=T))
hist.genus.bdcNAP<- beta.div.comp(hist.genus.abundnonNAP[,5:110], coef="S", quant=TRUE, save.abc=FALSE)
hist.genus.bdcNAP$part

#SPL
hist.genus.abundnonSPL<- as.data.frame(subset(hist.genus.abundnon, WSA_ECOREGION == "SPL", drop=T))
hist.genus.bdcSPL<- beta.div.comp(hist.genus.abundnonSPL[,5:110], coef="S", quant=TRUE, save.abc=FALSE)
hist.genus.bdcSPL$part

#TPL
hist.genus.abundnonTPL<- as.data.frame(subset(hist.genus.abundnon, WSA_ECOREGION == "TPL", drop=T))
hist.genus.bdcTPL<- beta.div.comp(hist.genus.abundnonTPL[,5:110], coef="S", quant=TRUE, save.abc=FALSE)
hist.genus.bdcTPL$part

#UMW
hist.genus.abundnonUMW<- as.data.frame(subset(hist.genus.abundnon, WSA_ECOREGION == "UMW", drop=T))
hist.genus.bdcUMW<- beta.div.comp(hist.genus.abundnonUMW[,5:110], coef="S", quant=TRUE, save.abc=FALSE)
hist.genus.bdcUMW$part

#WMT
hist.genus.abundnonWMT<- as.data.frame(subset(hist.genus.abundnon, WSA_ECOREGION == "WMT", drop=T))
hist.genus.bdcWMT<- beta.div.comp(hist.genus.abundnonWMT[,5:110], coef="S", quant=TRUE, save.abc=FALSE)
hist.genus.bdcWMT$part


##RAREFIED GENUS RICHNESS, ALPHA AND GAMMA DIVERISTY  

##Use quantitative data for surface and historical and subset by ecoregion. 

surf.genus.abundnon
surf.genus.abundnonCPL
surf.genus.abundnonNAP
surf.genus.abundnonSPL
surf.genus.abundnonTPL
surf.genus.abundnonUMW
surf.genus.abundnonWMT

surf.genus.abundnon
hist.genus.abundnonCPL
hist.genus.abundnonNAP
hist.genus.abundnonSPL
hist.genus.abundnonTPL
hist.genus.abundnonUMW
hist.genus.abundnonWMT

##Rarefied richness

#Surface- all
surf.Srar <- rarefy(surf.genus.abundnon[,5:101], min(rowSums(surf.genus.abundnon[,5:101]))) 
surf.Srar #String
surf.Srar<- as.data.frame(surf.Srar)
summary(surf.Srar) #can get mean 

#Surface- CPL
surf.Srar.CPL <- rarefy(surf.genus.abundnonCPL[,5:101], min(rowSums(surf.genus.abundnonCPL[,5:101]))) 
surf.Srar.CPL<- as.data.frame(surf.Srar.CPL)
summary(surf.Srar.CPL) 

#Surface- NAP
surf.Srar.NAP <- rarefy(surf.genus.abundnonNAP[,5:101], min(rowSums(surf.genus.abundnonNAP[,5:101]))) 
surf.Srar.NAP<- as.data.frame(surf.Srar.NAP)
summary(surf.Srar.NAP) 

#Surface- SPL
surf.Srar.SPL <- rarefy(surf.genus.abundnonSPL[,5:101], min(rowSums(surf.genus.abundnonSPL[,5:101]))) 
surf.Srar.SPL<- as.data.frame(surf.Srar.SPL)
summary(surf.Srar.SPL) 

#Surface- TPL
surf.Srar.TPL <- rarefy(surf.genus.abundnonTPL[,5:101], min(rowSums(surf.genus.abundnonTPL[,5:101]))) 
surf.Srar.TPL<- as.data.frame(surf.Srar.TPL)
summary(surf.Srar.TPL) 

#Surface- UMW
surf.Srar.UMW <- rarefy(surf.genus.abundnonUMW[,5:101], min(rowSums(surf.genus.abundnonUMW[,5:101]))) 
surf.Srar.UMW<- as.data.frame(surf.Srar.UMW)
summary(surf.Srar.UMW) 

#Surface- WMT
surf.Srar.WMT <- rarefy(surf.genus.abundnonWMT[,5:101], min(rowSums(surf.genus.abundnonWMT[,5:101]))) 
surf.Srar.WMT<- as.data.frame(surf.Srar.WMT)
summary(surf.Srar.WMT) 


#Historical- all
hist.Srar <- rarefy(hist.genus.abundnon[,5:110], min(rowSums(hist.genus.abundnon[,5:110]))) 
hist.Srar<- as.data.frame(hist.Srar)
summary(hist.Srar) 

#Historical- CPL
hist.Srar.CPL <- rarefy(hist.genus.abundnonCPL[,5:110], min(rowSums(hist.genus.abundnonCPL[,5:110]))) 
hist.Srar.CPL<- as.data.frame(hist.Srar.CPL)
summary(hist.Srar.CPL) 

#Historical- NAP
hist.Srar.NAP <- rarefy(hist.genus.abundnonNAP[,5:110], min(rowSums(hist.genus.abundnonNAP[,5:110]))) 
hist.Srar.NAP<- as.data.frame(hist.Srar.NAP)
summary(hist.Srar.NAP) 

#Historical- SPL
hist.Srar.SPL <- rarefy(hist.genus.abundnonSPL[,5:110], min(rowSums(hist.genus.abundnonSPL[,5:110]))) 
hist.Srar.SPL<- as.data.frame(hist.Srar.SPL)
summary(hist.Srar.SPL) 

#Historical- TPL
hist.Srar.TPL <- rarefy(hist.genus.abundnonTPL[,5:110], min(rowSums(hist.genus.abundnonTPL[,5:110]))) 
hist.Srar.TPL<- as.data.frame(hist.Srar.TPL)
summary(hist.Srar.TPL) 

#Historical- UMW
hist.Srar.UMW <- rarefy(hist.genus.abundnonUMW[,5:110], min(rowSums(hist.genus.abundnonUMW[,5:110]))) 
hist.Srar.UMW<- as.data.frame(hist.Srar.UMW)
summary(hist.Srar.UMW) 

#Historical- WMT
hist.Srar.WMT <- rarefy(hist.genus.abundnonWMT[,5:110], min(rowSums(hist.genus.abundnonWMT[,5:110]))) 
hist.Srar.WMT<- as.data.frame(hist.Srar.WMT)
summary(hist.Srar.WMT) 

#Check relationship between historical rarfied species richness and latitude and longitude (for all sites)

surf.Srar<- as.data.frame(cbind(surf.Srar, surf.genus.abundnon$LAT_DD, surf.genus.abundnon$LONG_DD))
colnames(surf.Srar)[1]<- 'Srar'
colnames(surf.Srar)[2]<- 'Lat'
colnames(surf.Srar)[3]<- 'Long'

surf.Srar.fit<- lm(Srar ~ Lat, data=surf.Srar) #Adj. R2 = -0.005224, P=0.72
summary(surf.Srar.fit)

surf.Srar.fit2<- lm(Srar ~ Long, data=surf.Srar) #Adj. R2 = 0.09, P<0.05
summary(surf.Srar.fit2)


hist.Srar<- as.data.frame(cbind(hist.Srar, hist.genus.abundnon$LAT_DD, hist.genus.abundnon$LONG_DD))
colnames(hist.Srar)[1]<- 'Srar'
colnames(hist.Srar)[2]<- 'Lat'
colnames(hist.Srar)[3]<- 'Long'

hist.Srar.fit<- lm(Srar ~ Lat, data=hist.Srar) #R2 = -0.0041, P=0.58
summary(hist.Srar.fit)

hist.Srar.fit2<- lm(Srar ~ Long, data=hist.Srar) #R2 = 0.09, P<0.05
summary(hist.Srar.fit2)

##Alpha diversity - Shannon + Simpson

#Surface- all
surf.Shannon<- diversity(surf.genus.abundnon[,5:101], index = "shannon")
summary(surf.Shannon)

surf.simpson<- diversity(surf.genus.abundnon[,5:101], index = "simpson")
summary(surf.simpson)

#Surface- CPL
surf.Shannon.CPL<- diversity(surf.genus.abundnonCPL[,5:101], index = "shannon")
summary(surf.Shannon.CPL)

surf.simpson.CPL<- diversity(surf.genus.abundnonCPL[,5:101], index = "simpson")
summary(surf.simpson.CPL)

#Surface- NAP
surf.Shannon.NAP<- diversity(surf.genus.abundnonNAP[,5:101], index = "shannon")
summary(surf.Shannon.NAP)

surf.simpson.NAP<- diversity(surf.genus.abundnonNAP[,5:101], index = "simpson")
summary(surf.simpson.NAP)

#Surface- SPL
surf.Shannon.SPL<- diversity(surf.genus.abundnonSPL[,5:101], index = "shannon")
summary(surf.Shannon.SPL)

surf.simpson.SPL<- diversity(surf.genus.abundnonSPL[,5:101], index = "simpson")
summary(surf.simpson.SPL)

#Surface- TPL
surf.Shannon.TPL<- diversity(surf.genus.abundnonTPL[,5:101], index = "shannon")
summary(surf.Shannon.TPL)

surf.simpson.TPL<- diversity(surf.genus.abundnonTPL[,5:101], index = "simpson")
summary(surf.simpson.TPL)

#Surface- UMW
surf.Shannon.UMW<- diversity(surf.genus.abundnonUMW[,5:101], index = "shannon")
summary(surf.Shannon.UMW)

surf.simpson.UMW<- diversity(surf.genus.abundnonUMW[,5:101], index = "simpson")
summary(surf.simpson.UMW)

#Surface- WMT
surf.Shannon.WMT<- diversity(surf.genus.abundnonWMT[,5:101], index = "shannon")
summary(surf.Shannon.WMT)

surf.simpson.WMT<- diversity(surf.genus.abundnonWMT[,5:101], index = "simpson")
summary(surf.simpson.WMT)


#Historical- all
hist.Shannon<- diversity(hist.genus.abundnon[,5:110], index = "shannon")
summary(hist.Shannon)

hist.simpson<- diversity(hist.genus.abundnon[,5:110], index = "simpson")
summary(hist.simpson)

#Historical- CPL
hist.Shannon.CPL<- diversity(hist.genus.abundnonCPL[,5:110], index = "shannon")
summary(hist.Shannon.CPL)

hist.simpson.CPL<- diversity(hist.genus.abundnonCPL[,5:110], index = "simpson")
summary(hist.simpson.CPL)

#Historical- NAP
hist.Shannon.NAP<- diversity(hist.genus.abundnonNAP[,5:110], index = "shannon")
summary(hist.Shannon.NAP)

hist.simpson.NAP<- diversity(hist.genus.abundnonNAP[,5:110], index = "simpson")
summary(hist.simpson.NAP)

#Historical- SPL
hist.Shannon.SPL<- diversity(hist.genus.abundnonSPL[,5:110], index = "shannon")
summary(hist.Shannon.SPL)

hist.simpson.SPL<- diversity(hist.genus.abundnonSPL[,5:110], index = "simpson")
summary(hist.simpson.SPL)

#Historical- TPL
hist.Shannon.TPL<- diversity(hist.genus.abundnonTPL[,5:110], index = "shannon")
summary(hist.Shannon.TPL)

hist.simpson.TPL<- diversity(hist.genus.abundnonTPL[,5:110], index = "simpson")
summary(hist.simpson.TPL)

#Historical- UMW
hist.Shannon.UMW<- diversity(hist.genus.abundnonUMW[,5:110], index = "shannon")
summary(hist.Shannon.UMW)

hist.simpson.UMW<- diversity(hist.genus.abundnonUMW[,5:110], index = "simpson")
summary(hist.simpson.UMW)

#Historical- WMT
hist.Shannon.WMT<- diversity(hist.genus.abundnonWMT[,5:110], index = "shannon")
summary(hist.Shannon.WMT)

hist.simpson.WMT<- diversity(hist.genus.abundnonWMT[,5:110], index = "simpson")
summary(hist.simpson.WMT)


##Gamma diversity 

#Surface- all
surf.sums<- as.data.frame(colSums(surf.genus.abundnon[,5:101]))
surf.sums<- t(surf.sums)
surf.gamma<- diversity(surf.sums[,1:97], index = "shannon")
summary(surf.gamma) 

#Surface- CPL
surf.sums.CPL<- as.data.frame(colSums(surf.genus.abundnonCPL[,5:101]))
surf.sums.CPL<- t(surf.sums.CPL)
surf.gamma.CPL<- diversity(surf.sums.CPL[,1:97], index = "shannon")
summary(surf.gamma.CPL) 

#Surface- NAP
surf.sums.NAP<- as.data.frame(colSums(surf.genus.abundnonNAP[,5:101]))
surf.sums.NAP<- t(surf.sums.NAP)
surf.gamma.NAP<- diversity(surf.sums.NAP[,1:97], index = "shannon")
summary(surf.gamma.NAP) 

#Surface- SPL
surf.sums.SPL<- as.data.frame(colSums(surf.genus.abundnonSPL[,5:101]))
surf.sums.SPL<- t(surf.sums.SPL)
surf.gamma.SPL<- diversity(surf.sums.SPL[,1:97], index = "shannon")
summary(surf.gamma.SPL) 

#Surface- TPL
surf.sums.TPL<- as.data.frame(colSums(surf.genus.abundnonTPL[,5:101]))
surf.sums.TPL<- t(surf.sums.TPL)
surf.gamma.TPL<- diversity(surf.sums.TPL[,1:97], index = "shannon")
summary(surf.gamma.TPL) 

#Surface- UMW
surf.sums.UMW<- as.data.frame(colSums(surf.genus.abundnonUMW[,5:101]))
surf.sums.UMW<- t(surf.sums.UMW)
surf.gamma.UMW<- diversity(surf.sums.UMW[,1:97], index = "shannon")
summary(surf.gamma.UMW) 

#Surface- WMT
surf.sums.WMT<- as.data.frame(colSums(surf.genus.abundnonWMT[,5:101]))
surf.sums.WMT<- t(surf.sums.WMT)
surf.gamma.WMT<- diversity(surf.sums.WMT[,1:97], index = "shannon")
summary(surf.gamma.WMT) 

#Historical- all
hist.sums<- as.data.frame(colSums(hist.genus.abundnon[,5:110]))
hist.sums<- t(hist.sums)
hist.gamma<- diversity(hist.sums[,1:106], index = "shannon")
summary(hist.gamma) 

#Historical- CPL
hist.sums.CPL<- as.data.frame(colSums(hist.genus.abundnonCPL[,5:110]))
hist.sums.CPL<- t(hist.sums.CPL)
hist.gamma.CPL<- diversity(hist.sums.CPL[,1:106], index = "shannon")
summary(hist.gamma.CPL) 

#Historical- NAP
hist.sums.NAP<- as.data.frame(colSums(hist.genus.abundnonNAP[,5:110]))
hist.sums.NAP<- t(hist.sums.NAP)
hist.gamma.NAP<- diversity(hist.sums.NAP[,1:106], index = "shannon")
summary(hist.gamma.NAP) 

#Historical- SPL
hist.sums.SPL<- as.data.frame(colSums(hist.genus.abundnonSPL[,5:110]))
hist.sums.SPL<- t(hist.sums.SPL)
hist.gamma.SPL<- diversity(hist.sums.SPL[,1:106], index = "shannon")
summary(hist.gamma.SPL) 

#Historical- TPL
hist.sums.TPL<- as.data.frame(colSums(hist.genus.abundnonTPL[,5:110]))
hist.sums.TPL<- t(hist.sums.TPL)
hist.gamma.TPL<- diversity(hist.sums.TPL[,1:106], index = "shannon")
summary(hist.gamma.TPL) 

#Historical- UMW
hist.sums.UMW<- as.data.frame(colSums(hist.genus.abundnonUMW[,5:110]))
hist.sums.UMW<- t(hist.sums.UMW)
hist.gamma.UMW<- diversity(hist.sums.UMW[,1:106], index = "shannon")
summary(hist.gamma.UMW) 

#Historical- WMT
hist.sums.WMT<- as.data.frame(colSums(hist.genus.abundnonWMT[,5:110]))
hist.sums.WMT<- t(hist.sums.WMT)
hist.gamma.WMT<- diversity(hist.sums.WMT[,1:106], index = "shannon")
summary(hist.gamma.WMT) 


#TEMPORAL BETA DIVERSITY (TBI) 
##Temporal BD between historical and surface using decompose.D2 (QUANTITATIVE DATA)
#Need non-transformed quantitative data with same species. 

#Matrices to stack: 
#surf.genus.abundnon #surface diatoms, quantitative, non-transformed, genus col 5:101
#hist.genus.abundnon #historical diatoms, quantitative, non-transformed, genus col 5:110

#Add a column with Surface or Historical to each matrix
#Surface
surf.genus.abundnon.label<- surf.genus.abundnon #create new dataframe that can use for manipulation
surf.genus.abundnon.label$TYPE<- rep("Surface", nrow(surf.genus.abundnon.label))

#Historical
hist.genus.abundnon.label<- hist.genus.abundnon
hist.genus.abundnon.label$TYPE<- rep("Historical", nrow(hist.genus.abundnon.label))

#Melt each so in long format: 
surf.long<- melt(surf.genus.abundnon.label, id.vars=c("SITE_ID", "WSA_ECOREGION", "LONG_DD", "LAT_DD", "TYPE"))
hist.long<- melt(hist.genus.abundnon.label, id.vars=c("SITE_ID", "WSA_ECOREGION", "LONG_DD", "LAT_DD", "TYPE"))

#Stack the long format dataframes on top of eachother: 
combo.long<- rbind(surf.long, hist.long)

#Cast into wideformat
combo.wide<- dcast(combo.long, SITE_ID + WSA_ECOREGION + LONG_DD + LAT_DD + TYPE ~ variable)
#Now for each site there is a surface and historical abundance measurement, common set of species 
#NAs for species not found in one of the matrices. 
#Replace NAs with zeros. 
combo.wide[is.na(combo.wide)]<-0

#Resubset into separate matrices
surf.genus.fortempBD<- as.data.frame(subset(combo.wide, TYPE == "Surface", drop=T)) #6:118 are genus columns
#Make rownames the SITE_ID
rownames(surf.genus.fortempBD)<- as.character(surf.genus.fortempBD[,1])

hist.genus.fortempBD<- as.data.frame(subset(combo.wide, TYPE == "Historical", drop=T)) #6: are genus columns
rownames(hist.genus.fortempBD)<- as.character(hist.genus.fortempBD[,1])


genus.TBI<- TBI(hist.genus.fortempBD[,6:118], surf.genus.fortempBD[,6:118],method="%difference", pa.tr=FALSE, nperm=99, permute.sp=1, BCD=TRUE, replace=FALSE, clock=FALSE)

#Looking at output
genus.TBI
genus.TBI.TBI<- as.data.frame(genus.TBI$TBI)
genus.TBI.pTBI<- as.data.frame(genus.TBI$p.TBI)
genus.TBI.adjpTBI<- as.data.frame(genus.TBI$p.adj)
genus.TBI.BCD<- as.data.frame(genus.TBI$BCD)



#Bind this quantitative output back with ecoregion data
genus.TBI.output<- as.data.frame(cbind(surf.genus.fortempBD[,1:4], genus.TBI.BCD, genus.TBI.pTBI, genus.TBI.adjpTBI))                             
colnames(genus.TBI.output)[5]<- 'Loss' #Genus loss from T1
colnames(genus.TBI.output)[6]<- 'Gain' #Genus gain in T2
colnames(genus.TBI.output)[7]<- 'Total_BD' #Dissimilarity 
colnames(genus.TBI.output)[8]<- 'p_TBI' #p-value
colnames(genus.TBI.output)[9]<- 'p_adj' #p-value adjusted for multiple testing

summary(genus.TBI.output)

TBI.output.summary<- as.data.frame(ddply(genus.TBI.output, "WSA_ECOREGION", summarize, 
                                         mean_Loss = mean(Loss, na.rm=T), 
                                         sd_Loss = sd(Loss, na.rm=T),
                                         mean_Gain = mean(Gain, na.rm=T),
                                         sd_Gain = sd(Gain, na.rm=T),
                                         mean_TotalBD = mean(Total_BD, na.rm=T),
                                         sd_TotalBD = sd(Total_BD, na.rm=T)))




##NOT DOING PART WITH SIGNIFICANT TBI. 


#URTS

##Surface LCBD and water quality
surf.LCBD.env<- as.data.frame(cbind(surf.LCBD.summary, waterqual.dat.red[,2:8], nla.surf$ORIGIN))
colnames(surf.LCBD.env)[14]<- 'ORIGIN'  #though only a few man-made, majority are natural. 

surf.fit2<- rpart(LCBD~ORIGIN + WSA_ECOREGION + LONG_DD + LAT_DD + SECMEAN + CHLA + PTL + NTL + MEAN_T + COND + PH_FIELD, method="anova", data=surf.LCBD.env)
printcp(surf.fit2) #Variables actually used in construction: Chla, Cond, Lat_DD, Long_DD, Mean temp, TN,Secchi
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
surf.LCBD.land<- as.data.frame(cbind(surf.LCBD.summary, nla.surf[,25:42]))

surf.fit3<-  rpart(LCBD~PCT_DEVELOPED_BSN + PCT_BARREN_BSN + PCT_FOREST_BSN + PCT_AGRIC_BSN + PCT_WETLAND_BSN, method="anova", data=surf.LCBD.land)
printcp(surf.fit3) #Variables actually used in construction: AGRIC, DEVELOPED, FOREST, WETLAND, BARREN
plot(surf.fit3)
summary(surf.fit3)
surf.fit3prune<-  prune(surf.fit3, cp=surf.fit3$cptable[which.min(surf.fit3$cptable[,"xerror"]),"CP"])
summary(surf.fit3prune)

par(mfrow=c(1,2)) 
rsq.rpart(surf.fit2prune)  # visualize cross-validation result
#R2 is 1-Rel Error
# R2      
1-surf.fit3$cptable[8,3]   #column 3, row 9 in cptable --> this gives you the adjusted R2 of the model. 


##Temporal BD and land cover
temporalBD.landuse<- as.data.frame(cbind(surf.genus.abundnon.label$SITE_ID, surf.genus.abundnon.label$WSA_ECOREGION, surf.genus.abundnon.label$LONG_DD, surf.genus.abundnon.label$LAT_DD, nla.surf$ORIGIN, nla.surf$LAKE_TYPE, nla.surf$PCT_DEVELOPED_BSN, nla.surf$PCT_FOREST_BSN, nla.surf$PCT_SHRUBLAND_BSN, nla.surf$PCT_GRASS_BSN, nla.surf$PCT_AGRIC_BSN, nla.surf$PCT_WETLAND_BSN, temporalBD.genus[,10:16]))
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

temporalland.fit1<- rpart(Total_BD~PCT_DEVELOPED_BSN + PCT_FOREST_BSN + PCT_SHRUBLAND_BSN + PCT_GRASS_BSN + PCT_AGRIC_BSN + PCT_WETLAND_BSN, method="anova", data=temporalBD.landuse) 
printcp(temporalland.fit1) #Variables actually used in construction: AGRIC, DEVELOPED, GRASS, FOREST, SHRUBLAND, WETLAND
summary(temporalland.fit1)
temporal.land1prune<-  prune(temporalland.fit1, cp=temporalland.fit1$cptable[which.min(temporalland.fit1$cptable[,"xerror"]),"CP"])
summary(temporal.land1prune)
#R2 is 1-Rel Error
# R2      
1-temporalland.fit1$cptable[10,3]  


##BASELGA ANALYSES

#Surface
#beta.pair (can do sorenson)
betapart.quant.sorenson<- beta.pair(surf.genus.abundnon[,5:101], index.family = "sorensen") #don't match abundance results with beta.div

#try pa
x<- decostand(surf.genus.abundnon[,5:101], "pa")
betapart.qual.sorenson<- beta.pair(x, index.family = "sorensen")
summary(betapart.qual.sorenson$beta.sor)


#Historical
y<- decostand(hist.genus.abundnon[,5:110], "pa")
betapart.qualhist.sorenson<- beta.pair(y, index.family = "sorensen")
summary(betapart.qualhist.sorenson$beta.sor)















##GENUS LEVEL FIGURES FOR RESUBMISSION
##CORRECTION 1- Re-do Figure 2(a) so that does not have double legend --> Decided to just make gray-scale.
#(From Script version V5)

all.states <- map_data("usa")
us.map<- ggplot()
us.map<- us.map + geom_polygon(data=all.states, aes(x=long, y=lat, group=group), colour="black", fill="grey92")

us.map
windows()
us.map.ecoregiong<- us.map + geom_point(data=nla.surf, aes(x=LONG_DD, y=LAT_DD, shape=WSA_ECOREGION), size=4)
us.map.ecoregiong<- us.map.ecoregiong + theme_bw() + labs(x= "Longitude (deg)", y="Latitude (deg)")
us.map.ecoregiong<- us.map.ecoregiong + theme(axis.text.x = element_text(colour="black", size=16))
us.map.ecoregiong<- us.map.ecoregiong + theme(axis.text.y = element_text(colour="black", size=16))
us.map.ecoregiong<- us.map.ecoregiong + theme(axis.title.x = element_text(size = rel(2), angle=00))
us.map.ecoregiong<- us.map.ecoregiong + theme(axis.title.y = element_text(size = rel(2), angle=90))
us.map.ecoregiong<- us.map.ecoregiong + scale_shape_manual(values = c("CPL" = 18, "NAP" = 16, "SPL" = 8, "TPL" = 17, "UMW" = 3, "WMT" = 15, "XER" = 7))
us.map.ecoregiong<- us.map.ecoregiong + theme(legend.title=element_text(size=16))
us.map.ecoregiong<- us.map.ecoregiong + theme(legend.text=element_text(size=16))
us.map.ecoregiong<- us.map.ecoregiong + annotate("text", x=-120, y=52, label="(a)", size=10)


#(B) Lake area 
lakearea.bp<- ggplot(nla.surf, aes(x=WSA_ECOREGION, y=LAKEAREA_KM2)) + geom_boxplot()
lakearea.bp<- lakearea.bp + theme_bw() + xlab("Ecoregion")
lakearea.bp<- lakearea.bp + theme_bw() + ylab(expression(paste("Area " (km^{2}))))
lakearea.bp<- lakearea.bp + theme(axis.text.x = element_text(colour="black", size=16))
lakearea.bp<- lakearea.bp + theme(axis.text.y = element_text(colour="black", size=16))
lakearea.bp<- lakearea.bp + theme(axis.title.x = element_text(size = rel(2), angle=00))
lakearea.bp<- lakearea.bp + theme(axis.title.y = element_text(size = rel(2), angle=90))
lakearea.bp<- lakearea.bp + scale_fill_brewer(type="qual", palette="Dark2", "Ecoregion (WSA)", breaks=c("CPL", "NAP", "NPL", "SPL", "TPL", "UMW", "WMT", "XER"), labels=c("Coastal Plains", "N. Appalachians", "N. Plains", "S. Plains", "Temperate Plains", "Upper MidWest", "W. Mountains", "Xeric"))  
lakearea.bp<- lakearea.bp + theme(legend.position="none")
lakearea.bp<- lakearea.bp + annotate("text", x=1, y=120, label="(b)", size=10)

#(C) Boxplots of Z_max
zmax.bp<- ggplot(nla.surf, aes(x=WSA_ECOREGION, y=Z_MAX)) + geom_boxplot()
zmax.bp<- zmax.bp + theme_bw() + labs(x= "Ecoregion", y="Maximum depth (m)")
zmax.bp<- zmax.bp + theme(axis.text.x = element_text(colour="black", size=16))
zmax.bp<- zmax.bp + theme(axis.text.y = element_text(colour="black", size=16))
zmax.bp<- zmax.bp + theme(axis.title.x = element_text(size = rel(2), angle=00))
zmax.bp<- zmax.bp + theme(axis.title.y = element_text(size = rel(2), angle=90))
zmax.bp<- zmax.bp + scale_fill_brewer(type="qual", palette="Dark2", "Ecoregion (WSA)", breaks=c("CPL", "NAP", "NPL", "SPL", "TPL", "UMW", "WMT", "XER"), labels=c("Coastal Plains", "N. Appalachians", "N. Plains", "S. Plains", "Temperate Plains", "Upper MidWest", "W. Mountains", "Xeric"))  
zmax.bp<- zmax.bp + theme(legend.position="none")
zmax.bp<- zmax.bp + annotate("text", x=1, y=60, label="(c)", size=10)

#(D) Boxplots of pH
ph.bp<- ggplot(waterqual.dat.red, aes(x=WSA_ECOREGION, y=PH_FIELD)) + geom_boxplot()
ph.bp<- ph.bp + theme_bw() + labs(x= "Ecoregion", y="pH")
ph.bp<- ph.bp + theme(axis.text.x = element_text(colour="black", size=16))
ph.bp<- ph.bp + theme(axis.text.y = element_text(colour="black", size=16))
ph.bp<- ph.bp + theme(axis.title.x = element_text(size = rel(2), angle=00))
ph.bp<- ph.bp + theme(axis.title.y = element_text(size = rel(2), angle=90))
ph.bp<- ph.bp + scale_fill_brewer(type="qual", palette="Dark2", "Ecoregion (WSA)", breaks=c("CPL", "NAP", "NPL", "SPL", "TPL", "UMW", "WMT", "XER"), labels=c("Coastal Plains", "N. Appalachians", "N. Plains", "S. Plains", "Temperate Plains", "Upper MidWest", "W. Mountains", "Xeric"))  
ph.bp<- ph.bp + theme(legend.position="none")
ph.bp<- ph.bp + annotate("text", x=1, y=10, label="(d)", size=10)

#(E) Boxplots of TP
tp.bp<- ggplot(waterqual.dat.red, aes(x=WSA_ECOREGION, y=PTL)) + geom_boxplot()
tp.bp<- tp.bp + theme_bw() + labs(x= "Ecoregion", y="Total phosphorus (ug/L)")
tp.bp<- tp.bp + theme(axis.text.x = element_text(colour="black", size=16))
tp.bp<- tp.bp + theme(axis.text.y = element_text(colour="black", size=16))
tp.bp<- tp.bp + theme(axis.title.x = element_text(size = rel(2), angle=00))
tp.bp<- tp.bp + theme(axis.title.y = element_text(size = rel(2), angle=90))
tp.bp<- tp.bp + scale_fill_brewer(type="qual", palette="Dark2", "Ecoregion (WSA)", breaks=c("CPL", "NAP", "NPL", "SPL", "TPL", "UMW", "WMT", "XER"), labels=c("Coastal Plains", "N. Appalachians", "N. Plains", "S. Plains", "Temperate Plains", "Upper MidWest", "W. Mountains", "Xeric"))  
tp.bp<- tp.bp + theme(legend.position="none")
tp.bp<- tp.bp + annotate("text", x=1, y=1500, label="(e)", size=10)

#Combo of B-E
windows()
be.panel<- grid.arrange(lakearea.bp, zmax.bp, ph.bp, tp.bp, nrow=2) #export 


#Figure 3
#GENUS LEVEL HISTORICAL
#Use percentage difference 
hist.beta.divB <- beta.div(hist.genus.abundnon[,5:110], method="percentagedifference", sqrt.D=FALSE, samp=TRUE, nperm=9999, save.D=FALSE, clock=FALSE)
summary(hist.beta.divB)
hist.beta.divB$SStotal_BDtotal 
hist.beta.divB$SCBD #Null 
hist.beta.divB$LCBD
hist.beta.divB$p.LCBD

#MUltiple testing correction
hist.p.adj<- p.adjust(hist.beta.divB$p.LCBD, method="holm", n=length(hist.beta.divB$p.LCBD)) #5 less than 0.05
summary(hist.p.adj<0.05)

#Make a dataframe with LCBD information for surface data using percentage difference
hist.LCBD.valuesB<- as.data.frame(hist.beta.divB$LCBD)
hist.LCBD.p<- as.data.frame(hist.beta.divB$p.LCBD)
hist.p.adj<- as.data.frame(hist.p.adj)
hist.LCBD.summary.adj<- as.data.frame(cbind(hist.genus.abundnon$SITE_ID, hist.genus.abundnon$WSA_ECOREGION, hist.genus.abundnon$LONG_DD, hist.genus.abundnon$LAT_DD, hist.LCBD.values, hist.LCBD.p, hist.p.adj))
colnames(hist.LCBD.summary.adj)[1]<- 'SITE_ID'
colnames(hist.LCBD.summary.adj)[2]<- 'WSA_ECOREGION'
colnames(hist.LCBD.summary.adj)[3]<- 'LONG_DD' 
colnames(hist.LCBD.summary.adj)[4]<- 'LAT_DD'
colnames(hist.LCBD.summary.adj)[5]<- 'LCBD'
colnames(hist.LCBD.summary.adj)[6]<- 'LCBD.P'
colnames(hist.LCBD.summary.adj)[7]<- 'LCBD.P.ADJ'
#verify, but looks like no significant LCBD after correcting for multiple testing. 

#GENUS LEVEL MODERN
#Use percentage difference 
surf.beta.divB <- beta.div(surf.genus.abundnon[,5:99], method="percentagedifference", sqrt.D=FALSE, samp=TRUE, nperm=9999, save.D=FALSE, clock=FALSE)
summary(surf.beta.divB)
surf.beta.divB$SStotal_BDtotal 
surf.beta.divB$SCBD 
surf.beta.divB$LCBD
surf.beta.divB$p.LCBD

#MUltiple testing correction
surf.p.adj<- p.adjust(surf.beta.divB$p.LCBD, method="holm", n=length(surf.beta.divB$p.LCBD)) # 1 less than 0.05
summary(surf.p.adj<0.05)

#Make a dataframe with LCBD information for surface data using percentage difference
surf.LCBD.valuesB<- as.data.frame(surf.beta.divB$LCBD)
surf.LCBD.p<- as.data.frame(surf.beta.divB$p.LCBD)
surf.p.adj<- as.data.frame(surf.p.adj)
surf.LCBD.summary.adj<- as.data.frame(cbind(surf.genus.abundnon$SITE_ID, surf.genus.abundnon$WSA_ECOREGION, surf.genus.abundnon$LONG_DD, surf.genus.abundnon$LAT_DD, surf.LCBD.values, surf.LCBD.p, surf.p.adj))
colnames(surf.LCBD.summary.adj)[1]<- 'SITE_ID'
colnames(surf.LCBD.summary.adj)[2]<- 'WSA_ECOREGION'
colnames(surf.LCBD.summary.adj)[3]<- 'LONG_DD' 
colnames(surf.LCBD.summary.adj)[4]<- 'LAT_DD'
colnames(surf.LCBD.summary.adj)[5]<- 'LCBD'
colnames(surf.LCBD.summary.adj)[6]<- 'LCBD.P'
colnames(surf.LCBD.summary.adj)[7]<- 'LCBD.P.ADJ'


#LCBD FIGURE
#GENUS LEVEL HISTORICAL LCBD- ADJUSTED P
windows()
us.map.histLCBD.adj<- us.map + geom_point(data=hist.LCBD.summary.adj, aes(x=LONG_DD, y=LAT_DD, size=LCBD, colour=LCBD.P.ADJ<0.05, shape=LCBD.P.ADJ<0.05)) #could change this to LCBD.P (but could show adjusted if keep this way)
us.map.histLCBD.adj<- us.map.histLCBD.adj + scale_size(name="LCBD") 
us.map.histLCBD.adj<- us.map.histLCBD.adj + scale_colour_manual(values=(c("blue", "red")))
us.map.histLCBD.adj<- us.map.histLCBD.adj + labs(x="Longitude (deg)", y="Latitude (deg)") 
us.map.histLCBD.adj<- us.map.histLCBD.adj + theme_bw()
us.map.histLCBD.adj<- us.map.histLCBD.adj + scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1))
us.map.histLCBD.adj<- us.map.histLCBD.adj + theme(axis.text.x = element_text(colour="black", size=16))
us.map.histLCBD.adj<- us.map.histLCBD.adj + theme(axis.text.y = element_text(colour="black", size=16))
us.map.histLCBD.adj<- us.map.histLCBD.adj + theme(axis.title.x = element_text(size = rel(2), angle=00))
us.map.histLCBD.adj<- us.map.histLCBD.adj + theme(axis.title.y = element_text(size = rel(2), angle=90))
us.map.histLCBD.adj<- us.map.histLCBD.adj + theme(legend.title=element_text(size=16))
us.map.histLCBD.adj<- us.map.histLCBD.adj + theme(legend.text=element_text(size=16))
us.map.histLCBD.adj<- us.map.histLCBD.adj + annotate("text", x=-110, y=51, label="(a) Historical LCBD", size=10) #for using in a panel with other LCBD map
us.map.histLCBD.adj<- us.map.histLCBD.adj + theme(legend.position="bottom")

windows()
us.map.histLCBD.adj


#GENUS LEVEL 2007 LCBD
us.map.surfLCBD.adj<- us.map + geom_point(data=surf.LCBD.summary.adj, aes(x=LONG_DD, y=LAT_DD, size=LCBD, colour=LCBD.P.ADJ<0.05, shape=LCBD.P.ADJ<0.05)) #could change this to LCBD.P (but could show adjusted if keep this way)
us.map.surfLCBD.adj<- us.map.surfLCBD.adj + scale_size(name="LCBD") 
us.map.surfLCBD.adj<- us.map.surfLCBD.adj + scale_colour_manual(values=(c("blue", "red")))
us.map.surfLCBD.adj<- us.map.surfLCBD.adj + labs(x="Longitude (deg)", y="Latitude (deg)") 
us.map.surfLCBD.adj<- us.map.surfLCBD.adj + theme_bw()
us.map.surfLCBD.adj<- us.map.surfLCBD.adj + scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1))
us.map.surfLCBD.adj<- us.map.surfLCBD.adj + theme(axis.text.x = element_text(colour="black", size=16))
us.map.surfLCBD.adj<- us.map.surfLCBD.adj + theme(axis.text.y = element_text(colour="black", size=16))
us.map.surfLCBD.adj<- us.map.surfLCBD.adj + theme(axis.title.x = element_text(size = rel(2), angle=00))
us.map.surfLCBD.adj<- us.map.surfLCBD.adj + theme(axis.title.y = element_text(size = rel(2), angle=90))
us.map.surfLCBD.adj<- us.map.surfLCBD.adj + theme(legend.title=element_text(size=16))
us.map.surfLCBD.adj<- us.map.surfLCBD.adj + theme(legend.text=element_text(size=16))
us.map.surfLCBD.adj<- us.map.surfLCBD.adj + annotate("text", x=-110, y=51, label="(b) Modern LCBD", size=10) #for using in a panel with other LCBD map
us.map.surfLCBD.adj<- us.map.surfLCBD.adj + theme(legend.position="bottom")

Windows()
us.map.surfLCBD.adj


windows()
fig3ab.adj<- grid.arrange(us.map.histLCBD.adj, us.map.surfLCBD.adj, nrow=1)


#GENUS LEVEL TBI - no sig
windows()
us.map.temp.genus<- us.map + geom_point(data=genus.TBI.output, aes(x=LONG_DD, y=LAT_DD, size=Total_BD))
us.map.temp.genus<- us.map.temp.genus + scale_size(name="TBI") 
us.map.temp.genus<- us.map.temp.genus + labs(x="Longitude (deg)", y="") 
us.map.temp.genus<- us.map.temp.genus + theme_bw()
us.map.temp.genus<- us.map.temp.genus + theme(axis.text.x = element_text(colour="black", size=16))
us.map.temp.genus<- us.map.temp.genus + theme(axis.text.y = element_text(colour="black", size=16))
us.map.temp.genus<- us.map.temp.genus + theme(axis.title.x = element_text(size = rel(2), angle=00))
us.map.temp.genus<- us.map.temp.genus + theme(axis.title.y = element_text(size = rel(2), angle=90))
us.map.temp.genus<- us.map.temp.genus + theme(legend.title=element_text(size=16))
us.map.temp.genus<- us.map.temp.genus + theme(legend.text=element_text(size=16))
us.map.temp.genus<- us.map.temp.genus + annotate("text", x=-118, y=51, label="(c) TBI", size=9) #for using in a panel with other LCBD maps
us.map.temp.genus<- us.map.temp.genus + theme(legend.position="bottom")


#Figure 4
plot(as.party(surf.fit2prune), tp_args = list(id = FALSE)) #Fig 4a

windows()
plot(as.party(surf.fit3prune), tp_args = list(id = FALSE))  #Fig 4b

plot(as.party(temporal.land1prune), tp_args = list(id = FALSE)) #Fig 5

develo#################COME BACK TO SPECIES LEVEL ANALYSES####################################
#SPECIES LEVEL- SURFACE
#Bind the non-transformed species columns to site, ecoregion and lat/long columns. 
surf.sp.abundnon<- as.data.frame(cbind(nla.surf.charles[,1], nla.surf.charles[,11], nla.surf.charles[,4], nla.surf.charles[,5], surf.sp.red[,4:930]))
colnames(surf.sp.abundnon)[1]<- 'SITE_ID'
colnames(surf.sp.abundnon)[2]<- 'WSA_ECOREGION'
colnames(surf.sp.abundnon)[3]<- 'LONG_DD'
colnames(surf.sp.abundnon)[4]<- 'LAT_DD'

#Use percentage difference 
surf.beta.sp <- beta.div(surf.sp.abundnon[,4:930], method="percentagedifference", sqrt.D=FALSE, samp=TRUE, nperm=999, save.D=FALSE, clock=FALSE)
summary(surf.beta.sp)
surf.beta.sp$SStotal_BDtotal #SS total =16.4 , BD total = 0.28
surf.beta.sp$SCBD 
surf.beta.sp$LCBD
surf.beta.sp$p.LCBD

#Make a dataframe with LCBD information for surface data using percentage difference
surf.LCBD.values.sp<- as.data.frame(surf.beta.sp$LCBD)
surf.LCBD.p.sp<- as.data.frame(surf.beta.sp$p.LCBD)
surf.LCBD.summary.sp<- as.data.frame(cbind(surf.sp.abundnon$SITE_ID, surf.sp.abundnon$WSA_ECOREGION, surf.sp.abundnon$LONG_DD, surf.sp.abundnon$LAT_DD, surf.LCBD.values.sp, surf.LCBD.p.sp))
colnames(surf.LCBD.summary.sp)[1]<- 'SITE_ID'
colnames(surf.LCBD.summary.sp)[2]<- 'WSA_ECOREGION'
colnames(surf.LCBD.summary.sp)[3]<- 'LONG_DD' 
colnames(surf.LCBD.summary.sp)[4]<- 'LAT_DD'
colnames(surf.LCBD.summary.sp)[5]<- 'LCBD'
colnames(surf.LCBD.summary.sp)[6]<- 'LCBD.P'

#Look at sites with significant (p<0.05)
surf.LCBD.summary.sp.sig<- as.data.frame(subset(surf.LCBD.summary.sp, LCBD.P<0.05)) #14 sites

#SPECIES LEVEL- HISTORICAL 
#####################################################################
