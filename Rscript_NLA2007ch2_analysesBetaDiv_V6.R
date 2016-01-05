##################################################################################################
####NLA 2007 CHAPTER 2 PROJECT: DIATOM BETA-DIVERSITY ACROSS THE US (historical/surface seds)    #
##################################################################################################

##In this script: Beta-diversity focused analyses and beta-diversity partitioning 
#Beta-diversity focused analyses 
#Previous script: Rscript_NLA2007ch2_analysesBetaDiv_V5.R
#R version: 3.1.2 (Pumpkin Helmet)

##Last update: January 5, 2016 
##Associated workspace: workspace_NLA2007ch2_analysesBetaDiv_V6.RData 
##Associated markdown: 
##Associated .txt of R script: Rscript_NLA2007ch2_analysesBetaDiv_V6.txt
##Github: ##NEED TO POST. 

##NOTE: This script is a duplicate of V5 with the same data and analyses- except for a reduction
#of sites because of issues in agreement between taxonomists. As such, sites for spatial 
#analyses reduced further so that dataset reduced to only include sites identified by
#taxnomists in Don Charles' group. 

##################################################################################################

##################################################################################################
####WORKING DIRECTORY                                                                            #
##################################################################################################
setwd("C:/Users/Winegardner/Documents/MCGILL/PhD chapters and projects/EPA National lakes Assessment/Chapter 2/Final_analyses_data")

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
library(corrplot) #Correlation matrices 

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


##################################################################################################
####DATA SET-UP                                                                                  #
##################################################################################################

##All final data files placed into working directory (above) 

####Master abiotic/site data####
##General data
nla.data2<- read.csv("nla2007_lakes_topbotsamples_UpdateNov2014_V2.csv") #nla2007_lakes_topbotsamples_UpdateNov2014_V2.csv
#Additional path: McGill/PhD Chapters and projects/EPA National lakes assessment/Chapter 2/ NLA data
#498 rows (249 sites x 2)
#Only cores greater than 30cm (249, 186 of those High confidence)

##Subsetting of data required for matching/merging dataframes

##In nla.data2 --> retain only high confidence cores
nla.data2.highC<- subset(nla.data2, HIGH_C_CORE == "YES", drop=T) #372 rows (186 top, 186 bottom)

#In nla.data2.highC --> subset into surface and historical sediments 
nla.highC.surf<- subset(nla.data2.highC, SED_TYPE == "Surface") #Surface/top sediments #186 sites, 63 variables
nla.highC.hist<- subset(nla.data2.highC, SED_TYPE == "Historical") #Historical/bottom sediments #186 sites, 63 variables

##For each of these subsets- add in column with the 2012 lake type re-classification (adding onto end of the dataframes)
new.classes<- read.csv("nla2007_comparelakeclasses.csv") #nla2007_comparelakeclasses.csv
#Additional path: McGill/PhD Chapters and projects/EPA National lakes assessment/Chapter 2

#Order new.classes by nla.highC.surf and .hist
new.classes<-new.classes[match(nla.highC.surf[,1], new.classes[,1]),]

#Define rownames for these datasets.
rownames(nla.highC.surf)<- as.character(nla.highC.surf[,1])
rownames(nla.highC.hist)<- as.character(nla.highC.hist[,1])

#Add 'LAKE_TYPE_RECLASSED_NAME' from new.classes 
nla.highC.surf<- as.data.frame(cbind(nla.highC.surf, new.classes$LAKE_TYPE_RECLASSED_NAME)) #now 64 variables
colnames(nla.highC.surf)[64]<- 'LAKE_TYPE_RECLASSED'

nla.highC.hist<- as.data.frame(cbind(nla.highC.hist, new.classes$LAKE_TYPE_RECLASSED_NAME)) #now 64 variables
colnames(nla.highC.hist)[64]<- 'LAKE_TYPE_RECLASSED'

##Remove 7 additional sites that were flagged from radiometric dating from both nla.highC.surf and .hist 
#Remove: NLA06608-2162, NLA06608-0753, NLA06608-0050, NLA06608-1414, NLA06608-0805, NLA06608-NELP-1041, 
#NLA06608-ELS:1E1-096

##Also remove 3 sites that are classified as coastal plains but in Eastern Seabord. 
#Remove: NLA06608-0198, NLA06608-0710, NLA06608-ELS:1D1-035

remove.list<- c("NLA06608-2162", "NLA06608-0753", "NLA06608-0050", "NLA06608-1414", "NLA06608-0805", "NLA06608-NELP-1041", "NLA06608-ELS:1E1-096", "NLA06608-0198", "NLA06608-0710", "NLA06608-ELS:1D1-035")

#Remove from nla.highC.surf- now have 64 columns (variables) and 176 sites 
nla.highC.surf<- nla.highC.surf[ !(rownames(nla.highC.surf) %in% remove.list),]

#Remove from nla.highC.hist- now have 64 columns (variables) and 176 sites 
nla.highC.hist<- nla.highC.hist[ !(rownames(nla.highC.hist) %in% remove.list),]

##DATA FOR n = 176 ANALYSES (TEMPORAL)
nla.highC.surf
nla.highC.hist

##Remove sites that were not identified by Don Charles' group
charles.sites<- read.csv("nla2007_CharlesGroup_sites.csv") #NOTE: Sub in if get additional information from Don Charles. 

#Keep only sites in nla.highC.surf and nla.highC.hist that are also found in charles.sites. 
rownames(charles.sites)<- as.character(charles.sites[,1])

#Surface
#Remove sites from nla.highC.surf that aren't in charles.sites
surface.selected<- (nla.highC.surf$SITE_ID %in% charles.sites$SITE_ID)
nla.red.surf<- nla.highC.surf[surface.selected,] 
#So if the list of sediments I have is accurate, only 35 of the screened sites I have been working with would have
#been done by Don Charles' group. --> so spatial analyses would proceed with only 35 sites. 

#Historical
#Remove sites from nla.highC.hist that aren't in charles.sites
historical.selected<- (nla.highC.hist$SITE_ID %in% charles.sites$SITE_ID)
nla.red.hist<- nla.highC.hist[historical.selected,] 
#So if the list of sediments I have is accurate, only 35 of the screened sites I have been working with would have
#been done by Don Charles' group. --> so spatial analyses would proceed with only 35 sites. 

##DATA FOR n = 35 ANALYSES (SPATIAL)
nla.red.surf
nla.red.hist


##Easiest way might be to add a column to the nla.highC.surf and .hist datasets to desginate Don Charles group or not, 
#then can just subset later on. 

##SURFACE
#Add a "CHARLES_GROUP" column to nla.red.surf with "Yes" 
nla.red.surf$CHARLES_GROUP<- "YES"

#join surface dataframes so that CHARLES_GROUP column is the last, then replace NA's with NO
nla.surf<- join(nla.highC.surf, nla.red.surf, by= NULL, type = "left", match = "all") #"left" means use the rows in the first dataframe
nla.surf$CHARLES_GROUP #shows NAs and "YES" --> need to replace NAs with "NO"

nla.surf$CHARLES_GROUP[is.na(nla.surf$CHARLES_GROUP)]<- "NO"
as.factor(nla.surf$CHARLES_GROUP)
summary(nla.surf$CHARLES_GROUP)

#So nla.surf --> all n = 176 sites
#Subset for Charles group
nla.surf.charles<- as.data.frame(subset(nla.surf, CHARLES_GROUP == "YES", drop=T)) #35 rows- good. 


##HIST
#Add a "CHARLES_GROUP" column to nla.red.hist with "Yes"
nla.red.hist$CHARLES_GROUP<- "YES"

#join surface dataframes so that CHARLES_GROUP column is the last, then replace NA's with NO
nla.hist<- join(nla.highC.hist, nla.red.hist, by= NULL, type = "left", match = "all") 
nla.hist$CHARLES_GROUP #shows NAs and "YES" --> need to replace NAs with "NO"

nla.hist$CHARLES_GROUP[is.na(nla.hist$CHARLES_GROUP)]<- "NO"
as.factor(nla.hist$CHARLES_GROUP)
summary(nla.hist$CHARLES_GROUP)

#So nla.hist --> all n = 176 sites
#Subset for Charles group
nla.hist.charles<- as.data.frame(subset(nla.hist, CHARLES_GROUP == "YES", drop=T)) #35 rows- good. 

####NOTE: 
##So moving forward- site/abiotic data are called: 
nla.surf
nla.hist
#these dataframes should be used to manipulate all additional files --> and for other analyses, can be subsetted to
#only include the Charles Group sites if use the CHARLES_GROUP column. 


####Add water quality data#### 
##Extract water quality data columns: 
extradat2<- read.csv("Beaulieu&Taranu2014_NLAdataset.csv") #Beaulieu&Taranu2014_NLAdataset.csv
#File path: C:\Users\Winegardner\Documents\MCGILL\PhD chapters and projects\EPA National lakes Assessment\Chapter 2
#All visit 1

#Keep only SECMEAN, CHLA, PTL, NTL, MEAN_T, COND, PH_FIELD, ELEV_PT and Silica
waterqual.dat<- as.data.frame(subset(extradat2, select=c("SITE_ID", "SECMEAN", "CHLA", "PTL", "NTL", "MEAN_T", "COND", "PH_FIELD", "ELEV_PT", "SIO2")))

#Remove sites from waterqual.dat that aren't in nla.highC.surf (using surface as template)
waterqual.dat.selected<- (waterqual.dat$SITE_ID %in% nla.surf$SITE_ID)
waterqual.dat.red<- waterqual.dat[waterqual.dat.selected,] #ok, same sites(CHECK), but in different order than in nla.highC.surf
#Need to order so that will be able to copy env distance into a column. 

#Re-order so can check that same- yes. 
waterqual.dat.red<- waterqual.dat.red[ order(match(waterqual.dat.red$SITE_ID, nla.surf$SITE_ID)), ]#ok

#Add WSA_ECOREGION
waterqual.dat.red<- as.data.frame(cbind(waterqual.dat.red, nla.surf$WSA_ECOREGION))
colnames(waterqual.dat.red)[11]<- 'WSA_ECOREGION'


####Presence-absence diatoms#### 
##Diatom data
nla.diat<- read.csv("combo.red30highC.diatoms_matchedpa_CLEAN.csv") #combo.red30highC.diatoms_matchedpa_CLEAN.csv
#File path: McGill/PhD Chapters and projects/EPA National lakes assessment/Chapter 2/ NLA data/ Files for analysesNov2014Rscript

#In nla.diat --> subset into surface and historical sediments 
nla.diat.surf<- subset(nla.diat, SED_TYPE == "Surface") #Surface/top sediments #186 
rownames(nla.diat.surf)<- as.character(nla.diat.surf[,1])
nla.diat.hist<- subset(nla.diat, SED_TYPE == "Historical") #Historical/bottom sediments #186
rownames(nla.diat.hist)<- as.character(nla.diat.hist[,1])

##Surface diatoms 
#In nla.diat.surf --> retain only sites (rows) found in nla.highC.surf (now down to 179)
#Remove the same sites as removed from nla.highC.surf
nla.diat.surf<- nla.diat.surf[ !(rownames(nla.diat.surf) %in% remove.list),]
#Sites in nla.diat.surf that are found in nla.highC.surf
surf.check<- (nla.diat.surf$SITE_ID %in% nla.surf$SITE_ID) #All match- now n=176 for surface diatoms. 

#Add MY_SAMPLE_ID column to the diatom file
#Add WSA_ECOREGION column to the diatom file
#Add LONG_DD
#Add LAT_DD
#Add CHARLES_GROUP
nla.diat.surf<- as.data.frame(cbind(nla.diat.surf, nla.surf[,2], nla.surf[,11], nla.surf[,4], nla.surf[,5], nla.surf[,65]))
summary(nla.diat.surf) #last 5 columns are now "MY_SAMPLE ID", "WSA_ECOREGION", "LONG_DD", "LAT_DD", "CHARLES_GROUP" (but labelled as column numbers)
#Rename those last two columns as "MY_SAMPLE_ID", "WSA_ECOREGION", "LONG_DD", "LAT_DD", "CHARLES_GROUP"
colnames(nla.diat.surf)[1215]<- 'MY_SAMPLE_ID'
colnames(nla.diat.surf)[1216]<- 'WSA_ECOREGION'
colnames(nla.diat.surf)[1217]<- 'LONG_DD'
colnames(nla.diat.surf)[1218]<- 'LAT_DD'
colnames(nla.diat.surf)[1219]<- 'CHARLES_GROUP'

#Ok- breakdown for new nla.diat.surf: 
#[,1]: SITE_ID
#[,2]: SED_TYPE
#[,3:1214]: diatom species columns
#[,1215]: MY_SAMPLE_ID (unique identifier that can match to abiotic file as well)
#[,1216]: WSA_ECOREGION 
#[,1217]: LONG_DD
#[,1218]: LAT_DD
#[,1219]: CHARLES_GROUP
#So will be able to subsample the diatom files based on ecoregion later on. 

##Historical diatoms 
#In nla.diat.hist --> retain only sites (rows) found in nla.highC.hist
#Remove the same sites as removed from nla.highC.hist
nla.diat.hist<- nla.diat.hist[ !(rownames(nla.diat.hist) %in% remove.list),]
#Sites in nla.diat.surf that are found in nla.highC.surf
hist.check<- (nla.diat.hist$SITE_ID %in% nla.highC.hist$SITE_ID) #All match, now historical sediments with n=176

#Add MY_SAMPLE_ID column to the diatom file
#Add WSA_ECOREGION column to the diatom file
#Add LONG_DD
#Add LAT_DD
#Add CHARLES_GROUP
nla.diat.hist<- as.data.frame(cbind(nla.diat.hist, nla.hist[,2], nla.hist[,11], nla.hist[,4], nla.hist[,5], nla.hist[,65]))
summary(nla.diat.hist) #last 4 columns are now "MY_SAMPLE ID", "WSA_ECOREGION", "LONG_DD", "LAT_DD", "CHARLES_GROUP" (but labelled as column numbers)
#Rename those last two columns as "MY_SAMPLE_ID", "WSA_ECOREGION", "LONG_DD", "LAT_DD", "CHARLES_GROUP"
colnames(nla.diat.hist)[1215]<- 'MY_SAMPLE_ID'
colnames(nla.diat.hist)[1216]<- 'WSA_ECOREGION'
colnames(nla.diat.hist)[1217]<- 'LONG_DD'
colnames(nla.diat.hist)[1218]<- 'LAT_DD'
colnames(nla.diat.hist)[1219]<- 'CHARLES_GROUP'

#Ok- breakdown for new nla.diat.hist: 
#[,1]: SITE_ID
#[,2]: SED_TYPE
#[,3:1214]: diatom species columns
#[,1215]: MY_SAMPLE_ID (unique identifier that can match to abiotic file as well)
#[,1216]: WSA_ECOREGION 
#[,1217]: LONG_DD
#[,1218]: LAT_DD
#[,1219]: CHARLES_GROUP
#So will be able to subsample the diatom files based on ecoregion later on. 


####Quantitative diatoms#### 
##Surface diatoms 
surf.diat.abund<- read.csv("topv1.red30highC.diatoms_abund_CLEAN.csv") #topv1.red30highC.diatoms_abund_CLEAN
#Additional path: C:\Users\Winegardner\Documents\MCGILL\PhD chapters and projects\EPA National lakes Assessment\Chapter 2\NLA data\Files for analysesNov2014R script
rownames(surf.diat.abund)<- as.character(surf.diat.abund$SITE_ID)

#Remove sites so that down to 176
surf.diat.abund<- surf.diat.abund[ !(rownames(surf.diat.abund) %in% remove.list),]
#176 sites
#927 species 

##Historical diatoms 
hist.diat.abund<- read.csv("botv1.red30highC.diatoms_abund_CLEAN.csv") #botv1.red30highC.diatoms_abund_CLEAN
#Additional path: C:\Users\Winegardner\Documents\MCGILL\PhD chapters and projects\EPA National lakes Assessment\Chapter 2\NLA data\Files for analysesNov2014R script
rownames(hist.diat.abund)<- as.character(hist.diat.abund$SITE_ID)

#Remove sites so that down to 176
hist.diat.abund<- hist.diat.abund[ !(rownames(hist.diat.abund) %in% remove.list),]
#176 sites
#1002 species 

##ADD CHARLES_GROUP LATER IF NEEDED. 


####Radiometric dating data#### 
##Data as of April 18, 2015 (still waiting on a few samples)
dates.dat<- read.csv("nla_dating_binary.csv") #nla_dating_binary.csv
#Additional path: C:\Users\Winegardner\Documents\MCGILL\PhD chapters and projects\EPA National lakes Assessment\Chapter 2\Radiometric dating
#1 means sample is thought to be pre-1850, 0 means sample NOT pre-1850

####Final data available####
nla.surf #Site info/land use data- surface sediments (high C sites only)   (can be subsetted by CHARLES-GROUP == YES to reduce)
nla.hist #Site info/land use data- historical sediments (high C sites only) (can be subsetted by CHARLES-GROUP == YES to reduce)
waterqual.dat.red #Water qulaity data (contemporary) for the 186 sites 
km.matrix #Spatial distances bwteen sites 
waterqual.dist #Water-quality (env) distances between sites 
dist.data #Dataframe with both spatial and env distances between sites
nla.diat.surf #Presence-absence surface diatoms (can be subsetted by CHARLES-GROUP == YES to reduce)
nla.diat.hist #Presence-absence historical diatoms (can be subsetted by CHARLES-GROUP == YES to reduce)
surf.diat.abund #Quantitative surface diatoms
hist.diat.abund #Quantitative historical diatoms
dates.dat #Dating data with  binary assignments for pre-1850 based on radiometric dating 


####NOTE TO CONTINUE: DATA SET-UP NOW READY TO MODIFY ANALYSES BASED ON REDUCED NUMBER OF SITES FROM CHARLES' GROUP.
#As re-do spatial beta diversity analyses, edit to subset by CHARLES_GROUP == "YES". 
#If get an updated list from Don Charles- re-run script with this new list, then continue. 
