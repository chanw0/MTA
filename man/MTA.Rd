\name{MTA}
\alias{MTA}
\title{The overall function in MTA framework}
\description{

This overall function performs three tasks: 1) capturing the common and major community level microbial dynamic trends for a given group of subjects and identifying the dominant taxa; 2) examining whether the microbial overall dynamic trends are significantly different or not in groups; 3) classifying a given subject based on its longitudinal microbial profiling.}


\usage{
MTA(Y1,Y_new=NULL,phy.tree=NULL,timevec=NULL, transf="None",
             M=NULL,proportion.explained=0.85,
             k=5,lambda1.set=c(0.01,0.1,1,10),
             lambda2.set=seq(0,1,0.001),lambda3.set=c(0),num.sam=10,alpha=0.05)
}

\arguments{

\item{Y1}{ Y1 can be an N*P*T numeric array containing compositional microbiome data for N subjects having P taxa across T time points, or a list, in which each element represents a group of subjects having P taxa across T time points. The number of groups is equal to or larger than 2. The list's name can be assigned as the group labels, otherwise, it is 1:length(Y1) by default. If Y1 is an N*P*T numeric array, MTA gives the microbial dynamics for Y1 and identifies the key taxa. If Y1 is a list, MTA provides the group comparison or classification.}

\item{Y_new}{
An optional N*P*T numeric array containing compositional microbiome data for N subjects having P taxa across T time points. Y_new represents the subjects who need to be classified based on the training set Y1. Default=NULL, no calculation for the classification part.
}

\item{phy.tree}{
  The phylogenetic tree (phylo-class) object for P taxa. Default=NULL, no Laplacian penalty dealing with the phylogenetical correlation among taxa. }

\item{timevec}{
An optional numeric vector of the continuous time point with length = T.
}

\item{transf}{
  A character value, which represents the transformation technique one wants to take. "None" represents MTA analyzes the microbial relative abundance directly; "Log", "ALR", and "CLR" represent the corresponding log, additive log-ratio transformation and centered log-ratio transformation respectively. Default="None".}

\item{M}{
The number of top common trends one wants to select. Default=NULL.
}

\item{proportion.explained}{
A numeric value used to determine the number of top common trends as an alternative. Default=0.85.
}

\item{k}{
A numeric value represents k-fold cross-validation to determine the tuning parameters. Default=5.
}


\item{lambda1.set}{
A numeric vector containing the candidate tunning parameters for the smoothing technique to deal with the smoothness of the extracted trends. Default=c(0.01,0.1,1,10).
}

\item{lambda2.set}{
A numeric vector containing the candidate tunning parameters for the Lasso penalty to deal with the high dimensionality of taxa. Default=seq(0,1,0.001).
}

\item{lambda3.set}{
A numeric vector containing the candidate tunning parameters for the Laplacian penalty to deal with the phylogenetical correlation. Default=c(0).
}

\item{num.sam}{
An integer value, the number of resampling for the group comparison. Default=10.
}


\item{alpha}{
An numeric value, the signficance level for the group comparison. Default=0.05.
}
}


\value{

When Y1 is an N*P*T numeric array, MTA performs the first task: capturing the common and major community level microbial dynamic trends and identifying the dominant taxa based on the Y1. A list which contains two elements:

  \describe{

  \item{Trend Plot}{The extracted common trends extracted from all subjects Y1}
  \item{Factor Score}{The estimated factor scores for all key taxa which contribute to the extracted common trends.}
}

When Y1 is a list involving two or more than two groups, but Y_new is NULL, MTA performs the second task:  examining whether the microbial overall dynamic trends are significantly different or not in groups and provides a list which contains comparing results between any two groups. The comparing results between any two groups contains five elements:

  \describe{

   \item{P value}{P value for the common trends representing the difference in the microbial dynamic pattern between these two groups}
  \item{Trends}{The extracted common trends among R resamplings}
  \item{Discover rate}{The discover rates for all taxa among R resamplings. This item is provided when the common trend is significant, otherwise, on this item.}
    \item{Factor Score}{The estimated factor scores for all taxa which contribute to the common trends. This item is provided when the common trend is significant, otherwise, on this item.}
    \item{Standard error}{The standard error of the estimated factor scores for all taxa among R resamplings. This item is provided when the common trend is significant, otherwise, on this item.}
}

When Y1 is a list involving two or more than two groups, and Y_new is an N*P*T numeric array, MTA performs the third task:  classifying the subjects in Y_new based on the training set Y1. The result contains an N*2 matrix, in which the first column is the subject's ID and the second column is its group label accordingly.
}

\author{
Chan Wang, Jiyuan Hu, Martin J. Blaser, Huilin Li.

Maintainer: Chan Wang <Chan.Wang@nyulangone.org>, Huilin Li <Huilin.Li@nyulangone.org>
}
\references{
Wang C, Hu J, Blaser M J, Li H (2021). Microbial Trend Analysis for Common Dynamic Trend, Group Comparison, and Classification in Longitudinal Microbiome Study.

}

\examples{

########### generation data
library(dirmult)
## sample size
N=20;
## time point
T=10;
##number of taxa
P=20;

beta0=abs(rnorm(P,10,20))
## key taxa
causal.ind=sample(1:P,3)

#### generating two groups
group1=array(NA, dim = c(N,P,T))
for(t in 1:T)  group1[, ,t] = rdirichlet(N,beta0)

group2=array(NA, dim = c(N,P,T))
for(t in 1:T){
  betaT=rep(0,P)
  betaT[causal.ind]=c(-1,2, -1)/2*min(beta0[causal.ind])*sin(t)
  group2[, ,t] = rdirichlet(N,(beta0+betaT))}


########## MTA captures the common trends shared by all subjects in group1
### extracting the common trend from a group of subjects with explained variation >=85%
MTA(Y1=group1)

### extracting the common trend from a group of subjects with 2 common trends
MTA(Y1=group1,M=2)


######## comparing between group1 and group2
MTA(list(group1,group2))


######## classifying subjects based on the training set
MTA(list(group1,group2),group2[1:5,,])

}

\keyword{MTA}
