#' Toy data used for illustration
#'
#' Toydata is a list with the length=3.
#'
#' @format A list with the length=3 with names Case, Control, and phy.tree:
#' \describe{
#'   \item{Case}{A 15 * 35 * 4 Array of relative abundances for 15 Cases}
#'   \item{Control}{A 15 * 35 * 4 Array of relative abundances for 15 Controls}
#'   \item{phy.tree}{The phylogenetic tree (phylo-class) object for 35 taxa}
#' }
#'
#' @examples
#'####### Toy data for illustration
#'library(ggpubr)
#'library(toOrdinal)
#'library(reshape2)
#'
#'Case=Toydata[[1]]
#'Control=Toydata[[2]]
#'phy.tree=Toydata[[3]]
#'
#' ### Group comparison
#'
#' aa=MTA(list(Control,Case),phy.tree=phy.tree)$`Comparison between group 1 and group 2`
#'
#'
#' aa1=aa[[4]]
#' aa2=aa[[5]]
#'
#'for(i in 1:nrow(aa2)) {
#'
#'  bb1=aa1[i,-1];bb2=aa2[i,-1]
#'  bb=order(bb1^2/sum(bb1^2),decreasing =TRUE)
#'  bb1[-c(bb[1:which(cumsum(cumsum(sort(bb1^2/sum(bb1^2),decreasing =TRUE))>0.99)==1)])]=0
#'
#'  bb1[((bb1+1.96*bb2)>=0 & (bb1-1.96*bb2)<=0)]=NA
#'
#'  aa1[i,-1]=bb1
#'}
#'effect1=aa1
#'colnames(effect1)=c("trend",unlist(dimnames(Control)[2]))
#'effect1=data.frame(effect1,check.names =FALSE)
#'effect1=melt(effect1,id.vars = 1,variable.name = "taxa",value.name = "score")
#'effect1$trend=paste(sapply(effect1$trend, toOrdinal), "common trend")
#'
#'effect1=na.omit(effect1)
#'
#'pp1=ggbarplot(effect1, x = "taxa", y = "score",
#'              ylab = "Estimated factor score", xlab="", facet.by = "trend",
#'              add = "mean")+
#'  facet_wrap(~trend,scales="free")+theme(axis.text.x = element_text(angle = 90, hjust = 1))
#'
#'######### trend figure
#'
#'trend1=aa[[2]]
#'
#'pp2=ggline(trend1, x = "time",y = "Escore",
#'           ylab = "Common trend", xlab="Week", facet.by = "trend",
#'           add = "mean")+
#'  facet_wrap(~trend,scales="free_y")
#'
#'
#'ggarrange(pp2,pp1,labels = c("(A)", "(B)"), ncol = 2)
#'
#'#### Classification
#'aa=MTA(list(Control,Case),Case, phy.tree=phy.tree)
#'
#'### Group label 1: the first group (Control); 2: the second group (Case)
#'print(aa)
#'
"Toydata"




