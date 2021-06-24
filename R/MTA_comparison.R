############# Group comparison

MTA_comparison=function(Controls, Cases,k,Laplacian.matrix,timevec,lambda1.set,
                        lambda2.set,lambda3.set,num.sam,alpha)
{
  N1=dim(Controls)[1];N2=dim(Cases)[1];
  N=min(N1,N2)

  T=dim(Controls)[3]
  if(is.null(timevec)) timevec=1:T

  BS = create.bspline.basis(timevec, norder=4)
  B = getbasismatrix(timevec, BS)  #basis function
  Omega = getbasispenalty(BS)

  cc.list=ff.list=list()
  res=cc.reference=NULL

  for(rr in 1:num.sam)
  {
    aa1=sample(1:N2,floor(N/5*4));aa2=sample(1:N1,floor(N/5*4))
    #aa1=sample(1:N2,N);aa2=sample(1:N1,N)

    Z=Cases[aa1,,]-Controls[aa2,,]

    Z1=Z
    ff=cc=pro=NULL
    for(i in 1:100)
    {
      AA1=MTA01(Z1,k,Laplacian.matrix,timevec,lambda1.set,lambda2.set,lambda3.set)

      cc=rbind(cc,AA1[[1]]);ff=cbind(ff,AA1[[2]])
      Z.pred=ff%*%cc%*%(t(B))

      pro=c(pro,dim(Z)[1]*sum(Z.pred^2)/sum(Z^2))

      tem=diff(pro)
      if(length(tem)!=0) if(tem[length(tem)]<=0) break;

      for(nn in 1:dim(Z1)[1]) Z1[nn,,]=Z1[nn,,]-Z.pred
    }

    cc=cc[1:length(tem),]; ff=ff[,1:length(tem)]

    if(is.null(dim(cc))) {cc=t(cc);ff=as.matrix(ff)}

    cc.list[[rr]]=cc;ff.list[[rr]]=ff

    aa1=sample(1:N1,floor(N1/2))
    aa2=setdiff(1:N1,aa1); aa3=min(length(aa1),length(aa2))

    Z=Controls[aa1[1:aa3],,]-Controls[aa2[1:aa3],,]
    AA2=MTA01(Z,k,Laplacian.matrix,timevec,lambda1.set,lambda2.set=seq(0,0.1,0.005),lambda3.set)
    ##lambda2.set=seq(0,0.2,0.001)

    cc.reference=rbind(cc.reference,AA2[[1]])
  }

  ff.tem=matrix(unlist(ff.list),ncol=dim(Z)[2],byrow = TRUE)

  if(min(apply(matrix(unlist(ff.list),nrow=num.sam,byrow = TRUE),1,function(x) max(abs(x))))!=0) {

  num.tpc=max(sapply(1:num.sam, function(x,cc.list) dim(cc.list[[x]])[1], cc.list=cc.list))

  sin.re=matrix(NA,nrow=num.sam,ncol=(num.tpc+1))
  for(i in 1:num.sam) { tem.va=rowSums((cc.list[[i]]%*%(t(B)))^2);
  sin.re[i,1:length(tem.va)]=tem.va;
  sin.re[i,(num.tpc+1)]=sum((cc.reference[i,]%*%(t(B)))^2)}

  pvalue.tpc=NULL
  for(j in 1:num.tpc)  if(length(na.omit(sin.re[,j]))>2) pvalue.tpc=c(pvalue.tpc,wilcox.test(sin.re[,j],sin.re[,(num.tpc+1)],alternative = c("greater"))$p.value)

  ########## Multiple comparison correct: Bonferroni correction
  #index.sig=which(pvalue.tpc<=alpha/length(pvalue.tpc))

  ### gatekeeping procedure
  pvalue.tpc=c(0,pvalue.tpc)
  for(j in 2:length(pvalue.tpc)) pvalue.tpc[j]=max(pvalue.tpc[j-1],pvalue.tpc[j])
  pvalue.tpc=pvalue.tpc[-1]
  index.sig=which(pvalue.tpc<=alpha)



  if (length(index.sig)!=0) {

    Alltrend=NULL
    for(j in index.sig) {
      common.trend=sapply(1:num.sam, function(x,cc.list) cc.list[[x]][j,]%*%(t(B)),cc.list=cc.list)
      sign.ind=which(sign(common.trend[1,])==-1)
      common.trend[,sign.ind]=common.trend[,sign.ind]*(-1)

      for(xx in sign.ind) ff.list[[xx]][,j]=ff.list[[xx]][,j]*(-1)

      common.trend=data.frame(common.trend,time=1:nrow(common.trend),trend=rep(j,nrow(common.trend)))
      Alltrend=rbind(Alltrend,common.trend) }

    plot.data=melt(Alltrend,id=c("time","trend"),variable.name ="factor",value.name = "Escore")

    plot.data$trend=sapply(plot.data$trend, toOrdinal)
    plot.data$trend=paste(plot.data$trend, "common trend")

    # qq=ggplot(plot.data, aes(x=time, y=Escore, group = factor))+
    #   geom_line() +geom_point()+theme_bw()+scale_x_continuous(breaks=unique(plot.data$time))+
    #   facet_wrap(~trend,scales="free_y")

    ########### Discovery rate

    AllDiscover=Effect=SE=NULL
    for(j in index.sig) {
      select.taxa=sapply(1:num.sam, function(x,ff.list) ff.list[[x]][,j],ff.list=ff.list)

      AllDiscover=rbind(AllDiscover, c(j,rowMeans(select.taxa!=0)))
      Effect=rbind(Effect,c(j,rowMeans(select.taxa)))
      SD=apply(select.taxa,1,sd); #print(SD)
      SE=rbind(SE,c(j,SD/sqrt(num.sam)))
    }

Allresults=list(pvalue.tpc,plot.data,AllDiscover,Effect,SE)

names(Allresults)=c("P value","Trends","Discover rate", "Factor score", "Standard error")

} else {Allresults=list(pvalue.tpc);
     names(Allresults)=c("P value")}} else Allresults="None of taxa were selected. The Lasso penalty is too strict, please select a more merciful range"

  return(Allresults) }
