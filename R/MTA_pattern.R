############Extracted the dynamic trends from a group of subjects

MTA_pattern=function(x,M,proportion.explained, k,Laplacian.matrix,timevec,lambda1.set,
               lambda2.set,lambda3.set)
{

  N=dim(x)[1]; P=dim(x)[2];T=dim(x)[3];

  if(is.null(timevec)) timevec=timevec=1:T

  BS = create.bspline.basis(timevec, norder=4)
  B = getbasismatrix(timevec, BS)  #basis function
  Omega = getbasispenalty(BS)


  if(is.null(M)) {

    ff=cc=NULL
    explained.variance=NULL
    xx=x
    for(i in 1:50)
    {

      predict.res=MTA01(xx,k,Laplacian.matrix,timevec,lambda1.set,lambda2.set,lambda3.set)
      c.updated=predict.res[[1]];f.updated=predict.res[[2]];

      cc=rbind(cc,c.updated);ff=cbind(ff,f.updated)
      predict.matrix=ff%*%cc%*%(t(B))

      explained.variance=sum(predict.matrix^2)*N/(sum(x^2))

      for(j in 1:N) xx[j,,]=xx[j,,]-predict.matrix

      if (explained.variance>proportion.explained) break;}} else {

        ff=cc=NULL
        xx=x
        for(i in 1:M)
        {

          predict.res=MTA01(xx,k,Laplacian.matrix,timevec,lambda1.set,lambda2.set,lambda3.set)
          c.updated=predict.res[[1]];f.updated=predict.res[[2]];

          cc=rbind(cc,c.updated);ff=cbind(ff,f.updated)
          predict.matrix=ff%*%cc%*%(t(B))

          for(j in 1:N) xx[j,,]=xx[j,,]-predict.matrix
        }}


  plot.data=data.frame(cbind(timevec,sapply(1:nrow(cc),function(x,cc,B){return(cc[x,]%*%t(B))},cc=cc,B=B)))
  colnames(plot.data)=c("time",paste(sapply(1:nrow(cc), toOrdinal), "common trend"))

  plot.data=melt(plot.data,id.vars = 1)

  pc.plot=ggplot(plot.data, aes(time, value))+geom_point()+
    # geom_smooth(se=FALSE)+
    geom_line()+theme_bw()+ylab("Microbial common trend")+
    facet_wrap(~variable,scales="free")

  AA=list(pc.plot,ff,cc)

  return(AA)

}
