
library("purrr")
library("parallel")
library("spatstat")
library("dplyr")
library("ggplot2")
library("ggthemes")
library("RColorBrewer")
library("beepr")
library("reshape2")
library("rlist")

##############Dont forget to specify the host number, i.e. sync with previous data

###############################################

newvar<-function(mydata){
  dftemp<-pluck(mydata)
  temp1 <- filter(dftemp, infected <= (900/4))#<------sync here
}

temp<-map(dflist,newvar)

temptest<-temp[c(1:3)]
###############################################


lineareglist<-list()
g<-1
linearpluck<-function(x){
  dftemp<-pluck(x)
  for (i in unique(dftemp$sim)){
    lmoutput<-lm(formula=log(infected)~time,data=filter(dftemp,sim==i))
    lmoutput_int<-as.numeric(exp(lmoutput$coefficients[1]))
    lmoutput_r<-as.numeric(lmoutput$coefficients[2])
    lmoutput_cof<-as.numeric(summary(lmoutput)$r.squared)
  dfnew<-data.frame(int=lmoutput_int, rsquared=lmoutput_cof, r=lmoutput_r, beta = sample(dftemp$beta,1),
             theta = sample(dftemp$theta,1),landscape=sample(dftemp$landscape,1),sim=i)
  lineareglist[[g]]<-dfnew
  print(g)
  g<-g+1
  }
  return(lineareglist)
}
r_lnreg<-function(i=NULL){


rdata<-map(temp,linearpluck)

}
cl <- makeCluster(mc <- getOption("cl.cores", 20))
clusterCall(cl,function() library("purrr"))
clusterCall(cl,function() library("dplyr"))
clusterExport(cl=cl, varlist=c("temp","lineareglist","linearpluck","g"),envir = environment())
par_r1<-parLapply(1,fun=r_lnreg,cl=cl)
stopCluster(cl)
rdataset <- do.call("rbind", par_r1)

save(rdataset,file="rdatasetbeta7theta7.Rda")

##################################################################################################

newvar2<-function(mydata){
  rdatasettemp<-(pluck(mydata))
  rdatasettemp1<- unlist(lapply(rdatasettemp,`[`, c("r")))
  rdatamean<-mean(rdatasettemp1)
  rdatasd<-sd(rdatasettemp1)
  rdataframesdmean<-data.frame(r_mean=rdatamean,rdatasd=rdatasd,
                               beta=rdatasettemp[[1]]["beta"],
                               theta=rdatasettemp[[1]]["theta"],
                               landscape=rdatasettemp[[1]]["landscape"])
}


r_lnreg1<-function(i=NULL){

rdatamean<-map(rdataset,newvar2)
rdatameantable<-do.call(rbind.data.frame,rdatamean)

}
cl <- makeCluster(mc <- getOption("cl.cores", 20))
clusterCall(cl,function() library("purrr"))
clusterCall(cl,function() library("dplyr"))
clusterExport(cl=cl, varlist=c("rdataset","lineareglist","linearpluck","newvar2"),envir = environment())
par_r1<-parLapply(1,fun=r_lnreg1,cl=cl)
stopCluster(cl)
rdataset1 <- do.call("rbind", par_r1)




#rdatamean<-rdataset%>%group_by(beta, theta,landscape)%>%summarise_at(vars(r),list(r_mean = mean))
#library(tidyverse)
#tempLongPreds <- pivot_longer(temp, cols=c(infected,predExp,predLog), names_to="estimate", values_to="inf")
ggplot(rdataset1)+geom_line(aes(x=theta,y=r_mean,group=as.factor(beta),colour=as.factor(beta)))
ggplotrexp1<-ggplot(rdataset1)+geom_line(aes(x=theta,y=r_mean,group=as.factor(beta),colour=as.factor(beta)))
