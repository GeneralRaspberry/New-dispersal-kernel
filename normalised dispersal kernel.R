dist<-seq(0,100,by=0.1)
b<-1
t<-100

mod<-1-(dist/t)^((exp(1)*b))
mod1<-mod*((exp(1)*b)+1)/(exp(1)*b*t)
plot(dist,mod1)#,ylim=c(0,1))
