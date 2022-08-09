library("spatstat")
library("dplyr")
library("ggplot2")
library("ggthemes")
library("RColorBrewer")
library("beepr")
## tau-leap Gillespie algorithm function
tauLeapG <- function(beta, # transmission rate
                     theta, # dispersal scale
                     b=1, # kernel shape parameter, 1 for exponential
                     sigma=0, # asymptomatic period, used for outputing the time series
                     q0=0, # starting incidence if ppp is without marks
                     q.end=1, # stoping condition 1: incidence lvl
                     t.end=Inf, # stoping condition 2: time after first simulated time step
                     area.host=10, # surface area occupied by one host
                     delta.t=1, # time step
                     ppp, # point pattern as a ppp object, optinally with marks 1/0 for infected/healthy
                     dist.mat=NULL){ # matrix distance if its computation is to be avoided here (for e.g. repeated calls)
  
  ## if the point pattern has no marks, generate some randomly that fits q0
  if (is.null(marks(ppp))){
    inf.start <- max(1, round(ppp$n * q0))
    marks(ppp) <- sample(c(rep(FALSE, ppp$n-inf.start), rep(TRUE, inf.start)))
  }
  
  ## compute distance matrix if not provided
  if (is.null(dist.mat)){ 
    ## add the kernel computation that can be added directly on the dist matrix to reduce comp time
    dist.mat <- 1-(pairdist(ppp)*theta)^(1/pairdist(ppp)*b)
    diag(dist.mat) <- NA
  }
  
  ## function that compute infection event probability, based on the dispersal kernel # constant part of the exponential power kernel
  infection <- function(infected, dist){
    inf <-  matrix(dist[infected,!infected],
                   nrow=sum(infected), byrow=FALSE)
    inf[is.na(inf)] <- 0
    inf
  }
  
  ## starting time
  time <- 0
  ## inititate the heavy dataframe that will aggregate all changes
  df.big <- data.frame(time=0, who=which(ppp$marks), t(ppp$marks))
  
  ## computation loop
  while (any(!ppp$marks) & time <= t.end & mean(ppp$marks) < q.end){
    ## infection event probaility
    events <- infection(ppp$marks, dist=dist.mat)
    ## random proisson draws
    new.infected <- which(!ppp$marks)[rpois(n=sum(!ppp$marks), lambda=apply(events, 2, sum) * delta.t) > 0]
    ## change marks of newly infected
    ppp$marks[new.infected] <- TRUE
    ## increment time
    time <- time + delta.t
    ## if some infection, increment the big dataframe
    if (length(new.infected) > 0){
      df.big <- rbind(df.big, data.frame(time=time, who=new.infected, t(ppp$marks)))
    }
    ## print a dot per new infection
    # cat(paste0(rep('.', length(new.infected)), collapse = '')) ## comment for quiet
  }
  
  ## make compact, time only, version of the big dataframe
  times.i <- unique(df.big[,1])
  times.d <- times.i + sigma
  times <- sort(unique(c(times.i, times.d)))
  infected <- sapply(times, FUN=function(t) sum(t >= df.big[,1]))
  detectable <- sapply(times, FUN=function(t) sum(t >= df.big[,1] + sigma))
  df.small <- data.frame(time=times, infected=infected, detectable=detectable)
  
  ## out put the simplified time series, and the big one
  list(df.small[df.small$time <= max(df.big$time),], df.big) 
} 




## meta parameters
delta.t <- 1 # time step (ALEX-THIS IS BIGGER THAN THE EXPERIMENT BELOW BECAUSE IT IS TAKING SO MUCH LONGER!)
iterations <- 10 # how many epidemic to simulate
hosts <- 1000 # number of hosts
dim <- 1000 # dimension of the landscape

## epidemic parameters
sigma <- 0 #this is the assymptomatic period, doesn't change yet

##The data I sent you, which is called data in R is the 1000 realisations of these parameters
theta1 <- 1000
theta<-1/theta1
b <- 1

infbegin<-1
randmod<-0
diseasename<-"Citrus canker"
years<-5

##################################add a timer##############################################################

ts<-proc.time()

###########################################################################################################
##Concatenating a list of metric values
##-----------------------------------------

datalist=list()

for(i in 1:iterations){  
  
  set.seed(seed=NULL)
  
landscape<- runifpoint(hosts, win = square(dim), nsim = 1)

  landscapedataframe<-as.data.frame(landscape)
  data <- data.frame(x=landscapedataframe$x, y=landscapedataframe$y, id=1:hosts)
  
  ## design a function that will be called
  
  
  set.seed(seed=NULL)
  marks(landscape)<- sample(c(rep(TRUE,infbegin), rep(FALSE, hosts-infbegin)))
  
  output <- tauLeapG(theta = theta, b = b,
                     sigma = sigma, delta.t = delta.t,
                     ppp = landscape)
  temp <- output[[2]][,1:2][order(output[[2]][,2]),]
  temp<-na.omit(temp)
 data<- data.frame(time=temp$time, who=temp$who, x=landscape$x[temp$who], y=landscape$y[temp$who]) ## what it exports will be concatenated in a list
 data$i<-i
 datalist[[i]]<-data
 print(paste0("simulation ",i," complete"))
}

big_data<-do.call(rbind,datalist)
##################################add a timer############################################################
proc.end<-proc.time()-ts
proc.end
beep()
###################################plot your data###########################################################
t2<- proc.time()

head(data)
data1<-data.frame(big_data)
times <- sort(unique(data1$time))
data_logistic <- 
  data1 %>% group_by(i)%>%
    do(data.frame(time=times, infected=sapply(times, function(x) sum(.$time <= x))))

## make a logistic df from this data

data_log<-data.frame(data_logistic)





############################################################################################################


temptimemax<-data_log%>%filter(infected<999)%>%filter(time==max(time))
temptimemax<-temptimemax[,"time"]
ggplot(data_log) + geom_line(aes(x=time, y=infected/hosts,group=i), size=.2,colour="gray70") +
theme_tufte()+xlim(0,max(times)) +
  ylab("Prevalence") +
  xlab("Time") 

proc.end2<-proc.time()-t2
proc.end2
beep()

###################################################colourplot###############################################



