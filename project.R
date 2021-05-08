
######################################
# SBP changing over time
######################################

library(dplyr )
library(sas7bdat)
library(msm)
library(ctmcmove)
baseline <- read.sas7bdat("/extraspace/ychen42/homework/data/baseline.sas7bdat")
all <-read.sas7bdat("/extraspace/ychen42/homework/data/all_visits.sas7bdat")


all$state <- ifelse(all$sbp<130,1,(ifelse(all$sbp<140,2,ifelse(all$sbp<160,3,ifelse(all$sbp<180,4,5)))))
all$interval_year<-all$interval/365
all$interval_year[all$cvisit==1] <- 0

data <- merge(baseline,all,by=c("SIID","PTID"),all=FALSE)
data$id <- data$SIID*10000+data$PTID


#data$state <- ifelse(data$sbp<130,1,(ifelse(data$sbp<140,2,ifelse(data$sbp<160,3,ifelse(data$sbp<180,4,5)))))

data$interval_month <- round(data$interval_year*365/12)
data$interval_week <- round(data$interval_month/4)

data <- data[order(data$id,data$interval_year),]

(statetable<- statetable.msm(state,id, data=data))

Q <- rbind ( c(0, 0.1, 0.1, 0.1,0.1,0.1),
             + c(0.1, 0,0.1,0.1,0.1,0.1),
             + c(0.1, 0.1, 0, 0.1,0.1,0.1),
             + c(0.1, 0.1, 0.1, 0,0.1,0.1),
             + c(0.1, 0.1, 0.1, 0.1,0,0.1),
             + c(0,0,0,0,0,0))

(Q.crude <- crudeinits.msm(state ~ interval_year, id, data=data,qmatrix=Q))

(msm.fit <- msm(state~interval_year,subject=id,data=data, qmatrix=Q,na.action=na.omit,exacttimes=TRUE))
(msm.fit_assign <- msm(state~interval_year,subject=id,data=data, qmatrix=Q,na.action=na.omit,exacttimes=TRUE,covariates=~assign2))
(msm.fit_cov <- msm(state~interval_year,subject=id,data=data, qmatrix=Q,na.action=na.omit,exacttimes=TRUE,covariates=~assign2+race+GENDER))
qmatrix.msm(msm.fit)
qmatrix.msm(msm.fit_cov)
qmatrix.msm(msm.fit_cov, covariates=list(assign2=0))
qmatrix.msm(msm.fit_cov, covariates=list(assign2=1))
pmatrix.msm(msm.fit, t=9)
pnext.msm(msm.fit)
pnext.msm(msm.fit_cov)
totlos.msm(msm.fit)
hazard.msm(msm.fit)

#plot(msm.fit,legend.pos=c(8, 1))
#(ematrix.msm(msm.fit, covariates="mean", ci=c("delta"),
#             cl=0.95, B=1000, cores=NULL))

#bootstrapping
q.list <- boot.msm(msm.fit, stat=function(x){qmatrix.msm(x)$estimates})
q.array <- array(unlist(q.list), dim=c(4,4,1000))
apply(q.array, c(1,2), sd)
apply(q.array, c(1,2), function(x)quantile(x, c(0.025, 0.975)))


#Model assessment
#Observed and expected prevalence
options(digits=3)
prevalence.msm(msm.fit, times=seq(0,20,2))
plot.prevalence.msm(msm.fit, mintime=0, maxtime=20)

#Pearson-type goodness-of-fit test, escaped
#### Next step
# Hidden markov chain
# Identify potential confounding


######################################
# High BP with death
######################################
outcomes <- read.sas7bdat("W:/Office/Documents/INVEST/SASData/outcomes.sas7bdat")

#start
start <- all[all$cvisit==1]
start$state <- ifelse(start$sbp<160,1,2)

#identify risky bp
hyper <- all[all$state>3,]
hyper$state <- 2

#post_cva
cva_raw <- outcomes[outcomes$cvaevent=1]
cva <- data.frame(SIID=cva_raw$SIID,PTID=cva_raw$PTID,state=3,interval_year=cva_raw$cvadur)

#post_mi
mi_raw <- outcomes[outcomes$mievent=1]
mi <- data.frame(SIID=mi_raw$SIID,PTID=mi_raw$PTID,state=4,interval_year=mi_raw$midur)

#death
death_raw <- outcomes[outcomes$deevent=1]
death <- data.frame(SIID=death_raw$SIID,PTID=death_raw$PTID,state=5,interval_year=death_raw$dedur)


comb <- rbind(mi,death)
comb <- comb[order(comb$interval_year, comb$state, comb$PTID, comb$SIID,decreasing=TRUE),]
# id var
#  4   2
#  3   5
#  2   3
#  2   1
#  1   4
#  1   2

# Keep only the first row for each duplicate of z$id; this row will have the
# largest value for z$var
z <- z[!duplicated(z$id),]

# Sort so it looks nice
z <- z[order(z$id, z$var),]
# id var
#  1   4
#  2   3
#  3   5
#  4   2


all$state <- ifelse(all$sbp<130,1,(ifelse(all$sbp<140,2,ifelse(all$sbp<160,3,ifelse(all$sbp<180,4,5)))))
all$interval_year<-all$interval/365
all$interval_year[all$cvisit==1] <- 0

visits <- data.frame(SIID=all$SIID,PTID=all$PTID,state=all$state,interval_year=all$interval_year)


sbp_death <- rbind(visits,death)


data <- merge(baseline,sbp_death,by=c("SIID","PTID"),all=FALSE)
data$id <- data$SIID*10000+data$PTID



#data$state <- ifelse(data$sbp<130,1,(ifelse(data$sbp<140,2,ifelse(data$sbp<160,3,ifelse(data$sbp<180,4,5)))))

data$interval_month <- round(data$interval_year*365/12)
data$interval_week <- round(data$interval_month/4)

data <- data[order(data$id,data$interval_year),]

(statetable<- statetable.msm(state,id, data=data))

Q <- rbind ( c(0, 0.1, 0.1, 0.1,0.1,0.1),
             + c(0.1, 0,0.1,0.1,0.1,0.1),
             + c(0.1, 0.1, 0, 0.1,0.1,0.1),
             + c(0.1, 0.1, 0.1, 0,0.1,0.1),
             + c(0.1, 0.1, 0.1, 0.1,0,0.1),
             + c(0,0,0,0,0,0))

(Q.crude <- crudeinits.msm(state ~ interval_year, id, data=data,qmatrix=Q))

(msm.fit <- msm(state~interval_year,subject=id,data=data, qmatrix=Q,na.action=na.omit,exacttimes=TRUE, deathexact = 6))
(msm.fit <- msm(state~interval_year,subject=id,data=data, qmatrix=Q,na.action=na.omit,exacttimes=TRUE,covariates=~assign2))
(msm.fit <- msm(state~interval_year,subject=id,data=data, qmatrix=Q,na.action=na.omit,exacttimes=TRUE,covariates=~assign2+race+GENDER))
qmatrix.msm(msm.fit)
qmatrix.msm(msm.fit, covariates=list(assign2=0))
qmatrix.msm(msm.fit, covariates=list(assign2=1))
pmatrix.msm(msm.fit, t=9)
pnext.msm(msm.fit)
totlos.msm(msm.fit)
hazard.msm(msm.fit)
plot(msm.fit,legend.pos=c(8, 1))
(ematrix.msm(msm.fit, covariates="mean", ci=c("delta"),
             cl=0.95, B=1000, cores=NULL))

#bootstrapping
q.list <- boot.msm(cav.msm, stat=function(x){qmatrix.msm(x)$estimates})
q.array <- array(unlist(q.list), dim=c(4,4,1000))
apply(q.array, c(1,2), sd)
apply(q.array, c(1,2), function(x)quantile(x, c(0.025, 0.975)))


#Model assessment
#Observed and expected prevalence
options(digits=3)
prevalence.msm(cav.msm, times=seq(0,20,2))
plot.prevalence.msm(cav.msm, mintime=0, maxtime=20)
