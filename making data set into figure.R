

data<-read.csv("qog_ei_ts_sept21.csv")
data
data[,c("cckp_temp")]
temp<-data[,c("cckp_temp")]
annual_avg_temp<-data[,c("cckp_temp")]
hist(annual_avg_temp)