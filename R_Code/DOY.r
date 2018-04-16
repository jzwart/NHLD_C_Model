#DOY.r
#Function to calculate DOY or dayfrac from date and time.
DOY<-function(dateTime){
     x<-as.POSIXlt(dateTime)
     doy<-(x$yday+(x$hour+x$min/60)/24)+1
     return(doy)
}
