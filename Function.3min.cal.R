# functions used in 3min calculation
# 20180129
# Zhen Li

library(mzR)
library(zoo)

chromatogram_fill<- function(hd,chromatogram_mat){
  average_level<- median(chromatogram_mat[,2])
  chromatogram_full<- cbind(hd$retentionTime,average_level)
  chromatogram_full[match(chromatogram_mat[,1],chromatogram_full[,1]),2]<- chromatogram_mat[,2]
  return(chromatogram_full)
}

peak_area_cal<- function(chromatogram_mat,RT){
  noise_level<- median(chromatogram_mat[,2])
  peak_area<- sapply(RT,function(x){
    cal_window<- chromatogram_mat
    peak_intensity<- chromatogram_mat[chromatogram_mat[,1]==x,2]
    index_out<- which(cal_window[,2]< 0.1*(peak_intensity-noise_level)+noise_level)
    index_peak<- which(cal_window[,1]==x)
    RT_window<- c(max(1,index_out[index_out<index_peak]),min(index_out[index_out>index_peak],dim(chromatogram_mat)[1]))
    cal_window<- cal_window[RT_window[1]:RT_window[2],]
    area<- area_cal(cal_window,noise_level)
    return(area)
  })
  return(peak_area)
}

area_cal<- function(intensity_mat,noise_level){
  n<- dim(intensity_mat)[1]-1
  ave_index<- sapply((1:n),function(x){return(intensity_mat[x,2]+intensity_mat[x+1,2]-2*noise_level)})/2
  area<- sum(diff(intensity_mat[,1])*ave_index)
  return(area)
}

peak_find <- function(x, y, width, span,SNR.th) {
  n <- length(y)
  y.smooth <- loess(y ~ x, span = span)$fitted
  y.max <- rollapply(zoo(y.smooth), 2*width+1, max, align="center")
  delta <- y.max - y.smooth[-c(1:width, n+1-1:width)]
  i.max <- which(delta <= 0) + width
  if (length(i.max)>0){
    SNR<- SNR_cal(x,y,i.max,width)
  SNR[is.na(SNR)]<- 1
  i.max<- i.max[SNR>max(SNR.th,max(SNR,na.rm=T)/3)]
  }
  return(list(RT=x[i.max], peak_index=i.max, y.hat=y.smooth))
}

#SNR calculation
SNR_cal<- function(x,y,peak_index,width){
  signal_center<- x[peak_index]
  signal_index<- sapply(peak_index,function(t){
    signal_window<- c(x[t]-width/2,x[t]+width/2)
    index<- which(x>= signal_window[1] & x<=signal_window[2])
    return(index)
  })
  signal_power<- sapply(peak_index,function(t){
    signal_window<- c(x[t]-width/2,x[t]+width/2)
    index<- which(x>= signal_window[1] & x<=signal_window[2])
    sig_power<- max(y[index])
    return(sig_power)
  })
  noise_index<- setdiff((1:length(x)),unlist(signal_index))
  noise_power<- mean(y[noise_index])
  SNR<- signal_power/noise_power
  return(SNR)
}
