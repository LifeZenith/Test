# functions used in 3min calculation
# 20180129
# Zhen Li

library(mzR)
library(zoo)
library(hash)

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

cal_3min<- function(filename,unique_mz,mz_window,RT_width,span,SNR.th){
  ms<- openMSfile(filename)
  hd<- header(ms)
  peaks<- peaks(ms)
  peaks<- lapply(peaks,function(x){
    return(t(x[x[,2]!=0,]))
  })
  length_peak<- unlist(lapply(peaks,function(x){return(dim(x)[2])}))
  peak_all<- matrix(unlist(peaks),ncol = 2,byrow = T)
  rt<- sapply(1:length(length_peak),function(x){return(rep(hd$retentionTime[x],length_peak[x]))})
  rt<- unlist(rt)
  peak_all<- cbind(peak_all,rt)
  peak_all<- peak_all[!is.na(peak_all[,1]),]
  colnames(peak_all)<- c("mz","intensity","RT")
  index<- floor(peak_all[,1])
  key<- unique(index)
  h<- hash()
  for (i in 1:length(key)){
    .set(h, keys = key[i], values = peak_all[index==key[i],])
  }
  tmp<- sapply(unique_mz,function(x){
    #tmp<- sapply(exact_mass[1:2],function(x){
    lowmz = x-mz_window
    highmz = x+mz_window
    #chromatogram_raw<- peak_all[peak_all[,1]>= lowmz & peak_all[,1] <= highmz,c(3,2)]
    if (floor(lowmz)!=floor(highmz)){
      peak_sub<- rbind(h[[as.character(floor(lowmz))]],h[[as.character(floor(highmz))]])
    } else {
      peak_sub<- h[[as.character(floor(x))]]
    }
    chromatogram_raw<- peak_sub[peak_sub[,1]>= lowmz & peak_sub[,1] <= highmz,c(3,2)]
    if (class(chromatogram_raw)=="numeric" | class(chromatogram_raw)=="NULL"){
      chromatogram_mat<- list("NULL")
      #return(chromatogram_mat)
    } else {
      if (dim(chromatogram_raw)[1]<3){
        chromatogram_mat<- list("NULL")
        #return(chromatogram_mat)
      } else {
        unique_rt<- unique(chromatogram_raw[,1])
        intensity<- sapply(unique_rt,function(x){
          return(sum(chromatogram_raw[chromatogram_raw[,1]==x,2]))
        })
        chromatogram_mat<- cbind(unique_rt,intensity)
        colnames(chromatogram_mat)<- c("RT","intensity")
        #average_level<- mean(chromatogram_mat[order(chromatogram_mat[,2])[(1+floor(dim(chromatogram_mat)[1]/10)):(dim(chromatogram_mat)[1]-floor(dim(chromatogram_mat)[1]/10))],2])
        #return(chromatogram_mat)
      }
    }
    if (class(chromatogram_mat)=="list"){
      mass<- NA
      RT_calculated<- NA
      peak_area<- NA
    } else {
      chromatogram_full<- chromatogram_fill(hd,chromatogram_mat)
      chromatogram_peak<- peak_find(chromatogram_full[,1],chromatogram_full[,2],width = RT_width, span = span,SNR.th = SNR.th)
      #plot chromatogram
      #plot_chromatogram(paste("Chromatogram",cpd[k,2],sub_folder_list[j],"png",sep = "."),chromatogram_mat[,1],chromatogram_mat[,2],RT_width, span = span, RT_peak = chromatogram_peak)
      #plot_chromatogram("test.png",chromatogram_full[,1],chromatogram_full[,2],RT_width, span = span, RT_peak = chromatogram_peak)
      #plot_chromatogram(paste(folder_list,cpd_group[j],x,"png",sep = "."),chromatogram_full[,1],chromatogram_full[,2],RT_width, span = span, RT_peak = chromatogram_peak)
      RT_calculated<- chromatogram_peak$RT
      if (length(RT_calculated)>0){
        peak_area<- peak_area_cal(chromatogram_mat,RT_calculated)
        mass<- rep(x,length(RT_calculated))
      } else {
        mass<- NA
        RT_calculated<- NA
        peak_area<- NA
      }
    }
    res<- cbind(mass,RT_calculated,peak_area)
    return(res)
  })
  return(tmp)
}
