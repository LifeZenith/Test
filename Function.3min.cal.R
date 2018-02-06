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

# reconfig the index to hash table
# input: target cpd Dr. Dai provided, first column cpd ID, second column mass, third column RT
# output: hash table, cpd with the same mz in one entry, key is mass
index_reconfig<- function(target_cpd){
  index<- unique(target_cpd[,2])
  h<- hash()
  for (i in 1:length(index)){
    sub_RT<- target_cpd[target_cpd[,2]==index[i],c(1,3)]
    if (dim(sub_RT)[1]>1){
      sub_RT<- sub_RT[order(sub_RT[,2]),]
      RT_dif<- diff(sub_RT[,2])
      k=1
      group_RT<- vector(length = dim(sub_RT)[1])
      for (j in 1:(dim(sub_RT)[1]-1)){
        group_RT[j] = k
        if (RT_dif[j]>10){
          k=k+1
        }
      }
      group_RT[dim(sub_RT)[1]]<- k
      RT_combine<- sapply(unique(group_RT),function(x){
        temp<- sub_RT[group_RT==x,]
        cpd_combine<- Reduce(paste,as.list(temp[,1]))
        lowRT<- min(temp[,2])-5
        highRT<- max(temp[,2])+5
        return(c(cpd_combine,lowRT,highRT))
        })
      RT_combine<- data.frame(t(RT_combine))
      RT_combine[,1]<- as.character(RT_combine[,1])
      RT_combine[,2:3]<- apply(RT_combine[,2:3],2,as.numeric)
    } else {
      RT_combine<- cbind(sub_RT[,1],sub_RT[,2]-5,sub_RT[,2]+5)
      RT_combine<- data.frame(RT_combine)
      RT_combine[,1]<- as.character(RT_combine[,1])
      RT_combine[,2:3]<- apply(RT_combine[,2:3],2,as.numeric)
    }
    .set(h,keys = index[i], values = RT_combine)
  }
  return(h)
}


RT_match<- function(index_mat,RT_calculated,peak_area,mass){
  index_map<- sapply(RT_calculated,function(x){
    left_judge<- x>index_mat[,2]
    right_judge<- x<index_mat[,3]
    if (length(which(left_judge & right_judge))==0){
      return(NA)
    } else{
      return(which(left_judge & right_judge))
    }
    })
  index_not_mapped<- which(is.na(index_map))
  res_not_mapped<- rbind(rep(mass,length(RT_calculated)),RT_calculated,peak_area)[,index_not_mapped]
  res_mapped<- sapply(1:dim(index_mat)[1],function(x){
    index_tmp<- which(index_map==x)
    if(length(index_tmp)==0){
      return(c(NA,NA))
    }
    if (length(index_tmp)==1){
      return(c(RT_calculated[index_tmp],peak_area[index_tmp]))
    }
    if (length(index_tmp)>1){
      index_max<- which.max(peak_area[index_tmp])
      index_tmp<- index_tmp[index_max]
      return(c(RT_calculated[index_tmp],peak_area[index_tmp]))
    }
    })
  return(list(res_mapped,res_not_mapped))
}


# main function of 3min calculation
# input: 
# filename - filename of mzXML file
# unique_mass - unique mass of targeted cpds, can be used to find corresponding RT in index hash table
# mz_window - search window size of mz, by default 0.01
# RT_width, span - parameters used in peak finding, by default RT_width = 20, span =0.05
# SNR.th - Signal to noise ratio threshold for the peaks, by default 3
# hash_RT - hash table format of the targeted cpds
# polarity - polarity of the experiment mode, -1: negative, 1: positive
cal_3min<- function(filename,unique_mass,mz_window = 0.01,RT_width = 20,span = 0.05,SNR.th = 3,hash_RT,polarity){
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
  tmp<- sapply(unique_mass,function(x){
    #tmp<- sapply(exact_mass[1:2],function(x){
    lowmz = x-mz_window+polarity*1.007276
    highmz = x+mz_window+polarity*1.007276
    #chromatogram_raw<- peak_all[peak_all[,1]>= lowmz & peak_all[,1] <= highmz,c(3,2)]
    if (floor(lowmz)!=floor(highmz)){
      peak_sub<- rbind(h[[as.character(floor(lowmz))]],h[[as.character(floor(highmz))]])
    } else {
      peak_sub<- h[[as.character(floor(lowmz))]]
    }
    if (class(peak_sub)=="numeric"){
      chromatogram_raw<- c(0,0)
    } else{
      chromatogram_raw<- peak_sub[peak_sub[,1]>= lowmz & peak_sub[,1] <= highmz,c(3,2)]
    }
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
    index_cpd<- hash_RT[[as.character(x)]]
    RT_map_res<- RT_match(index_cpd,RT_calculated,peak_area,x)

    #res<- cbind(mass,RT_calculated,peak_area)
    return(RT_map_res)
  })
  expr_mapped<- matrix(unlist(tmp[1,]),ncol = 2, byrow = T)
  expr_not_mapped<- matrix(unlist(tmp[2,]),ncol = 3, byrow = T)
  expr_not_mapped<- expr_not_mapped[!is.na(expr_not_mapped[,2]),]
  return(list(expr_mapped,expr_not_mapped))
}

# using multiple mz window to calculate
cal_3min_multi_mz_window<- function(filename,unique_mass,mz_window = 0.01,RT_width = 20,span = 0.05,SNR.th = 3,hash_RT,polarity){
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
  for (i in 1:length(mz_window)){
    tmp<- sapply(unique_mass,function(x){
    #tmp<- sapply(exact_mass[1:2],function(x){
    lowmz = x-mz_window[i]+polarity*1.007276
    highmz = x+mz_window[i]+polarity*1.007276
    #chromatogram_raw<- peak_all[peak_all[,1]>= lowmz & peak_all[,1] <= highmz,c(3,2)]
    if (floor(lowmz)!=floor(highmz)){
      peak_sub<- rbind(h[[as.character(floor(lowmz))]],h[[as.character(floor(highmz))]])
    } else {
      peak_sub<- h[[as.character(floor(lowmz))]]
    }
    if (class(peak_sub)=="numeric"){
      chromatogram_raw<- c(0,0)
    } else{
      chromatogram_raw<- peak_sub[peak_sub[,1]>= lowmz & peak_sub[,1] <= highmz,c(3,2)]
    }
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
    index_cpd<- hash_RT[[as.character(x)]]
    RT_map_res<- RT_match(index_cpd,RT_calculated,peak_area,x)

    #res<- cbind(mass,RT_calculated,peak_area)
    return(RT_map_res)
  })
    expr_mapped<- matrix(unlist(tmp[1,]),ncol = 2, byrow = T)
    colnames(expr_mapped)<- paste(c("RT","Area"),mz_window[i],sep = "_")
    if (i==1){
      res_mat<- expr_mapped
    } else{
      res_mat<- cbind(res_mat,expr_mapped)
    }
  }
  return(res_mat)
}

Cal_by_directory<- function(directory,mz_window,row_head,hash_RT,unique_mass,polarity,RT_width){
  filelist<- dir(directory)
  filelist<- filelist[grepl("mzXML",filelist)]
  filelist<- filelist[!grepl("Blank",filelist)]
  res_mapped_list_tmp<- list()
  #res_not_mapped_list_tmp<- list()
  for (j in 1:length(filelist)){
    print(j)
    res_3min<- cal_3min_multi_mz_window(paste(directory,filelist[j],sep = "/"),unique_mass = unique_mass,mz_window = mz_window,hash_RT = hash_RT,polarity = -1)
    colnames(res_3min)<- paste(colnames(res_3min),substr(filelist[j],1,nchar(filelist)-6),sep = "_")
    res_mapped_list_tmp[[j]]<- res_3min[[1]]
  }
  res_mapped<- Reduce(cbind,res_mapped_list_tmp)
  res_mapped<- data.frame(res_mapped)
  #colnames(res_mapped)[seq(1,(2*length(filelist)),2)] <- paste("RT",substr(filelist,1,nchar(filelist)-6),sep = "_")
  #colnames(res_mapped)[seq(2,(2*length(filelist)),2)] <- paste("Area",substr(filelist,1,nchar(filelist)-6),sep = "_")
  res_mapped<- cbind(row_head,res_mapped)
  
  #res_unmapped<- Reduce(rbind,res_not_mapped_list_tmp)
  #unmapped_length<- unlist(lapply(res_not_mapped_list_tmp,function(x){return(dim(x)[1])}))
  #file_tag<- sapply(1:length(unmapped_length),function(x){return(rep(filelist[x],unmapped_length[x]))})
  #res_unmapped<- cbind(file_tag,res_unmapped)
  #colnames(res_unmapped)<- c("sample","mass","RT","area")
  return(res_mapped)
  #return(list(res_mapped,res_unmapped))
}
