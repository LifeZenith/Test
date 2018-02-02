# analyze 3min 5V5 data
setwd("/home/zhen/Metabolite_test/3min/data")
source("/home/zhen/code/Function.3min.cal.R")
wd<- paste("/home/zhen/Metabolite_test/3min/data/20180125_serum_5V5_",c("Neg","Pos"),sep = "")
res_mapped_big_list<- list()
res_not_mapped_big_list<- list()
res_mapped_list_tmp<- list()
res_not_mapped_list_tmp <- list()
target_cpd<- read.csv("/home/zhen/Metabolite_test/3min_sample_test_plasma_You_20180108/identified_merged_sort_index.csv")
target_cpd<- target_cpd[,2:4]
hash_RT<- index_reconfig(target_cpd)
unique_mass<- unique(target_cpd[,2])
row_head<- matrix(unlist(lapply(as.list(hash_RT),t)),ncol = 3, byrow = T)
colnames(row_head)<- c("ID","lowRT","highRT")

for (i in 1:2){
  setwd(wd[i])
  filelist<- dir()
  filelist<- filelist[grepl("mzXML",filelist)]
  for (j in 1:length(filelist)){
    print(j)
    res_3min<- cal_3min(filelist[j],unique_mass = unique_mass,hash_RT = hash_RT,polarity = -1)
    res_mapped_list_tmp[[j]]<- res_3min[[1]]
    res_not_mapped_list_tmp[[j]] <- res_3min[[2]]
  }
  res_mapped_big_list[[i]]<- res_mapped_list_tmp
  res_not_mapped_big_list[[i]]<- res_not_mapped_list_tmp
}

res_mapped_neg<- res_mapped_big_list[[1]]
res_mapped_neg<- Reduce(cbind,res_mapped_neg)
res_mapped_neg<- data.frame(res_mapped_neg)
colnames(res_mapped_neg)[seq(1,394,2)] <- paste("RT",substr(filelist,1,nchar(filelist)-6),sep = "_")
colnames(res_mapped_neg)[seq(2,394,2)] <- paste("Area",substr(filelist,1,nchar(filelist)-6),sep = "_")
res_mapped_neg<- cbind(row_head,res_mapped_neg)

res_mapped_pos<- res_mapped_big_list[[2]]
res_mapped_pos<- Reduce(cbind,res_mapped_pos)
res_mapped_pos<- data.frame(res_mapped_pos)
colnames(res_mapped_pos)[seq(1,394,2)] <- paste("RT",substr(filelist,1,nchar(filelist)-6),sep = "_")
colnames(res_mapped_pos)[seq(2,394,2)] <- paste("Area",substr(filelist,1,nchar(filelist)-6),sep = "_")
res_mapped_pos<- cbind(row_head,res_mapped_pos)

res_unmapped_neg<- res_not_mapped_big_list[[1]]
res_unmapped_neg<- Reduce(rbind,res_unmapped_neg)
colnames(res_unmapped_neg)<- c("mass","RT","area")

res_unmapped_pos<- res_not_mapped_big_list[[2]]
res_unmapped_pos<- Reduce(rbind,res_unmapped_pos)
colnames(res_unmapped_pos)<- c("mass","RT","area")
save(res_mapped_neg,res_unmapped_neg,file = "5V5.negative.result.Rdata")
save(res_mapped_pos,res_unmapped_pos,file = "5V5.positive.result.Rdata")


#test
#debug
filename<- filelist[31]
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
#peak_all<- peak_all[!is.na(peak_all[,1]),]
colnames(peak_all)<- c("mz","intensity","RT")
index<- floor(peak_all[,1])
key<- unique(index)
h<- hash()
for (i in 1:length(key)){
  .set(h, keys = key[i], values = peak_all[index==key[i],])
}

hash_RT<- index_reconfig(target_cpd1)
unique_mass<- unique(target_cpd1[,2])
mz_window = 0.01
RT_width = 20
span = 0.05
SNR.th = 3
tmp<- sapply(unique_mass[322:323],function(x){
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
  
  #dim(chromatogram_raw)
  if (class(chromatogram_raw)=="numeric" | class(chromatogram_raw)=="NULL"){
    chromatogram_mat<- list("NULL")
    #return(chromatogram_mat)
  } else {
    if (dim(chromatogram_raw)[1]<5){
      chromatogram_mat<- list("NULL")
      #return(chromatogram_mat)
    } else {
      unique_rt<- unique(chromatogram_raw[,1])
      intensity<- sapply(unique_rt,function(x){
        return(sum(chromatogram_raw[chromatogram_raw[,1]==x,2]))
      })
      chromatogram_mat<- cbind(unique_rt,intensity)
      #plot(chromatogram_mat[,1],chromatogram_mat[,2])
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
