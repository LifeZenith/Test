# using different RT window to calculate 600 standard data
# Date: 20180203
# Author: Zhen
setwd("/home/storage/storage5/MS_data_storage/MS_data_Tianjin/compound_discoverer_Standard_library_build_TJ/3min")
target_cpd<- read.csv("/home//storage/storage5/MS_data_storage/MS_data_US/3min_sample_test_plasma_You_20180108/identified_merged_sort_index.csv",as.is = T)
target_cpd<- target_cpd[,2:4]
hash_RT<- index_reconfig(target_cpd)
unique_mass<- unique(target_cpd[,2])
test<- sapply(unique_mass,function(x){
  return(t(hash_RT[[as.character(x)]]))
})
row_head<- t(Reduce(cbind,test))
#test1<- matrix(unlist(test),ncol = 3, byrow = T)

#row_head<- matrix(unlist(lapply(as.list(hash_RT),t)),ncol = 3, byrow = T)
colnames(row_head)<- c("ID","lowRT","highRT")

Cal_by_directory<- function(directory,mz_window,row_head,hash_RT,unique_mass,polarity,RT_width){
  filelist<- dir(directory)
  filelist<- filelist[grepl("mzXML",filelist)]
  for (j in 1:length(filelist)){
    print(j)
    res_3min<- cal_3min(paste(directory,filelist[j],sep = "/"),unique_mass = unique_mass,mz_window = mz_window,hash_RT = hash_RT,polarity = -1)
    res_mapped_list_tmp[[j]]<- res_3min[[1]]
    res_not_mapped_list_tmp[[j]] <- res_3min[[2]]
  }
  res_mapped<- Reduce(cbind,res_mapped_list_tmp)
  res_mapped<- data.frame(res_mapped)
  colnames(res_mapped)[seq(1,(2*length(filelist)),2)] <- paste("RT",substr(filelist,1,nchar(filelist)-6),sep = "_")
  colnames(res_mapped)[seq(2,(2*length(filelist)),2)] <- paste("Area",substr(filelist,1,nchar(filelist)-6),sep = "_")
  res_mapped<- cbind(row_head,res_mapped)
  
  res_unmapped<- Reduce(rbind,res_not_mapped_list_tmp)
  unmapped_length<- unlist(lapply(res_not_mapped_list_tmp,function(x){return(dim(x)[1])}))
  file_tag<- sapply(1:length(unmapped_length),function(x){return(rep(filelist[x],unmapped_length[x]))})
  res_unmapped<- cbind(file_tag,res_unmapped)
  colnames(res_unmapped)<- c("sample","mass","RT","area")
  
  return(list(res_mapped,res_unmapped))
}


folder_list<- dir()
folder_list<- folder_list[!grepl(".png",folder_list)]
folder_list<- folder_list[!grepl(".xlsx",folder_list)]
folder_list<- folder_list<- paste(getwd(),folder_list,sep = "/")
res_mapped_neg_list<- list()
res_mapped_pos_list<- list()



# mz_window 0.01
for (i in 1:length(folder_list)){
  tmp_neg<- Cal_by_directory(paste(folder_list[i],"neg_organic",sep = "/"),mz_window = 0.01, row_head = row_head, hash_RT = hash_RT, unique_mass = unique_mass, polarity = -1)
  tmp_pos<- Cal_by_directory(paste(folder_list[i],"pos_organic",sep = "/"),mz_window = 0.01, row_head = row_head, hash_RT = hash_RT, unique_mass = unique_mass, polarity = 1)
  res_mapped_neg_list[[i]]<- tmp_neg[[1]]
  res_mapped_pos_list[[i]]<- tmp_pos[[1]]
}
