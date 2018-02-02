# using new 3 min data to analyze DBS, serum and plasma 10*3*3
setwd("/home/storage/storage5/MS_data_storage/MS_data_Tianjin/3min/10x3x3")
source("/home/zhen/code/Function.3min.cal.R")

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

Cal_one_day<- function(directory){
  filelist<- dir(directory)
  filelist<- filelist[grepl("mzXML",filelist)]
  for (j in 1:length(filelist)){
    print(j)
    res_3min<- cal_3min(paste(directory,filelist[j],sep = "/"),unique_mass = unique_mass,hash_RT = hash_RT,polarity = -1)
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
folder_list<- folder_list[!grepl(".xlsx",folder_list)]
folder_list<- paste(getwd(),folder_list,sep = "/")

expr_20180126_neg<- Cal_one_day(paste(folder_list[1],"neg",sep = "/"))
expr_20180126_pos<- Cal_one_day(paste(folder_list[1],"pos",sep = "/"))

expr_20180129_neg<- Cal_one_day(paste(folder_list[2],"neg",sep = "/"))
expr_20180129_pos<- Cal_one_day(paste(folder_list[2],"pos",sep = "/"))

expr_20180130_neg<- Cal_one_day(paste(folder_list[3],"neg",sep = "/"))
expr_20180130_pos<- Cal_one_day(paste(folder_list[3],"pos",sep = "/"))
