# using different RT window to calculate 600 standard data
# Date: 20180203
# Author: Zhen
setwd("/home/storage/storage5/MS_data_storage/MS_data_Tianjin/compound_discoverer_Standard_library_build_TJ/3min")
target_cpd<- read.csv("/home/zhen/Metabolite_test/3min_sample_test_plasma_You_20180108/identified_merged_sort_index.csv")
target_cpd<- target_cpd[,2:4]
hash_RT<- index_reconfig(target_cpd)
unique_mass<- unique(target_cpd[,2])
row_head<- matrix(unlist(lapply(as.list(hash_RT),t)),ncol = 3, byrow = T)
colnames(row_head)<- c("ID","lowRT","highRT")

folder_list<- dir()
folder_list<- folder_list[!grepl(".png",folder_list)]
folder_list<- folder_list[!grepl(".xlsx",folder_list)]
folder_list<- folder_list<- paste(getwd(),folder_list,sep = "/")
res_mapped_neg_list<- list()
res_mapped_pos_list<- list()

for (i in 1:length(folder_list)){
  tmp_neg<- Cal_one_day(paste(folder_list[i],"neg_organic",sep = "/"),row_head = row_head, hash_RT = hash_RT, unique_mass = unique_mass, polarity = -1)
  tmp_pos<- Cal_one_day(paste(folder_list[i],"pos_organic",sep = "/"),row_head = row_head, hash_RT = hash_RT, unique_mass = unique_mass, polarity = 1)
  res_mapped_neg_list[[i]]<- tmp_neg[[1]]
  res_mapped_pos_list[[i]]<- tmp_pos[[1]]
}
