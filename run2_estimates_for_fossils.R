# 2016 03 14 I.Zliobaite
# code for public acompanying manuscript 
# estimates for fossil records
# the regression formulas are handcoded in this script, so if you change the way regression is estimated in the previous script, then the formulas need to be modified manually in this script

input_file <- 'data/Turkana_47_for_code.csv'
input_file_knn <- 'data/data_Africa_modern.csv'
input_file_bins <- 'data/Turkana_47_for_code_bins.csv'
input_file_members <- 'data/Turkana_47_for_code_members.csv'

out_file <- 'data/Turkana_47_with_estimates.csv'
out_file_bins <- 'data/Turkana_47_with_estimates_bins.csv'
out_file_members <- 'data/Turkana_47_with_estimates_members.csv'

param_k <- 15 #how many neighbors for knn

data_all <- read.csv(input_file, header = TRUE,sep = ',')
data_knn <- read.csv(input_file_knn, header = TRUE,sep = ',')
data_bins <- read.csv(input_file_bins, header = TRUE,sep = ',')
data_members <- read.csv(input_file_members, header = TRUE,sep = ',')

#add mean HYP and mean LOP, add estimates

un_comlocalities <- unique(data_all[,'ComLoc'])
data_comlocalities <- c()
for (sk in 1:length(un_comlocalities)){
  #print(sk)
  comloc_now <- un_comlocalities[sk]
  ind <- which(data_all[,'ComLoc']==comloc_now)
  data_now <- data_all[ind,]
  ind_unique_species <- which(data_now[,'unique_species_used']==1)
  # mean HYP and LOP over unique species (like in modern day models)
  meanHYP <- round(mean(data_now[ind_unique_species,'HYP']),digits = 2)
  meanLOP <- round(mean(data_now[ind_unique_species,'LOP']),digits = 2)
  MATreg <- round(27.8 - 1.1*meanHYP - 1.2*meanLOP,digits = 1)
  MAPreg <- round(1251.9 - 460.9*meanHYP + 2237.1*meanLOP - 823.7*meanHYP*meanLOP)
  if (MAPreg < 0){MAPreg <- 0}
  #knn estimation
  ref_now <- c(meanHYP,meanLOP)
  names(ref_now) <- c('HYP','LOP')
  dd <- dist(rbind(ref_now,data_knn[,c('HYP','LOP')]))
  dd <- dd[1:dim(data_knn)[1]] #take only first column of the triangle
  ord <- order(dd)
  MATknn <- round(mean(data_knn[ord[1:(param_k)],'MAT']),digits = 1)
  MAPknn <- round(mean(data_knn[ord[1:(param_k)],'MAP']))
  data_comlocalities <- rbind(data_comlocalities,cbind(as.vector(comloc_now),meanHYP,meanLOP,MAPreg,MATreg,MAPknn,MATknn,length(ind_unique_species),length(ind)))
}
colnames(data_comlocalities) <- c('ComLoc','meanHYP','meanLOP','MAPreg','MATreg','MAPknn','MATknn','no_species','no_specimen')

data_out <- data_all
data_out[,c('meanHYP','meanLOP','MAPreg','MATreg','MAPknn','MATknn')] <- NA

for (sk in 1:dim(data_comlocalities)[1]){
  comloc_now <- data_comlocalities[sk,'ComLoc']
  ind <- which(data_out[,'ComLoc'] == comloc_now)
  data_out[ind,'meanHYP'] <- data_comlocalities[sk,'meanHYP']
  data_out[ind,'meanLOP'] <- data_comlocalities[sk,'meanLOP']
  data_out[ind,'MAPreg'] <- data_comlocalities[sk,'MAPreg']
  data_out[ind,'MATreg'] <- data_comlocalities[sk,'MATreg']
  data_out[ind,'MAPknn'] <- data_comlocalities[sk,'MAPknn']
  data_out[ind,'MATknn'] <- data_comlocalities[sk,'MATknn']
}
write.table(data_out, file = out_file, row.names = FALSE, col.names = TRUE, sep = ',', quote = FALSE)


# same for bins
un_bins <- unique(data_bins[,'ComLoc'])
data_comlocalities_bins <- c()
for (sk in 1:length(un_bins)){
  #print(sk)
  bin_now <- un_bins[sk]
  ind <- which(data_bins[,'ComLoc']==bin_now)
  data_now <- data_bins[ind,]
  ind_unique_species <- which(data_now[,'unique_species_used']==1)
  # mean HYP and LOP over unique species (like in modern day models)
  meanHYP <- round(mean(data_now[ind_unique_species,'HYP']),digits = 2)
  meanLOP <- round(mean(data_now[ind_unique_species,'LOP']),digits = 2)
  MATreg <- round(27.8 - 1.1*meanHYP - 1.2*meanLOP,digits = 1)
  MAPreg <- round(1251.9 - 460.9*meanHYP + 2237.1*meanLOP - 823.7*meanHYP*meanLOP)
  if (MAPreg < 0){MAPreg <- 0}
  #knn estimation
  ref_now <- c(meanHYP,meanLOP)
  names(ref_now) <- c('HYP','LOP')
  dd <- dist(rbind(ref_now,data_knn[,c('HYP','LOP')]))
  dd <- dd[1:dim(data_knn)[1]] #take only first column of the triangle
  ord <- order(dd)
  MATknn <- round(mean(data_knn[ord[1:(param_k)],'MAT']),digits = 1)
  MAPknn <- round(mean(data_knn[ord[1:(param_k)],'MAP']))
  data_comlocalities_bins <- rbind(data_comlocalities_bins,cbind(as.vector(bin_now),meanHYP,meanLOP,MAPreg,MATreg,MAPknn,MATknn,length(ind_unique_species),length(ind)))
}
colnames(data_comlocalities_bins) <- c('ComLoc','meanHYP','meanLOP','MAPreg','MATreg','MAPknn','MATknn','no_species','no_specimen')

data_out_bins <- data_bins
data_out_bins[,c('meanHYP','meanLOP','MAPreg','MATreg','MAPknn','MATknn')] <- NA

for (sk in 1:dim(data_comlocalities_bins)[1]){
  bin_now <- data_comlocalities_bins[sk,'ComLoc']
  ind <- which(data_out_bins[,'ComLoc'] == bin_now)
  data_out_bins[ind,'meanHYP'] <- data_comlocalities_bins[sk,'meanHYP']
  data_out_bins[ind,'meanLOP'] <- data_comlocalities_bins[sk,'meanLOP']
  data_out_bins[ind,'MAPreg'] <- data_comlocalities_bins[sk,'MAPreg']
  data_out_bins[ind,'MATreg'] <- data_comlocalities_bins[sk,'MATreg']
  data_out_bins[ind,'MAPknn'] <- data_comlocalities_bins[sk,'MAPknn']
  data_out_bins[ind,'MATknn'] <- data_comlocalities_bins[sk,'MATknn']
}
write.table(data_out_bins, file = out_file_bins, row.names = FALSE, col.names = TRUE, sep = ',', quote = FALSE)


# same for bins East and West separately (needed for Fig 7)
un_members <- unique(data_members[,'ComLoc'])
data_comlocalities_members <- c()
for (sk in 1:length(un_members)){
  #print(sk)
  bin_now <- un_members[sk]
  ind <- which(data_members[,'ComLoc']==bin_now)
  data_now <- data_members[ind,]
  ind_unique_species <- which(data_now[,'unique_species_used']==1)
  # mean HYP and LOP over unique species (like in modern day models)
  meanHYP <- round(mean(data_now[ind_unique_species,'HYP']),digits = 2)
  meanLOP <- round(mean(data_now[ind_unique_species,'LOP']),digits = 2)
  MATreg <- round(27.8 - 1.1*meanHYP - 1.2*meanLOP,digits = 1)
  MAPreg <- round(1251.9 - 460.9*meanHYP + 2237.1*meanLOP - 823.7*meanHYP*meanLOP)
  if (MAPreg < 0){MAPreg <- 0}
  #knn estimation
  ref_now <- c(meanHYP,meanLOP)
  names(ref_now) <- c('HYP','LOP')
  dd <- dist(rbind(ref_now,data_knn[,c('HYP','LOP')]))
  dd <- dd[1:dim(data_knn)[1]] #take only first column of the triangle
  ord <- order(dd)
  MATknn <- round(mean(data_knn[ord[1:(param_k)],'MAT']),digits = 1)
  MAPknn <- round(mean(data_knn[ord[1:(param_k)],'MAP']))
  data_comlocalities_members <- rbind(data_comlocalities_members,cbind(as.vector(bin_now),meanHYP,meanLOP,MAPreg,MATreg,MAPknn,MATknn,length(ind_unique_species),length(ind)))
}
colnames(data_comlocalities_members) <- c('ComLoc','meanHYP','meanLOP','MAPreg','MATreg','MAPknn','MATknn','no_species','no_specimen')

data_out_members <- data_members
data_out_members[,c('meanHYP','meanLOP','MAPreg','MATreg','MAPknn','MATknn')] <- NA

for (sk in 1:dim(data_comlocalities_members)[1]){
  member_now <- data_comlocalities_members[sk,'ComLoc']
  ind <- which(data_out_members[,'ComLoc'] == member_now)
  data_out_members[ind,'meanHYP'] <- data_comlocalities_members[sk,'meanHYP']
  data_out_members[ind,'meanLOP'] <- data_comlocalities_members[sk,'meanLOP']
  data_out_members[ind,'MAPreg'] <- data_comlocalities_members[sk,'MAPreg']
  data_out_members[ind,'MATreg'] <- data_comlocalities_members[sk,'MATreg']
  data_out_members[ind,'MAPknn'] <- data_comlocalities_members[sk,'MAPknn']
  data_out_members[ind,'MATknn'] <- data_comlocalities_members[sk,'MATknn']
}
write.table(data_out_members, file = out_file_members, row.names = FALSE, col.names = TRUE, sep = ',', quote = FALSE)