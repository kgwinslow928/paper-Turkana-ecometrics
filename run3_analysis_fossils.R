# 2016 03 16 I.Zliobaite
# code for public acompanying manuscript 
# analysis plots

input_file <- 'data/Turkana_47_with_estimates.csv'
input_file_bins <- 'data/Turkana_47_with_estimates_bins.csv'
input_file_members <- 'data/Turkana_47_with_estimates_members.csv'
input_file_carnivors <- 'data/Turkana_47_for_code_carnivors.csv'

#for plotting
ht <- 5 #plot height
wd <- 6 #plot width
param_span <- 0.75 #for loess
param_deg <- 1 #for loess

data_all <- read.csv(input_file, header = TRUE,sep = ',')
data_bins <- read.csv(input_file_bins, header = TRUE,sep = ',')
data_members <- read.csv(input_file_members, header = TRUE,sep = ',')
data_carnivors <- read.csv(input_file_carnivors, header = TRUE,sep = ',')
data_KBS <- data_all

# filtering age >0.9'
ind <- which(data_all[,'MidMemberAge']>=0.9)
data_all <- data_all[ind,]
ind <- which(data_bins[,'TurkanaTimeBin']>=0.9)
data_bins <- data_bins[ind,]
ind <- which(data_members[,'MidMemberAge']>=0.9)
data_members <- data_members[ind,]
ind <- which(data_carnivors[,'MidMemberAge']>=0.9)
data_carnivors <- data_carnivors[ind,]


ind <- !duplicated(data_all$ComLoc)
data_sites <- data_all[ind,]
ind <- !duplicated(data_bins$ComLoc)
data_sites_bins <- data_bins[ind,]
ind <- !duplicated(data_members$ComLoc)
data_sites_members <- data_members[ind,]

#main plots
plot_predictions <- function(plot_name_now,data_sites,fet_plot,fet_age,param_span,param_deg,ht,wd,title_suffix){
if (substr(fet_plot,1,3)=='MAP'){
  ymin <- -100
  ymax <- 1700
  ylb <- 'MAP estimate'
}else{
  if (substr(fet_plot,1,3)=='MAT'){
    ymin <- 18
    ymax <- 30  
    ylb <- 'MAT estimate'
  }
}
if (substr(fet_plot,4,6)=='reg'){
  mn = paste('regression',title_suffix)
}else{
  if (substr(fet_plot,4,6)=='knn'){
    mn = paste('k-nearest neighbour model',title_suffix)
  }
}
mid_point <- data_sites[,fet_age]
target <- data_sites[,fet_plot]
data_smoo <- as.data.frame(cbind(mid_point,target))
data_plot = data.frame(mid_point = seq(0,10,0.1),target = NA)
fit.loess = loess(target ~ mid_point,data_smoo,family='gaussian', span=param_span, degree=param_deg)
pred = cbind(predict(fit.loess, data_plot),
      predict(fit.loess, data_plot)+predict(fit.loess, data_plot, se=TRUE)$se.fit*qnorm(1-.05/2),
      predict(fit.loess, data_plot)-predict(fit.loess, data_plot, se=TRUE)$se.fit*qnorm(1-.05/2))  
pdf(plot_name_now,width = wd,height = ht)
plot(NA,NA,xlim = c(0,8),ylim = c(ymin,ymax),xlab = 'age',ylab = ylb, main = mn)
ind3 <- which(!is.na(pred[,3])) #this removes early and late points which are out of age range
ind2 <- which(!is.na(pred[,2])) #should be the same values as above
polygon(c(data_plot[ind3,1],rev(data_plot[ind2,1])),c(pred[ind2,2],rev(pred[ind3,3])),col="goldenrod1",border = NA)    
points(mid_point,data_sites[,fet_plot],pch = 20,col = 'black',cex = 0.8)
lines(data_plot[,1],pred[,1],col = 'darkorange2',lwd = 4)
dev.off()
}
  
plot_predictions('figures/fig4a_est_MAPreg.pdf',data_sites,'MAPreg','MidMemberAge',param_span,param_deg,ht,wd,'by ComLoc')
plot_predictions('figures/fig4c_est_MAPknn.pdf',data_sites,'MAPknn','MidMemberAge',param_span,param_deg,ht,wd,'by ComLoc')
plot_predictions('figures/fig3a_est_MATreg.pdf',data_sites,'MATreg','MidMemberAge',param_span,param_deg,ht,wd,'by ComLoc')
plot_predictions('figures/fig3c_est_MATknn.pdf',data_sites,'MATknn','MidMemberAge',param_span,param_deg,ht,wd,'by ComLoc')
plot_predictions('figures/fig4b_est_MAPreg_bins.pdf',data_sites_bins,'MAPreg','TurkanaTimeBin',param_span,param_deg,ht,wd,'by Bin')
plot_predictions('figures/fig4d_est_MAPknn_bins.pdf',data_sites_bins,'MAPknn','TurkanaTimeBin',param_span,param_deg,ht,wd,'by Bin')
plot_predictions('figures/fig3b_est_MATreg_bins.pdf',data_sites_bins,'MATreg','TurkanaTimeBin',param_span,param_deg,ht,wd,'by Bin')
plot_predictions('figures/fig3d_est_MATknn_bins.pdf',data_sites_bins,'MATknn','TurkanaTimeBin',param_span,param_deg,ht,wd,'by Bin')

# East West
plot_predictions_EW <- function(plot_name_now,data_sites,fet_plot,fet_age,param_span,param_deg,ht,wd,title_suffix){
  ymin <- -100
  ymax <- 1700
  ylb <- 'MAP estimate'
  if (substr(fet_plot,4,6)=='reg'){
    mn = paste('regression',title_suffix)
  }else{
    if (substr(fet_plot,4,6)=='knn'){
      mn = paste('k-nearest neighbour model',title_suffix)
    }
  }
  mid_point <- data_sites[,fet_age]
  target <- data_sites[,fet_plot]
  data_smoo <- as.data.frame(cbind(mid_point,target))
  data_plot = data.frame(mid_point = seq(0,10,0.1),target = NA)
  indE <- which(data_sites[,'EastWest']=='E')
  indW <- which(data_sites[,'EastWest']=='W')
  fit.loess_E = loess(target ~ mid_point,data_smoo[indE,],family='gaussian', span=param_span, degree=param_deg)
  fit.loess_W = loess(target ~ mid_point,data_smoo[indW,],family='gaussian', span=param_span, degree=param_deg)
  pred_E <- predict(fit.loess_E, data_plot)
  pred_W <- predict(fit.loess_W, data_plot)  
  data_plot_E <- data_plot
  data_plot_W <- data_plot
  pdf(plot_name_now,width = wd,height = ht)
  plot(NA,NA,xlim = c(0,8),ylim = c(ymin,ymax),xlab = 'age',ylab = ylb,main = mn)
  points(mid_point[indE],data_sites[indE,fet_plot],pch = 20,col = 'deepskyblue',cex = 0.8)
  points(mid_point[indW],data_sites[indW,fet_plot],pch = 20,col = 'darkorange',cex = 0.8)  
  lines(data_plot_E[,1],pred_E,col = 'deepskyblue3',lwd = 4)
  lines(data_plot_W[,1],pred_W,col = 'darkorange3',lwd = 4)
  leg_hei <- ymin + (ymax-ymin)/5
  legend(6,leg_hei, c("East", "West"), col = c('deepskyblue3','darkorange3'),lwd = c(4,4),cex = 0.8, bty = 'n')  
  dev.off()
}

plot_predictions_EW('figures/fig6a_EW_MAPreg.pdf',data_sites,'MAPreg','MidMemberAge',param_span,param_deg,ht,wd,'by ComLoc')
plot_predictions_EW('figures/fig6c_EW_MAPknn.pdf',data_sites,'MAPknn','MidMemberAge',param_span,param_deg,ht,wd,'by ComLoc')
plot_predictions_EW('figures/fig6b_EW_MAPreg_members.pdf',data_sites_members,'MAPreg','MidMemberAge',param_span,param_deg,ht,wd,'by Member')
plot_predictions_EW('figures/fig6d_EW_MAPknn_members.pdf',data_sites_members,'MAPknn','MidMemberAge',param_span,param_deg,ht,wd,'by Member')


# KBS zoom

ind <- which(data_KBS[,'MidMemberAge']>=0.9)
data_KBS <- data_KBS[ind,]
ind <- which(data_KBS[,'EastWest']=='E')
data_KBS <- data_KBS[ind,]
#filter out large members
un_mem <-  as.vector(unique(data_KBS[,'Member']))
ind_large <- c()
for (sk11 in 1:length(un_mem)){
  ind <- which(data_KBS[,'Member']==un_mem[sk11])
  un_ComLocs <- unique(data_KBS[ind,'ComLoc'])
  if (length(un_ComLocs)>=7){
    ind_large <- c(ind_large,ind)
  }
}
data_KBS <- data_KBS[ind_large,]
ind <- !duplicated(data_KBS$ComLoc)
data_KBS <- data_KBS[ind,]


plot_predictions_KBS <- function(plot_name_now,data_sites,fet_plot,fet_age,param_span,param_deg,ht,wd){
  if (substr(fet_plot,1,3)=='MAP'){
    ymin <- -100
    ymax <- 1700
    ylb <- 'MAP estimate'
    txlev <- 100
    txdif <- 0.075
  }else{
    if (substr(fet_plot,1,3)=='MAT'){
      ymin <- 18
      ymax <- 30  
      ylb <- 'MAT estimate'
      txlev <- 20
      txdif <- 0
    }
  }
  mn = 'regression by ComLoc'
  
  mid_point <- data_sites[,fet_age]
  target <- data_sites[,fet_plot]
  data_smoo <- as.data.frame(cbind(mid_point,target))
  data_plot = data.frame(mid_point = seq(0,10,0.1),target = NA)
  indE <- which(data_sites[,'EastWest']=='E')
  un_mmp_E <- unique(mid_point[indE])
  un_mmp_E <- un_mmp_E[order(un_mmp_E)]
  data_plot_E = data.frame(mid_point = un_mmp_E,target = NA)
  pred_E <- c()
  pred_E_ci <- c()
  for (sk8 in 1:length(un_mmp_E))
  {
    ind_temp <- which(mid_point[indE]==un_mmp_E[sk8])
    mea <- mean(target[indE[ind_temp]])
    sda <- 1.96*sd(target[indE[ind_temp]])/sqrt(length(target[indE[ind_temp]]))
    pred_E <- rbind(pred_E,mea)
    pred_E_ci <- rbind(pred_E_ci,sda)
  }
  pdf(plot_name_now,width = wd,height = ht)
  plot(NA,NA,xlim = c(1,3.5),ylim = c(ymin,ymax),xlab = 'age',ylab = ylb,main = mn)
  points(mid_point[indE],data_sites[indE,fet_plot],pch = 20,col = 'deepskyblue',cex = 0.9)
  lines(data_plot_E[,1],pred_E,col = 'deepskyblue3',lwd = 3)
  mem_txt <- c()
  for (sk12 in 1:dim(data_plot_E)[1])
  {
    ind_temp <- which(data_sites[,fet_age] == data_plot_E[sk12,1])
    mem <- as.vector(unique(data_sites[ind_temp,'Member']))
    if (length(mem)>1){print('!!problem with members in plotting')}
    mem_txt <- c(mem_txt,mem)
  }
  text(data_plot_E[,1]-txdif,txlev,mem_txt,cex = 0.7,srt=90, col='darkorange')
  # hack: we draw arrows but with very special "arrowheads" - for error bars
  arrows(data_plot_E[,1], pred_E-pred_E_ci, data_plot_E[,1], pred_E+pred_E_ci, length=0.05, angle=90, code=3,col = 'black', lwd = 1.5)
  dev.off()
}

plot_predictions_KBS('figures/fig9b_KBS_MAP.pdf',data_KBS,'MAPreg','MidMemberAge',param_span,param_deg,4,4)
plot_predictions_KBS('figures/fig9a_KBS_MAT.pdf',data_KBS,'MATreg','MidMemberAge',param_span,param_deg,4,4)


# species vs. specimen 
n <- dim(data_all)[1]
un_sites <- unique(data_all[,'ComLoc'])
data_species <- c()
param_smoothing_spec <- 4
for (sk in 1:length(un_sites))
{
  site_now <- un_sites[sk]
  ind <- which(data_all[,'ComLoc']==site_now)
  data_species <- rbind(data_species,cbind(as.vector(site_now),sum(data_all[ind,'unique_species_used']),length(ind),as.vector(data_all[ind[1],'Member'])))
}
colnames(data_species) <- c('ComLoc','species','specimen','Member')
mypalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
pdf('figures/fig7_species_specimen.pdf', height = 5.5, width = 5)
ind_KBS <- which(data_species[,'Member'] == 'KBS')
ind_UB <- which(data_species[,'Member'] == 'Upper Burgi')
ind_TB <- which(data_species[,'Member'] == 'Tulu Bor')
ind_OK <- which(data_species[,'Member'] == 'Okote')
#KBS
plot(as.numeric(data_species[ind_KBS,'species']),as.numeric(data_species[ind_KBS,'specimen']),xlab = 'no. species',ylab = 'no. specimen',log="y",col = mypalette[7], main = 'By ComLoc',lwd = 2) 
lines(smooth.spline(as.numeric(data_species[ind_KBS,'species']),as.numeric(data_species[ind_KBS,'specimen']), df = param_smoothing_spec), lty = 1, col = mypalette[7],lwd = 3)
#OK
points(as.numeric(data_species[ind_OK,'species']),as.numeric(data_species[ind_OK,'specimen']),col = mypalette[2], lwd = 2)
lines(smooth.spline(as.numeric(data_species[ind_OK,'species']),as.numeric(data_species[ind_OK,'specimen']), df = param_smoothing_spec), lty = 1, col = mypalette[2],lwd = 3)
#UB
points(as.numeric(data_species[ind_UB,'species']),as.numeric(data_species[ind_UB,'specimen']),col = mypalette[4], lwd = 2)
lines(smooth.spline(as.numeric(data_species[ind_UB,'species']),as.numeric(data_species[ind_UB,'specimen']), df = param_smoothing_spec), lty = 1, col = mypalette[4],lwd = 3)
#TB
points(as.numeric(data_species[ind_TB,'species']),as.numeric(data_species[ind_TB,'specimen']),col = mypalette[3], lwd = 2)
lines(smooth.spline(as.numeric(data_species[ind_TB,'species']),as.numeric(data_species[ind_TB,'specimen']), df = param_smoothing_spec), lty = 1, col = mypalette[3],lwd = 3)
legend('bottomright', c('Okote',"KBS",'Upper Burgi',"Tulu Bor"), col = c(mypalette[2],mypalette[7],mypalette[4],mypalette[3]),pch = c(1,1,1,1), lwd = c(2,2,2,2), cex = 0.8, bty = "n")
dev.off()


# species identified proportion plots
  make_data_prop <- function(data_all,ind)
  {
    Member <- as.vector(data_all[ind,'Member'])
    OrderFamily <- as.vector(data_all[ind,'Order'])
    Family <- as.vector(data_all[ind,'Family'])
    ind <- which(OrderFamily == 'Artiodactyla')
    OrderFamily[ind] <- Family[ind]
    data_proportions <- as.data.frame(cbind(Member,OrderFamily))
    data_proportions <- transform(data_proportions, Member = factor(Member, levels = c('Okote','KBS','Upper Burgi','Tulu Bor')))
    #manual cleaning
    ind_temp <- which(!is.na(data_proportions[,'Member']))
    data_proportions <- data_proportions[ind_temp,]
    ind_temp <- which(data_proportions[,'OrderFamily']!='Camelidae')
    data_proportions <- data_proportions[ind_temp,]
    ind_temp <- which(data_proportions[,'OrderFamily']!='Hippopotamidae/Suidae')
    data_proportions <- data_proportions[ind_temp,]
    ind_temp <- which(data_proportions[,'OrderFamily']!='Rodentia')
    data_proportions <- data_proportions[ind_temp,]
    ind_temp <- which(data_proportions[,'OrderFamily']!='Tubulidentata/Carnivora')
    data_proportions <- data_proportions[ind_temp,]
    ind_temp <- which(data_proportions[,'OrderFamily']!='Tubulidentata')
    data_proportions <- data_proportions[ind_temp,]
    ind_temp <- which(data_proportions[,'OrderFamily']!='Lagomorpha')
    data_proportions <- data_proportions[ind_temp,]
    #
    return(data_proportions)
  }
  ind <- which(!is.na(data_carnivors[,'ComLoc']))
  data_proportions <- make_data_prop(data_carnivors,ind)
  library(ggplot2)
  data_proportions$OrderFamily <- factor(data_proportions$OrderFamily, levels = rev(levels(data_proportions$OrderFamily)))
  plot_now <- ggplot(data_proportions,aes(x = Member,fill = OrderFamily,order = -as.numeric(OrderFamily))) + geom_bar(position = "fill") + scale_fill_manual(values=c("#CC79A7","#009E73","#0072B2","#56B4E9","#999999","#F0E442","#D55E00","#E69F00")) + ggtitle("By specimen")
  ggsave(plot_now, file='figures/fig9c_specimen.pdf',width=5, height=4)
  ind <- which(!is.na(data_carnivors[,'ComLoc']))
  ind2 <- which(data_carnivors[,'unique_species_used']==1)
  ind <- intersect(ind,ind2)
  data_proportions <- make_data_prop(data_carnivors,ind)
  data_proportions$OrderFamily <- factor(data_proportions$OrderFamily, levels = rev(levels(data_proportions$OrderFamily)))
  plot_now <- ggplot(data_proportions,aes(x = Member,fill = OrderFamily,order = -as.numeric(OrderFamily))) + geom_bar(position = "fill") + scale_fill_manual(values=c("#CC79A7","#009E73","#0072B2","#56B4E9","#999999","#F0E442","#D55E00","#E69F00")) + ggtitle("By identified species") 
  ggsave(plot_now, file='figures/fig9c_species.pdf',width=5, height=4)