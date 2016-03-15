# 2016 03 10 I.Zliobaite
# code for public acompanying manuscript 
# fitting modern regression models

input_file <- 'data/data_Africa_modern.csv'
plot_name_MAP <- 'out/fig_MAP.pdf'
plot_name_MAT <- 'out/fig_MAT.pdf'
out_data_syn <- 'out/data_syn.csv'

param_step <- 0.01 # step for generating synthetic visualization data
param_k <- 15 #how many neighbors

do_syn <- TRUE #write synthetic data and estimates

data_all <- read.csv(input_file, header = TRUE,sep = ',')
  
# function compute R2
compute_R2 <- function(ttrue,ppred)
{
  R2 <- 1 - sum((ttrue - ppred)^2)/sum((ttrue - mean(ttrue))^2)
  return(R2)
}

#recover coefficients for PLS
recover_coefficients <- function(fit_pls,comp,data_sites)
{
  model_now <- coef(fit_pls,ncomp = comp)
  n <- length(model_now)
  coefs <- model_now
  data <-data_sites*0
  data <- data[1,]
  for (sk in 1:n)
  {
    fet_now <- rownames(model_now)[sk]
    data1 <- data
    data1[fet_now] <- 1
    pred <- predict(fit_pls, newdata = data1)[,,comp] - predict(fit_pls, newdata = data)[,,comp]
    model_now[sk] <- pred
  }
  model_now <- c(predict(fit_pls, newdata = data)[,,comp],model_now)
  return(model_now)
}

# regression for precipitation
fml_MAP <- as.formula('MAP ~ HYP + LOP + HYP*LOP')
fit_MAP <- glm(fml_MAP,data=data_all) 
cf_MAP <- coefficients(fit_MAP)
R2_MAP <- compute_R2(data_all[,'MAP'],predict(fit_MAP))
print('Model for MAP')
print(cf_MAP)

pdf(plot_name_MAP,width=5,height=5)
plot(predict(fit_MAP),data_all[,'MAP'],xlab='estimated',ylab='true',type='p',main = paste('R2 =',round(R2_MAP,digits=2)))
dev.off()

#regression for temperature
library(pls)
comp = 1 #parameter for PLS regression, how many projected components to consider
fml_MAT <- as.formula('MAT ~ HYP + LOP')
fit_MAT <- mvr(fml_MAT, comp, method = 'svdpc', data = data_all,scale = TRUE)
cf_MAT <- recover_coefficients(fit_MAT,comp,data_all)
print('Model for MAT')
print(cf_MAT)

R2_MAT <- compute_R2(data_all[,'MAT'],predict(fit_MAT))
pdf(plot_name_MAT,width=5,height=5)
plot(predict(fit_MAT),data_all[,'MAT'],xlab='estimated',ylab='true',type='p',main = paste('R2 =',round(R2_MAT,digits=2)))
dev.off()


# visualize models
generate_syn_data <- function(step)
{
  data <- c()
  seqHYP <- seq(1,3,step)
  seqLOP <- seq(0,2,step)
  for (sk1 in 1:length(seqHYP))
  {
    for (sk2 in 1:length(seqLOP))
    {
      data <- rbind(data,c(seqHYP[sk1],seqLOP[sk2]))
    }
  }
  data <- cbind(data,data[,1]*data[,2])
  ind <- which((data[,1]-data[,2])>=0)
  data <- data[ind,]
  ind <- which((data[,1]-data[,2]-2)<=0)
  data <- data[ind,]
  colnames(data) <- c('HYP','LOP','HYPLOP')
  return(data)
}

make_predictions <- function(data_syn,cf)
{
  pred <- cf[1]
  for (sk in 2:length(cf))
  {
    pred <- pred + data_syn[,sk-1]*cf[sk]
  }
  return(pred)
}

make_predictions_knn <- function(data_syn,data_ref,k,fet_target)
{
  pred <- c()
  for (sk in 1:dim(data_syn)[1])
  {
    dd <- dist(rbind(data_syn[sk,1:2],data_ref[,c('HYP','LOP')]))
    dd <- dd[1:dim(data_ref)[1]]
    ord <- order(dd)
    #print(data_syn[sk,1:2])
    #print(data_ref[ord[1:k],fet_target])
    #print(dd[ord[1:kk]])
    pred <- rbind(pred,mean(data_ref[ord[1:k],fet_target]))
  }
  return(pred)
}

if (do_syn){
  data_syn <- generate_syn_data(param_step)
  reg_MAP <- make_predictions(data_syn,cf_MAP)
  reg_MAT <- make_predictions(data_syn,cf_MAT)
  knn_MAP <- as.vector(make_predictions_knn(data_syn,data_all,param_k,'MAP'))
  knn_MAT <- as.vector(make_predictions_knn(data_syn,data_all,param_k,'MAT'))
  data_syn <- cbind(data_syn,reg_MAP,reg_MAT,knn_MAP,knn_MAT)
  ind <- which(data_syn[,'reg_MAP']<0)
  data_syn[ind,'reg_MAP'] <- 0
  min_T <- min(c(data_syn[,'reg_MAT'],data_syn[,'knn_MAT']))
  max_T <- max(c(data_syn[,'reg_MAT'],data_syn[,'knn_MAT']))
  min_P <- min(c(data_syn[,'reg_MAP'],data_syn[,'knn_MAP']))
  max_P <- max(c(data_syn[,'reg_MAP'],data_syn[,'knn_MAP']))
  data_syn <- rbind(c(1.1,0,0,min_P,min_T,min_P,min_T),data_syn)
  data_syn <- rbind(c(1.2,0,0,max_P,max_T,max_P,max_T),data_syn)
  write.table(data_syn, file = out_data_syn, row.names = FALSE, col.names = TRUE, sep = ',', quote = FALSE)
}

data_syn <- read.csv(out_data_syn, header = TRUE,sep = ',')

#plotting

scatter_fill <- function (x, y, z,pale,xlim=c(min(x),max(x)),ylim=c(min(y),max(y)),zlim=c(min(z),max(z)),
                          nlevels = 20, plot.title, plot.axes, 
                          key.title, key.axes, asp = NA, xaxs = "i", 
                          yaxs = "i", las = 1, 
                          axes = TRUE, frame.plot = axes, ...) 
{
  mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
  on.exit(par(par.orig))
  w <- (3 + mar.orig[2L]) * par("csi") * 2.54
  layout(matrix(c(2, 1), ncol = 2L), widths = c(1, lcm(w)))
  par(las = las)
  mar <- mar.orig
  mar[4L] <- mar[2L]
  mar[2L] <- 1
  par(mar = mar)
  
  # choose colors to interpolate
  levels <- seq(zlim[1],zlim[2],length.out = nlevels)
  if (pale == 'T'){
    #col <- colorRampPalette(c('red','yellow',"dark green"))(nlevels)    
    col <- colorRampPalette(c('black','darkblue','cornflowerblue','khaki1','goldenrod1'))(nlevels)    
  }else{
    if (pale == 'P'){
      col <- colorRampPalette(c('red','yellow','dark green','black'))(nlevels)      
    }else{
      col <- colorRampPalette(c('black','gainsboro'))(nlevels)      
    }
  }
  
  colz <- col[cut(z,nlevels)]  
  #   
  plot.new()
  plot.window(xlim = c(0, 1), ylim = range(levels), xaxs = "i", yaxs = "i")
  
  rect(0, levels[-length(levels)], 1, levels[-1L],col=col,border=col) 
  if (missing(key.axes)) {if (axes){axis(4)}}
  else key.axes
  box()
  if (!missing(key.title)) 
    key.title
  mar <- mar.orig
  mar[4L] <- 1
  par(mar = mar)
  
  # points
  plot(x,y,type = "n",xaxt='n',yaxt='n',xlab="",ylab="",xlim=xlim,ylim=ylim,bty="n")
  points(x,y,col = colz,xaxt='n',yaxt='n',xlab="",ylab="",bty="n",...)
  
  ## options to make mapping more customizable
  
  if (missing(plot.axes)) {
    if (axes) {
      title(main = "", xlab = "", ylab = "")
      Axis(x, side = 1)
      Axis(y, side = 2)
    }
  }
  else plot.axes
  if (frame.plot) 
    box()
  if (missing(plot.title)) 
    title(...)
  else plot.title
  invisible()
}

#pdf('figures/fig10a_map_MAT.pdf')
#scatter_fill(data_all[,'lon'],data_all[,'lat'],data_all[,'MAT'],'T',nlevels=20,pch=".",cex=8,main = 'MAT',xlab = 'Longitude',ylab = 'Latitude')
#dev.off()

#pdf('figures/fig10b_map_MAP.pdf')
#scatter_fill(data_all[,'lon'],data_all[,'lat'],data_all[,'MAP'],'P',nlevels=20,pch=".",cex=8,main = 'MAP',xlab = 'Longitude',ylab = 'Latitude')
#dev.off()

#pdf('figures/fig10d_map_HYP.pdf')
#scatter_fill(data_all[,'lon'],data_all[,'lat'],data_all[,'HYP'],'k',nlevels=20,pch=".",cex=8,main = 'HYP',xlab = 'Longitude',ylab = 'Latitude')
#dev.off()

#pdf('figures/fig10c_map_LOP.pdf')
#scatter_fill(data_all[,'lon'],data_all[,'lat'],data_all[,'LOP'],'k',nlevels=20,pch=".",cex=8,main = 'LOP',xlab = 'Longitude',ylab = 'Latitude')
#dev.off()

pdf('figures/fig1a_reg_MAT.pdf')
scatter_fill(data_syn[,'HYP'],data_syn[,'LOP'],data_syn[,'reg_MAT'],'T',nlevels=40,pch=".",cex=8,main = 'regression MAT',xlab = 'HYP',ylab = 'LOP')
dev.off()

#pdf('figures/fig1b_reg_MAP.pdf')
#scatter_fill(data_syn[,'HYP'],data_syn[,'LOP'],data_syn[,'reg_MAP'],'P',nlevels=40,pch=".",cex=8,main = 'regression MAP',xlab = 'HYP',ylab = 'LOP')
#dev.off()
print('kar')
pdf('figures/fig1c_knn_MAT.pdf')
scatter_fill(data_syn[,'HYP'],data_syn[,'LOP'],data_syn[,'knn_MAT'],'T',nlevels=40,pch=".",cex=8,main = 'knn MAT',xlab = 'HYP',ylab = 'LOP')
dev.off()

#pdf('figures/fig1d_knn_MAP.pdf')
#scatter_fill(data_syn[,'HYP'],data_syn[,'LOP'],data_syn[,'knn_MAP'],'P',nlevels=40,pch=".",cex=8,main = 'knn MAP',xlab = 'HYP',ylab = 'LOP')
#dev.off()