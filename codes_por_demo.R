
# useful R functions used for STM/SFC 
rm(list=ls())
require(maptools)
library(maps)
library(fields)
library(lattice)
library(ggplot2)
library(ape)
library(R.matlab)
library(rgdal)
library(class)
library(e1071)
library(classInt)
library(maptools)
library(spdep)
library(RColorBrewer)
require(rgeos) #for gIntersection
require(Hmisc) #xYplot Cbind
windowsFonts(Times=windowsFont("Times New Roman"))
trellis.par.set(grid.pars = list(fontfamily = "Times"))
lattice.options(default.theme = canonical.theme(color = FALSE))  



#++++++++++++++++++++++++ read data
if(TRUE){
  dirs <- "./NUTS III"        #tmp <- as(shape, "data.frame")
  shape <- readShapePoly(paste(dirs, "/NUTS_3.shp", sep=""))
  snam <- as.character(shape$NUTs_NUT3)
  n <- length(shape$NUTs_NUT3)
  W = read.table('neighbor.txt',head=F); W <- W - diag(n)
  coords <- coordinates(shape)
  
  dat <- read.csv('new input data 1991 2009.csv')
  # Freeman-Tukey transformation
  fac = 1e3
  y <- sqrt(fac*dat[,'births']/dat[,'women']) + sqrt(fac*(dat[,'births']+1)/dat[,'women'])  
  dat$y <- y
  
  dat0 <- readOGR(dsn="./NUTS III", "NUTS_3") 
  shape2 <- spTransform(dat0, CRS("+proj=longlat +datum=WGS84"))
  
  mat <- dat[,c('Year','regions','age.groups','y','f_all','births','women')]
  mat$Year <- mat$Year-1990
  mat <- mat[order(mat$age.groups),]
  
  grpnam <- paste0(c('15-19','20-24','25-29','30-34','35-39','40-44','45-50'),' yrs')
  mat$grp <- factor(grpnam[mat$age.groups])
  mat$times <- mat$Year+1990
  
  J <- length(grpnam); T <- nrow(mat)/J/n; 
  timetick <- unique(mat$times)[seq(1,T,len=4)]
  
  my.padding <- list(strip = .8, top.padding = 1, main.key.padding = 0, key.axis.padding = 0,
                     axis.ylab.padding = 0, ylab.key.padding  = 0, right.padding = 0,
                     left.padding = 0, key.sub.padding = 0, bottom.padding = 0,
                     axis.right = 0, axis.top = 0, key.right = 0)
  
  my.scales <- list(alternating = F, cex = 1,
                    x = list(at = timetick, tck = c(.2, 0), rot= c(54, 0), cex=.7,labels = timetick), y = list(cex=.75, tck=c(.2,0)))
  
  # customize the color for spatial mapping
  perm <- c(21,12,6,18,28,16,15,5,9,8,26,24,14,25,27,10,11,13,22,20,17,23,2,4,1,7,19,3) 
  # index in shape -> index in data
  IND <- order(perm) # index in data -> index in shape
  cvcols <- topo.colors(n)[rank(coords[,2])]
}  #read data


if(TRUE){
  #++++++++++++++ get the boundary for the two clusters in Castro et al. (2015)
  findSharedBoundary <- function(s1, s2){
    mat1 <- shape2@polygons[[which(snam==s1)]]@Polygons[[1]]@coords
    poly1 <- SpatialPolygons(list(Polygons(list(Polygon(coords=mat1)), ID=c("a") )))
    mat2 <- shape2@polygons[[which(snam==s2)]]@Polygons[[1]]@coords
    poly2 <- SpatialPolygons(list(Polygons(list(Polygon(coords=mat2)), ID=c("a") )))
    a <- gIntersection(poly1, poly2)
    # plot(shape2, axes=T)
    # lines(y=mat1[,2], x=mat1[,1], col='blue', lwd=2, lty=2)
    # lines(y=mat2[,2], x=mat2[,1], col='red', lwd=2)
    # lines(a, col='green2', lwd=2)
    return(a)
  }
  
  n1 <- length(s1s <- c('OESTE','LEZÍRIA DO TEJO','ALTO ALENTEJO'))
  n2 <- length(s2s <- c('PINHAL LITORAL','MÉDIO TEJO','PINHAL INTERIOR SUL','BEIRA INTERIOR SUL'))
  bds <- as.list(rep(NA, n1*n2))
  for(i in 1:n1) for(j in 1:n2)  bds[[(i-1)*n2+j]] <- findSharedBoundary(s1s[i], s2s[j])
  
  # plot(shape2, axes=T)
  # for(i in 1:length(bds)) lines(bds[[i]], col='red', lwd=3)
  
  # check completeness
  # for(i in 1:28) lines(shape2@polygons[[which(snam==snam[i])]]@Polygons[[1]]@coords, col='green')
  
  # seems the polygon for OESTE (bds[[1]]) was problematic. Need manual correction
  mat1 <- shape2@polygons[[which(snam=='PINHAL LITORAL')]]@Polygons[[1]]@coords
  inds <- which(mat1[,1]< -8.759377 & mat1[,2]<39.75)
  # lines(y=mat1[inds,2], x=mat1[inds,1], col='red', lwd=3)
  bds[[1]] <- mat1[inds,]
  
  coords2 <- coordinates(shape2)
  # palette(topo.colors(n)); par(mar=c(1.5,1.5,0,0)+.1, mgp = c(1.2,.3,0), 
  #                              cex.main=1, cex.lab=.1, tcl=-0.2, cex.axis=.1, srt=0)
  # l1 <- list("SpatialPolygonsRescale", layout.north.arrow(type=1), offset = c(-6.5,37.2), scale=.5)
  # l2 = list("sp.text", coords2, shape$NUTs_NUT3, cex=.5, srt=5)
  # spplot(shape2, axes=T, col.regions=adjustcolor(cvcols,alpha.f=0.7), colorkey=F, 
  #        scales=list(draw=TRUE, cex=.7), xlim=c(-10,-6), sp.layout=list(l1,l2), col="gray70")
  # dev.off()
  # # note: need to crop the PDF. Hard to control the margin for spplot pdf output. 
  
  # # without using spplot. Just use simple plot
  # version to indicate the regions as legend
  # jpeg(paste('./tex/STmapnew.jpg',sep=''),quality=100, 
  #      height=4, width=4, units='in', pointsize=9, res=600)
  pdf(paste0('./tex/STmapnew.pdf'), height=4, width=4.4, pointsize=11, paper='special')
  layout(mat=matrix(1:2,nrow=1), widths = c(.63,.37))
  par(oma=c(0,0,0,0))
  par(mar=c(1.4,1.4,0,0)+.2,mgp=c(1,.3,0), tck=-0.01, 
      cex.axis=1.1, cex.lab=1, cex.main=1)
  plot(shape2, axes=T, col=adjustcolor(cvcols,alpha.f=0.6), border='gray60')
  # add boundary for the two clusters
  lines(y=bds[[1]][,2], x=bds[[1]][,1], col='red', lwd=3)
  for(i in 2:length(bds)) lines(bds[[i]], col='red', lwd=3)
  text(coords2, labels=1:length(snam), srt=0, cex=1)
  
  par(mar=c(0,0,0,0), xpd=T)
  plot(1, type = "n", axes=FALSE, xlab="", ylab="")
  legend('left',legend=paste(1:length(snam),snam),cex=.7,
         inset=c(-.07,0), title='NUTS III region', bty='n')
  par(xpd=F)
  dev.off()
  
}  #STmapnew.jpg


if(TRUE){
  #++++++++++++++++++++++++ for the figure STdata
  # xyplot(y~times|grp, data=mat, layout=c(7,1))
  FIG <- xyplot(y ~ times|grp, group=regions, data=mat, 
                xlab='', ylab='', ylim=c(0, 28), layout=c(7,1), #Trasformed Fertility Rates
                panel = function(x,y,...) { 
                  panel.refline(v=timetick, h=seq(0,25,by=5), lty=2, col='gray60') 
                  panel.xyplot(x,y,..., type='b', pch=16, lty=1, cex=.5, col=cvcols[perm]) #'gray50'
                  #panel.lines(x0*100, meancuves[,useperm[which(Ms==mean(y))]], lwd=2, col='red')
                },
                scales = my.scales, par.settings = list(layout.heights = my.padding, layout.widths = my.padding), strip = strip.custom(par.strip.text = list(cex=.8, font=3))
  )
  
  # postscript(file=paste('./tex/STdata.ps',sep=''),pointsize=12,width=7.2,height=2.2,horizontal=F,paper='special')
  # jpeg(paste0('./tex/STdata.jpg'),quality=100, 
  #      height=2.2,width=7.2,units='in', pointsize=28, res=600)
  pdf(paste0('./tex/STdata.pdf'), height=2.2, width=7.2, pointsize=11, paper='special')
  print(FIG)
  dev.off()
  
}  #STdata.jpg


if(TRUE){
  # ================================ map Portugal data
  # -------------------------------   results mapping
  #dat <- readOGR(".", "NUTS_3")
  ##dat1 <- spTransform(dat, CRS("+init=epsg:4267"))  #27700
  ##ID <- formatC(dat1$ID, width = 2, flag = "0")
  ##dat <- spChFIDs(dat1, ID)
  #bluepal <- colorRampPalette(c("azure1", "steelblue4"))
  #perm <- c(21,12,6,18,28,16,15,5,9,8,26,24,14,25,27,10,11,13,22,20,17,23,2,4,1,7,19,3)
  #IND <- order(perm)
  ## dirs <- "./NUTS III"        #tmp <- as(shape, "data.frame")
  ##shape <- readShapePoly(paste(dirs, "/NUTS_3.shp", sep=""))
  ##shape$NUTs_NUT3[perm]
  #labs <- readMat('./labsC.mat')$labsC
  #labs <- as.numeric(labs)
  ## palette(topo.colors(length(labs)))
  #dat$labs <- labs[IND]
  #spplot(dat, c("labs"), at=c(1:max(labs)), col.regions = bluepal(16),names.attr=c("Cluster"))
  
  # use my own function
  perm <- c(21,12,6,18,28,16,15,5,9,8,26,24,14,25,27,10,11,13,22,20,17,23,2,4,1,7,19,3)
  IND <- order(perm)
  dirs <- "./NUTS III"        #tmp <- as(shape, "data.frame")
  shape <- readShapePoly(paste(dirs, "/NUTS_3.shp", sep=""))
  shape$NUTs_NUT3[perm]# for figures in the manuscript
  coords <- coordinates(shape)
  tmp <- readMat('./SS.mat')
  tmpC <- tmp$labsC
  plot(tmp$SensitivityC, tmp$SpecificityC, type='p', pch=16, xlim=c(0.5,1), ylim=c(0.5,1))
  dev.off()
  
  
  ###################################################   produce ccAll.jpg
  
  #jpeg(paste0('./tex/ccAll.jpg'),quality=100, height=2000,width=6000, pointsize=25, res=300)
  pdf(paste0('./tex/ccAll.pdf'), height=2.2, width=6.6, pointsize=11, paper='special')
  facs <- c(4,4,5,4,5,5,2)
  agenam <- c('15-19','20-24','25-29','30-34','35-39','40-44','45-49')
  
  par(mfrow=c(1,7), mar=c(0,0,0,0), srt=0)  #   srt=15
  
  ### for loop
  for(nvar in 1:7){
    
    tmp <- readMat('dat_portugal_st.mat')
    # Y <- tmp$Y[1:16,,nvar]
    Y <- tmp$Bir[1:16,,nvar]/tmp$E[1:16,,nvar]   # fertility rates
    ntime <- nrow(Y); N <- ncol(Y)
    tmp <- readMat(paste('./useparas',nvar,'.mat', sep=''))
    tmpC <- tmp$labsC
    d <- max(tmpC)
    palette(rainbow(d))
    cols <- 1:d
    plot(shape, col=cols[tmpC[IND]], border='gray60')  #, main=agenam[nvar]
    
    for(r in 1:d){
      if(r==2 & nvar==2){ ccs <- coords[tmpC[IND]==r, ]; ccs <- ccs[1,] }
      else if(r==2 & nvar==3){ ccs <- coords[tmpC[IND]==r, ]; ccs <- ccs[1,] }
      else if(r==7 & nvar==1){ ccs <- coords[tmpC[IND]==r, ]; ccs <- ccs[1,] }
      else if(r==1 & nvar==6){ ccs <- coords[tmpC[IND]==r, ]; ccs <- ccs[1,] }
      else{
        if(sum(tmpC[IND]==r)>1) ccs <- colMeans(coords[tmpC[IND]==r, ])
        if(sum(tmpC[IND]==r)==1) ccs <- coords[tmpC[IND]==r, ]
      }
      text(x=ccs[1], y=ccs[2], labels=r, col='black', cex=1.3, font=2)#paste('cluster',r)
    }
    legend('top', legend=agenam[nvar], bty='n', cex=1.5, inset=c(0,-.03))
  } #################################################################
  
  dev.off()
  
  # check inverse Freeman-Tukey transformation
  #x1 <- seq(1e-2,,10,len=1e2)
  #y1 <- sqrt(x1) + sqrt(x1+1)
  #plot(x1,y1,type='l')
  #y2 <- seq(1,60,len=1e2)
  #x2 <- 0.5*((y2^4+1)/(2*y2^2)-1)
  #points(x2,y2,pch=16,col='red')
  
}  #ccAll.jpg


if(TRUE){
  ###################################################   produce clusterbeta.jpg
  facs <- c(4,4,5,4,5,5,2)
  agenam <- c('15-19','20-24','25-29','30-34','35-39','40-44','45-49')
  Means <- list(7)
  getFR <- function(x){ n <- nrow(x); return((x[n,]-x[1,])/x[1,]) }
  
  for(nvar in 1:7){
    
    tmp <- readMat('dat_portugal_st.mat')
    # Y <- tmp$Y[1:16,,nvar]
    Y <- tmp$Bir[1:16,,nvar]/tmp$E[1:16,,nvar]   # fertility rates
    ntime <- nrow(Y); N <- ncol(Y)
    tmp <- readMat(paste0('./useparas',nvar,'.mat', sep=''))
    tmpC <- tmp$labsC
    max(tmpC)
    palette(rainbow(max(tmpC)))
    mes <- matrix(tmp$mat2[,1], nrow=ntime)
    lbs <- matrix(tmp$mat2[,2], nrow=ntime)
    ubs <- matrix(tmp$mat2[,3], nrow=ntime)
    Means[[nvar]] <- mes
    mylab <- paste('VR = ',round(getFR(mes),4)*100,'%',sep='')
    
    dat <- matrix(Y, length(Y), 1); dat <- as.data.frame(dat); names(dat) <- "y"
    dat$grp <- as.factor(paste(agenam[nvar],': cluster ',rep(tmpC, each=ntime), sep=''))
    #dat <- cbind(dat, ym)
    dat$x <- rep(1:ntime, N)+1990
    dat$sites <- rep(1:N,each=ntime)
    
    timetick <- unique(dat$x)[seq(1,ntime,len=4)]
    my.padding <- list(strip = .8, top.padding = 0.4, main.key.padding = 0, key.axis.padding = 0,
                       axis.ylab.padding = 0, ylab.key.padding  = 0, right.padding = 0,
                       left.padding = 0, key.sub.padding = 0, bottom.padding = 0,
                       axis.right = 0, axis.top = 0, key.right = 0)
    
    my.scales <- list(alternating = F, cex = 1,
                      x = list(at = timetick, tck = c(.2, 0), rot= c(54, 0), 
                               cex=.7,labels = timetick), y = list(cex=.75, tck=c(.2,0)))
    
    Ms <- as.numeric(aggregate(dat$y,list(dat$grp),mean)[,2])
    
    # loc1 <- 8
    loc2 <- max(dat$y)*1.01
    widths <- 5.5201; heights <- 3.3667
    
    layout1 <- c(max(tmpC),1)
    if(nvar==1) layout1 <- c(4,2)                     
    FIG <- 
      xyplot(y ~ x|grp, group=sites, data=dat, 
             #xlab='Year', ylab='Trasformed Fertility Rates', 
             xlab='', ylab='',
             layout=layout1,  #subscripts = TRUE,
             scales = my.scales, 
             par.settings = list(layout.heights = my.padding, 
                                 layout.widths = my.padding), 
             strip = strip.custom(par.strip.text = list(cex=.8, font=3)), 
             panel = function(x,y,...) { 
               panel.refline(v=timetick, h=myhs, lty=2, col='gray70') 
               panel.xyplot(x,y,..., type='b', pch=16, lty=1, cex=.5, col=tmpC) #'gray50'   
               r <- which(Ms==mean(y))
               panel.lines(1990+(1:ntime),mes[,r], lwd=2, lty=1, col='black')
               panel.lines(1990+(1:ntime),lbs[,r], lwd=1, lty=1, col='black')
               panel.lines(1990+(1:ntime),ubs[,r], lwd=1, lty=1, col='black')
               #panel.lines(x0*100, meancuves[,useperm[which(Ms==mean(y))]], lwd=2, col='red')
               panel.text(mean(x),loc2,labels=mylab[r],cex=.75,col='black')
             }
      )
    
    # jpeg(paste('./tex/clusterbeta',nvar,'.jpg',sep=''),quality=100, height=heights,width=widths, 
    #      pointsize=35, res=600, units='in')
    # postscript(file=paste('../../POR/tex/clusterbeta',nvar,'.ps',sep=''),pointsize=12,width=widths,height=heights,horizontal=F,paper='special')    
    pdf(paste0('./tex/clusterbeta',nvar,'.pdf'), height=heights, width=widths, 
        pointsize=11, paper='special')
    lattice.options(default.theme = canonical.theme(color = FALSE))  
    print(FIG, default.theme = canonical.theme(color = FALSE))
    dev.off()
    
  } #################################################################
  
  # calculating the variation rate
  lapply(Means, getFR)
  
  # save it as Excel
  nam <- mmat <- frs <- numeric()
  for(nvar in 1:length(Means)){
    mmat <- cbind(mmat, Means[[nvar]])
    nam <- c(nam, paste(paste('group ',nvar,',',sep=''),'custer',1:ncol(Means[[nvar]])))
    frs <- c(frs, getFR(Means[[nvar]]))
  }
  mmat <- rbind(mmat, as.numeric(frs))
  mmat <- as.data.frame(mmat); names(mmat) <- nam 
  write.csv(mmat,file='meancurve.csv',row.names=F)
  
}  #clusterbeta*.jpg


if(TRUE){ #pred.jpg 
  tmp <-  readMat('./paper3/Rmat.mat')
  dat <- tmp$Rmat 
  indj0 <- tmp$indj0 
  indt0 <- tmp$indt0 
  
  for(j in 1:J){
    
    mat <- cbind(indt0, dat); mat <- mat[indj0==j,]
    N <- nrow(mat)
    mat2 <- mat[,1:2]; mat2 <- cbind(mat2, NA, NA)
    mat2 <- rbind(mat2, mat[,c(1,3:5)])
    mat2 <- rbind(mat2, mat[,c(1,6:8)])
    mat2 <- rbind(mat2, mat[,c(1,9:11)])
    mat2 <- as.data.frame(mat2)
    names(mat2) <- c('year','means','lower','upper')
    mat2$group <- rep(levs <- c('Observed','1','2','3'), each=N)
    mat2$group <- factor(mat2$group, levels=levs) 
    mat2$site <- rep(rep(shape$NUTs_NUT3[perm], 3), 4)
    mat2$year <- mat2$year + rep(c(0,.1,-.1,.2), each=N)+1990
    
    lattice.options(default.theme = canonical.theme(color = FALSE))
    my.padding <- list(strip = .7, top.padding = 0, main.key.padding = 0, 
                       key.axis.padding = 0, axis.ylab.padding = 0, ylab.key.padding  = 0, 
                       right.padding = 0, left.padding = 0, key.sub.padding = 0, 
                       bottom.padding = 0, axis.right = 0, key.right = 0)
    
    
    FIG <- xYplot(Cbind(means, lower, upper) ~ year|site, groups=group, data=mat2, 
                  ylab='', xlab='', pch = c(16,18,18,18), cex = c(1.2,1,1,1)*.8, 
                  col=c('black','red','green','blue'), layout=c(6,5), 
                  ylim=range(mat2[,c('means','lower','upper')])+c(-1,1)*.05, 
                  strip = strip.custom(par.strip.text = list(cex=.5, font=1))
    )
    
    FIG <- update(
      FIG, 
      scales = list(alternating=F, cex=.8, abbreviate=F, 
                    y=list(cex=.7, tck=c(.3,0)), 
                    x = list(at=c(1990+(17:19)), cex=.7, tck=c(.3, 0), rot=c(26,0)) ),  
      par.settings=list(layout.heights=my.padding, layout.widths=my.padding), #relation = "free",
      key = list(text=list(c("Observed", "Model (1)", "Model (2)","Model (3)")), 
                 border=F, corner=c(1.02, 1), cex=.75,
                 lines=list(lwd=c(2,1,1,1), col=c('black','red','green','blue'), lty=c(1,1,1,1)),
                 points=list(cex=c(1.2,1,1,1), pch=c(16,18,18,18),col=c('black','red','green','blue'))   
      )
    )
    
    # postscript(paste0('./fig',j,'.eps'), width = 7.5, height = 6, 
    #            onefile=T, pointsize=12, horizontal=F,paper='special')
    fac <- 1.3
    # jpeg(paste0('./tex/pred',j,'.jpg'), quality=100, height=3.12*fac, width=5.2*fac, 
    #      pointsize=1, res=600, units='in')  #, family='Times'
    pdf(paste0('./tex/pred',j,'.pdf'), height=3.12*fac, width=5.2*fac, 
        pointsize=11, paper='special')
    print(FIG)
    dev.off()
    
  }
}  #pred.jpg


#++++++++++++++++++++++++ additional

if(TRUE){
  #############################   visualize simulated data, save for saTScan
  ## for simulation data
  sobj <- readMat('simuY1sig2_06.mat')  #sfc_demo.m, ID=1, pa.sigma2=0.6
  tmp <- sobj$simuY #simulated data T-by-N matrix
  nsite <- ncol(tmp);   ntime <- nrow(tmp)
  t0 <- 1991:2009
  sdat <- data.frame(s=rep(1:nsite, each=ntime), t=t0[rep(1:ntime, nsite)], y=as.vector(tmp))
  write.csv(file='sdat_portugal_st.csv', sdat, row.names=F) #for SaTScan
  
  # prepare the Non-Euclidian Neighbors File for SaTScan
  W0 <- matrix(NA,n,1+max(rowSums(W)))
  for(i in 1:n){
    inds <- which(W[i,]==1); W0[i,1:(length(inds)+1)] <- c(i, inds)
  }
  write.table(W0, file='W0.txt', na='', row.names=FALSE, col.names=FALSE)
  
  # in SaTScan: choose Type of Analysis = Space-Time, Probability Model = Normal
  # Scan For Area with = High or Low values. Column Output Format: check first four
  
  # run SaTScan and load results to summarize in R
  tmp <- read.csv2('satOut-ID1Sigma2_06_scanHighNLow.txt')
  for(i in 1:nrow(tmp)) if(tmp[i,]=='CLUSTERS DETECTED') i1 <- i
  for(i in 1:nrow(tmp)) if(tmp[i,]=='PARAMETER SETTINGS') i2 <- i-2
  k <- i1+1;  labS <- rep(NA, n)
  while(k < i2){
    r <- as.numeric(substr(tmp[k,], 1,1))
    indr <- (substr(tmp[k,], gregexpr(':', tmp[k,])[[1]][1]+1, 1000L))
    indr <- as.numeric(strsplit(indr,split=", ",fixed=TRUE)[[1]])
    labS[indr] <- r
    k <- k+9  #fixed gap (#lines) between clusters in the SaTScan output
  }
  print(labS) 
  which(is.na(labS))  #1  2  9 10  #some regions didn't receive any cluster labels
  
  # visualize simulated data
  pdf(paste0('./tex/figSimu.pdf'), height=3, width=6.2, pointsize=11, paper='special')
  par(mfrow=c(1,4), mar=c(2,2,1,0)+.7,mgp=c(1.2,.3,0), tck=-0.02, 
      cex.axis=.9, cex.lab=1.1, cex.main=1.2, font.main=2)
  sig2vec <- c(.2,.6,1)
  t0 <- 1991:2009
  cols <- c('#e41a1c','#377eb8','#4daf4a','#984ea3')
  colsL <- adjustcolor(cols, alpha.f=.5)
  for(i in 1:length(sig2vec)){
    sig2 <- sig2vec[i];  ss <- as.character(sig2*10); ss <- paste0('0',ss) #if(nchar(ss)==1) 
    # sobj <- readMat(paste0('simuY1sig2_',ss,'.mat'))  #sfc_demo.m, ID=1, pa.sigma2=0.6
    if(i==1) sobj <- readMat('simuY2sig2_02.mat') 
    if(i==2) sobj <- readMat('simuY1sig2_06.mat')
    if(i==3) sobj <- readMat('simuY3sig2_010.mat')
    if(i==1){
      labs <- as.numeric(sobj$labs[IND])
      centers <- as.numeric(sobj$centers)
      indc <- rep(NA,n); indc[centers] <- 1:4; indc <- indc[IND]; ind1 <- which(!is.na(indc))
      plot(shape2, axes=F, col=colsL[labs], border='gray60')  #gray60
      text(x=coords2[ind1,1], y=coords2[ind1,2], 
           labels=paste0('C',indc[ind1]), col='black', cex=1.3, font=2)
    }
    par(xaxt='n')
    matplot(sobj$simuY, type='n', xlab='Year', ylab='Transformed fertility rate',
            main=substitute(paste(sigma^2,'=',sig2), list(sig2=sig2)), ylim=c(10,24))
    d <- max(sobj$labs)
    for(r in c(4,2,3,1)){  #bring C3, C1 to foreground
      indr <- which(sobj$labs==r)
      matplot(sobj$simuY[,indr], col=colsL[r], lwd=1, lty=1, type='l', add=T)
      lines(sobj$meanCurve[,r], lwd=2, col=cols[r])
    }
    legend('topright', legend=paste0('C',1:d), lwd=2, col=cols, bty='n', inset=c(.02,0))
    par(xaxt='s'); i0 <- seq(1,length(t0),by=5); axis(side=1, at=i0, labels=t0[i0])
  }
  dev.off()
  
}

# not run







