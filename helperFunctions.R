##Author: Dr. Jiri Pyörälä (University of Helsinki, Finnish Geospatial Institute), 2022
##jiri.pyorala@helsinki.fi
##Research data and methods used in Pehkonen et al. (Manuscript) "How knot structure and stem geometry reflect ring properties in Norway spruce?"
#https://osf.io/984tk/S

organizeLogs <- function(df, vrbl){

  # Stem & knots ------------------------------------------------------------

  if(nrow(df) > 5){
    
    H = max(df$H[!is.na(df$H)])
    DBH = max(df$DBH[!is.na(df$DBH)])
    
    ##STEMCURVE
    d1 = data.frame("TreeID"=character(),"cl"=character(),"d"=numeric(),"h"=numeric(), "relh" = numeric(),
                    "H"=numeric(), "DBH"=numeric(),
                    "vrbl"=numeric())
    
    df3 = subset(df, df$StemSection==1)
    
    #loop through butt logs <- select other logs from same tree
    for(t in levels(as.factor(df3$TreeID))){
      d2=subset(df, df$TreeID==t)
      dbh1<-mean(d2[,"DBH"], na.rm=T)
      h1<-mean(d2[,"H"], na.rm=T)
      diam0<-d2[1,"buttD"]/10
      d20<- mean(d2[,"topD"]/10)
      
      d1 = rbind(d1, data.frame("TreeID"=t, 
                                "cl"=d2[1,"StemSection"],
                                "d"=d2[1,"buttD"]/10, 
                                "h"=d2[1,"buttH"]/100,
                                "relh"=(d2[1,"buttH"]/100)/h1,
                                "reld"=(d2[1,"buttD"]/10)/d20,
                                "H" = h1, "DBH" = dbh1,
                                "vrbl" = d2[1,vrbl]))
      
      #loop through all logs of that tree
      for(h in seq(1,nrow(d2))){
        
        d1 = rbind(d1, data.frame("TreeID"=t,"cl"=d2[h,"StemSection"],"d"=d2[h,"topD"]/10, "h"=d2[h,"topH"]/100, 
                                  "relh" = (d2[h,"topH"]/100)/h1, "reld"=(d2[h,"topD"]/10)/d20,
                                  "H" = h1, "DBH" = dbh1,
                                  "vrbl" = d2[h,vrbl]))
      }
      
      d1 = rbind(d1, data.frame("TreeID"=t,"cl"="0", "d"=0, "h"=h1, "relh"=1, "reld"=0,
                                "H" = h1, "DBH" = dbh1,
                                "vrbl" = 0))
    }
    
    d1[d1<0] <- 0
    d1[d1==Inf] <- 0
    d1[d1==0] <- 0.0001

  }
  
  return(d1)

}


rings2logs <- function(ringFile, stemFile) {
  
  rings = read.delim(ringFile, sep="\t", header=T)
  logs = read.delim(stemFile)
  
  rings = rings[rings$TreeID %in% logs$TreeID,]
  
  rings$ForestType = NA
  rings$DBH = NA
  rings$H = NA
  rings$buttH = NA
  rings$topH = NA
  rings$buttD = NA
  rings$topD = NA
  #rings$VolumeMeasuredXray = NA
  
  for(s in seq(1,nrow(rings))){
    stand = rings[s,"Stand"]
    tree = rings[s,"TreeID"]
    log = rings[s,"StemSection"]
    
    data = subset(logs, (logs$TreeID==tree & logs$StemSection==log))
    if(nrow(data) > 0){
      rings$ForestType[s] = data[1,"ForestType"]
      rings$DBH[s] = data[1,"DBH"]
      rings$H[s] = data[1,"H"]
      rings$buttH[s] = data[1,"buttH"]
      rings$topH[s] = data[1,"topH"]
      rings$buttD[s] = data[1,"buttD"]
      rings$topD[s] = data[1,"topD"]
    }
    else{
      data=subset(logs, (logs$TreeID==tree & logs$StemSection==log-1))
      rings$ForestType[s] = data[1,"ForestType"]
      rings$DBH[s] = data[1,"DBH"]
      rings$H[s] = data[1,"H"]
      rings$buttH[s] = data[1,"topH"]
      rings$topH[s] = data[1,"H"]*100
      rings$buttD[s] = data[1,"topD"]
      rings$topD[s] = 0
    }
  }
  
  
  return(rings)
  
}

knots2Rings <- function(ringdata, group1, stemdata, knotGradients){
  
  library(lme4)
  library(nlme)
  library(reshape2)
  library(ggplot2)
  library(timeDate)
  library(polynom)
  
  ##define stem and knottiness variables to be calculated to each ring
  ##See Table 2 in the manuscript
  ringdata$TreeAge=0  #Tree age each year, a
  ringdata$h=0 # Sample height, m
  ringdata$DBHr=0 #Diameter-at-breast-height each year, cm
  ringdata$Hr=0 #Tree height each year, m
  #ringdata$DBH #already in data
  #ringdata$H #already in data
  ringdata$D=0 #Final stem diameter at given height, cm
  ringdata$relD=0 #Final stem diameter at given height relative to DBH
  ringdata$relH= #Sample height relative to tree height each year
  ringdata$Vr=0 #volume of the stem each year
  ringdata$Vr_Hr_h=0 #volume of the stem above the location each year
  ringdata$SL0=0 #length of last shoot each year
  ringdata$SL1=0 #length of previous years shoot
  ringdata$KIr=0 #knot index just above the location
  ringdata$KSIr=0 #knot size index just above the location
  ringdata$KImax=0 #maximum knot index (annual, above h)
  ringdata$KImean=0 #mean knot index (annual, above h)
  ringdata$KSImax=0 #maximum knot size index (annual, above h)
  ringdata$KSImean=0 #mean knot size index (annual, above h)
  ringdata$hCrKI=0 #height of the maximum knot index
  ringdata$hCrKSI=0 #height of the maximum knot size index
  ringdata$relhCrKI=0 #Height of the maximum knot index, relative to tree height each year
  ringdata$relhCrKSI=0 #Height of the maximum knot size index, relative to tree height each year
  ringdata$CrlKI=0 #"crown length": Distance from the maximum knot index to tree-top each year
  ringdata$CrlKSI=0 #"crown length": Distance from the maximum knot size index to tree-top each year
  ringdata$relCrlKI=0 #"crown length" relative to tree height
  ringdata$relCrlKSI=0 #"crown length" relative to tree height
  ringdata$dCrKI=0 #distance from the maximum knot index : negative=below crown, positive=within crown
  ringdata$dCrKSI=0 #distance from the maximum knot size index : negative=below crown, positive=within crown
  ringdata$reldCrKI=0 #distance from crown base relative to tree height each year: negative=below crown, positive=within crown
  ringdata$reldCrKSI=0 #distance from crown base relative to tree height each year: negative=below crown, positive=within crown
  
  for(f in levels(factor(ringdata[,group1]))){
    
    knotGradient1 = knotGradients[[1]][as.character(f)]
    knotGradient2 = knotGradients[[2]][as.character(f)]
    knotGradient3 = knotGradients[[3]][as.character(f)]
    
    df1 = organizeLogs(df = stemdata[stemdata$Stand==f,], vrbl = vrbl)
    stemModel=lm(reld ~ poly(relh, 4), data=df1)
    polyStem=polynomial(as.list(stemModel$coefficients))
    
    H = max(df1$H[!is.na(df1$H)])
    DBH = max(df1$DBH[!is.na(df1$DBH)])
    D20 = DBH/predict(stemModel, newdata = data.frame(relh=1.3/H))
    D0 = predict(stemModel, newdata = data.frame(relh=0)) * D20
    
    ###Predict stem & knot features to rings:
    for(ring in 1:nrow(ringdata)){
      if(ringdata[ring, "Stand"]==f&ringdata[ring,"R"]>0){
        
        year = ringdata[ring,"Year"]
        TreeAge = year - min(ringdata[ringdata$TreeID == ringdata[ring, "TreeID"], "Year"])
        h = ringdata[ring,"buttH"]/100
        d=max(ringdata[c(ringdata$TreeID == ringdata[ring, "TreeID"]&ringdata$Year == year),"R"], na.rm=T)*0.2 #stem butt diameter at given year, cm
        if(d>D0){d<-D0}
        d20r=d/predict(stemModel, newdata=data.frame(relh=0))
        
        Hr = H * (d/D0) #tree height at given year, assume fixed H:D ratio
        if(Hr>H){Hr<-H}
        
        
        DBHr = predict(stemModel, newdata=data.frame(relh=(1.3/Hr))) * d20r #dbh at given year
        Rh = predict(stemModel, newdata=data.frame(relh=h/H)) * D20 #final stem diameter at height h
        
        knot1 = c()
        knot2 = c()
        knot3 = c()
        for(i in 1:(Hr*10)-1){
          knot1 = append(knot1, knotGradient1[[1]][i, (predict(stemModel, newdata=data.frame(relh=(i/10)/Hr)) * d20r)*10])
          knot2 = append(knot2, knotGradient2[[1]][i, (predict(stemModel, newdata=data.frame(relh=(i/10)/Hr)) * d20r)*10])
          knot3 = append(knot3, knotGradient3[[1]][i, (predict(stemModel, newdata=data.frame(relh=(i/10)/Hr)) * d20r)*10])
          }#Annual knot features of the ring year
        knot1 = knot1[!is.na(knot1)]
        knot2 = knot2[!is.na(knot2)]
        knot3 = knot3[!is.na(knot3)]
        
        ###ADD FEATURES TO RING DATA
        ringdata$Year[ring]=year
        ringdata$TreeAge[ring] = TreeAge
        ringdata$h[ring]=h
        ringdata$DBHr[ring]=DBHr
        ringdata$Hr[ring]=Hr
        ringdata$D[ring]=Rh
        ringdata$relD[ring] = Rh / DBHr
        ringdata$relH[ring] = h/Hr 
        ringdata$Vr[ring] = (d/100)^2 * Hr * (pi/4) * integral(polyStem, limits=c(0,1))^2
        ringdata$Vr_Hr_h[ring] = (d/100)^2 * Hr * (pi/4) * integral(polyStem, limits=c(h/Hr,1))^2 #(pi * r1^2 * (Hr-h))/3
        
        if(length(knot1)>1){
          ringdata$SL0[ring] = mean(knot3[length(knot3)]) #length of last shoot
          ringdata$SL1[ring] = mean(knot3[length(knot3)-1]) #length of previous years shoot
          ringdata$KIr[ring] <- knot1[h+1] #knot index just above the location
          ringdata$KSIr[ring] <- knot2[h+1] #knot size index just above the location
          ringdata$KImax[ring] <- max(knot1[(h+1):length(knot1)],na.rm=T) #maximum knot index (annual, above h)
          ringdata$KImean[ring] <- mean(knot1[(h+1):length(knot1)], na.rm=T) #mean knot index (annual, above h)
          ringdata$KSImax[ring] <- max(knot2[(h+1):length(knot2)],na.rm=T) #maximum knot size index (annual, above h)
          ringdata$KSImean[ring] <- mean(knot2[(h+1):length(knot2)], na.rm=T) #mean knot size index (annual, above h)
          ringdata$hCrKI[ring] <- min(which(knot1 == max(knot1, na.rm=T)))/10 #height of the maximum knot index, m
          ringdata$hCrKSI[ring] <- min(which(knot2 == max(knot2, na.rm=T)))/10 #height of the maximum knot size index, m
          ringdata$relhCrKI[ring] <- h/ringdata[ring,"hCrKI"] - 1 #relative position within crown: negative=below crown, positive=within crown
          ringdata$relhCrKSI[ring] <- h/ringdata[ring,"hCrKSI"] - 1 #relative position within crown: negative=below crown, positive=within crown
          ringdata$CrlKI[ring] <- Hr-ringdata[ring,"hCrKI"]
          ringdata$CrlKSI[ring] <- Hr-ringdata[ring,"hCrKSI"]
          ringdata$relCrlKI[ring] <- ringdata[ring,"CrlKI"]/Hr
          ringdata$relCrlKSI[ring] <- ringdata[ring,"CrlKSI"]/Hr
          ringdata$dCrKI[ring] <- h - ringdata[ring,"hCrKI"] #distance from crown base: negative=below crown, positive=within crown
          ringdata$dCrKSI[ring] <- h - ringdata[ring,"hCrKSI"] #distance from crown base: negative=below crown, positive=within crown
          ringdata$reldCrKI[ring] <- ringdata[ring,"dCrKI"]/Hr #distance from crown base: negative=below crown, positive=within crown
          ringdata$reldCrKSI[ring] <- ringdata[ring,"dCrKSI"]/Hr #distance from crown base: negative=below crown, positive=within crown
          }
        } 
      }
    }
    return(ringdata)
  }
