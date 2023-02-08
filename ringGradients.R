##Author: Jiri Pyörälä (University of Helsinki, Finnish Geospatial Institute), 2022
##Part of methods used in Pehkonen et al. (Manuscript) "How knot structure and stem geometry reflect ring properties in Picea abies H. Karst?"


##input: modeling data with ring properties and stem and knottiness features, and ring model
##loop stands
 ##predict ring properties to stems

ringGradients <- function(dfx, model_data, group1){
  
  library(nlme)
  library(reshape2)
  library(ggplot2)

  results = data.frame(group1=character(), vrbl = character(), "R2"=double(), "Sigma"=double(), "DF"=double())

  vrbl1 = na.exclude(model_data[,"CA"])
  vrbl1 = vrbl1[vrbl1 != Inf]
  vrbl1[vrbl1<0] <- 0
  cols1 = seq(quantile(vrbl1, probs=0.001), quantile(vrbl1, probs=0.95), length=48)
  
  vrbl2 = na.exclude(model_data[,"RA"])
  vrbl2 = vrbl2[vrbl2 != Inf]
  vrbl2[vrbl2<0] <- 0
  cols2 = seq(quantile(vrbl2, probs=0.0001), quantile(vrbl2, probs=0.999), length=48)
  
  vrbl3 = na.exclude(model_data[,"LWP"])
  vrbl3 = vrbl3[vrbl3 != Inf]
  vrbl3[vrbl3<0] <- 0
  cols3 = seq(quantile(vrbl3, probs=0.0001), quantile(vrbl3, probs=0.999), length=48)
  
  vrbl4 = na.exclude(model_data[,"RD"])
  vrbl4 = vrbl4[vrbl4 != Inf]
  vrbl4[vrbl4<0] <- 0
  cols4 = seq(quantile(vrbl4, probs=0.0001), quantile(vrbl4, probs=0.999), length=48)

  for(f in levels(factor(dfx[,group1]))){

    df2 = subset(dfx,dfx[,group1]==f)
    
    CA_MM = lm(formula(CA_MM1), data=model_data[model_data$Stand==f,])
    RA_MM = lm(formula(RA_MM1), data=model_data[model_data$Stand==f,])
    LWP_MM = lm(formula(LWP_MM1), data=model_data[model_data$Stand==f,])
    RD_MM = lm(formula(RD_MM1), data=model_data[model_data$Stand==f,])
    
    if(nrow(df2) > 5){
    
      H = max(df2$H[!is.na(df2$H)])
      DBH = max(df2$DBH[!is.na(df2$DBH)])
      
      ##STEMCURVE
      d1 = data.frame("TreeID"=character(),"d"=numeric(),
                      "h"=numeric(), "relh" = numeric(),
                      "H"=numeric(), "DBH"=numeric())
      
      ##Loop through butt logs in x-ray data 
      df3 = subset(df2, df2$StemSection==1)
        
      #loop through trees of that log type <- select other logs from same tree
      for(t in levels(as.factor(df3$TreeID))){
        d2=subset(dfx, dfx$TreeID==t)
        dbh1<-mean(d2[,"DBH"])*10
        h1<-mean(d2[,"H"])
        diam0<-d2[1,"DiameterMeasuredBottomAverage"]
        d20<- mean(d2[,"DiameterSorted"]/10)
        
        d1 = rbind(d1, data.frame("TreeID"=t, 
                                  "cl"=d2[1,"StemSection"],
                                  "d"=d2[1,"DiameterMeasuredBottomAverage"]/10, 
                                  "h"=d2[1,"buttH"]/100,
                                  "relh"=(d2[1,"buttH"]/100)/h1,
                                  "reld"=(d2[1,"DiameterMeasuredBottomAverage"]/10)/d20,
                                  "H" = h1, "DBH" = dbh1))
        
        #loop through all logs of that tree
        for(h in seq(1,nrow(d2))){
          
          d1 = rbind(d1, data.frame("TreeID"=t,"cl"=d2[h,"StemSection"],"d"=d2[h,"DiameterSorted"]/10, "h"=d2[h,"topH"]/100, 
                                    "relh" = (d2[h,"topH"]/100)/h1, "reld"=(d2[h,"DiameterSorted"]/10)/d20,
                                    "H" = h1, "DBH" = dbh1))
        }
        
        d1 = rbind(d1, data.frame("TreeID"=t,"cl"="0", "d"=0, "h"=h1, "relh"=1, "reld"=0,
                                  "H" = h1, "DBH" = dbh1))
      }
      
      d1[d1<0] <- 0
      d1[d1==Inf] <- 0
      d1[d1==0] <- 0.001
      stemModel=lm(reld ~ poly(relh, 4), data=d1)
      #summary(lmStem)
      D20 = DBH/predict(stemModel, newdata = data.frame(relh=1.3/H))
      D0 = predict(stemModel, newdata = data.frame(relh=0)) * D20
      
      ##RING PROPERTIES
      d1 = subset(model_data, model_data[,group1]==f)
      
      #TreeAge = max(d1$TreeAge)
      
      #loop trees
      if(nrow(d1)>5){
        
        #lmRing = lm(ringModel, data = d1)
        #sm=summary(lmRing)
        #results = rbind(results, data.frame(group1=f, vrbl = ringVrbl, "R2"=sm$r.squared, "Sigma"=sm$sigma, "DF"=sm$df[2]))
        
        CAGradient = matrix(data=NA_real_, nrow = 350, ncol = 650, dimnames = list(h = seq(1,350), d = seq(1,650)))
        RAGradient = matrix(data=NA_real_, nrow = 350, ncol = 650, dimnames = list(h = seq(1,350), d = seq(1,650)))
        LWPGradient = matrix(data=NA_real_, nrow = 350, ncol = 650, dimnames = list(h = seq(1,350), d = seq(1,650)))
        RDGradient = matrix(data=NA_real_, nrow = 350, ncol = 650, dimnames = list(h = seq(1,350), d = seq(1,650)))
      
        ###Predict wood properties to d x h matrix:
        for(d in seq(1,650)){ #loop through d in "stem", 1 mm intervals

          for(h in seq(1,350)){ #loop through heights
            
            h=h/10
            
            Hr = H * (d/(D0*10)) #tree height at given year, assume fixed H:D ratio
  
            d20 = d/predict(stemModel, newdata=data.frame(relh=0)) #diam at 20-% height, mm
            #DBHr = predict(stemModel, newdata=data.frame(relh=1.3/Hr)) * d20
            
            r = predict(stemModel, newdata=data.frame(relh=h/Hr))*d20 /2 #r at h, mm
            Rh = predict(stemModel, newdata=data.frame(relh=h/H)) * D20 * 10 #final stem diameter at height h, mm
            relH = h/Hr
            #relD = r/(Rh/2)

            if(c(r*2<=Rh&&r>2)){
              d_pred = data.frame("Stand"=f, "h"=h, "H"=H, "R"=r, "DBH"=DBH, "Hr"=Hr, "relH"=relH)
              
              CA = exp(predict(CA_MM, newdata = d_pred))
              d_pred$CA = CA
              
              RA = exp(predict(RA_MM, newdata = d_pred))
              d_pred$RA = RA
              
              LWP = exp(predict(LWP_MM, newdata = d_pred))
              d_pred$LWP = LWP
              
              RD = exp(predict(RD_MM, newdata = d_pred))
              
              h=h*10
              r=r*2
              CAGradient[h,r] <- CA
              RAGradient[h,r] <- RA
              LWPGradient[h,r] <- LWP
              RDGradient[h,r] <- RD
              
              if(is.na(CAGradient[h,r-1])&&r-1>0){
                CAGradient[h,r-1] <- CA
                RAGradient[h,r-1] <- RA
                LWPGradient[h,r-1] <- LWP
                RDGradient[h,r-1] <- RD}
              if(is.na(CAGradient[h,r-2])&&r-2>0){
                CAGradient[h,r-2] <- CA
                RAGradient[h,r-2] <- RA
                LWPGradient[h,r-2] <- LWP
                RDGradient[h,r-2] <- RD}
              
            }
          }
        }
      }
        
        
      CAGradient[CAGradient<0]<-0
      
      CAGradient[CAGradient<=quantile(vrbl1, probs=0.001)]<-quantile(vrbl1, probs=0.001)
      CAGradient[CAGradient>quantile(vrbl1, probs=0.95)] <- quantile(vrbl1, probs=0.95)
      RAGradient[RAGradient<quantile(vrbl2, probs=0.0001)]<-quantile(vrbl2, probs=0.0001)
      RAGradient[RAGradient>quantile(vrbl2, probs=0.999)] <- quantile(vrbl2, probs=0.999)
      LWPGradient[LWPGradient<quantile(vrbl3, probs=0.0001)]<-quantile(vrbl3, probs=0.0001)
      LWPGradient[LWPGradient>quantile(vrbl3, probs=0.999)] <- quantile(vrbl3, probs=0.999)
      RDGradient[RDGradient<quantile(vrbl4, probs=0.0001)]<-quantile(vrbl4, probs=0.0001)
      RDGradient[RDGradient>quantile(vrbl4, probs=0.999)] <- quantile(vrbl4, probs=0.999)

       
      
      ringTable=as.data.frame(as.table(CAGradient))
      ringVrbl = "CA"
      names(ringTable)[names(ringTable) == "Freq"] <- "ringVrbl"
      fig1 = ggplot(color="transparent") +
        ggtitle(paste(f, ringVrbl)) +
        theme(text=element_text(size=28), axis.text = element_text(size=14),
              axis.title=element_text(size=12,face="bold"),
              panel.background = element_rect(fill = "transparent",colour = NA),
              plot.background = element_rect(fill = "transparent",colour = NA),
              legend.position = "None") +
        coord_cartesian(xlim=c(0,650), ylim=c(0,350)) +
        scale_x_continuous(name="Stem diameter (mm)", breaks=seq(0,650, by=100), limits=c(0,650, by=100)) +
        scale_y_continuous(name="Height (dm)", breaks=seq(0,350, by=50), limits=c(0,350))  +
        scale_fill_gradientn(name= ringVrbl, colors=c("yellow", "lightgreen", "darkgreen", "purple", "darkblue", "black"), values=cols1, limits=c(min(cols1), max(cols1)), rescaler=function(x,...) x, na.value="white", guide="colorbar") +
        scale_color_gradientn(name= ringVrbl, colors=c("yellow", "lightgreen", "darkgreen", "purple", "darkblue", "black"), values=cols1, limits=c(min(cols1), max(cols1)), rescaler=function(x,...) x, na.value="white", guide="colorbar") +
        geom_raster(data=ringTable, aes(x=as.numeric(d), y=as.numeric(h), fill=ringVrbl))

      
      ringTable=as.data.frame(as.table(RAGradient))
      ringVrbl = "RA"
      names(ringTable)[names(ringTable) == "Freq"] <- "ringVrbl"
      fig2 = ggplot(color="transparent") +
        ggtitle(paste(f, ringVrbl)) +
        theme(text=element_text(size=28), axis.text = element_text(size=14),
              axis.title=element_text(size=12,face="bold"),
              panel.background = element_rect(fill = "transparent",colour = NA),
              plot.background = element_rect(fill = "transparent",colour = NA),
              legend.position = "None") +
        coord_cartesian(xlim=c(0,650), ylim=c(0,350)) +
        scale_x_continuous(name="Stem diameter (mm)", breaks=seq(0,650, by=100), limits=c(0,650, by=100)) +
        scale_y_continuous(name="Height (dm)", breaks=seq(0,350, by=50), limits=c(0,350))  +
        scale_fill_gradientn(name= ringVrbl, colors=c("yellow", "lightgreen", "darkgreen", "purple", "darkblue", "black"), values=cols2, limits=c(min(cols2), max(cols2)), rescaler=function(x,...) x, na.value="white", guide="colorbar") +
        scale_color_gradientn(name= ringVrbl, colors=c("yellow", "lightgreen", "darkgreen", "purple", "darkblue", "black"), values=cols2, limits=c(min(cols2), max(cols2)), rescaler=function(x,...) x, na.value="white", guide="colorbar") +
        geom_raster(data=ringTable, aes(x=as.numeric(d), y=as.numeric(h), fill=ringVrbl))
      
      ringTable=as.data.frame(as.table(LWPGradient))
      ringVrbl = "LWP"
      names(ringTable)[names(ringTable) == "Freq"] <- "ringVrbl"
      fig3 = ggplot(color="transparent") +
        ggtitle(paste(f, ringVrbl)) +
        theme(text=element_text(size=28), axis.text = element_text(size=14),
              axis.title=element_text(size=12,face="bold"),
              panel.background = element_rect(fill = "transparent",colour = NA),
              plot.background = element_rect(fill = "transparent",colour = NA),
              legend.position = "None") +
        coord_cartesian(xlim=c(0,650), ylim=c(0,350)) +
        scale_x_continuous(name="Stem diameter (mm)", breaks=seq(0,650, by=100), limits=c(0,650, by=100)) +
        scale_y_continuous(name="Height (dm)", breaks=seq(0,350, by=50), limits=c(0,350))  +
        scale_fill_gradientn(name= ringVrbl, colors=c("yellow", "lightgreen", "darkgreen", "purple", "darkblue", "black"), values=cols3, limits=c(min(cols3), max(cols3)), rescaler=function(x,...) x, na.value="white", guide="colorbar") +
        scale_color_gradientn(name= ringVrbl, colors=c("yellow", "lightgreen", "darkgreen", "purple", "darkblue", "black"), values=cols3, limits=c(min(cols3), max(cols3)), rescaler=function(x,...) x, na.value="white", guide="colorbar") +
        geom_raster(data=ringTable, aes(x=as.numeric(d), y=as.numeric(h), fill=ringVrbl))
      
      ringTable=as.data.frame(as.table(RDGradient))
      ringVrbl = "RD"
      names(ringTable)[names(ringTable) == "Freq"] <- "ringVrbl"
      fig4 = ggplot(color="transparent") +
        ggtitle(paste(f, ringVrbl)) +
        theme(text=element_text(size=28), axis.text = element_text(size=14),
              axis.title=element_text(size=12,face="bold"),
              panel.background = element_rect(fill = "transparent",colour = NA),
              plot.background = element_rect(fill = "transparent",colour = NA),
              legend.position = "None") +
        coord_cartesian(xlim=c(0,650), ylim=c(0,350)) +
        scale_x_continuous(name="Stem diameter (mm)", breaks=seq(0,650, by=100), limits=c(0,650, by=100)) +
        scale_y_continuous(name="Height (dm)", breaks=seq(0,350, by=50), limits=c(0,350))  +
        scale_fill_gradientn(name= ringVrbl, colors=c("yellow", "lightgreen", "darkgreen", "purple", "darkblue", "black"), values=cols4, limits=c(min(cols4), max(cols4)), rescaler=function(x,...) x, na.value="white", guide="colorbar") +
        scale_color_gradientn(name= ringVrbl, colors=c("yellow", "lightgreen", "darkgreen", "purple", "darkblue", "black"), values=cols4, limits=c(min(cols4), max(cols4)), rescaler=function(x,...) x, na.value="white", guide="colorbar") +
        geom_raster(data=ringTable, aes(x=as.numeric(d), y=as.numeric(h), fill=ringVrbl))
      
      show(fig1)
      show(fig2)
      show(fig3)
      show(fig4)
      
      # ggsave(filename = paste0("ringFigures/",f,"_CA",".png"), plot = fig1, width = 6, height=7.5, device = png())
      # dev.off()
      # ggsave(filename = paste0("ringFigures/",f,"_RA",".png"), plot = fig2, width = 6, height=7.5, device = png())
      # dev.off()
      # ggsave(filename = paste0("ringFigures/",f,"_LWP",".png"), plot = fig3, width = 6, height=7.5, device = png())
      # dev.off()
      # ggsave(filename = paste0("ringFigures/",f,"_RD",".png"), plot = fig4, width = 6, height=7.5, device = png())
      # dev.off()
    }
  }
  
  return(results)
}



