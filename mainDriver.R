##Authors: Jiri Pyörälä (University of Helsinki, Finnish Geospatial Institute) and Mika Pehkonen (University of Helsinki, Finnish Geospatial Institute), 2022
##jiri.pyorala@helsinki.fi
##Research data and methods used in Pehkonen et al. (Manuscript) "How knot structure and stem geometry reflect ring properties in Norway spruce?"
#https://osf.io/984tk/S
library(rsq)
library(nlme)
library(ggplot2)
source("helperFunctions.R")
source("woodGradient.R")
source("gradientFigure.R")
source("mixedModels.R")

stemdata = read.delim("stem_data.txt", header=TRUE) # read X-ray data
ringdata = rings2logs(ringFile="ring_data.txt", stemFile="stem_data.txt") # match measured rings from ring-sample trees to corresponding x-ray data

stemdata = stemdata[!stemdata$TreeID %in% ringdata$Tree,]  #remove ring-sample trees from the modeling data

##CREATE OUTPUT VECTORS FOR MODEL RESULTS AND PREDICTIONS (GRADIENTS) FOR DESIRED VARIABLES (3) AND STANDS (14)
knotModel_results = data.frame("Stand"=character(), "vrbl" = character(), "R2"=double(), "Sigma"=double(), "DF"=double())
knot_gradients = list(vector("list", 14),
                     vector("list", 14),
                     vector("list", 14))

#INITIATE KNOT MODEL STRUCTURES
knotModel1 = function(H, DBH, h, d) (h/H) + log(d/DBH)
knotModel2 = function(H, DBH, h, d) (h/H) + log(d/DBH)
knotModel3 = function(H, DBH, h, d) (H*DBH) + (h/H) + log(d/DBH) 

#DEFINE COLOR FOR VISUALIZATION
KI = na.exclude(stemdata[,"KnotIndex"])
KI = KI[KI != Inf]
KIcols = seq(0, 0.125, length=48) #%/100

KSI = na.exclude(stemdata[,"KnotSizeIndex"])
KSI = KSI[KSI != Inf]
KSIcols = seq(0, 2.5, length=48) #cm3

WTW = na.exclude(stemdata[,"KnotClusterDistanceAverage"])
WTW = WTW[WTW != Inf]
WTWcols = seq(0, 600, length=48) #mm


#LOOP THROUGH STANDS AND TARGET VARIABLES: CREATE DxH GRADIENTS FOR EACH STAND AND VARIABLE 
for(f in levels(as.factor(stemdata$Stand))){
  n=0
   for(vrbl in c("KnotIndex","KnotSizeIndex","KnotClusterDistanceAverage")){
     n=n+1
     if(n==1){cols=KIcols}else if(n==2){cols=KSIcols}else{cols=WTWcols}
    #SUBSET DATA FOR THE SPECIFIC VARIABLE AND STAND, AND ORGANIZE INTO STEMS TO ENABLE FUNCTIONIZING BY H x D
    modeldata1 = organizeLogs(df = stemdata[stemdata$Stand==f,], vrbl = vrbl)

    #REPARAMETERIZE THE STEM AND KNOT MODELS
    stemModel=lm(reld ~ poly(relh, 4), data=modeldata1)
    model=as.formula(paste(vrbl,
                           paste(as.character(body(paste0("knotModel",n))), collapse=("+")),
                           sep="~"))
    names(modeldata1)[names(modeldata1)=="vrbl"] <- vrbl

    lmModel = lm(model, data=modeldata1)
    sm=summary(lmModel)
    knotModel_results = rbind(knotModel_results,
                              data.frame("Stand"=f, "vrbl" = vrbl, "R2"=sm$r.squared, "Sigma"=sm$sigma, "DF"=sm$df[2]))

    ##CALCULATE THE GRADIENT TO A MATRIX
    gradient <- woodGradient(DBH=max(modeldata1$DBH, na.rm=T), H=max(modeldata1$H, na.rm=T),
                                             stem=stemModel, model=lmModel,
                                             vrbl=vrbl, groupName=f)
    #SAVE THE GRADIENT
    knot_gradients[[n]][[f]] <- gradient
    
    #LITTLE CLEAN-UP AND VISUALIZATION
    gradient[gradient>max(cols)] <- max(cols)
    gradientFigure(gradient, f, 
                   vrbl=vrbl, 
                   cols = cols, 
                   view=F, 
                   save=F, folder="Figures/")
  
    }
  }


#LIST OF KNOT FEATURES
expl = c("Year", #Year             
         "CA", #cambial age of the ring, a        
         "RA", #area of the ring, mm^2         
         "LWP", #latewood percentage of the ring, %/100        
         "RD", #mean density of the ring, kg/m3/100        
         "R", #stem radius at the ring, or "distance-from-the-pith", mm      
         "TreeAge", #total tree age at given year, a        
         "h", #absolute sample height, m          
         "DBHr", #DBH at given year, cm         
         "Hr", #tree height at given year, m        
         "H",              
         "DBH",              
         "D",              
         "relD", #stem diameter (R*2) relative to DBH at given year (DBHr)    
         "relH", #sample height relative to tree height at given year (Hr)    
         "Vr", #volume of the stem at given year, m3      
         "Vr_Hr_h", #volume of the stem above the location at given year, m3   
         "SL0", #length of last shoot at given year, cm      
         "SL1", #length of previous years shoot, cm        
         "KIr", #knot index just above the location        
         "KSIr",#, #knot size index just above the location       
         "KImax", #maximum knot index (annual, above h)        
         "KImean", #mean knot index (annual, above h)        
         "KSImax", #maximum knot size index (annual, above h)       
         "KSImean", #mean knot size index (annual, above h)       
         "hCrKI", #height of the maximum knot index, m       
         "hCrKSI", #height of the maximum knot size index, m      
         "relhCrKI", #height of the maximum knot index, m       
         "relhCrKSI", #height of the maximum knot size index, m      
         "CrlKI", #distance from the maximum knot index height to tree top, m   
         "CrlKSI", #distance from the maximum knot size index height to tree top, m  
         "relCrlKI", #distance from the maximum knot index height to tree top, m   
         "relCrlKSI", #distance from the maximum knot size index height to tree top, m  
         "dCrKI", #distance from crown base (negative=below crown, positive=within crown), m     
         "dCrKSI", #distance from crown base (negative=below crown, positive=within crown), m     
         "reldCrKI", #distance from crown base (negative=below crown, positive=within crown), m     
         "reldCrKSI") #distance from crown base (negative=below crown, positive=within crown), m     


ringdata = ringdata[!is.na(ringdata$buttH),]

###CALCULATE KNOT FEATURES (AS ABOVE) TO RINGS
modeldata2 = knots2Rings(ringdata, group1="Stand", stemdata, knot_gradients)
#write.table(modeldata2, "modeldata2.txt",  sep=" ", quote = F, append = F, row.names = F)
#modeldata2 = read.delim("modeldata2.txt", header=T, sep=" ")

###FILTER DATA
modeldata2_filt = subset(modeldata2, c(modeldata2$relH<=1&modeldata2$Hr>0))
modeldata2_filt = modeldata2_filt[complete.cases(modeldata2_filt),]
modeldata2_filt$RD=modeldata2_filt$RD*1000
modeldata2_filt$LWP=modeldata2_filt$LWP*100
#write.table(modeldata2_filt, "modeldata2_filtered.txt",  sep=" ", quote = F, append = F, row.names = F)
#modeldata2_filt = read.delim("modeldata2_filtered.txt", header=T, sep=" ")

######TEST STAND-LEVEL SIMPLE MIXED MODELS OF THE RESPONSE RING VARIABLES
resp= c("CA", "RA", "LWP","RD")
simpleMixedModel_results <- simpleMixedModels(resp, expl, model_data=modeldata2_filt)
#write.table(simpleMixedModel_results,"simpleMixedModel_results.txt",  sep=" ", quote = F, append = F, row.names = F)
#simpleMixedModel_results = read.delim("simpleMixedModel_results.txt", header=T, sep=" ")

multipleMixedModel_results <- multipleMixedModels(model_data=modeldata2_filt)
#write.table(multipleMixedModel_results,"multipleMixedModel_results.txt",  sep=" ", quote = F, append = F, row.names = F)
#multipleMixedModel_results = read.delim("multipleMixedModel_results.txt", header=T, sep=" ")

save.image("/projappl/project_2003078/jpyorala/Scripts/ResearchData_Pehkonenetal23/.RData")

#INITIATE RING MODEL STRUCTURES: SEE mixedModels.R: multipleMixedModels()
ringModel1 = function(H, DBH, h, d) 0 + relH + log(R)
ringModel2 = function(H, DBH, h, d) relH + poly(log(R),2)*log(CA)
ringModel3 = function(H, DBH, h, d) relH + log(RA) + log(R)*log(CA)
ringModel4 = function(H, DBH, h, d) relH + poly(log(R),2) + log(CA) + log(RA) + log(LWP)

CA_model=as.formula(paste(paste0("log(CA)"),
                          paste(as.character(body(paste0("ringModel",1))), collapse=("+")), 
                          sep="~"))
RA_model=as.formula(paste(paste0("log(RA)"),
                          paste(as.character(body(paste0("ringModel",2))), collapse=("+")), 
                          sep="~"))
LWP_model=as.formula(paste(paste0("log(LWP)"),
                          paste(as.character(body(paste0("ringModel",3))), collapse=("+")), 
                          sep="~"))
RD_model=as.formula(paste(paste0("log(RD)"),
                          paste(as.character(body(paste0("ringModel",4))), collapse=("+")), 
                          sep="~"))


#names(ringdata)[names(ringdata)=="Tree"] <- "TreeID"
ring_gradients = list(vector("list", 14),
                 vector("list", 14),
                 vector("list", 14),
                 vector("list", 14))

ringModel_results = data.frame("Stand"=character(), "vrbl" = character(), "R2"=double(), "Sigma"=double(), "DF"=double())

#DEFINE COLORS FOR THE VISUALIZATION
vrbl1 = na.exclude(modeldata2_filt[,"CA"])
vrbl1 = vrbl1[vrbl1 != Inf]
cols1 = seq(0, 120, length=48) #years

vrbl2 = na.exclude(modeldata2_filt[,"RA"])
vrbl2 = vrbl2[vrbl2 != Inf]
cols2 = seq(0, 4200, length=48) #mm2

vrbl3 = na.exclude(modeldata2_filt[,"LWP"])
vrbl3 = vrbl3[vrbl3 != Inf]
cols3 = seq(5, 60, length=48) #%

vrbl4 = na.exclude(modeldata2_filt[,"RD"])
vrbl4 = vrbl4[vrbl4 != Inf]
cols4 = seq(250, 650, length=48) #kg/m3

#LOOP THROUGH STANDS: VISUALIZE DxH GRADIENTS FOR EACH STAND AND VARIABLE 
for(f in levels(as.factor(modeldata2_filt$Stand))){

  #SUBSET DATA FOR THE SPECIFIC STAND, AND REPARAMETERIZE STEM MODEL
  modeldata3 = organizeLogs(df = stemdata[stemdata$Stand==f,], vrbl = "Stand")
  stemModel=lm(reld ~ poly(relh, 4), data=modeldata3)

  CA_lmModel = lm(CA_model, data=modeldata2_filt[modeldata2_filt$Stand==f,])
  sm=summary(CA_lmModel)
  ringModel_results = rbind(ringModel_results,
                            data.frame("Stand"=f, "vrbl" = "CA", "R2"=sm$r.squared, "Sigma"=sm$sigma, "DF"=sm$df[2]))

  ##CALCULATE THE GRADIENT TO A MATRIX
  gradient <- woodGradient(DBH=max(modeldata3$DBH, na.rm=T),
                                           H=max(modeldata3$H, na.rm=T),
                                           stem=stemModel,
                                           model=CA_lmModel,
                                           vrbl="CA",
                                           groupName=f,
                                           logt=TRUE)
  #SAVE THE GRADIENT
  ring_gradients[[1]][[f]] <- gradient
  
  #LITTLE CLEAN-UP BEFORE THE VISUALIZATION
  gradient[gradient<0]<-0
  gradient[gradient>120]<-120
  
  gradientFigure(gradient, f, 
                 vrbl="CA", 
                 cols = cols1, 
                 view=T, 
                 save=T, 
                 folder="Figures/")
  
  
  RA_lmModel = lm(RA_model, data=modeldata2_filt[modeldata2_filt$Stand==f,])
  sm=summary(RA_lmModel)
  ringModel_results = rbind(ringModel_results,
                            data.frame("Stand"=f, "vrbl" = "RA", "R2"=sm$r.squared, "Sigma"=sm$sigma, "DF"=sm$df[2]))

  gradient <- woodGradient(DBH=max(modeldata3$DBH, na.rm=T),
                                           H=max(modeldata3$H, na.rm=T),
                                           stem=stemModel,
                                           model=RA_lmModel,
                                           vrbl="RA",
                                           groupName=f,
                                           logt=TRUE,
                                           CA_grad = ring_gradients[[1]][[f]])
  
  ring_gradients[[2]][[f]]  <- gradient
  
  gradient[gradient<0]<-0
  gradient[gradient>4200]<-4200
  
  gradientFigure(gradient, f, 
                 vrbl="RA", 
                 cols = cols2, 
                 view=T, 
                 save=T, 
                 folder="Figures/")
  
  
  LWP_lmModel = lm(LWP_model, data=modeldata2_filt[modeldata2_filt$Stand==f,])
  sm<-summary(LWP_lmModel)
  ringModel_results = rbind(ringModel_results,
                            data.frame("Stand"=f, "vrbl" = "LWP", "R2"=sm$r.squared, "Sigma"=sm$sigma, "DF"=sm$df[2]))

  gradient <- woodGradient(DBH=max(modeldata3$DBH, na.rm=T),
                                            H=max(modeldata3$H, na.rm=T),
                                            stem=stemModel,
                                            model=LWP_lmModel,
                                            vrbl="LWP",
                                            groupName=f,
                                            logt=TRUE,
                                            CA_grad = ring_gradients[[1]][[f]],
                                            RA_grad = ring_gradients[[2]][[f]])
  
  ring_gradients[[3]][[f]] <- gradient
  
  gradient[gradient<5]<-5
  gradient[gradient>60]<-60

  gradientFigure(gradient, f, 
                 vrbl="LWP", 
                 cols = cols3, 
                 view=T, 
                 save=T, 
                 folder="Figures/")
  
  
  RD_lmModel = lm(RD_model, data=modeldata2_filt[modeldata2_filt$Stand==f,])
  sm<-summary(RD_lmModel)
  ringModel_results = rbind(ringModel_results,
                            data.frame("Stand"=f, "vrbl" = "RD", "R2"=sm$r.squared, "Sigma"=sm$sigma, "DF"=sm$df[2]))

  gradient <- woodGradient(DBH=max(modeldata3$DBH, na.rm=T),
                                           H=max(modeldata3$H, na.rm=T),
                                           stem=stemModel,
                                           model=RD_lmModel,
                                           vrbl="RD",
                                           groupName=f,
                                           logt=TRUE,
                                           CA_grad = ring_gradients[[1]][[f]],
                                           RA_grad = ring_gradients[[2]][[f]],
                                           LWP_grad = ring_gradients[[3]][[f]])

  ring_gradients[[4]][[f]] <- gradient
  
  gradient[gradient<250]<-250
  gradient[gradient>650]<-650
  
  gradientFigure(gradient, f, 
                 vrbl="RD", 
                 cols = cols4, 
                 view=T, 
                 save=T, 
                 folder="Figures/")
  
}

save.image(".RData")

