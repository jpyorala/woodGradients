##Authors: Jiri Pyörälä (University of Helsinki, Finnish Geospatial Institute) and Mika Pehkonen (University of Helsinki, Finnish Geospatial Institute), 2022
##Methods used in Pehkonen et al. (Manuscript) "How knot structure and stem geometry reflect ring properties in Norway spruce?"
library(nlme)
library(rsq)
source("ringGradients.R")

model_data = read.delim("model_data.txt", header=T, sep="\t")
stem_data = read.delim("stem_data.txt", header=TRUE) # read X-ray data

###LIST OF EXPLANATORY VARIABLES
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



###RESPONSE RING VARIABLES
resp= c("CA", "RA", "LWP","RD")

results = data.frame("Resp"=character(), "Expl"=character(), 
                     "R2_fixed"=double(), "R2_random"=double(), 
                     "int"=double(), "p_int"=double(), 
                     "coef"=double(), "p_coef"=double(),
                     "RandIntSD"=double(), "RandResSD"=double(),
                     "RMSE"=double(), "RMSE%"=double())

###TEST STAND-LEVEL SIMPLE MIXED MODELS
for(y in resp){
  
  for(x in expl){
    
    if(y!=x){
    
      MM = lme(formula(paste(y,x,sep="~")), random= ~1|Stand, data=model_data, 
               control=lmeControl(opt="optim", maxIter=100),na.action=na.omit)
      s=summary(MM)
      a=anova(MM)
      r2=rsq.lmm(MM)
      obs=model_data[c(!is.na(model_data[,y])&!is.na(model_data[,x])),y] #observed values
      fit=fitted(MM) #fitted values of model (takaisin muunnos -> huomioi variance)
      res=resid(MM) #residuals
      
      # 
      results=rbind(results,
             data.frame("Resp"=y, "Expl"=x,
             "R2_fixed"=r2$fixed, "R2_random"=r2$random,
             "int"=s$coefficients$fixed[1], "p_int"=a$p[1],
             "coef"=s$coefficients$fixed[2], "p_coef"=a$p[2],
             "RandIntSD"=sd(s$coefficients$random$Stand), "ResSD"=s$sigma,
             "RMSE"=sqrt((sum((fit-obs)^2))/length(obs)),
             "RMSE%"=(sqrt((sum((fit-obs)^2))/length(obs)))/mean(obs)*100))
      }
    }
  }


#write.table(results,"simpleMixedModel_results.txt",  sep=" ", quote = F, append = F, row.names = F)

###FINAL STAND-LEVEL MULTIPLE MIXED MODELS

##CAMBIAL AGE##
CA_MM1 = lme(log(CA) ~ 0 + relH + log(R), 
             random= ~relH + log(R)|Stand, data=model_data)
summary(CA_MM1)
anova(CA_MM1)
rsq.lmm(CA_MM1)
plot(model_data$CA, exp(fitted(CA_MM1)), col=model_data$Stand)

CA_MM2 = lme(log(CA) ~ log(Vr_Hr_h) + dCrKI + KImax, 
             random= ~ log(Vr_Hr_h) + dCrKI + KImax |Stand, data=model_data, 
             control=lmeControl(maxIter=1000, msMaxIter=1000,opt="optim", tolerance=0.001),
             na.action=na.omit)
summary(CA_MM2)
anova(CA_MM2)
rsq.lmm(CA_MM2)
plot(model_data$CA, exp(fitted(CA_MM2)), col=model_data$Stand)


##RING AREA##
RA_MM1 = lme(log(RA) ~ relH + poly(log(R),2)*log(CA),
             random= ~relH + poly(log(R),2)*log(CA)|Stand, data=model_data, 
             control=lmeControl(maxIter=1000, msMaxIter=1000,opt="optim", tolerance=0.001))
summary(RA_MM1)
anova(RA_MM1)
rsq.lmm(RA_MM1)
plot(model_data$RA, exp(fitted(RA_MM1)), col=model_data$Stand)
abline(0,1)

plot(model_data$R, exp(fitted(RA_MM1)), col=model_data$Stand)

RA_MM2 = lme(log(RA) ~log(Vr_Hr_h) + reldCrKSI + KImax, 
             random= ~log(Vr_Hr_h) + reldCrKSI + KImax|Stand, data=model_data,
             control=lmeControl(maxIter=1000, msMaxIter=1000,opt="optim", tolerance=0.001))
summary(RA_MM2)
anova(RA_MM2)
rsq.lmm(RA_MM2)
plot(model_data$RA, exp(fitted(RA_MM2)), col=model_data$Stand)
abline(0,1)


##LATEWOOD PERCENTAGE##
LWP_MM1 = lme(log(LWP) ~ relH + log(RA) + log(R)*log(CA), random= ~relH + log(RA) + log(R)*log(CA)|Stand, data=model_data,
              control=lmeControl(maxIter=1000, msMaxIter=1000,  tolerance=0.001, opt="optim"))
summary(LWP_MM1)
anova(LWP_MM1)
rsq.lmm(LWP_MM1)
plot(model_data$LWP, exp(fitted(LWP_MM1)), col=model_data$Stand)
abline(0,1)

plot(model_data$relH, exp(fitted(LWP_MM1)), col=model_data$Stand)


LWP_MM2 = lme(log(LWP) ~ log(Vr)+poly(dCrKI,2),
              random= ~log(Vr)+poly(dCrKI,2)|Stand, data=model_data,
              control=lmeControl(maxIter=1000, msMaxIter=1000,  tolerance=0.001, opt="optim"))
summary(LWP_MM2)
anova(LWP_MM2)
rsq.lmm(LWP_MM2)
plot(model_data$LWP, exp(fitted(LWP_MM2)), col=model_data$Stand)
abline(0,1)


##RING DENSITY##
RD_MM1 = lme(log(RD) ~ relH + poly(log(R),2) + log(CA) + log(RA) + log(LWP), 
             random= ~relH + poly(log(R),2) + log(CA) + log(RA) + log(LWP)|Stand, data=model_data,
             control=lmeControl(maxIter=1000, msMaxIter=1000,  tolerance=0.001, opt="optim"))
summary(RD_MM1)
anova(RD_MM1)
rsq.lmm(RD_MM1)
plot(model_data$RD, exp(fitted(RD_MM1)), col=model_data$Stand)
abline(0,1)

plot(model_data$relH, exp(fitted(RD_MM1)), col=model_data$Stand)


RD_MM2 = lme(log(RD) ~ log(Vr) + dCrKI + relCrlKI, 
             random= ~ log(Vr) + dCrKI + relCrlKI|Stand, data=model_data,
             control=lmeControl(maxIter=1000, msMaxIter=1000,  tolerance=0.001, opt="optim"),
             na.action=na.omit)
summary(RD_MM2)
anova(RD_MM2)
rsq.lmm(RD_MM2)
plot(model_data$RD, exp(fitted(RD_MM2)), col=model_data$Stand)
abline(0,1)


##VISUALIZE THE REFERENCE RING PROPERTY MODELS TO STEMS PER STAND
ringGradients(stem_data, model_data, group1="Stand")
