##Author: Dr. Jiri Pyörälä (University of Helsinki, Finnish Geospatial Institute), 2023
##jiri.pyorala@helsinki.fi
##Research data and methods used in Pehkonen et al. (Manuscript) "How knot structure and stem geometry reflect ring properties in Norway spruce?"
#https://osf.io/984tk/S

simpleMixedModels <- function(resp, expl, model_data){
  
  simpleMixedModel_results = data.frame("Resp"=character(), "Expl"=character(), 
                       "R2_fixed"=double(), "R2_random"=double(), 
                       "int"=double(), "p_int"=double(), 
                       "coef"=double(), "p_coef"=double(),
                       "RandIntSD"=double(), "RandResSD"=double(),
                       "RMSE"=double(), "RMSE%"=double())
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
        simpleMixedModel_results=rbind(simpleMixedModel_results,
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
  return(simpleMixedModel_results)
}



###FINAL STAND-LEVEL MULTIPLE MIXED MODELS

multipleMixedModels <- function(model_data){
  multipleMixedModel_results=data.frame("Resp"=character(), "Expl"=character(), 
                     "R2_fixed"=double(), "R2_random"=double(), 
                     "int"=double(), "int_SE"=double(), 
                     "coef1"=double(), "coef2"=double(), 
                     "coef3"=double(), "coef4"=double(),
                     "coef5"=double(), "coef6"=double(),
                     "RandIntSD"=double(), "RandResSD"=double(),
                     "RMSE"=double(), "RMSE%"=double())
  ##CAMBIAL AGE##
  CA_MM1 = lme(log(CA) ~ 0 + relH + log(R), 
               random= ~relH + log(R)|Stand, data=model_data)
  s=summary(CA_MM1)
  a=anova(CA_MM1)
  r2=rsq.lmm(CA_MM1)
  obs=model_data[!is.na(model_data[,"CA"]),"CA"] #observed values
  fit=exp(fitted(CA_MM1))
  res=exp(resid(CA_MM1))
  
  plot(obs, fit, col=model_data$Stand, main="RD pred2")
  abline(0,1)
  plot(fit, res, col=model_data$Stand, main="RD res2")
  abline(h=0)
  
  print("CA model1")
  print(s)
  print(r2)
  
  multipleMixedModel_results = rbind(multipleMixedModel_results,
                                     data.frame("Resp"="CA", "Expl"="model1",
                                                "R2_fixed"=r2$fixed, "R2_random"=r2$random,
                                                "int"=s$coefficients$fixed[1],
                                                "coef1"=s$coefficients$fixed[2],"coef2"=s$coefficients$fixed[3],
                                                "coef3"=s$coefficients$fixed[4],"coef4"=s$coefficients$fixed[5],
                                                "coef5"=s$coefficients$fixed[6],"coef6"=s$coefficients$fixed[7],
                                                
                                                "RandIntSD"=sd(s$coefficients$random$Stand), "ResSD"=s$sigma,
                                                "RMSE"=sqrt((sum((fit-obs)^2))/length(obs)),
                                                "RMSE%"=(sqrt((sum((fit-obs)^2))/length(obs)))/mean(obs)*100))
                                     
  
  CA_MM2 = lme(log(CA) ~ log(Vr_Hr_h) + dCrKI + KImax, 
               random= ~ log(Vr_Hr_h) + dCrKI + KImax |Stand, data=model_data, 
               control=lmeControl(maxIter=1000, msMaxIter=1000,opt="optim", tolerance=0.001),
               na.action=na.omit)
  
  s=summary(CA_MM2)
  a=anova(CA_MM2)
  r2=rsq.lmm(CA_MM2)
  obs=model_data[!is.na(model_data[,"CA"]),"CA"] #observed values
  fit=exp(fitted(CA_MM2))
  res=exp(resid(CA_MM2))
  
  plot(obs, fit, col=model_data$Stand, main="RD pred2")
  abline(0,1)
  plot(fit, res, col=model_data$Stand, main="RD res2")
  abline(h=0)
  
  print("CA model2")
  print(s)
  print(r2)
  
  multipleMixedModel_results = rbind(multipleMixedModel_results,
                                     data.frame("Resp"="CA", "Expl"="model2",
                                                "R2_fixed"=r2$fixed, "R2_random"=r2$random,
                                                "int"=s$coefficients$fixed[1],
                                                "coef1"=s$coefficients$fixed[2],"coef2"=s$coefficients$fixed[3],
                                                "coef3"=s$coefficients$fixed[4],"coef4"=s$coefficients$fixed[5],
                                                "coef5"=s$coefficients$fixed[6],"coef6"=s$coefficients$fixed[7],
                                                
                                                "RandIntSD"=sd(s$coefficients$random$Stand), "ResSD"=s$sigma,
                                                "RMSE"=sqrt((sum((fit-obs)^2))/length(obs)),
                                                "RMSE%"=(sqrt((sum((fit-obs)^2))/length(obs)))/mean(obs)*100))
  
  
  ##RING AREA##
  RA_MM1 = lme(log(RA) ~ relH + poly(log(R),2)*log(CA),
               random= ~relH + poly(log(R),2)*log(CA)|Stand, data=model_data, 
               control=lmeControl(maxIter=1000, msMaxIter=1000,opt="optim", tolerance=0.001))
  s=summary(RA_MM1)
  a=anova(RA_MM1)
  r2=rsq.lmm(RA_MM1)
  obs=model_data[!is.na(model_data[,"RA"]),"RA"] #observed values
  fit=exp(fitted(RA_MM1))
  res=exp(resid(RA_MM1))
  
  plot(obs, fit, col=model_data$Stand, main="RD pred2")
  abline(0,1)
  plot(fit, res, col=model_data$Stand, main="RD res2")
  abline(h=0)
  
  print("RA model1")
  print(s)
  print(r2)
  
  multipleMixedModel_results = rbind(multipleMixedModel_results,
                                     data.frame("Resp"="RA", "Expl"="model1",
                                                "R2_fixed"=r2$fixed, "R2_random"=r2$random,
                                                "int"=s$coefficients$fixed[1],
                                                "coef1"=s$coefficients$fixed[2],"coef2"=s$coefficients$fixed[3],
                                                "coef3"=s$coefficients$fixed[4],"coef4"=s$coefficients$fixed[5],
                                                "coef5"=s$coefficients$fixed[6],"coef6"=s$coefficients$fixed[7],
                                                
                                                "RandIntSD"=sd(s$coefficients$random$Stand), "ResSD"=s$sigma,
                                                "RMSE"=sqrt((sum((fit-obs)^2))/length(obs)),
                                                "RMSE%"=(sqrt((sum((fit-obs)^2))/length(obs)))/mean(obs)*100))
  
  RA_MM2 = lme(log(RA) ~log(Vr_Hr_h) + reldCrKSI + KImax, 
               random= ~log(Vr_Hr_h) + reldCrKSI + KImax|Stand, data=model_data,
               control=lmeControl(maxIter=1000, msMaxIter=1000,opt="optim", tolerance=0.001))
  s=summary(RA_MM2)
  a=anova(RA_MM2)
  r2=rsq.lmm(RA_MM2)
  obs=model_data[!is.na(model_data[,"RA"]),"RA"] #observed values
  fit=exp(fitted(RA_MM2))
  res=exp(resid(RA_MM2))
  
  plot(obs, fit, col=model_data$Stand, main="RD pred2")
  abline(0,1)
  plot(fit, res, col=model_data$Stand, main="RD res2")
  abline(h=0)
  
  print("RA model2")
  print(s)
  print(r2)
  
  multipleMixedModel_results = rbind(multipleMixedModel_results,
                                     data.frame("Resp"="RA", "Expl"="model2",
                                                "R2_fixed"=r2$fixed, "R2_random"=r2$random,
                                                "int"=s$coefficients$fixed[1],
                                                "coef1"=s$coefficients$fixed[2],"coef2"=s$coefficients$fixed[3],
                                                "coef3"=s$coefficients$fixed[4],"coef4"=s$coefficients$fixed[5],
                                                "coef5"=s$coefficients$fixed[6],"coef6"=s$coefficients$fixed[7],
                                                
                                                "RandIntSD"=sd(s$coefficients$random$Stand), "ResSD"=s$sigma,
                                                "RMSE"=sqrt((sum((fit-obs)^2))/length(obs)),
                                                "RMSE%"=(sqrt((sum((fit-obs)^2))/length(obs)))/mean(obs)*100))
  
  
  ##LATEWOOD PERCENTAGE##
  LWP_MM1 = lme(log(LWP) ~ relH + log(RA) + log(R)*log(CA), random= ~relH + log(RA) + log(R)*log(CA)|Stand, data=model_data,
                control=lmeControl(maxIter=1000, msMaxIter=1000,  tolerance=0.001, opt="optim"))
  s=summary(LWP_MM1)
  a=anova(LWP_MM1)
  r2=rsq.lmm(LWP_MM1)
  obs=model_data[!is.na(model_data[,"LWP"]),"LWP"] #observed values
  fit=exp(fitted(LWP_MM1))
  res=exp(resid(LWP_MM1))
  
  plot(obs, fit, col=model_data$Stand, main="RD pred2")
  abline(0,1)
  plot(fit, res, col=model_data$Stand, main="RD res2")
  abline(h=0)
  
  print("LWP model1")
  print(s)
  print(r2)
  
  multipleMixedModel_results = rbind(multipleMixedModel_results,
                                     data.frame("Resp"="LWP", "Expl"="model1",
                                                "R2_fixed"=r2$fixed, "R2_random"=r2$random,
                                                "int"=s$coefficients$fixed[1],
                                                "coef1"=s$coefficients$fixed[2],"coef2"=s$coefficients$fixed[3],
                                                "coef3"=s$coefficients$fixed[4],"coef4"=s$coefficients$fixed[5],
                                                "coef5"=s$coefficients$fixed[6],"coef6"=s$coefficients$fixed[7],
                                                
                                                "RandIntSD"=sd(s$coefficients$random$Stand), "ResSD"=s$sigma,
                                                "RMSE"=sqrt((sum((fit-obs)^2))/length(obs)),
                                                "RMSE%"=(sqrt((sum((fit-obs)^2))/length(obs)))/mean(obs)*100))
  
  
  LWP_MM2 = lme(log(LWP) ~ log(Vr)+poly(dCrKI,2),
                random= ~log(Vr)+poly(dCrKI,2)|Stand, data=model_data,
                control=lmeControl(maxIter=1000, msMaxIter=1000,  tolerance=0.001, opt="optim"))
  s=summary(LWP_MM2)
  a=anova(LWP_MM2)
  r2=rsq.lmm(LWP_MM2)
  obs=model_data[!is.na(model_data[,"LWP"]),"LWP"] #observed values
  fit=exp(fitted(LWP_MM2))
  res=exp(resid(LWP_MM2))
  
  plot(obs, fit, col=model_data$Stand, main="RD pred2")
  abline(0,1)
  plot(fit, res, col=model_data$Stand, main="RD res2")
  abline(h=0)
  
  print("LWP model2")
  print(s)
  print(r2)
  
  multipleMixedModel_results = rbind(multipleMixedModel_results,
                                     data.frame("Resp"="LWP", "Expl"="model2",
                                                "R2_fixed"=r2$fixed, "R2_random"=r2$random,
                                                "int"=s$coefficients$fixed[1],
                                                "coef1"=s$coefficients$fixed[2],"coef2"=s$coefficients$fixed[3],
                                                "coef3"=s$coefficients$fixed[4],"coef4"=s$coefficients$fixed[5],
                                                "coef5"=s$coefficients$fixed[6],"coef6"=s$coefficients$fixed[7],
                                                
                                                "RandIntSD"=sd(s$coefficients$random$Stand), "ResSD"=s$sigma,
                                                "RMSE"=sqrt((sum((fit-obs)^2))/length(obs)),
                                                "RMSE%"=(sqrt((sum((fit-obs)^2))/length(obs)))/mean(obs)*100))
  
  
  ##RING DENSITY##
  RD_MM1 = lme(log(RD) ~ relH + poly(log(R),2) + log(CA) + log(RA) + log(LWP), 
               random= ~relH + poly(log(R),2) + log(CA) + log(RA) + log(LWP)|Stand, data=model_data,
               control=lmeControl(maxIter=1000, msMaxIter=1000,  tolerance=0.001, opt="optim"))
  s=summary(RD_MM1)
  a=anova(RD_MM1)
  r2=rsq.lmm(RD_MM1)
  obs=model_data[!is.na(model_data[,"RD"]),"RD"] #observed values
  fit=exp(fitted(RD_MM1))
  res=resid(RD_MM1)
  
  print("RD model1")
  print(s)
  print(r2)
  
  plot(obs, fit, col=model_data$Stand, main="RD pred1")
  abline(0,1)
  plot(fit, res, col=model_data$Stand, main="RD res1")
  abline(h=0)
  multipleMixedModel_results = rbind(multipleMixedModel_results,
                                     data.frame("Resp"="RD", "Expl"="model1",
                                                "R2_fixed"=r2$fixed, "R2_random"=r2$random,
                                                "int"=s$coefficients$fixed[1],
                                                "coef1"=s$coefficients$fixed[2],"coef2"=s$coefficients$fixed[3],
                                                "coef3"=s$coefficients$fixed[4],"coef4"=s$coefficients$fixed[5],
                                                "coef5"=s$coefficients$fixed[6],"coef6"=s$coefficients$fixed[7],
                                                
                                                "RandIntSD"=sd(s$coefficients$random$Stand), "ResSD"=s$sigma,
                                                "RMSE"=sqrt((sum((fit-obs)^2))/length(obs)),
                                                "RMSE%"=(sqrt((sum((fit-obs)^2))/length(obs)))/mean(obs)*100))
  
  
  RD_MM2 = lme(log(RD) ~ log(Vr) + dCrKI + relCrlKI, 
               random= ~ log(Vr) + dCrKI + relCrlKI|Stand, data=model_data,
               control=lmeControl(maxIter=1000, msMaxIter=1000,  tolerance=0.001, opt="optim"),
               na.action=na.omit)
  s=summary(RD_MM2)
  a=anova(RD_MM2)
  r2=rsq.lmm(RD_MM2)
  obs=model_data[!is.na(model_data[,"RD"]),"RD"] #observed values
  fit=exp(fitted(RD_MM2))
  res=resid(RD_MM2)
  
  plot(obs, fit, col=model_data$Stand, main="RD pred2")
  abline(0,1)
  plot(fit, res, col=model_data$Stand, main="RD res2")
  abline(h=0)
  
  print("RD model2")
  print(s)
  print(r2)
  
  multipleMixedModel_results = rbind(multipleMixedModel_results,
                                     data.frame("Resp"="RD", "Expl"="model2",
                                                "R2_fixed"=r2$fixed, "R2_random"=r2$random,
                                                "int"=s$coefficients$fixed[1],
                                                "coef1"=s$coefficients$fixed[2],"coef2"=s$coefficients$fixed[3],
                                                "coef3"=s$coefficients$fixed[4],"coef4"=s$coefficients$fixed[5],
                                                "coef5"=s$coefficients$fixed[6],"coef6"=s$coefficients$fixed[7],
                                                
                                                "RandIntSD"=sd(s$coefficients$random$Stand), "ResSD"=s$sigma,
                                                "RMSE"=sqrt((sum((fit-obs)^2))/length(obs)),
                                                "RMSE%"=(sqrt((sum((fit-obs)^2))/length(obs)))/mean(obs)*100))
  
  return(multipleMixedModel_results)

}
