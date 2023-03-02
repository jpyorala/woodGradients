##Author: Jiri Pyörälä (University of Helsinki, Finnish Geospatial Institute), 2022
##jiri.pyorala@helsinki.fi
##Parts of research data and methods used in Pehkonen et al. (Manuscript) "How knot structure and stem geometry reflect ring properties in Norway spruce?"
#https://osf.io/984tk/S

##input: modeling data with ring properties and stem and knottiness features, and ring model
##loop stands
##predict ring properties to stems

woodGradient <- function(DBH, H, stem, model, vrbl, groupName, logt=FALSE, CA_grad=NULL, RA_grad=NULL, LWP_grad=NULL){
  
  library(nlme)
  library(ggplot2)
  
  D20 = DBH/predict(stem, newdata = data.frame(relh=1.3/H))
  D0 = predict(stem, newdata = data.frame(relh=0)) * D20
    
  gradient = matrix(NA, nrow = H*10, ncol = D0*12, dimnames = list(h = seq(1,H*10), d = seq(1,D0*12)))
  
  ###Predict wood properties to d x h matrix:
  for(d in seq(1,D0*12)){
    
    for(h in seq(1, nrow(gradient)-1)){ #loop through heights
      
      h=h/10
      Hr = H * (d/(D0*10)) #tree height at given year, assume fixed H:D ratio, m
      d20 = d/predict(stem, newdata=data.frame(relh=0)) #diam at 20-% height, mm
      #DBHr = predict(stemModel, newdata=data.frame(relh=1.3/Hr)) * d20
      
      Rh = predict(stem, newdata=data.frame(relh=h/H)) * D20 #final stem diameter at height h, mm
      relH = h/Hr
      #relD = r/(Rh/2)
      
      if(d/10<=Rh & d>2 & Rh >= 0){
        
        if(!is.null(CA_grad)){CA=CA_grad[h*10,d]} else{CA=NA}
        if(!is.null(RA_grad)){RA=RA_grad[h*10,d]} else{RA=NA}
        if(!is.null(LWP_grad)){LWP=LWP_grad[h*10,d]} else{LWP=NA}
        
        d_pred = data.frame("Stand"=f, "h"=h, "H"=H, "d"=d, "R"=d/2, "DBH"=DBH, "Hr"=Hr, "relH"=relH, "CA"=CA, "RA"=RA, "LWP"=LWP)
        
        if(logt==TRUE){pred1 = exp(predict(model, newdata = d_pred))}else{pred1 = predict(model, newdata = d_pred)}
        
        h=h*10
        gradient[h,d] <- pred1

      }
    }
  }
  return(gradient)
}



