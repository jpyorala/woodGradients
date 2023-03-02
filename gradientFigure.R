
gradientFigure <- function(gradient, groupName, vrbl, cols, view=TRUE, save=FALSE, folder=getwd()){

  table1=as.data.frame(as.table(gradient))
  names(table1)[names(table1) == "Freq"] <- "vrbl"
  
  fig1 = ggplot(color="transparent") +
    ggtitle(paste(groupName, vrbl)) +
    theme(text=element_text(size=28), axis.text = element_text(size=14),
          axis.title=element_text(size=12,face="bold"),
          panel.background = element_rect(fill = "transparent",colour = NA),
          plot.background = element_rect(fill = "transparent",colour = NA),
          legend.position = "None") +
    coord_cartesian(xlim=c(0,650), ylim=c(0,350)) +
    scale_x_continuous(name="Stem diameter (mm)", breaks=seq(0,650, by=100), limits=c(0,650, by=100)) +
    scale_y_continuous(name="Height (dm)", breaks=seq(0,350, by=50), limits=c(0,350))  +
    scale_fill_gradientn(name= vrbl, colors=c("yellow", "lightgreen", "darkgreen", "purple", "darkblue", "black"), values=cols, limits=c(min(cols), max(cols)), rescaler=function(x,...) x, na.value="white", guide="colorbar") +
    scale_color_gradientn(name= vrbl, colors=c("yellow", "lightgreen", "darkgreen", "purple", "darkblue", "black"), values=cols, limits=c(min(cols), max(cols)), rescaler=function(x,...) x, na.value="white", guide="colorbar") +
    geom_raster(data=table1, aes(x=as.numeric(d), y=as.numeric(h), fill=vrbl))
  
  if(view==TRUE){
    show(fig1)}
  
  if(save==TRUE){
    ggsave(filename = paste0(folder,groupName,"_",vrbl,".png"), plot = fig1, width = 6, height=7.5, device = png())
    dev.off()
  }
}
