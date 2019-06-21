#written by noah friedman
#a template for R scripts for plotting

library(ggplot2)
library(grid)
require(cowplot)
library(egg)
library(dplyr)
library(data.table); setDTthreads(6)


plot_signature_waterfall <- function(df, title, includeNMut_Mb_Info=TRUE){
  
  plottingLevels <- c('Signature.AGE', 'Signature.APOBEC', 'Signature.HRD', 'Signature.SMOKING',
                      'Signature.5', 'Signature.MMR', 'Signature.UV', 'Signature.POLE',
                      'Signature.TMZ', 'Signature.POLE_plus_MMR','other')
  barColorPalette = c(
    "#00DFFF", "#FF0000", "#FF1493", "#FFA500", 
    "#FFB6C1", "#267574", "#FFF600", "#ADFF2F",
    "#2A52BE","#551A8B","#D3D3D3"
  )
  
  plt <- ggplot(df)+
    
    #The first bar of the predominant signature in the positive direction
    geom_bar(aes(x = reorder(Tumor_Sample_Barcode, orderingVal), y=signatureOfInterestMagnitude, 
                 fill = factor(signatureOfInterestName, levels=plottingLevels)), stat="identity")+
    
    #TODO FIX THE THIRD BAR SO IT 
    #THIRD BAR (plot it first so it is covered by the second bar)
    geom_bar(aes(x = reorder(Tumor_Sample_Barcode, orderingVal), y=-secondPredominantSigMagnitude - thirdPredominantSigMagnitude, 
                 fill = factor(thirdPredominantSigName, levels=plottingLevels)), stat = "identity")+ #color the lower signature column by which signature it is                                                             
    
    #The second bar of the second predominant signature in the negative direction
    geom_bar(aes(x = reorder(Tumor_Sample_Barcode, orderingVal), y=-secondPredominantSigMagnitude, 
                 fill = factor(secondPredominantSigName, levels=plottingLevels)), stat = "identity")+ #color the lower signature column by which signature it is                                                                                
    
    scale_fill_manual(values=barColorPalette, drop = FALSE)+
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size=10))+
    ylab('Signature Fraction')+
    ylim(-1,1)+
    guides(fill=guide_legend(title="Signature"))+
    ggtitle(title)
  if(includeNMut_Mb_Info == TRUE){
    #INCLUDE dots to mark mutation burden if needed
    plt <- plt + geom_point(aes(x = reorder(Tumor_Sample_Barcode, orderingVal), y=0, size=Nmut_Mb), alpha=0.75)
    return(plt)
    }
  else{
    return(plt)
  }
}

#insert the path to your signatures tsv file here
df <- read.table('~/Desktop/mnt/ifs/work/taylorlab/friedman/myAdjustedDataFiles/signatureWaterfallPlot.csv',sep = '\t', header=TRUE)

plt <- plot_signature_waterfall(df, title='Test title',
                                includeNMut_Mb_Info=TRUE) #set to false to just return the waterfall plot

ggsave('~/Desktop/plot.pdf', plot=plt,  width = 15, height = 10, units = c("in"))





