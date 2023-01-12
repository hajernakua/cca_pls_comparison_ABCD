Barplot_LV <- function(Data, column2plot = 1, threshold = 0,
                          title = NULL, sort = TRUE,
                          color4pos = "#612975", color4neg = "darkolivegreen",
                          font.size = 4, horizontal = FALSE,
                          ylim.min = NULL, ylim.max = NULL,
                          xaxis.textsize = 10,
                          xaxis.lwd = 1){
  if (sort){
    val2plot <- as.matrix(sort(setNames(Data[,column2plot], rownames(Data))))
  }else{
    val2plot <- Data[,column2plot]
  }
  
  if (is.null(ylim.min)){
    ylim.min <- min(min(Data[,column2plot])-0.15*diff(range(Data[,column2plot])), 0)
  }
  if (is.null(ylim.max)){
    ylim.max <- max(max(Data[,column2plot])+0.15*diff(range(Data[,column2plot])), 0)
  }
  PrettyBarPlot2(val2plot[,1],
                 threshold = threshold,
                 color4bar = ifelse(val2plot[,1] > 0, color4pos, color4neg), 
                 color4ns = ifelse(val2plot[,1] > 0, lighten(color4pos), lighten(color4neg)),
                 horizontal = horizontal,
                 font.size = font.size,
                 main = title
  )+
    ylim(ylim.min,ylim.max) +
    theme(axis.text.x = element_text(color = "black", size = xaxis.textsize),
          axis.line.x = element_line(size = xaxis.lwd))
  
}
