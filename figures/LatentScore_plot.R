LxLyplot <- function(Lx, Ly, 
                           column2plot.Lx = 1,
                           column2plot.Ly = 1,
                           title = NULL,
                           Name4X = "X",
                           Name4Y = "Y",
                           textsize = 10,
                           col.line = "navy",
                           col.points = "royalblue3",
                           line.width = 1.5,
                           alpha.points = 0.2){
  scores2plot <- cbind(as.matrix(Lx)[,column2plot.Lx], as.matrix(Ly)[,column2plot.Ly])
  lxplot <- createFactorMap(scores2plot,
                            col.points = col.points,
                            pch = 16,
                            display.labels = FALSE,
                            col.background = NULL,
                            col.axes = "gray60",
                            alpha.axes = 0.5,
                            alpha.points = alpha.points,
                            title = title)
  
  lxplot$zeMap +
    xlab(paste0("Behavior Latent Scores")) +
    ylab(paste0("Brain Latent Scores")) + 
    theme(axis.title = element_text(size=textsize), 
          axis.text.x = element_text(size=textsize),
          axis.text.y = element_text(size=textsize)) + 
    geom_smooth(method=lm, se=TRUE, color = col.line, lwd = line.width)
  
}
