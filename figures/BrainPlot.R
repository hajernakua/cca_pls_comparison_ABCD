Brainplot <- function(Data, column2plot = 1,
                            palette = "RdBu", 
                            limits = c(-0.6, 0.6), 
                            values = c(0, 0.47, 0.5, 0.53, 1),
                            title = NULL,
                            brain.position = position_brain(side ~ hemi)){
  dk_brain_test <- c("bankssts", "caudal anterior cingulate", "caudal middle frontal", "cuneus", "entorhinal",
                     "fusiform", "inferior parietal", "inferior temporal", "isthmus cingulate", 
                     "lateral occipital", "lateral orbitofrontal", "lingual", "medial orbitofrontal", 
                     "middle temporal", "parahippocampal", "paracentral", "pars opercularis", "pars orbitalis", 
                     "pars triangularis","pericalcarine", "postcentral", "posterior cingulate", "precentral",
                     "precuneus", "rostral anterior cingulate", "rostral middle frontal", "superior frontal",
                     "superior parietal", "superior temporal", "supramarginal", "frontal pole", "temporal pole", 
                     "transverse temporal", "insula")
  
  Data <- as.data.frame(Data)
  col2plot.name <- colnames(Data)[column2plot]
  
  Data4brain <- dk %>% 
    as_tibble() %>% 
    left_join(tibble(Data, region = rep(dk_brain_test,2)))
  
  Data4brain <- as.data.frame(Data4brain)

  Data4brain %>% as.data.frame %>%
    ggplot() +
    geom_brain(mapping = aes(fill = Data4brain[, col2plot.name]),
               atlas = dk,
               position = brain.position
    ) +
    scale_fill_distiller(name = col2plot.name, palette = palette, limits = limits, values = values) +
    ggtitle(title) + 
    theme(axis.text.y.left = element_blank(), 
          axis.text.x.bottom = element_blank()) + 
    theme_brain(text.family = "Arial")
  
}
