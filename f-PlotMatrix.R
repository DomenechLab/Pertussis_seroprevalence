#######################################################################################################
# Function to visualize a contact matrix
#######################################################################################################

PlotMatrix <- function(M_in, plot_title = "") {
  
  #browser()
  if(!is.null(colnames(M_in))) colnames(M_in) <- NULL
  M_in <- as.matrix(M_in)
  
  M_long <- reshape2::melt(M_in)
  
  pl <- ggplot(data = M_long, mapping = aes(x = Var1, y = Var2, fill = value)) + 
    geom_tile() + 
    #scale_fill_gradient(low = "white", high = "blue", trans = "identity") +
    scale_fill_viridis(option = "rocket", trans = "sqrt") + 
    theme_minimal() + 
    theme(legend.position = "right", panel.grid = element_blank()) + 
    labs(x = "Age 1", y = "Age 2", fill = "Contact rate", title = plot_title)
  print(pl)
}