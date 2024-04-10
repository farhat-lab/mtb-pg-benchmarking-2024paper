#' plot_gain_loss
#' 
#' @importFrom rlang .data
#'
#' @description Plots a expected gain and loss events onto the branches of the phylogeny.
#'
#' @param fit a fitted pangenome model output by running 'panstripe'
#' @param tip_label whether or not to show tree tip labels (default=TRUE)
#' @param text_size adjusts the size of text in the plot
#' @param color_pallete the colour pallete to use. A number between 1 & 9. See 'scale_colour_brewer' for more details
#'
#' @return a plot of the phylogeny coloured by the inferred total gene gain/loss events per branch 
#'
#' @examples
#'
#' sim <- simulate_pan(rate=1e-3)
#' fA <- panstripe(sim$pa, sim$tree, nboot=0)
#' plot_gain_loss_v2(fA, color_pallete=7)
#'
#' @export
plot_gain_loss_v2 <- function(fit,
                           change_type = "all",  # New parameter to choose the trait
                           tip_label=TRUE,
                           text_size=14,
                           color_pallete=7,
                           aspect_ratio=1.5){  # Added parameter for aspect ratio
  
  #check inputs
  #if (class(fit) != 'panfit') stop('fit is not of class `panfit`!')
  #validate_panfit(fit)
  
  # Choose the trait based on trait_type
  if (change_type == "all") {
    trait_to_plot <- fit$datafull$acc
  } else if (change_type == "loss") {
    trait_to_plot <- fit$datafull$acc_loss
  } else if (change_type == "gain") {
    trait_to_plot <- fit$datafull$acc_gain
  } else {
    stop("Invalid change_type Options are 'all', 'gain', or 'loss'.")
  }
  
  
  stopifnot(all(fit$data$core==fit$tree$edge.length))
  

  # Choose the trait based on trait_type
  if (change_type == "all") {
      gt <- dplyr::full_join(ggtree::fortify(fit$tree), 
                            data.frame(node = fit$tree$edge[,2],
                                       trait = fit$datafull$acc), 
                            by = 'node')
      legend_label = 'Total genes\ngained & lost'  

  } else if (change_type == "loss") {

      gt <- dplyr::full_join(ggtree::fortify(fit$tree), 
                            data.frame(node = fit$tree$edge[,2],
                                       trait = fit$datafull$acc_loss), 
                            by = 'node')
      legend_label = 'Total genes\nlost'  

  } else if (change_type == "gain") {

      gt <- dplyr::full_join(ggtree::fortify(fit$tree), 
                            data.frame(node = fit$tree$edge[,2],
                                       trait = fit$datafull$acc_gain), 
                            by = 'node')
      
      legend_label = 'Total genes\ngained'  

  } else {
    stop("Invalid change_type Options are 'all', 'gain', or 'loss'.")
  }
  




  gg <- ggtree::ggtree(gt, ggplot2::aes(color=.data$trait), size=1) +
    ggplot2::labs(colour=legend_label) +
    ggplot2::scale_color_binned(type = 'viridis') +
    ggplot2::theme(aspect.ratio = aspect_ratio)  # Set the aspect ratio
  
  if (tip_label){
    gg <- gg + ggtree::geom_tiplab(align=TRUE, colour='black')
  }
    
  return(gg)
}