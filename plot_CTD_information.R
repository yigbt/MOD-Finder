
#' This function plots the Top40 genes that have changed expression profiles 
#' due to a certain chemical accroding to CTDbase.org
#' @param ctd_chem_object The CTD object created with ctd_chem_query().
#' @param compound The chemical that's being analyzed.
#' @return The ggplot object capturing the plot.
plot_chemical_gene_interaction <- function( ctd_chem, compound){

  ## get table of gene interactions
  df <- data.frame( get_table( ctd_chem, index_name = "gene interactions"))
  s <- unlist( strsplit( df$Interaction.Actions, "\\|"))
  dfs <- data.frame( interaction = s, stringsAsFactors = TRUE)
  summary( dfs)
  
  
  ## Plan: extract all references from df that correspond to increased and decreased expression
  ## combine rows based on their interaction type (decreased or increased) and their gene symbol
  ## because sometimes the mRNA of a gene is mentioned and sometimes the protein is referenced etc...
  ## for simplicity let's put them together
  ## then extract the top40 or whatever genes, based on their cumulative reference count and plot a horizontal bar chart.
  inc <- df[ grep( "increases\\^expression", df$Interaction.Actions), c(10, 4, 8)]
  dec <- df[ grep( "decreases\\^expression", df$Interaction.Actions), c(10, 4, 8)]
  
  ## summarize the reference counts
  inc2 <- inc %>% dplyr::group_by( Gene.Symbol) %>% dplyr::summarise( increase.Ref = sum( Reference.Count)) %>% dplyr::arrange( desc( increase.Ref))
  dec2 <- dec %>% dplyr::group_by( Gene.Symbol) %>% dplyr::summarise( decrease.Ref = sum( Reference.Count)) %>% dplyr::arrange( desc( decrease.Ref))
  genes <- unique( c( inc2$Gene.Symbol, dec2$Gene.Symbol))
  
  
  # combine both tibble to a single data.frame
  plotme <- data.frame( Gene.Symbol = genes, increase.Ref = rep( 0, length( genes)), decrease.Ref = rep( 0, length( genes)), abs = rep(0, length( genes)))
  plotme$increase.Ref[ match( inc2$Gene.Symbol, plotme$Gene.Symbol)] <- log( inc2$increase.Ref, base = 2)
  plotme$decrease.Ref[ match( dec2$Gene.Symbol, plotme$Gene.Symbol)] <- log( dec2$decrease.Ref, base = 2)
  plotme$abs <- pmax( plotme$increase.Ref, plotme$decrease.Ref)
  plotme <- plotme %>% dplyr::arrange( desc(abs)) %>% slice(1:40)
  plotme$decrease.Ref <- -1 * plotme$decrease.Ref
  plotme.df <- melt( plotme[, c(1,2,3)], value.name = "Reference.Count", variable.name = "Direction")
  
  max <- max( c( inc2$increase.Ref, dec2$decrease.Ref))
  
  if( max > 1000 ){
    
    logbreak <- log( c(1, 10, 100, 1000, 10000), base =2)
    breaks <- c( -1 * rev( logbreak), logbreak[2:5])
    labels <- c( "10000", "1000", "100", "10", "0", "10", "100", "1000", "10000")
    limit <- 14
    
    
  } else if( max > 100){
    
    logbreak <- log( c(1, 10, 100, 1000), base =2)
    breaks <- c( -1 * rev( logbreak), logbreak[2:4])
    labels <- c( "1000", "100", "10", "0", "10", "100", "1000") 
    limit <- 10
    
  } else{
    
    logbreak <- log( c(1, 10, 100), base =2)
    breaks <- c( -1 * rev( logbreak), logbreak[2:3])
    labels <- c( "100", "10", "0", "10", "100") 
    limit <- 7
    
  }
  
  plt <- ggplot( data = plotme.df, aes( x = Gene.Symbol, y = Reference.Count, fill = Direction)) + geom_bar(stat="identity") + coord_flip()
  plt <- plt + scale_y_continuous( expand = c(0,0), limits = c( -1*limit, limit), breaks=breaks, labels=labels)
  plt <- plt + scale_fill_discrete(name = "Change in Expression", labels = c("Increased", "Decreased"))
  plt <- plt + ggtitle( paste0("Changes in gene expression triggered by ", compound))

  return( plt)
  
}
