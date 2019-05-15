
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
  aff <- df[ grep( "affects\\^expression", df$Interaction.Actions), c( 10, 4, 8)]
  
  ## summarize the reference counts
  inc2 <- inc %>% dplyr::group_by( Gene.Symbol) %>% dplyr::summarise( increase.Ref = sum( Reference.Count)) %>% dplyr::arrange( desc( increase.Ref))
  dec2 <- dec %>% dplyr::group_by( Gene.Symbol) %>% dplyr::summarise( decrease.Ref = sum( Reference.Count)) %>% dplyr::arrange( desc( decrease.Ref))
  aff2 <- aff %>% dplyr::group_by( Gene.Symbol) %>% dplyr::summarise( affects.Ref = sum( Reference.Count)) %>% dplyr::arrange( desc( affects.Ref))
  genes <- unique( c( inc2$Gene.Symbol, dec2$Gene.Symbol))
  
  
  ## combine both tibbles to a single data.frame
  ## change the reference counts to log scale and transform the axis accordingly
  max <- max( c( inc2$increase.Ref, dec2$decrease.Ref))
  logAxis <- get_log_labels( max, mirror = TRUE)
  
  plotme <- data.frame( Gene.Symbol = genes, increase.Ref = rep( 0, length( genes)), decrease.Ref = rep( 0, length( genes)), abs = rep(0, length( genes)))
  
  if( max > 10){
    inc2$increase.Ref[ which( inc2$increase.Ref == 1)] <- sqrt( 2)
    dec2$decrease.Ref[ which( dec2$decrease.Ref == 1)] <- sqrt( 2)
    plotme$increase.Ref[ match( inc2$Gene.Symbol, plotme$Gene.Symbol)] <- log( inc2$increase.Ref, base = 2)
    plotme$decrease.Ref[ match( dec2$Gene.Symbol, plotme$Gene.Symbol)] <- log( dec2$decrease.Ref, base = 2)
  } else{
    plotme$increase.Ref[ match( inc2$Gene.Symbol, plotme$Gene.Symbol)] <- inc2$increase.Ref
    plotme$decrease.Ref[ match( dec2$Gene.Symbol, plotme$Gene.Symbol)] <- dec2$decrease.Ref
  }
  plotme$abs <- pmax( plotme$increase.Ref, plotme$decrease.Ref)
  plotme <- plotme %>% dplyr::arrange( desc(abs)) %>% slice(1:40)
  plotme$decrease.Ref <- -1 * plotme$decrease.Ref
  plotme.df <- melt( plotme[, c(1,2,3)], value.name = "Reference.Count", variable.name = "Direction")
  plotme.df$Direction <- ifelse( plotme.df$Direction == "increase.Ref", "Increased", "Decreased")
  
  
  plt <- ggplot( data = plotme.df, aes( x = Gene.Symbol, y = Reference.Count, fill = Direction)) + geom_bar(stat="identity") + coord_flip()
  plt <- plt + scale_y_continuous( expand = c(0,0), limits = c( -1*logAxis$limit, logAxis$limit), breaks=logAxis$breaks, labels=logAxis$labels)
  plt <- plt + scale_fill_manual( values = c("skyblue", "salmon", "#999999"), name = "Attribute", labels = c("Increased", "Decreased", "Affects"), drop = FALSE, limits = c( "Increased", "Decreased", "Affects"))
  plt <- plt + xlab( "Gene Symbol") + ylab( "Reference Counts (Decreased/Increased)")
  plt <- plt + ggtitle( paste0("Changes in gene expression triggered by ", compound))

  
  ## create a separated plot for affects on gene expression
  ## first, shrink the dataframe to those genes plotted in the above plot
  aff_plotme <- data.frame( Gene.Symbol = plotme$Gene.Symbol, affects.Ref = rep(0, length( plotme$Gene.Symbol)), affects.RefLog = rep( 0, length( plotme$Gene.Symbol)))
  m <- match( aff_plotme$Gene.Symbol, aff2$Gene.Symbol)
  n <- is.na(m)
  m <- m[ !n]
  aff_plotme$affects.RefLog[ !n] <- log( aff2$affects.Ref[m], base = 2)
  aff_plotme$affects.Ref[ !n] <- aff2$affects.Ref[m]
  aff_plotme$Direction <- rep( "Affects", nrow( aff_plotme))
  
  ## adapt the axis
  max <- max( aff_plotme$affects.Ref)
  
  if( max > 10){
    logAxis <- get_log_labels( max)
    
    aff_plt <- ggplot( data = aff_plotme, aes( x = Gene.Symbol, y = affects.RefLog, fill = Direction))
    
  } else{
    logAxis <- list( labels = seq(0,10), breaks = seq( 0, 10), limit = 10)  
    aff_plt <- ggplot( data = aff_plotme, aes( x = Gene.Symbol, y = affects.Ref, fill = Direction))
    
  }
  
  aff_plt <- aff_plt + geom_bar( stat = "identity") + coord_flip()
  aff_plt <- aff_plt + scale_y_continuous( expand = c(0,0), limits = c( 0, logAxis$limit), breaks=logAxis$breaks, labels=logAxis$labels)
  aff_plt <- aff_plt + scale_fill_manual( values = c("skyblue", "salmon", "#999999"), name = "Attribute", labels = c("Increased", "Decreased", "Affects"), drop = FALSE, limits = c("Increased", "Decreased", "Affects"))
  aff_plt <- aff_plt + xlab( "") + ylab( "Reference Counts (Affects)")
  aff_plt <- aff_plt + ggtitle( "")
  
  final_plt <- ggarrange( plt, aff_plt, common.legend = TRUE, widths = c(5,2), legend = "bottom")

  
  return( final_plt)
  
}


#' This function plots a distribution of interaction actions that are known to be
#' caused by a certain chemical accroding to CTDbase.org
#' @param ctd_chem_object The CTD object created with ctd_chem_query().
#' @param compound The chemical that's being analyzed.
#' @param cas The CAS-RN of the chemical compound.
#' @return The ggplot object capturing the bar plot.
plot_interaction_actions <- function( ctd_chem, compound, cas){
  
  ctd_genes_tab <- get_table( ctd_chem, index_name = "gene interactions")
  if( !is.null(cas)){
    ctd_genes_tab <- as.data.frame( ctd_genes_tab[ which( ctd_genes_tab$CAS.RN == cas), c(1,2,3,4,5,7)])
  }else{
    ctd_genes_tab <- as.data.frame( ctd_genes_tab[ , c(1,2,3,4,5,7)])
  }

  ## How many different interactions are there
  ia <- ctd_genes_tab$Interaction.Actions
  iau <- unique(ia)
  ia_count <- iau %>% length()

  splitted <- unlist(lapply(iau, function(x) {
    str_split(string=x, pattern='\\|', simplify=FALSE)
  }))
  
  IA <- data.frame(InteractionActions=splitted) %>%
    separate(InteractionActions, sep="\\^",
             into=c("Attribute", "Term"),
             remove=FALSE, fill="right")
  IAtermcount <- IA %>% select(Term) %>% group_by(Term) %>% count()
  IAtermattrcount <- IA %>% select(Term, Attribute) %>% group_by(Attribute, Term) %>% count()
  IAdistr <- left_join(IAtermattrcount, IAtermcount, by="Term") %>% arrange(desc(n.y), Term)
  factorlevels <- IAtermcount %>% arrange(desc(n)) %>% pull(Term) 
  IAdistr$term <- factor(IAdistr$Term,
                         levels=factorlevels)
  
  
  termdict <- unique(IA$term)
  attrdict <- unique(IA$attribute)
  
  map <- do.call(rbind, lapply(iau, function(l) {
    #l <- iau[i]
    ls <- unlist(str_split(string=l, pattern='\\|', simplify=FALSE))
    data.frame("Interaction.Actions"=l, "str"=ls) %>%
      separate(str, sep="\\^",
               into=c("Attribute", "Term"),
               remove=FALSE, fill="right")
  }))
  
  data <- left_join(ctd_genes_tab, map)
  data2 <- data %>% select("Gene.Symbol", "Term", "Attribute") %>% unique()
  
  data2_termcount <- data2 %>% select(Term) %>% group_by(Term) %>% count()
  data2_termattrcount <- data2 %>% select(Term, Attribute) %>%
    group_by(Attribute, Term) %>% count()
  data2distr <- left_join(data2_termattrcount, data2_termcount, by="Term") %>%
    arrange(desc(n.y), Term)
  factorlevels <- data2_termcount %>% arrange(desc(n)) %>% pull(Term)
  data2distr$term <- factor(data2distr$Term, levels=factorlevels)


  ## create the bar chart diagram
  pBars <- ggplot(data=data2distr, aes(x=term, y=n.x, fill=Attribute))
  pBars <- pBars + geom_bar(stat="identity")#, position=position_fill())
  pBars <- pBars + scale_fill_manual(values=c("#999999", "skyblue", "salmon"), name="Attribute", labels=c("Affects", "Increased", "Decreased"))
  pBars <- pBars + theme(axis.text.x = element_text(angle=45, hjust=1))
  pBars <- pBars + ylab( "Reference Counts") + xlab( "Interaction Term")
  pBars <- pBars + ggtitle( paste0( "Distribution of attributes and terms of InteractionAction column of Chemicals-gene association caused by ", compound))

  return( pBars)
  
}


#' This function plots a summary of diseases that are known to be
#' linked with a certain chemical accroding to CTDbase.org.
#' @param ctd_chem_object The CTD object created with ctd_chem_query().
#' @param compound The chemical that's being analyzed.
#' @param cas The CAS-RN of the chemical compound.
#' @return The ggplot object capturing the bar plot.
plot_diseases <- function( ctd_chem, compound, cas){

  disease <- as.data.frame( get_table( ctd_chem, index_name = "diseases"))
  if( !is.null( cas)){
    disease <- disease[ which( disease$CAS.RN == cas), ]
  }
  
  ## use the top40 diseases for plotting, based on reference count and inference score
  rfcount <- disease[, c(5,8,9)] %>% arrange( desc( Reference.Count)) %>% slice( 1:40) %>% select( Disease.ID)
  infscore <- disease[, c(5,8,9)] %>% arrange( desc( Inference.Score)) %>% slice( 1:40) %>% select( Disease.ID)
  
  df <- disease[ match( unique( c( rfcount$Disease.ID, infscore$Disease.ID)), disease$Disease.ID), c(2,3,4,5,6,8,9)] %>% arrange( desc( Reference.Count))
  df$status <- as.factor( ifelse( df$Direct.Evidence == "", "Inferred", "Curated"))
  
  ## change all Reference Counts of 1 to sqrt( 2)
  ## to avoid 0 in the log-transformed scale
  df$Reference.Count[ which( df$Reference.Count == 1)] <- sqrt(2)
  
  max <- max( df$Reference.Count)

  logAxis <- get_log_labels( max)

  plt <- ggplot( data = df, aes( y = log(Reference.Count, base = 2), x = reorder( Disease.Name, Reference.Count), fill = Inference.Score)) + geom_bar( stat="identity") + coord_flip()

  ## check if the status has empty entries, if that's not the case
  ## apply a facet grid to devide diseses into curated and inferred
  if( summary( df$status)[1] != 0 && summary( df$status)[2] != 0){
    plt <- plt + facet_grid( status ~ ., scales = "free_y", space = "free_y")
  }
  
  plt <- plt + scale_y_continuous( expand = c(0,0), limits = c( 0, logAxis$limit), breaks=logAxis$breaks, labels=logAxis$labels)
  plt <- plt + scale_fill_distiller( palette = "YlOrRd", name = "Inference Score")
  plt <- plt + xlab( "Disease Name") + ylab( "Reference Counts")
  plt <- plt + ggtitle( paste0("Diseases that are asscociated with ", compound, " or its descendants"))

  return( plt)
  
}


#' This function plots a summary of enriched pathways that are known to be
#' linked with a certain chemical accroding to CTDbase.org.
#' @param ctd_chem_object The CTD object created with ctd_chem_query().
#' @param compound The chemical that's being analyzed.
#' @param cas The CAS-RN of the chemical compound.
#' @return The ggplot object capturing the bar plot.
plot_pathways <- function( ctd_chem, compound, cas){
  
  path <- as.data.frame( get_table( ctd_chem, index_name = "kegg pathways"))
  
  df <- path %>% arrange(Corrected.P.value) %>% slice( 1:40) %>% select( Pathway, Corrected.P.value, Annotated.Genes.Quantity, Genome.Frequency)
  df$Genes <- as.integer( gsub( "/.*", "", df$Genome.Frequency))
  df$Frequence <- df$Annotated.Genes.Quantity / df$Genes
  df$Label <- paste0( df$Pathway, " (FDR: ", df$Corrected.P.value, ")")

  plt <- ggplot( data = df, aes( y = Frequence, x= reorder( Label, -Corrected.P.value), fill = Genes)) + geom_bar( stat = "identity") + coord_flip()
  plt <- plt + ylim( c(0,1)) + xlab( "KEGG Pathway") + ylab( "Ratio( enriched genes / annotated genes)")
  plt <- plt + scale_fill_distiller( palette = "YlOrRd", name = "Annotated Genes")
  plt <- plt + ggtitle( "Pathway enrichment analysis based on KEGG pathways")
  
  return( plt)
  
  
}


#' Design the axis labels for the plots
#' depending on the actual plot, axis label can be mirrored,
#' axis limit vary between 10, 100, 1000, and over 1000
#' all axis with a maximal value above 10 are created in logscale.
#' @param max Specify the maximal value that should be plottedt on this axis.
#' @param mirror Optional, if TRUE, axis labels and breaks will be mirrored at 0.
#' @return List containing the axis labels and breaks. 
get_log_labels <- function( max, mirror = FALSE){
  
  root <- sqrt(2)
  
  if( mirror == TRUE ){
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
      
    } else if( max > 10){
      
      logbreak <- log( c(1, 10, 100), base =2)
      breaks <- c( -1 * rev( logbreak), logbreak[2:3])
      labels <- c( "100", "10", "0", "10", "100") 
      limit <- 7
      
    }else{
      
      labels <- c( seq( 10, 0), seq( 1, 10))
      breaks <- c( seq( -10, 10))
      limit <- 10
      
    }
    
    
  } else {
    
    
    if( max > 1000 ){
      
      breaks <- log( c(1, root, 2, 4, 6, 8, 10, 20, 40, 60, 80, 100, 200, 400, 600, 800, 1000, 2000, 4000, 6000, 8000, 10000), base =2)
      labels <- c( "0", "1", "2", "4", "6", "8", "10", "20", "40", "60", "80", "100", "200", "400", "600", "800", "1000", "2000", "4000", "6000", "8000", "10000")
      limit <- 14
      
      
    } else if( max > 100){
      
      breaks <- log( c(1, root, 2, 4, 6, 8, 10, 20, 40, 60, 80, 100, 200, 400, 600, 800, 1000), base =2)
      labels <- c( "0", "1", "2", "4", "6", "8", "10", "20", "40", "60", "80", "100", "200", "400", "600", "800", "1000") 
      limit <- 10
      
    } else{
      
      breaks <- log( c(1, root, 2, 4, 6, 8, 10, 20, 40, 60, 80, 100), base =2)
      labels <- c( "0", "1", "2", "4", "6", "8", "10", "20", "40", "60", "80", "100") 
      limit <- 7
      
    }
  }
  
  return( logAxis <- list( labels = labels, breaks = breaks, limit = limit))
  
}
