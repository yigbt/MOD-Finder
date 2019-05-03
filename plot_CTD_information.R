
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
  plt <- plt + scale_fill_manual( values = c("skyblue", "salmon"), name = "Attribute", labels = c("Increased", "Decreased"))
  plt <- plt + xlab( "Gene Symbol") + ylab( "Reference Counts")
  plt <- plt + ggtitle( paste0("Changes in gene expression triggered by ", compound))

  return( plt)
  
}


#' This function plots a distribution of interaction actions that are known to be
#' caused by a certain chemical accroding to CTDbase.org
#' @param ctd_chem_object The CTD object created with ctd_chem_query().
#' @param compound The chemical that's being analyzed.
#' @param cas The CAS-RN of the chemical compound.
#' @return The ggplot object capturing the bar plot.
plot_interaction_actions <- function( ctd_chem, compound, cas){
  
  ctd_genes_tab <- get_table( ctd_chem, index_name = "gene interactions")
  ctd_genes_tab <- as.data.frame( ctd_genes_tab[ which( ctd_genes_tab$CAS.RN == cas), c(1,2,3,4,5,7)])

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

  ## create the pie chart diagrams 
  pPies <- ggplot(data=data2distr, aes(x="", y=n.x, fill=Attribute))
  pPies <- pPies + geom_bar(stat="identity", position=position_fill())
  pPies <- pPies + geom_text(aes(label=n.x), position=position_fill(vjust=0.5))
  pPies <- pPies + coord_polar(theta="y")
  pPies <- pPies + facet_wrap(~ Term)
  pPies <- pPies + scale_fill_manual(values=c("#999999", "skyblue", "salmon"))
#  pPies <- pPies + theme_minimal()
  pPies <- pPies + theme(legend.position='top')
  pPies <- pPies + guides(fill=guide_legend(nrow=2,byrow=TRUE))
  pPies <- pPies + ggtitle("Distribution of attributes and terms",
                           subtitle="- of InteractionAction column of Chemicals-gene association -")
  
  ## create the bar chart diagram
  pBars <- ggplot(data=data2distr, aes(x=term, y=n.x, fill=Attribute))
  pBars <- pBars + geom_bar(stat="identity")#, position=position_fill())
#  pBars <- pBars + scale_fill_discrete(name = "Attribute", labels = c("Affects", "Increased", "Decreased"))
  pBars <- pBars + scale_fill_manual(values=c("#999999", "skyblue", "salmon"), name="Attribute", labels=c("Affects", "Increased", "Decreased"))
  #pBars <- pBars + geom_text(aes(label=n.x), vjust=0.5, size=3.5)
#  pBars <- pBars + blank_theme_ticks
  pBars <- pBars + theme(axis.text.x = element_text(angle=45, hjust=1))
  pBars <- pBars + ylab( "Reference Counts")
#  pBars <- pBars + theme(legend.position='top')
  pBars <- pBars + ggtitle( paste0( "Distribution of attributes and terms of InteractionAction column of Chemicals-gene association caused by ", compound))

  return( pBars)
  
}