#' Return the description of the stacked bar plot that visualizes
#' the groups of InteractionActions that are linked with the user given compound input.
#' @param chemical The chemical that triggered the InteractionActions
#' @return Text that explains the bar plot structure
get_generalInformationDescription <- function( chemical){
  
  string <- paste0( "In the figure below, all mentioned 'InteractionActions' (IAs) are clustered that are collected by the CTDbase with regard to ", 
                    chemical,
                    " IAs group interactions into classes and can be combined such that multiple 
                    IAs can be used to describe chemical-gene relations. The different IAs contain one or more interaction terms,
                    separated by '|' character, from a dictionary of 79 entries. Each of these terms are mostly prefixed by an 
                    attribute such as 'increases', 'decreases', 'affects', yielding IAs in the form of 
                    'increases^expression | affects^activity'. All IAs that are connected with chemical-gene interactions triggered by ",
                    chemical,
                    " are displayed in the following plot, where  the number of references of each IA is shown in a stacked bar plot.")
  
  return( string)
  
}


#' Return the description of the sideways bar plot that displays
#' the top40 genes and their respective number of references.
#' @param chemical The chemical that was used as search term
#' @return Text that explains the structure of the plot
get_chemicalGeneDescription <- function( chemical){
  
  string <-"With a focus on the effects triggered by the exposure, the number of references that mention either a deacrease
           or increase in expression changes where plotted for the top40 genes (highest number of references). Publications 
           that merely mention an 'affect' rather than specifying the direction are addtionally visualized at the right side of the figure."
  
  return( string)
  
}



#' Return the description of the bar plot illustrating the most
#' referenced diseases that were linked with a certain compound.
#' @param chemical The chemical that was used as search term
#' @return Text that explains the structure of the plot
get_diseasesDescription <- function( chemical){
  
  string <- paste0("The following figure summarizes those diseases that have the largest reference counts and/or inference score associating a connection between ", chemical,
            " and the disease. Each of the CTDbase collected associations is either curated or inferred (via a curated gene interaction)",
            " The unique list of both top40 diseases with the most reference counts AND the top40 list of diseases with the highest inference score,
            is displayed as bar plot. The bar length indicates the actual number of references while the color scheme indicates the inference score.
            In brief, the inference score reflects the degree of similarity between CTD chemical–gene–disease networks and a similar scale-free random network. 
            The higher the score, the more likely the inference network has atypical connectivity. For more information about inference scores, please see: 
            King et al, 2012, PLoS One, 'Ranking Transitive Chemical-Disease Inferences Using Local Network Topology in the Comparative Toxicogenomics Database.'")
  
  return( string)
  
}



#' Return the description of the bar plot illustrating the most
#' enriched KEGG pathways.
#' @param chemical The chemical that was used as search term
#' @return Text that explains the structure of the plot
get_pathwaysDescription <- function( chemical){
  
  string <- paste0( "Finally, the top40 KEGG and REACTOME represented pathways, based on their adjusted P-Values are shown. The length of each bar 
                    reflects the ration between enriched genes and the total number of annotated genes in this respective pathway. Based on the CTDbase,
                    a pathway is considered 'enriched' if the proportion of genes annotated to it in a test set is significantly larger than the 
                    proportion of all genes annotated to it in the genome.", 
                    " The color scale of each bar reflects the amount of annotated genes in each pathway. The actual adjusted P-Value is stated in the
                    pathway label on the y axis.")
  
  return( string)
  
}