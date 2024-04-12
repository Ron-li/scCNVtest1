#' Create offtarget and ontarget reads
#'
#' This function uses the terminal to find overlap between fragments file and
#' peaks. Bedtools and Tabix need to be installed in order for this
#' preprocessing to occur
#'
#' @param x number
#' @return the x^2
#' @export
testfunc <- function(x){
  return(x^2)
}


#' Convert Genomic Tiles to a dataframe
#'
#' This function converts a GRanges object to a dataframe
#' with appropriate row names
#'
#' @param tiles GRanges Genomic Tiles Object
#' @return A dataframe with genomic features as row names
tileTodf <- function(tiles){
  df <- as.data.frame(tiles)[,c(1:4)]
  row_names <- do.call(paste, c(df[,c(1:3)], sep="-"))
  rownames(df) <- row_names
  return(df)
}
