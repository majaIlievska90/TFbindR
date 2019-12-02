#' Occurence of motifs in the gene promoter sequence.
#'
#' A dataset containing 3 matrices. Each matrix has the motifs on  
#' the rows and genes on the columns. The entries of the matrix are intergers
#' repsenting the number of times the motif occurs in the promoter
#' sequence of the gene given in the column. 
#' 
#'
#' @format An object with 3 matrices. 
#' \describe{
#'   \item{plus}{The occurence of motifs in the sense (plus) strand of the promoter sequence}
#'   \item{minus}{The occurence of motifs in the antisense (minus) strand of the promoter sequence}
#'  \item{both}{The occurence of motifs in the plus and minus strand summed together}
#' }
#' 
"Strands"