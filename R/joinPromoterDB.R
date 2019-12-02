#' Join promoterDB objects.
#' 
#' Combing together promoter.DB objects into one promoter.DB.
#' 
#' %% ~~ If necessary, more details than the description above ~~
#' 
#' @param DB1 prmoterDB object
#' @param DB2 promoterDB object
#' @param out.DB the output file for the combined promoterDB object.
#' @return A promoterDB object - list with 3 elements, \code{plus} - a count
#' matrix of the occurence of motifs in the "plus" strand, \code{minus} - a
#' count matrix of the occurence of motifs in the "minus" strand and
#' \code{both} - a count matrix of the occurence of motifs in the "minus" and
#' "plus" strands together.
#' @note %% ~~further notes~~
#' @author %% ~~who you are~~
#' @seealso %% ~~objects to See Also as \code{\link{help}}, ~~~
#' @references %% ~put references to the literature/web site here ~
#' @examples
#' 
#' 
#' @export joinPromoterDB
joinPromoterDB <- function(DB1, DB2, out.DB) {
	# Joins 2 promoterDB objects together. 
	# Args: 
	#	DB1: promoterDB object
	#	DB2: promoterDB object 
	#	out.DB: output file for the resulting object.
	# Returns: 
    load(DB1)
    Strands1 = Strands
    load(DB2)
    Strands2 = Strands
    Strands = list()
    for (i in 1:3) {
        Strands[[i]] = cbind(Strands1[[i]], Strands2[[i]])
    }
    names(Strands) = c("minus", "plus", "both")
    save(Strands, file = out.DB)
}
