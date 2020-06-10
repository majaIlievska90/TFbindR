#' makePromoterDBfromRE
#' 
#' Identify and count the occurence of sequence motifs, described by regular
#' expressions, in a set of promoter sequences. Both, sense and antisense
#' strands are scanned for the occurence of each motif and the number of
#' matches is reported.
#' 
#' %% ~~ If necessary, more details than the description above ~~
#' 
#' @param Motifs a character vector of the motifs to be matched. A motif can be
#' defined with the full sequences or with special characters.
#' @param promoter.fasta Fasta file of the promoter sequences.
#' @param promoterDB Output file where the results will be saved.
#' @return A list names Strands with 3 elements, \item{plus}{A count matrix of
#' the occurence of motifs in the "minus" strand.} \item{minus}{A count matrix
#' of the occurence of motifs in the "plus" strand.} \item{both}{A count matrix
#' of the occurence of motifs in the "minus" and "plus" strands together.}
#' @note %% ~~further notes~~
#' @author %% ~~who you are~~
#' @seealso %% ~~objects to See Also as \code{\link{help}}, ~~~
#' @references %% ~put references to the literature/web site here ~
#' @examples
#' 
#' # Read a file of sample sequence motifs.
#' motifs.dir = system.file("extdata", "REcomplete.txt", package="TFbindR", mustWork = TRUE)
#' Motifs=read.delim(file = motifs.dir ,header=TRUE, sep="\t", as.is=T)
#' 
#' # Remove duplicated motifs.
#' Motifs=Motifs[which(!duplicated(Motifs$Motifs)),]
#' 
#' # Read a sample fasta file. 
#' promoter.dir = system.file("extdata","TAIR10_upstream_1000_translation_start_20101028.fa",package="TFbindR")
#' 
#' # Specify the output directory.
#' output.dir="."
#' 
#' # Run the function to create the DB and save it into R file type .rda which is automatically compressed
#' makePromoterDBfromRE(Motifs = Motifs$Motif,promoter.fasta = promoter.dir, promoterDB=file.path(output.dir,"TAIR10_500bp_upstream.RE.db.rda"))
#' 
#' @export makePromoterDBfromRE
makePromoterDBfromRE <- function(Motifs, promoter.fasta, promoterDB = "promoter.DB") {
	# Counts the occurence of motifs, as defined by regular expressions, in promoter regions. 
	
	# Args: 
	#	Motifs: A character vector of sequence motifs. The length of a seuence motif a motif should not exceed the length of the promoter sequence.  
	#	promoter.fasta: A fasta file of the promoter sequences. 
	#	promoterDB: The path, and the name of the object where the result will be stored. 
	#
	# Returns: 
	#	A list of 3 elements, holding the count matrices for plus, minus and both strands together correspondingly. 
    library(seqinr)
    Nmotif = length(Motifs)
    # Read the forward promoter sequences.
    Pr = read.fasta(promoter.fasta, seqtype = "DNA", as.string = TRUE)
    Allgenes = names(Pr)
    
    # Read the reverse complement of the promoter sequences. 
    Pr2 = vector("list", length(Allgenes))
    for (i in 1:length(Allgenes)) {
        Pr2[i] = c2s(rev(comp(unlist(getSequence(Pr[i])))))
    }
    names(Pr2) = names(Pr)
    
    # Create the plus and minus strand matrices, populate the elements with the number of times the motif is found in the promoter sequence.  
    MinusStr = matrix(0, length(Allgenes), Nmotif, dimnames = list(Allgenes, Motifs))
    PlusStr = matrix(0, length(Allgenes), Nmotif, dimnames = list(Allgenes, Motifs))
    for (i in 1:Nmotif) {
        t1 = gregexpr(Motifs[i], Pr, ignore.case = T)
        MinusStr[, i] = sapply(t1, function(x) {
            if (any(x != -1)) {
                return(length(x))
            } else {
                return(0)
            }
        })
        t2 = gregexpr(Motifs[i], Pr2, ignore.case = T)
        PlusStr[, i] = sapply(t2, function(x) {
            if (any(x != -1)) {
                return(length(x))
            } else {
                return(0)
            }
        })
        message(i, "/", Nmotif)
    }
    # Sum the minus and plus into the Both strands count matrix.
    Both = MinusStr + PlusStr
    rownames(Both) = Allgenes
    colnames(Both) = Motifs
    MinusStr = MinusStr[, which(colSums(MinusStr > 0) > 0)]
    PlusStr = PlusStr[, which(colSums(PlusStr > 0) > 0)]
    Both = Both[, which(colSums(Both > 0) > 0)]
    Strands = vector("list", 3)
    Strands[[1]] = MinusStr
    Strands[[2]] = PlusStr
    Strands[[3]] = Both
    names(Strands) = c("minus", "plus", "both")
    save(Strands, file = promoterDB)
    return(1)
}
