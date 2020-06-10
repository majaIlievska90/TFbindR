#' Count the Occurrence of Sequence Motifs from Fimo File.
#' 
#' The function reads the input sequences from a FASTA file and identifies the
#' genes. Next, it splits the strands from the FIMO file into "+" and "-", and
#' counts the occurence of gene and motif combination for each strand
#' separately. The counts are reported in a matrix.
#' 
#' %% ~~ If necessary, more details than the description above ~~
#' 
#' @param fimo A fimo file containing the motifs of interest.
#' @param promoter.fasta Fasta file of the promoter sequences.
#' @param promoterDB The path, and the name of the object where results will be
#' saved.
#' @return A list named Strands with 3 elements, \item{plus}{A count matrix of
#' the occurence of motifs in the "plus" strand.} \item{minus}{A count matrix
#' of the occurence of motifs in the "minus" strand.} \item{both}{A count
#' matrix of the occurence of motifs in the "minus" and "plus" strands
#' together.}
#' @note %% ~~further notes~~
#' @author %% Maja Ilievska
#' @seealso %% ~~objects to See Also as \code{\link{help}}, ~~~
#' @references Grant, Charles E., Timothy L. Bailey, and William Stafford
#' Noble. "FIMO: scanning for occurrences of a given motif." Bioinformatics
#' 27.7 (2011): 1017-1018.
#' @examples
#' 
#' # Read a sample fimo file. The data is downloaded from the Arabidopsis Cistrome database.
#' fimo.dir = base::system.file("extdata", "JASPAR_cistrome.fimo.txt", package="TFbindR", mustWork = TRUE)
#' 
#' # Read a sample fasta file. 
#' promoter.dir = base::system.file("extdata","TAIR10_upstream_1000_translation_start_20101028.fa",package="TFbindR")
#' 
#' # Specify the output directory.
#' output.dir="."
#' 
#' # Run the function to create the DB and save it into R file type .rda which is automatically compressed. 
#' makePromoterDBfromFimo(fimo = fimo.dir, promoter.fasta=promoter.dir, promoterDB=file.path(output.dir,"TAIR10_500bp_upstream.fimo.db.rda"))
#' 
#' @export makePromoterDBfromFimo
makePromoterDBfromFimo <- function(fimo, promoter.fasta, promoterDB = "promoter.DB") {
    # Counts the occurence of motifs in a FIMO file. 
    # Args:
    #	fimo: a FIMO file, as given by the motif-based sequence analysis tool from MEME Suite. The same promoter sequences should be used in the FIMO analysis, as provided in the promoter.fasta file.		
    #	promoter.fasta: A fasta file of the promoter sequences. 
    #	promoterDB: The path, and the name of the object where results will be saved.
    # Returns: 
    #	#	A list of 3 elements, holding the count matrices for plus, minus and both strands together correspondingly. 
    library(seqinr)
    # Read the forward promoter sequences.
    Pr = read.fasta(promoter.fasta, seqtype = "DNA", as.string = TRUE)
    Allgenes = names(Pr)
    
    # Read the fimo file.
    D = read.delim(fimo, as.is = T)
    if (names(D)[2] == "motif_alt_id") {
        D[, 1] = apply(D[, 1:2], 1, paste, collapse = "_")
    }
    Motifs = unique(D[, 1])
    Nmotif = length(Motifs)
    
    # Split the plus and minus strands.
    D.split = split(D, D$strand)
    if (sum(names(D.split) %in% c("-", "+")) != 2) {
        message("Odd strand number ", length(D.split), " Check input file!")
        return(-1)
    }
    
    # Identify the target gene, stored in the column "sequence_name" in FIMO object. 
    seq.idx = grep("sequence", names(D))
    Strands = vector("list", 3)
    # For both strands, count the occurences of each gene-motif pair. 
    for (i in 1:2) {
        DD = D.split[[i]]
        out = matrix(0, length(Allgenes), Nmotif, dimnames = list(Allgenes, Motifs))
        genemotif = sort(apply(DD[, c(1, seq.idx)], 1, paste, collapse = ";"))
        motifcounts = rle(genemotif)
        for (j in 1:length(motifcounts$values)) {
            loc = strsplit(motifcounts$values[j], ";")[[1]]
            out[loc[2], loc[1]] = motifcounts$lengths[j]
        }
        Strands[[i]] = out
    }
    Strands[[3]] = Strands[[1]] + Strands[[2]]
    names(Strands) = c("minus", "plus", "both")
    for (i in 1:length(Strands)) Strands[[i]] = Strands[[i]][, which(colSums(Strands[[i]] > 
        0) > 0)]
    
    save(Strands, file = promoterDB)
    return(1)
}
