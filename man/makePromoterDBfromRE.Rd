\name{makePromoterDBfromRE}
\alias{makePromoterDBfromRE}
\title{makePromoterDBfromRE
}
\description{
Identify and count the occurence of sequence motifs, described by regular expressions, in a set of promoter sequences. Both, sense and antisense strands are scanned for the occurence of each motif and the number of matches is reported. 
}
\usage{
makePromoterDBfromRE(Motifs, promoter.fasta, promoterDB = "promoter.DB")
}
\arguments{
  \item{Motifs}{
a character vector of the motifs to be matched. A motif can be defined with the full sequences or with special characters. 
}
  \item{promoter.fasta}{
Fasta file of the promoter sequences.
}
  \item{promoterDB}{
Output file where the results will be saved.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list names Strands with 3 elements, \item{plus}{A count matrix of the occurence of motifs in the "minus" strand.} \item{minus}{A count matrix of the occurence of motifs in the "plus" strand.} \item{both}{A count matrix of the occurence of motifs in the "minus" and "plus" strands together.}
}

\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# Read a file of sample sequence motifs.
motifs.dir = system.file("extdata", "AGRIS_PLACE_motifs_short_mb2.txt", package="TFbindR", mustWork = TRUE)
Motifs=read.delim(file = motifs.dir ,header=TRUE, sep="\t", as.is=T)

# Remove duplicated motifs.
Motifs=Motifs[which(!duplicated(Motifs$Motif)),]

# Read a sample fasta file. 
promoter.dir = system.file("extdata","TAIR10_upstream_1000_translation_start_20101028.fa",package="TFbindR")

# Specify the output directory.
output.dir="."

# Run the function to create the DB and save it into R file type .rda which is automatically compressed
makePromoterDBfromRE(Motifs = Motifs$Motif,promoter.fasta = promoter.dir, promoterDB=file.path(output.dir,"TAIR10_500bp_upstream.RE.db.rda"))
}