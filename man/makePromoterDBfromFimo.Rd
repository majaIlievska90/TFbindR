\name{makePromoterDBfromFimo}
\alias{makePromoterDBfromFimo}
\title{
Count the Occurence of Sequence Motifs from Fimo File.
}
\description{
The function reads the input sequences from a FASTA file and identifies the genes. Next, it splits the strands from the FIMO file into "+" and "-", and counts the occurence of gene and motif combination for each strand separately. The counts are reported in a matrix. 
}
\usage{
makePromoterDBfromFimo(fimo, promoter.fasta, promoterDB = "promoter.DB")
}
\arguments{
  \item{fimo}{
A fimo file containing the motifs of interest.  
}
  \item{promoter.fasta}{
Fasta file of the promoter sequences.
}
  \item{promoterDB}{
The path, and the name of the object where results will be saved.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list named Strands with 3 elements, \item{plus}{A count matrix of the occurence of motifs in the "plus" strand.} \item{minus}{A count matrix of the occurence of motifs in the "minus" strand.} \item{both}{A count matrix of the occurence of motifs in the "minus" and "plus" strands together.}
}
\references{
Grant, Charles E., Timothy L. Bailey, and William Stafford Noble. "FIMO: scanning for occurrences of a given motif." Bioinformatics 27.7 (2011): 1017-1018.
}
\author{
%%  Maja Ilievska
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# Read a sample fimo file. The data is downloaded from the Arabidopsis Cistrome database.
fimo.dir = system.file("extdata", "cistrome.fimo.txt", package="TFbindR", mustWork = TRUE)

# Read a sample fasta file. 
promoter.dir = system.file("extdata","TAIR10_upstream_1000_translation_start_20101028.fa",package="TFbindR")

# Specify the output directory.
output.dir="."

# Run the function to create the DB and save it into R file type .rda which is automatically compressed. 
makePromoterDBfromFimo(fimo = fimo.dir, promoter.fasta=promoter.dir, promoterDB=file.path(output.dir,"TAIR10_500bp_upstream.fimo.db.rda"))
}