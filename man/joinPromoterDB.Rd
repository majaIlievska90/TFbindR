\name{joinPromoterDB}
\alias{joinPromoterDB}
\title{
Join promoterDB objects.

}
\description{
Combing together promoter.DB objects into one promoter.DB.  
}
\usage{
joinPromoterDB(DB1, DB2, out.DB)
}
\arguments{
  \item{DB1}{
prmoterDB object}
  \item{DB2}{
promoterDB object}
  \item{out.DB}{
the output file for the combined promoterDB object.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A promoterDB object - list with 3 elements, \code{plus} - a count matrix of the occurence of motifs in the "plus" strand, \code{minus} - a count matrix of the occurence of motifs in the "minus" strand and \code{both} - a count matrix of the occurence of motifs in the "minus" and "plus" strands together.}
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
}
