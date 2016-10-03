\name{group.center}
\alias{group.center}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Convenience function for calculating group-centered covariates.}
\description{Creates a within-study (or within-cluster) version of the 
covariate in question.}
\usage{group.center(var, grp)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{var}{The covariate to be group centered.}
  \item{grp}{A vector corresponding to the group identification.}
}
\value{A column or vector containing the group.centered covariate.}


\examples{

# Load data
data(corrdat) 

# Create a group centered covariate 
males_c <- group.center(corrdat$males, corrdat$studyid)

}
\keyword{robumeta}

