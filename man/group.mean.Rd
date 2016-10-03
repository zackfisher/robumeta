\name{group.mean}
\alias{group.mean}
\title{Convenience function for calculating group-mean covariates.}
\description{Creates a between-study (or between-cluster) version of the 
covariate in question.}
\usage{group.mean(var, grp)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{var}{The covariate cotaining the values to be group averaged.}
  \item{grp}{The group from which the average should be calculated.}
}
\value{A column or vector containing the group.mean covariate.}

\examples{

# Load data
data(corrdat)

# Create a group mean covariate 
age_m <- group.mean(corrdat$age, corrdat$studynum)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{robumeta}

