\name{sensitivity}
\alias{sensitivity.robu}
\alias{sensitivity}
\title{Sensitivity Analysis for Correlated Effects RVE}
\description{
  \code{sensitivity} is used to assess the impact of differing rho values on
  the correlated effects meta-regression model.
}
\usage{sensitivity(x)}
\arguments{
  \item{x}{
    A dataframe containing values of rho, tau squared, coefficient 
    estimates, and standard errors.
  }
}



\references{

  Hedges, L.V., Tipton, E., Johnson, M.C. (2010) Robust variance estimation in 
  meta-regression with dependent effect size estimates. \emph{Research Synthesis 
  Methods}. \bold{1}(1): 39--65. Erratum in \bold{1}(2): 164--165. DOI: 
  10.1002/jrsm.5

  Tipton, E. (in press) Small sample adjustments for robust variance estimation 
  with meta-regression. \emph{Psychological Methods}.
}

\examples{

# Correlated Effects Model
CorrMod   <-  robu(formula = effectsize ~ followup + males + binge + college, 
                   data = corrdat, studynum = studyid, var.eff.size = var, 
                   rho = .8, modelweights = "CORR", small = FALSE)

sensitivity(CorrMod) # Output sensitivity
}

\keyword{robu}