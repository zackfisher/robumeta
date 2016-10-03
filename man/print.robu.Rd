\name{print.robu}
\alias{print.robu}
\title{Outputs Model Information}
\description{Prints relevant information from robu function. }

\usage{
  \method{print}{robu}(x, digits = 3, ...)
}

\arguments{

  \item{x}{
    Object from robu class.
  }
  
  \item{digits}{
  Controls the number of digits to print when printing numeric values.
  }
  \item{...}{
    Additional arguments to be passed to the fitting function.
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

# Load data
data(hierdat)

### Small-Sample Corrections - Hierarchical Dependence Model
HierMod  <-  robu(formula = effectsize ~ binge + followup + sreport
                   + age, data = hierdat, studynum = studyid, 
                   var.eff.size = var, modelweights = "HIER", small = FALSE)

print(HierMod) # Output results
}

\keyword{robu}
