\name{robu}
\alias{robu}
\alias{CORR}
\alias{HIER}
\alias{USER}

\title{
Fitting Robust Variance Meta-Regression Models
}

\description{
  \code{robu} is used to meta-regression models using robust variance
  estimation (RVE) methods. \code{robu} can be used to estimate correlated and 
  hierarchical effects models using the original (Hedges, Tipton and Johnson, 
  2010) and small-sample corrected (Tipton, 2013) RVE methods. In addition, 
  \code{robu} contains options for fitting these models using user-specified 
  weighting schemes (see the Appendix of Tipton (2013) for a discussion of non-
  efficient weights in RVE).
}

\usage{
robu(formula, data, studynum, var.eff.size, userweights, modelweights = 
     c("CORR", "HIER"), rho = 0.8, small = TRUE,...)
}
\arguments{

  \item{formula}{ 
    An object of class \code{"formula"}. A typical meta-regression formula will 
    look similar to \code{y ~ x1 + x2...}, where \code{y} is a vector of effect 
    sizes and \code{x1 + x2...} are (optional) user-specified covariates. An 
    intercept only model can be specified with \code{y ~ 1} and the intercept 
    can be ommitted as follows \code{y ~ -1 +...}.
  }
  
  \item{data}{
    A data frame, list or environment or an object coercible by as.data.frame to 
    a data frame.
  }
  
  \item{studynum}{
    A vector of study numbers to be used in model fitting. \code{studynum} must 
    be a numeric or factor variable that uniquely identifies each study. 
  }
  
  \item{var.eff.size}{
    A vector of user-calculated effect-size variances.
  }
  
  \item{rho}{
    User-specified within-study effect-size correlation used to fit correlated 
    (\code{modelweights = "CORR"}) effects meta-regression models. The value of 
    \code{rho} must be between 0 and 1. The default value for \code{rho} is 0.8. 
    \code{rho} is not specified for hierarchical (\code{modelweights = "HIER"})
    effects models. 
  }
  
  \item{modelweights}{
    User-specified model weighting scheme.  The two two avialable options are 
    \code{modelweights = "CORR"} and \code{modelweights = "HIER"}.  The default 
    is \code{"CORR"}. See Hedges, Tipton and Johnson (2010) and Tipton (2013)
    for extended explanations of each weighting scheme.
  }
  
  \item{userweights}{
      A vector of user-specified weights if non-efficient weights are of 
      interest. Users interested in non-efficient weights should see the 
      Appendix of Tipton (2013) for a discussion of the role of non-efficient 
      weights in RVE).
  }
  
  \item{small}{
   \code{small = TRUE} is used to fit the meta-regression models with the small-
   sample corrections for both the residuals and degrees of freedom, as detailed 
   in Tipton (2013). Users wishing to use the original RVE estimator must 
   specify \code{small = FALSE} as the corrected estimator is the default 
   option.
   }
   \item{...}{
   Additional arguments to be passed to the fitting function.
   }
}  

\value{

  \item{output}{
    A data frame containing some combination of the robust coefficient names and 
    values, standard errors, t-test value, confidence intervals, degrees of 
    freedom and statistical significance.
  }

  \item{n}{The number of studies in the sample \code{n}}.

  \item{k}{The number of effect sizes in the sample \code{k}}.

  \item{k descriptives}{the minimum \code{min.k}, mean \code{mean.k}, median 
    \code{median .k}, and maximum \code{max.k} number of effect sizes per study.
  }

  \item{tau.sq.}{
    \code{tau.sq} is the between study variance component in the correlated 
    effects meta-regression model and the between-cluster variance component in 
    the hierarchical effects model. \code{tau.sq} is calculated using the 
    method-of-moments estimator provided in Hedges, Tipton, and Johnson (2010). 
    For the correlated effects model the method-of-moments estimar depends on 
    the user-specified value of rho.
  }
  
  \item{omega.sq.}{
    \code{omega.sq} is the between-studies-within-cluster variance component 
    for the hierarchical effects meta-regression model. \code{omega.sq} is 
    calculated using the method-of-moments estimator provided in Hedges, 
    Tipton, and Johnson (2010) erratum.
  }

  \item{I.2}{
    \code{I.2} is a test statistics used to quantify the amount of variability 
    in effect size estimates due to effect size heterogeneity as opposed to
    random variation.  
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

# Small-Sample Corrections - Hierarchical Dependence Model
HierModSm <-  robu(formula = effectsize ~ binge + followup + sreport
                   + age, data = hierdat, studynum = studyid, 
                   var.eff.size = var, modelweights = "HIER", small = TRUE)

print(HierModSm) # Output results
}

\keyword{robu}