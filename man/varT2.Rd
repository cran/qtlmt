\name{varT2}
\Rdversion{0.1}
\alias{varT2}
\title{
Group variables
}
\description{
Group variables via variable selection such that the grouped variables are optimal for multivariate analysis.
}
\usage{
varT2(x, grp, equalVar=T)
}
\arguments{
  \item{x}{
data matrix/frame. Columns are variables to select from.
}
  \item{grp}{
grouping indicator of observations.
}
  \item{equalVar}{
whether assume the same variance-covariance in two groups.
}
}

\value{
Hotelling's T^2 test statistic.
}

\examples{
data(etrait)
varT2(traits,mdat[,42])
}

