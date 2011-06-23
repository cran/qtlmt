\name{sureEst}
\Rdversion{0.1}
\alias{sureEst}
\title{
SURE model parameter estimation
}
\description{
Estimate parameters in a SURE model.
}
\usage{
sureEst(y,x,v,sigma,iter=250,tol=1e-8)
}
\arguments{
  \item{y}{
n by p matrix, whose columns are dependent variables.
}
  \item{x}{
n by m matrix, whose columns are predictor variables to select from.
}
  \item{v}{
list; v[[j]] indicates which x's in the model for y[,j].
}
  \item{sigma}{
residual variance-covarance matrix (if given).
}
  \item{iter}{
maximum number of iterations in a numerical process to estimate model parameters.
}
  \item{tol}{
convergence tolerance.
}
}

\value{
a list with the following components:
  \item{loglik:}{log-likelihood of the model}
  \item{b:}{estimates of model coefficients}
  \item{sigma:}{estimated residual variance-covariance}
}

\seealso{
\code{\link{sureStep}}
}
\examples{
data(etrait)
x<- as.matrix(3/2-mdat)
y<- as.matrix(traits)[,1:3]
v<- list()
upper<- list()
for(k in 1:ncol(y)){
   v[[k]]<- numeric(0)
   upper[[k]]<- 1:ncol(x)
}
\dontrun{
g<- sureStep(y, x, v=v, upper=upper, k=19, direction="both",
   iter=250, max.terms=250, steps=2000, tol=1e-6)
ge<- sureEst(y, x, v=g$v, iter=250, tol=1e-8)
}
}

