\name{sureStep}
\Rdversion{0.1}
\alias{sureAdd1}
\alias{sureDrop1}
\alias{sureStep}
\alias{sureStep2}
\title{
SURE model selection
}
\description{
Select a SURE model via model selection.
}
\usage{
sureAdd1(obj, y, x, range=NULL, iter=250, tol=1e-8, ext=FALSE)

sureDrop1(obj, y, x, range=NULL, iter=250, tol=1e-8, ext=FALSE)

sureStep(y, x, v, lower, upper, k, direction=c("both","backward",
   "forward"), iter=250, max.terms=200, steps=1000, tol=1e-8)

sureStep2(obj, y, x, cv, direction=c("both","backward","forward"),
   range=NULL, iter=250, steps=1000, tol=1e-8, ext=FALSE)
}
\arguments{
  \item{obj}{
initial model in model search; can be an object of \code{\link{sureEst}}.
}
  \item{y}{
n by p matrix, whose columns are dependent variables.
}
  \item{x}{
n by m matrix, whose columns are predictor variables to select from.
}
  \item{cv}{
critical value used in the likelihood ratio test to determine adding/dropping a term.
}
  \item{v}{
list; v[[j]] indicates which x's to start with for y[,j].
}
  \item{range}{
list; range[[j]] indicates which x's (all if \code{NULL}) correspond to which y[,j].
}
  \item{lower}{
list in the format of \code{v}; lower scope of the model.
}
  \item{upper}{
list in the format of \code{v}; upper scope of the model.
}
  \item{k}{
penalty, 0 if missing or <0.
}
  \item{direction}{
forward selection, backward elimination or stepwise.
}
  \item{iter}{
maximum number of iterations in a numerical process to estimate model parameters.
}
  \item{max.terms}{
maximum number of terms in the final model.
}
  \item{steps}{
maximum number of search steps.
}
  \item{tol}{
convergence tolerance.
}
  \item{ext}{
whether to perform an extensive search.
}
}
\value{
a list with the following components:
  \item{loglik:}{log-likelihood of the model}
  \item{b:}{estimates of model coefficients}
  \item{sigma:}{estimates of residual variance-covariance}
  \item{v:}{variables selected in the model}
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
g1<- sureStep(y, x, v=v, upper=upper, k=19)
g2<- sureStep2(g1, y, x, cv=1/2, ext=FALSE)
}
}

