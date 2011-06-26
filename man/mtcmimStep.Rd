\name{mtcmimStep}
\Rdversion{0.1}
\alias{mtcmimAdd1}
\alias{mtcmimDrop1}
\alias{mtcmimStep}
\alias{mtcmimStep2}
\title{
MTCMIM model selection
}
\description{
Model selection for multiple-trait composite multiple-interval mapping.
}
\usage{
mtcmimAdd1(obj, y, x, xid, mpos, mdat, pp=1, len=1, type=1,
   iter=2500, tol=1e-8, ext=FALSE)

mtcmimDrop1(obj, y, x, xid, mpos, mdat, pp=1, len=1, type=1,
   iter=2500, tol=1e-8, ext=FALSE)

mtcmimStep(obj, y, x, xid, mpos, mdat, cv=0,
   direction=c("both", "backward", "forward"), pp=1, len=1,
   type=1, iter=2500, tol=1e-8, ext=FALSE)
}
\arguments{
  \item{obj}{
object of class \code{\link{mtcmim}}.
}
  \item{y}{
n by p matrix, whose columns are dependent variables.
}
  \item{x}{
covariates; n by m numerical matrix.
}
  \item{xid}{
list of length p; xid[[j]] specifies columns of x as covariates for y[,j] .
}
  \item{mpos}{
data frame (id=marker index, ch=chromosome id, m=marker index on the chromosome, dist=genetic position in cM on the chromosome).  Chromosome id should be an integer.
}
  \item{mdat}{
matrix of n rows; marker genotypes (1 or 0). columns should correspond to markers in the order.
}
  \item{pp}{
mapping population: BC-1, RIL-selfing-2, RIL-brother-sister-mating-3.
}
  \item{len}{
step length in search.
}
  \item{type}{
1 if traits can have the different sets of covariates and QTL, 2 if all have the same set of covariates and QTL.
}
  \item{ext}{
whether to perform an extensive search.
}
  \item{cv}{
critical value used in the likelihood ratio test to determine adding/dropping a QTL.
}
  \item{direction}{
forward selection, backward elimination or both directions.
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
  \item{loglik:}{log-likelihood of the final model}
  \item{a:}{covariate effects}
  \item{b:}{QTL effects}
  \item{sigma:}{residual variance-covariance}
  \item{qtl:}{QTL for each trait}
  \item{dists:}{QTL locations}
}

\note{
Currently, not able to include epistatic effects.
}

\seealso{
\code{\link{mtcmim}}
}
\examples{
data(etrait)
y<- traits[,1:5]
qtl<- vector("list",5); qtl[[1]]<- c(1,2)
eps<- data.frame(y=1,q1=1,q2=2)
dists<- dists[c(4,11),]
x<- mdat - 3/2
\dontrun{
o<- mtcmim(y,mpos,mdat,dists=dists,
   qtl=qtl,eps=eps,win=5,range=-1,pp=2,len=1)
of<- mtcmimAdd1(o,y=y,mpos=mpos,mdat=mdat,pp=2,len=3)
os<- mtcmimStep(of,y=y,mpos=mpos,mdat=mdat,cv=25,
   direction="both",pp=2,len=3)
}
}

