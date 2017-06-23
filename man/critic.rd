\name{critic}
\Rdversion{1.1}
\alias{critic}
\title{ Pierre - Goldwin correlogram }
\description{
An integro delayed correlogram to find critical periods for a biological phenomenon driven by a climatic factor
}
\usage{
    critic(t, Y, fac = NULL, dinf = 10, durinf = 2, dsup = 90, dursup = 90, 
    nperm = 0, nboot = 0, period = 365, 
    dt = 1, seriesName = "year", grtype = "image",
    roll = FALSE, alpha = 0.05, ps.print = FALSE)
}
\arguments{
  \item{t}{
vector : The climatic time series. In this version, must be annual and sampled dayly. Its length must be a multiple of 365. February 29 must be discarded.
}
  \item{Y}{
vector : the observations to regress. Must be of the same length as the number of years in t. One observation per year if fac==NULL (the default). If fac is not null, there may be several observations per year. See \code{fac} and \code{details}
}
  \item{fac}{
factor grouping the observations per year. Its levels number must be equal to the number of years in \code{t}
}
  \item{dinf}{
integer : the number of the day taken as first beginning period to scan in the year
}
  \item{durinf}{
numeric : the number of days taken as lower span of the periods to scan in the year
}
  \item{dsup}{
numeric : the number of the day taken as last beginning period to scan in the year
}
  \item{dursup}{
numeric : the number of days taken as largen span of the periods to scan in the year
}
  \item{nperm}{
numeric : number of random permutations
}
  \item{nboot}{
numeric : number of bootstrap subsamples
}
  \item{period}{
numeric : Number of time units per period. Default = 365 (days in a year)
}
  \item{dt}{
numeric : value of the time increment for integration. Default = 1
}
  \item{seriesName}{
string : name of the replicates of the time series. Default = "year"
}
  \item{grtype}{
type of map to draw. grtype may take the values "image","contour","filledcontour","persp". These codes call the correspondig R base functions.
}
  \item{roll}{
logical : only used if grtype=="persp" in what case the perspective plot rotates slowly to show all aspects of the perspective.
}
  \item{alpha}{
numeric: significance level for the tests. Default=0.05
}
  \item{ps.print}{
logical: Pseudovalues of the bootstrap must be printed (TRUE) or not (FALSE). Default = FALSE 
}
}
\details{
For each replication (by default: year) calculates the sums of the time series t, begining at a time i varying from dinf to dsup, and ending a time varying from i+durinf to i+dursup. Then correlates these sums to the vector Y of independent observations. The result is the map rho[i,j] giving  the correlation between Y and the corresponding sum of j elements (duration) after the time i. The significant level where the map can be cut is obtained by random permutations the number of which is defined by nperm. The confidence interval of the maximum correlation, as well as its bivariate confidence interval, are obtained by optional bootstrap. If nperm = 0 (default), no permutation is done. If nboot = 0, no bootstrap is done.    
}
\value{
z : a matrix containing the correlation coefficients of Y with the sum of j days  
}
\author{
Jean-Sebastien Pierre
}
\references{Pierre, J. S., Guillome, M. and Querrien, M. T. 1986. A Statistical and Graphic Method for Seeking in Which Periods of the Year Are the Animal Populations Peculiarly Sensitive to a Given Weather Component (Critical Periods of Time) - Application to the Case of Cereal Aphids (in french) - Acta Oecologica-Oecologia Generalis 7: 365-380}

\seealso{
\code{\link{image}},  \code{\link{contour}},  \code{\link{filled.contour}},  \code{\link{persp}} for graphical representations of the correlogram. 
}
\examples{
data(sit,time3)
critic(t=time3,Y=sit,dinf=50,dsup=90,durinf=20,dursup=50)
}
\keyword{ ts }
