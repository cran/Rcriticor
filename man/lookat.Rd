\name{lookat}
\alias{lookat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Looking at the correlogram}
\description{
A function to inspect a point of the correlogram, to plot a scatterplot of dependent variable against the sum of the independent variable, and return the corresponding regression.   
}
\usage{
lookat(u, Y)
}
\arguments{
  \item{u}{
An object of class \code{criticor}
}
  \item{Y}{
numeric vector : the dependent variable
}
}
\details{
This functions requires the display of a correlogram plot on the current device, and that it has been previously saved as an object of class criticor. It then calls the function \code{locator}. The user has then to click on the correlogram at the point he wishes to inspect. A plot is then displayed on a new device, showing the scatterplot of the summed independent variable for the beginning of the series (Y axis on the correlogram) and the duration of the sum (X axis on the correlogram). Several point may be tried, only the results for the last one are returned by the function (see Value). The process stops when a left click outside the correlogram.
}
\value{
An object of class lm
}
\references{
Pierre, J. S., Guillome, M. and Querrien, M. T. 1986. A Statistical and Graphic Method for Seeking in Which Periods of the Year Are 
the Animal Populations Peculiarly Sensitive to a Given Weather Component (Critical Periods of Time) - Application to the Case of Cereal Aphids. 
- Acta Oecologica-Oecologia Generalis 7: 365-380. (in french, english summary)
}
\author{
Jean-Sebastien Pierre ;  jean-sebastien.pierre@univ-rennes1.fr
}
\section{Warning}{

	1. The scatterplot is currently recovered by the correlogram, whose device is reactivated after drawing \\
           the scatterplot. You have to move the correlogram apart with the mouse to see the scatterplot.\\ 
        2. lookat works only on correlograms drawn with the options image or contour The option filledcontour \\
           (calling filled.contour) gives erroneous results because its particular scale which is not accessible.\\
           It does not work at all with the option persp.
}
\seealso{
\code{\link{critic}},\code{\link{locator}},\code{\link{image}},\code{\link{contour}},\code{\link{filled.contour}},\code{\link{persp}} for graphical representations. 
}
\examples{
data(sit,time3)
U=critic(t=time3,Y=sit,dinf=1,dsup=90,durinf=2,dursup=100)
\dontrun{lookat(U,sit)}
}
\keyword{ ts }
