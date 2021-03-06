\name{sit2}
\alias{sit2}
\docType{data}
\title{
Sitobion with replicates
}
\description{
Peak densities of the cereal aphid Sitobion avenae, from 1974 to 1980, three replicates per year. Intendid to be used with
critic, together with fac
}
\usage{data(sit2)}
\format{
  The format is:
 num [1:21] 0.34 -0.175 0.147 2.977 1.477 ...
}
\source{
Pierre, J. S., Guillome, M. and Querrien, M. T. 1986. A Statistical and Graphic Method for Seeking in Which Periods of the Year Are 
the Animal Populations Peculiarly Sensitive to a Given Weather Component (Critical Periods of Time) - Application to the Case of Cereal Aphids. 
- Acta Oecologica-Oecologia Generalis 7: 365-380. (in french, english summary)
}
\examples{
data(sit2,fac)
plot(fac,sit2)
## maybe str(sit2) ; plot(sit2) ...
}
\keyword{datasets}
