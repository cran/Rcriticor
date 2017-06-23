\name{time3}
\Rdversion{1.1}
\alias{time3}
\docType{data}
\title{
mean daily temperature data in Rennes (France)
}
\description{
mean daily temperature data in Rennes (France) from 1975 to 1981. february 29 of bissextile years discarded
}
\usage{data(time3)}
\format{
  The format is:
 num [1:2555] 1 2.3 1.7 6.1 1.9 5.2 8.9 10.5 8.4 6.3 ...
}
\source{
PIERRE J.S., GUILLOME M., QUERRIEN M.T., 1986. Une methode statistique et graphique de recherche des periodes de l'annee ou les populations animales sont particulierement sensibles a une composante donnee du climat (periodes critiques). Application au cas des pucerons des cereales. Acta Oecologica, Oecol. Gener., 7, 365-380.
}
\examples{
data(time3)
summary(time(3))
}

\keyword{datasets}
