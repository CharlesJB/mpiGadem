 
 
\name{readPWMfile}
\alias{readPWMfile}
\title{Read Transfac File} 

\description{
This function is use to read standard Transfac type file.} 
\usage{
  readPWMfile(file)
}    
\arguments{
	\item{file}{Transfac file's name.}
	}
\details{This function is designed to read standard Transfac type file. For more information about the format, please refere to \url{http://mcast.sdsc.edu/doc/transfac-format.html}}
\author{Arnaud Droit <\email{arnaud.droit@ircm.qc.ca}>}
\value{A list of matrix.}
\examples{
#####Database and Scores#####
path <- system.file("extdata","jaspar2009.txt",package="rGADEM")
jaspar <- readPWMfile(path)

}
\keyword{misc}
 
