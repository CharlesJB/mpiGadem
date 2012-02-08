\name{gadem-class}
\Rdversion{1.1}
\docType{class}
\alias{gadem}
\alias{gadem-class}
\alias{[,gadem-method}
\alias{[[,gadem-method}
\alias{dim,gadem-method}
\alias{consensus,gadem-method}
\alias{startPos,gadem-method}
\alias{endPos,gadem-method}
\alias{parameters,gadem-method}
\alias{consensus}
\alias{nOccurrences}
\alias{nOccurrences,gadem-method}
\alias{getPWM}
\alias{getPWM,gadem-method}
\alias{length,gadem-method}
\alias{startPos}
\alias{endPos}
\alias{summary,list-method}
\alias{nMotifs,gadem-method}
\alias{names,gadem-method}
\alias{names<-,gadem-method}
\alias{plot}
\alias{plot,gadem-method}
\alias{plot,motif-method}
\alias{getPWM,motif-method}

\title{Class "gadem" }
\description{This object contains all gadem output information.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("gadem", ...)}.
}
\section{Slots}{
	 \describe{
    \item{motifList}{List of input PWM. }
	\item{parameters}{List of rGADEM parameters. }
  }
}

\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "gadem")}: subset gadem object.}
	\item{[[]}{\code{signature(x = "gadem")}: subset gadem object.}
	\item{nMotifs}{\code{signature(x = "gadem")}: Number of motifs identified}
	\item{names}{\code{signature(x = "gadem")}: Assign motifs names.}
	\item{dim}{\code{signature(x = "gadem")}: Number of sequences identified for each motifs.}
	\item{consensus}{\code{signature(x = "gadem")}:Sequence of consensus motifs.}
	\item{nOccurrences}{\code{signature(x = "gadem")}:View of PWMs.}
	\item{startPos}{\code{signature(x = "gadem")}:Start position for each sequences.}
	\item{endPos}{\code{signature(x = "gadem")}:End position for each sequences.}
	\item{getPWM}{\code{signature(x = "gadem")}:End position for each sequences.}
	 }
}


\author{Arnaud Droit \email{arnaud.droit@crchuq.ulaval.ca}}

\seealso{
	 \code{\linkS4class{motif}} , \code{\linkS4class{align}}, \code{\linkS4class{parameters}}
}
\examples{
showClass("gadem")
}
\keyword{classes}
