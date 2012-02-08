\name{motif-class}
\Rdversion{1.1}
\docType{class}
\alias{motif}
\alias{motif-class}
\alias{nMotifs}



\title{Class "motif" }
\description{This object contains contains PWM, motif consensus, motif length and all aligned sequences for a specific motif}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("motif_gadem", ...)}.
}
\section{Slots}{
	 \describe{
    \item{pwm}{:PWM results. }
	\item{consensus}{:Sequences consensus. }
	\item{alignList}{:List of sequences alignment. }
	\item{name}{:Name of sequences. }
  }
}

\author{Arnaud Droit \email{arnaud.droit@crchuq.ulaval.ca}}

\seealso{
	 \code{\linkS4class{gadem}} , \code{\linkS4class{align}}, \code{\linkS4class{parameters}}
}
\examples{
showClass("gadem")
}
\keyword{classes}