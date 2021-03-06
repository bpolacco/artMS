% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MSstats_functions.R
\name{artmsSILACtoLong}
\alias{artmsSILACtoLong}
\title{Convert the SILAC evidence file to MSstats format}
\usage{
artmsSILACtoLong(evidence_file, output = NULL, verbose = TRUE)
}
\arguments{
\item{evidence_file}{(char) Text filepath to the evidence file}

\item{output}{(char) Text filepath of the output name. If NULL it does not
write the output}

\item{verbose}{(logical) \code{TRUE} (default) shows function messages}
}
\value{
(data.frame) with SILAC data processed for MSstats (and output file)
}
\description{
Converting the evidence file from a SILAC search to a format
compatible with MSstats. It basically modifies the Raw.files adding the
Heavy and Light label
}
\examples{
\dontrun{
evidence2silac <- artmsSILACtoLong(evidence_file = "silac.evicence.txt",
                                   output = "silac-evidence.txt")
}
}
\keyword{convert,}
\keyword{evidence}
\keyword{silac,}
