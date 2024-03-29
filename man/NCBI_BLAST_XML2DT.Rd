% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convert_NCBI_BLAST_XML_results_into_datatable.R
\name{NCBI_BLAST_XML2DT}
\alias{NCBI_BLAST_XML2DT}
\title{Converts NCBI BLAST XML result file into an R data.table.}
\usage{
NCBI_BLAST_XML2DT(xml_file)
}
\arguments{
\item{xml_file}{A \code{character} specifying the path to a XML file
containing NCBI BLAST results. Such XML file can either be:
\itemize{
 \item{downloaded directly from the NCBI BLAST result page.}
 \item{retrieved from the result directory you specified when using \link{submit_NCBI_BLAST}.}
 \item{generated using blastSeq() function from the R package hoardeR.}
}}
}
\value{
A \code{data.table} of the NCBI BLAST results.
}
\description{
Converts NCBI BLAST XML result file into an R data.table.
}
\examples{
NCBI_BLAST_XML2DT(xml_file = "path/to/my/ncbi_blast_result.xml")
}
\references{
\itemize{
 \item{\href{https://www.ncbi.nlm.nih.gov/BLAST/tutorial/Altschul-1.html}{The Statistics of Sequence Similarity Scores.}}
 \item{\href{https://cran.r-project.org/web/packages/hoardeR/index.html}{hoardeR: Collect and Retrieve Annotation Data for Various Genomic Data Using Different Webservices.}}
 \item{\href{https://cran.r-project.org/web/packages/xml2/index.html}{xml2: Parse XML.}}
}
}
\author{
Yoann Pageaud.
}
