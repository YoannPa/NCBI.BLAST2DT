
#' Submits a list of DNA sequences to NCBI BLAST for alignment to a sequence
#' database.
#' 
#' @param seq.list  A named \code{list} of DNA sequences stored as character
#'                  strings to be BLASTed on NCBI servers databases.
#' @param res.dir   A \code{character} string specifying the directory where
#'                  NCBI BLAST submission results will be stored. Each BLAST
#'                  result is stored in a separate folder within the given
#'                  directory. Blast hits for each submission are available as a
#'                  XML file in each matching folder.
#' @param delay.req Mandatory seconds between BLAST request submissions as an
#'                  \code{integer} (Default: delay.req = 10).
#' @param email     Mandatory E-mail adress for NCBI BLAST request, specified as
#'                  a \code{character} string.
#' @param db        A \code{character} string specifying an NCBI genome or
#'                  reference sequence set database on which the submitted
#'                  sequences will be BLASTed (Default using GRCh38 (hg38)
#'                  genome assembly database:
#'                  db = "genomic/9606/GCF_000001405.39").
#'                  For more supported databases see Details.
#' @details
#' Other supported NCBI genome databases include:
#' \itemize{
#'  \item{Human GRCh37 (hg19): \code{db = "genomic/9606/GCF_000001405.25"}}
#'  \item{Mus musculus MGSCv37 (mm9): \code{db = "genomic/10090/GCF_000001635.18"}}
#'  \item{Mus musculus GRCm38.p6 (mm10): \code{db = "genomic/10090/GCF_000001635.26"}}
#'  \item{Mus musculus GRCm39 (mm39): \code{db = "genomic/10090/GCF_000001635.27"}}
#'  \item{Drosophila melanogaster (fruit fly) Release 6 plus ISO1 MT: \code{db = "genomic/7227/GCF_000001215.4"}}
#'  \item{Danio rerio (zebrafish) GRCz11: \code{db = "genomic/7955/GCF_000002035.6"}}
#' }
#' You can try other db strings based on the following database name
#' structure:\cr
#' \code{db = "genomic/{taxonomy ID}/{RefSeq GCF assembly accession ID}"}\cr\cr
#' Once you start the submission, some logs will appear in the console.\cr
#' The unique submission ID will be displayed on one line such as 
#' "Run G01J99FG013 : 00:00:01".\cr If you wish to track your submission
#' directly from the NCBI BLAST web interface, you can go
#' \href{https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=GetSaved&RECENT_RESULTS=on}{here} 
#' and paste your submission ID (in this case: G01J99FG013) in the "Request ID" 
#' field.\cr If your submission is done processing, it will give you access to
#' all results, and other file formats to save them.
#' @return A \code{type} object returned description.
#' @author Yoann Pageaud.
#' @export
#' @examples
#' #Create example list of sequences you wish to BLAST using NCBI BLAST API.
#' ls.seq <- list(
#'   "7qtel" = "CCCTAACACTGTTAGGGTTATTATGTTGACTGTTCTCATTGCTGTCTTAG",
#'   "1ptel" = "GATCCTTGAAGCGCCCCCAAGGGCATCTTCTCAAAGTTGGATGTGTGCAT",
#'   "17qtel" = "CCCTAACCCTAAACCCTAGCCCTAGCCCTAGCCCTAGCCCTAGCCCTAGC")
#' #Submit the list of sequences to NCBI BLAST using the GRCh38/hg38 genome assembly database.
#' submit_NCBI_BLAST(
#'   seq.list = ls.seq, res.dir = "~/result_directory",
#'   email = "myemailadress@dkfz.de", db = "genomic/9606/GCF_000001405.39")
#' @references

submit_NCBI_BLAST <- function(
  seq.list, res.dir, delay.req = 10, email,
  db = "genomic/9606/GCF_000001405.39"){
  if(!file.exists(res.dir)){ dir.create(res.dir) }
  res <- lapply(X = seq_along(seq.list), FUN = function(i){
    cat(paste0("BLASTing ", names(seq.list)[i],":\n"))
    dist.dir <- file.path(res.dir, names(seq.list)[i])
    if(!file.exists(dist.dir)){
      blast_request <- hoardeR::blastSeq(
        seq = seq.list[[i]], delay_req = delay.req, email = email,
        database = db, logFolder = dist.dir, xmlFolder = dist.dir)
    } else { blast_request <- "Already processed" }
    cat(paste0("Done.\n"))
    return(blast_request)
  })
}