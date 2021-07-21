
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
#' @param email     Mandatory E-mail adress for NCBI BLAST request, specifiedas
#'                  as a \code{character} string.
#' @param db        A \code{character} string specifying an NCBI genome or
#'                  reference sequence set database on which the submitted
#'                  sequences will be BLASTed (Default using GRCh37/hg19 genome
#'                  assembly database: db = "genomic/9606/GCF_000001405.25").
#' @return A \code{type} object returned description.
#' @author Yoann Pageaud.
#' @export
#' @examples
#' #Create example list of sequences you wish to BLAST using NCBI BLAST API.
#' ls.seq <- list(
#'   "7qtel" = "CCCTAACACTGTTAGGGTTATTATGTTGACTGTTCTCATTGCTGTCTTAG",
#'   "1ptel" = "GATCCTTGAAGCGCCCCCAAGGGCATCTTCTCAAAGTTGGATGTGTGCAT",
#'   "17qtel" = "CCCTAACCCTAAACCCTAGCCCTAGCCCTAGCCCTAGCCCTAGCCCTAGC")
#' #Submit the list of sequences to NCBI BLAST using the GRCh37/hg19 genome assembly database.
#' submit_NCBI_BLAST(
#'   seq.list = ls.seq, res.dir = "~/result_directory",
#'   email = "myemailadress@dkfz.de", db = "genomic/9606/GCF_000001405.25")
#' @references

submit_NCBI_BLAST <- function(
  seq.list, res.dir, delay.req = 10, email,
  db = "genomic/9606/GCF_000001405.25"){
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