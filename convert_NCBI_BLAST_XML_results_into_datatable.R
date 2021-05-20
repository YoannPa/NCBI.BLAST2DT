
#' Converts NCBI BLAST XML result file into an R data.table.
#' 
#' @param xml_file A \code{character} specifying the path to a XML file
#'                containing NCBI BLAST results. Such XML file can either be
#'                downloaded directly from the NCBI BLAST result page, or
#'                automatically generated using blastSeq() function from the
#'                R package hoardeR.
#' @return A \code{data.table} of the NCBI BLAST results.
#' @author Yoann Pageaud.
#' @export
#' @examples
#' NCBI_BLAST_XML2DT(xml_file = "path/to/my/ncbi_blast_result.xml")
#' @references
#' \itemize{
#'  \item{\href{https://www.ncbi.nlm.nih.gov/BLAST/tutorial/Altschul-1.html}{The Statistics of Sequence Similarity Scores.}}
#'  \item{\href{https://cran.r-project.org/web/packages/hoardeR/index.html}{hoardeR: Collect and Retrieve Annotation Data for Various Genomic Data Using Different Webservices.}}
#'  \item{\href{https://cran.r-project.org/web/packages/xml2/index.html}{xml2: Parse XML.}}
#' }

NCBI_BLAST_XML2DT <- function(xml_file){
  #Load NCBI BLAST XML results
  blast_xml <- xml2::read_xml(xml_file)
  #Get BLAST hits data
  hits_dt <- data.table::data.table(
    "hits_num" = as.integer(
      xml2::xml_text(xml2::xml_find_all(x = blast_xml, xpath = "//Hit_num"))),
    "hits_id" = xml2::xml_text(
      xml2::xml_find_all(x = blast_xml, xpath = "//Hit_id")),
    "hits_def" = xml2::xml_text(
      xml2::xml_find_all(x = blast_xml, xpath = "//Hit_def")),
    "hits_accession" = xml2::xml_text(
      xml2::xml_find_all(x = blast_xml, xpath = "//Hit_accession")),
    "hits_len" = as.integer(
      xml2::xml_text(xml2::xml_find_all(x = blast_xml, xpath = "//Hit_len"))))
  #Get BLAST hit samples data
  hsps_dt <- data.table::data.table(
    "subject_num" = as.integer(
      xml2::xml_text(xml2::xml_find_all(x = blast_xml, xpath = "//Hsp_num"))),
    "bit-score" = as.numeric(xml2::xml_text(
      xml2::xml_find_all(x = blast_xml, xpath = "//Hsp_bit-score"))),
    "score" = as.integer(
      xml2::xml_text(xml2::xml_find_all(x = blast_xml, xpath = "//Hsp_score"))),
    "E-value" = as.numeric(xml2::xml_text(
      xml2::xml_find_all(x = blast_xml, xpath = "//Hsp_evalue"))),
    "query_start" = as.integer(xml2::xml_text(
      xml2::xml_find_all(x = blast_xml, xpath = "//Hsp_query-from"))),
    "query_end" = as.integer(xml2::xml_text(
      xml2::xml_find_all(x = blast_xml, xpath = "//Hsp_query-to"))),
    "subject_start" = as.integer(xml2::xml_text(
      xml2::xml_find_all(x = blast_xml, xpath = "//Hsp_hit-from"))),
    "subject_end" = as.integer(xml2::xml_text(
      xml2::xml_find_all(x = blast_xml, xpath = "//Hsp_hit-to"))),
    "query_frame" = as.integer(xml2::xml_text(
      xml2::xml_find_all(x = blast_xml, xpath = "//Hsp_query-frame"))),
    "subject_frame" = as.integer(xml2::xml_text(
      xml2::xml_find_all(x = blast_xml, xpath = "//Hsp_hit-frame"))),
    "N_identities" = as.integer(xml2::xml_text(
      xml2::xml_find_all(x = blast_xml, xpath = "//Hsp_identity"))),
    "N_positives" = as.integer(xml2::xml_text(
      xml2::xml_find_all(x = blast_xml, xpath = "//Hsp_positive"))),
    "N_gaps" = as.integer(xml2::xml_text(
      xml2::xml_find_all(x = blast_xml, xpath = "//Hsp_gaps"))),
    "alignment_length" = as.integer(xml2::xml_text(
      xml2::xml_find_all(x = blast_xml, xpath = "//Hsp_align-len"))),
    "query_alignments" = xml2::xml_text(
      xml2::xml_find_all(x = blast_xml, xpath = "//Hsp_qseq")),
    "subject_alignments" = xml2::xml_text(
      xml2::xml_find_all(x = blast_xml, xpath = "//Hsp_hseq")),
    "midline_alignments" = xml2::xml_text(
      xml2::xml_find_all(x = blast_xml, xpath = "//Hsp_midline")))
  
  #Add empty 'hits_num' column when no alignment found in the XML file
  if(nrow(hsps_dt) == 0){ hsps_dt[, hits_num := integer()] } else {
    #Add 'hits_num'
    hsps_dt[, hits_num := cumsum(c(1, diff(as.integer(xml2::xml_text(
      xml2::xml_find_all(x = blast_xml, xpath = "//Hsp_num")))) != 1))]
  }
  #Create full BLAST data.table
  BLAST.dt <- data.table::merge.data.table(x = hits_dt, y = hsps_dt,
                                           by = "hits_num", all.y = TRUE)
  #Return result
  return(BLAST.dt)
}

#' Aggregates NCBI BLAST results from multiple XML files into an R data.table.
#' 
#' @param dir.to.xmls A \code{character} specifying a directory to be searched
#'                    recursively for NCBI BLAST XML results.
#' @param ncores      An \code{integer} specifying the number of cores or
#'                    threads to be used for parallel processing.
#' @return A \code{data.table} aggregating all NCBI BLAST results from all XML
#'         files found.
#' @author Yoann Pageaud.
#' @export
#' @examples
#' #Aggregate my NCBI BLAST results using 2 cores into an R data.table 
#' aggregate_NCBI_BLAST_XMLs2DT(dir.to.xmls = "my_dir/with/NCBI_BLAST_results/",
#'   ncores = 2)
#' @references
#' \itemize{
#'  \item{\href{https://www.ncbi.nlm.nih.gov/BLAST/tutorial/Altschul-1.html}{The Statistics of Sequence Similarity Scores.}}
#'  \item{\href{https://cran.r-project.org/web/packages/hoardeR/index.html}{hoardeR: Collect and Retrieve Annotation Data for Various Genomic Data Using Different Webservices.}}
#'  \item{\href{https://cran.r-project.org/web/packages/xml2/index.html}{xml2: Parse XML.}}
#' }

aggregate_NCBI_BLAST_XMLs2DT <- function(dir.to.xmls, ncores = 1){
  #List all XML files in a directory
  ls.xmls <- list.files(path = dir.to.xmls, include.dirs = FALSE,
                        recursive = TRUE, pattern = ".xml")
  ls.path <- file.path(dir.to.xmls, ls.xmls)
  names(ls.path) <- ls.xmls
  #Convert all XMLs found into data.tables in parallel
  ls.dt <- parallel::mclapply(
    X = ls.path, mc.cores = ncores, FUN = NCBI_BLAST_XML2DT)
  bind.dt <- data.table::rbindlist(l = ls.dt, idcol = "file_name")
  return(bind.dt)
}