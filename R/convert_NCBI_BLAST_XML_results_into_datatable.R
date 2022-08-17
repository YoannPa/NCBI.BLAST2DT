
#' Converts NCBI BLAST XML result file into an R data.table.
#' 
#' @param xml_file A \code{character} specifying the path to a XML file
#'                containing NCBI BLAST results. Such XML file can either be:
#'                \itemize{
#'                 \item{downloaded directly from the NCBI BLAST result page.}
#'                 \item{retrieved from the result directory you specified when using \link{submit_NCBI_BLAST}.}
#'                 \item{generated using blastSeq() function from the R package hoardeR.}
#'                }
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
  #Get BLAST query data
  query_dt <- data.table::data.table(
    "query_id" = xml2::xml_text(
      xml2::xml_find_all(x = blast_xml, xpath = "//Iteration_query-ID")),
    "query-def" = xml2::xml_text(
      xml2::xml_find_all(x = blast_xml, xpath = "//Iteration_query-def")),
    "query-len" = as.integer(
      xml2::xml_text(xml2::xml_find_all(
        x = blast_xml, xpath = "//Iteration_query-len"))))
  #Get all hits and map them to all queries
  iter_hits <- xml2::xml_find_all(x = blast_xml, xpath = "//Iteration_hits")
  ls_query_hits <- lapply(X = seq_along(iter_hits), FUN = function(i){
    if(xml2::xml_text(xml2::read_xml(as.character(iter_hits[[i]]))) == "\n"){
      data.table::data.table("query_id" = query_dt[i]$query_id, "hits_num" = NA)
    } else {
      data.table::data.table(
        "query_id" = query_dt[i]$query_id,
        "hits_num" = as.integer(xml2::xml_text(xml2::xml_find_all(
          xml2::read_xml(as.character(iter_hits[[i]])), xpath = "//Hit_num"))))
    }
  })
  dt_query_hits <- data.table::rbindlist(l = ls_query_hits)
  rm(ls_query_hits)
  #Merge query_dt & dt_query_hits
  query_dt <- data.table::merge.data.table(
    x = query_dt, y = dt_query_hits, by = "query_id")
  
  #Get BLAST hits data
  ls_hits <- lapply(X = seq_along(iter_hits), FUN = function(i){
    iter_hits_xml <- xml2::read_xml(as.character(iter_hits[[i]]))
    if(xml2::xml_text(iter_hits_xml) == "\n"){
      data.table::data.table(
        "query_id" = query_dt[i]$query_id, "hits_num" = NA, "hits_id" = NA,
        "hits_def" = NA, "hits_accession" = NA, "hits_len" = NA)
    } else {
      data.table::data.table(
        "query_id" = query_dt[i]$query_id,
        "hits_num" = as.integer(xml2::xml_text(xml2::xml_find_all(
          x = xml2::read_xml(as.character(iter_hits[[i]])),
          xpath = "//Hit_num"))),
        "hits_id" = xml2::xml_text(xml2::xml_find_all(
          x = xml2::read_xml(as.character(iter_hits[[i]])),
          xpath = "//Hit_id")),
        "hits_def" = xml2::xml_text(xml2::xml_find_all(
          x = xml2::read_xml(as.character(iter_hits[[i]])),
          xpath = "//Hit_def")),
        "hits_accession" = xml2::xml_text(xml2::xml_find_all(
          x = xml2::read_xml(as.character(iter_hits[[i]])),
          xpath = "//Hit_accession")),
        "hits_len" = as.integer(xml2::xml_text(xml2::xml_find_all(
          x = xml2::read_xml(as.character(iter_hits[[i]])),
          xpath = "//Hit_len"))))
    }
  })
  hits_dt <- data.table::rbindlist(l = ls_hits)
  rm(ls_hits)
  # hits_dt <- data.table::data.table(
  #   "hits_num" = as.integer(
  #     xml2::xml_text(xml2::xml_find_all(x = blast_xml, xpath = "//Hit_num"))),
  #   "hits_id" = xml2::xml_text(
  #     xml2::xml_find_all(x = blast_xml, xpath = "//Hit_id")),
  #   "hits_def" = xml2::xml_text(
  #     xml2::xml_find_all(x = blast_xml, xpath = "//Hit_def")),
  #   "hits_accession" = xml2::xml_text(
  #     xml2::xml_find_all(x = blast_xml, xpath = "//Hit_accession")),
  #   "hits_len" = as.integer(
  #     xml2::xml_text(xml2::xml_find_all(x = blast_xml, xpath = "//Hit_len"))))
  
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
  #Merge hits and hsps
  dt_hits_hsps <- data.table::merge.data.table(x = hits_dt, y = hsps_dt,
                                               by = "hits_num", all.x = TRUE)
  #Create new query_hits key in query_dt & hits_dt
  query_dt[, query_hits := paste(query_id, hits_num, sep = "-")]
  dt_hits_hsps[, query_hits := paste(query_id, hits_num, sep = "-")]
  #Merge query_dt & hits_dt
  dt_query_hits_hsps <- data.table::merge.data.table(
    x = query_dt, y = dt_hits_hsps, by = "query_hits")
  #Create full BLAST data.table
  data.table::setnames(
    x = dt_query_hits_hsps, old = "query_id.x", new = "query_id")
  data.table::setnames(
    x = dt_query_hits_hsps, old = "hits_num.x", new = "hits_num")
  BLAST.dt <- dt_query_hits_hsps[, -c("query_id.y", "hits_num.y"), ]
  # BLAST.dt <- data.table::merge.data.table(x = hits_dt, y = hsps_dt,
  #                                          by = "hits_num", all.y = TRUE)
  
  # rm intermediate tables
  rm(
    dt_hits_hsps, dt_query_hits, dt_query_hits_hsps, hits_dt, hsps_dt,
    iter_hits, query_dt)
  #Return result
  return(BLAST.dt)
}

#' Aggregates NCBI BLAST results from multiple XML files into an R data.table.
#' 
#' @param dir.to.xmls A \code{character} specifying a directory to be searched
#'                    recursively for NCBI BLAST XML results (For more
#'                    information about XML results from NCBI BLAST see
#'                    \link{NCBI_BLAST_XML2DT}).
#' @param seq.names   A \code{character} vector to specify sequence names from
#'                    which you wish to load BLAST results in the R data.table.
#'                    Useful if your result folder contains results from other
#'                    BLAST submissions. If NULL, all results in 'dir.to.xmls'
#'                    will be loaded into the data.table.
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

aggregate_NCBI_BLAST_XMLs2DT <- function(dir.to.xmls, seq.names = NULL, ncores = 1){
  #List all XML files in a directory
  ls.xmls <- list.files(path = dir.to.xmls, include.dirs = FALSE,
                        recursive = TRUE, pattern = ".xml")
  if(!is.null(seq.names)){ #load only results for given sequence names 
    file.names <- paste0(seq.names, "/1.xml")
    ls.xmls <- ls.xmls[ls.xmls %in% file.names]
  }
  ls.path <- file.path(dir.to.xmls, ls.xmls)
  names(ls.path) <- ls.xmls
  #Convert all XMLs found into data.tables in parallel
  ls.dt <- parallel::mclapply(
    X = ls.path, mc.cores = ncores, FUN = NCBI_BLAST_XML2DT)
  bind.dt <- data.table::rbindlist(l = ls.dt, idcol = "file_name")
  return(bind.dt)
}