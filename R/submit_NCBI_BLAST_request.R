
#' Prepares a set of GenBank Accession IDs for NCBI BLAST submission
#' 
#' @param GBaccess.bed A \code{data.table} or a \code{data.frame} in a BED-like
#'                     format containing columns as following:
#'                     \itemize{
#'                      \item{column 1: The Genbank accession IDs}
#'                      \item{column 2: The start position of the sequence to
#'                      keep}
#'                      \item{column 3: The end position of sequence to keep}
#'                     }
#' @param ncores       An \code{integer} specifying the number of cores or
#'                     threads to be used for parallel processing.
#' @return A \code{list} of sequences based on information provided in GBaccess.
#' @author Yoann Pageaud.
#' @export
#' @examples
#' #Create example 1 row data.frame
#' df <- data.frame("GBaccessID" = "AC073318", "Start" = 71401, "End" = 120576)
#' #Prepare sequence for submission to NCBI BLAST
#' seq.list <- prepare.gb.access(GBaccess.bed = df)
#' # seq.list can now be used by submit_NCBI_BLAST
#' @references
#' \href{https://academic.oup.com/bioinformatics/article/35/3/526/5055127}{Paradis E. & Schliep K. 2019. ape 5.0: an environment for modern phylogenetics and evolutionary analyses in R. Bioinformatics 35: 526–528. doi:10.1093/bioinformatics/bty633. HAL: ird-01920132.}

prepare.gb.access <- function(GBaccess.bed, ncores = 1){
  #Get Genbank accession IDs
  if(data.table::is.data.table(GBaccess.bed)){
    res <- ape::read.GenBank(access.nb = GBaccess.bed[[1]], as.character = TRUE)
  } else if(is.data.frame(GBaccess.bed)){
    res <- ape::read.GenBank(access.nb = GBaccess.bed[,1], as.character = TRUE)
  }
  res <- parallel::mclapply(
    X = res, mc.cores = ncores, FUN = paste, collapse = "")
  #Extract positions in sequences
  seq.list <- parallel::mclapply(
    X = seq_along(res), mc.cores = ncores, FUN = function(i){ substr(
      x = res[[i]], start = GBaccess.bed[i, 2], stop = GBaccess.bed[i, 3])})
  #Create sequence names
  if(data.table::is.data.table(GBaccess.bed)){
    names(seq.list) <- paste(names(res), paste(
      GBaccess.bed[[2]], GBaccess.bed[[3]], sep = "-"), sep = ":")
  } else if(is.data.frame(GBaccess.bed)){
    names(seq.list) <- paste(names(res), paste(
      GBaccess.bed[, 2], GBaccess.bed[, 3], sep = "-"), sep = ":")
  }
  return(seq.list)
}

#' Splits a data.frame of Genbank IDs and intervals into a data.table of smaller
#' intervals.
#' 
#' @param df A \code{data.frame} or a \code{data.table} in a BED-like format
#'           following the same specifications as 'GBaccess.bed' in
#'           \link{prepare.gb.access}.
#' @param by An \code{integer} to specify the length of the output DNA sequences
#'           you wish your queries to be splitted into, in base pair
#'           (Default: by = 7000).
#' @return A \code{data.table} of the Genbank IDs with intervals of the
#'         specified length. You can then subset the intervals of interest
#'         before submitting the data.table to \link{prepare.gb.access} or
#'         \link{get.NCBI.BLAST2DT}.
#' @author Yoann Pageaud.
#' @export
#' @examples
#' #Using an example data.frame of 1 Genbank ID
#' split.queries.df(
#'   df = data.frame("GB.access" = "AC073318", "Start" = 71401, "End" = 120576),
#'   by = 7025)

split.queries.df <- function(df, by = 7000){
  ls.dt <- lapply(X = seq(nrow(df)), FUN = function(i){
    #Compute breaks
    brks <- seq(from = df[i, 2], to = df[i, 3], by = by)
    brks <- brks[2:(length(brks)-1)]
    #Create intervals based on breaks
    dt <- data.table::data.table(
      c(df[i, 2]:df[i, 3]), findInterval(c(df[i, 2]:df[i, 3]), brks))
    #Create data.table of coordinates based on the intervals created
    do.call(rbind, by.default(dt, dt$V2, function(g){ data.table::data.table(
      GB.access = df[i, 1], Start = min(g$V1), End = max(g$V1))}))
  })
  #Rbind all data.frames into a data.table
  dt.res <- data.table::rbindlist(l = ls.dt)
  return(dt.res)
}


#' Splits a set of BLAST queries into a list of smaller DNA sequences.
#' 
#' @param x      Supports 2 possible formats:
#'               \itemize{
#'                \item{A \code{data.frame} in a BED-like format following the
#'                same specifications as 'GBaccess.bed' in
#'                \link{prepare.gb.access}.}
#'                \item{A \code{list} of DNA sequences.}
#'               }
#' @param by     An \code{integer} to specify the length of the output DNA
#'               sequences you wish your queries to be splitted into, in base
#'               pair (Default: by = 7000).
#' @param ncores An \code{integer} specifying the number of cores or threads to
#'               be used for parallel processing.
#' @return A \code{list} of DNA sequences of the requested length.
#' @author Yoann Pageaud.
#' @export
#' @examples
#' #Using an example data.frame of 1 Genbank ID
#' ls.seq <- split.queries(
#'   x = data.frame("GB.access" = "AC073318", "Start" = 71401, "End" = 120576),
#'   by = 7025)
#' #Using an example list of DNA sequences
#' ls.seq <- list(
#'   "7qtel" = "CCCTAACACTGTTAGGGTTATTATGTTGACTGTTCTCATTGCTGTCTTAG",
#'   "1ptel" = "GATCCTTGAAGCGCCCCCAAGGGCATCTTCTCAAAGTTGGATGTGTGCAT",
#'   "17qtel" = "CCCTAACCCTAAACCCTAGCCCTAGCCCTAGCCCTAGCCCTAGCCCTAGC")
#' ls.seq <- split.queries(x = ls.seq, by = 10)

split.queries <- function(x = x, by = 7000, ncores = 1){
  if(is.data.frame(x)){
    dt.res <- split.queries.df(df = x, by = by)
    #Convert to list of sequences
    ls.seq <- prepare.gb.access(GBaccess.bed = dt.res, ncores = ncores)
  } else if(is.list(x)){
    ls.seq <- lapply(X = x, FUN = function(i){
      brks <- seq(from = 1, to = nchar(i), by = by)
      brks <- brks[2:(length(brks)-1)]
      brks
      dt <- data.table::data.table(
        c(1:nchar(i)), findInterval(c(1:nchar(i)), brks))
      by.seq <- by.default(dt, dt$V2, function(g){ 
        substr(x = i, start = min(g$V1), stop = max(g$V1))
      })
      by.names <- by.default(dt, dt$V2, function(g){ 
        paste(min(g$V1), max(g$V1), sep = "-")
      })
      sub.ls <- as.list(by.seq)
      names(sub.ls) <- c(by.names)
      sub.ls
    })
    #Flatten output list
    ls.seq <- lapply(X = unlist(ls.seq), FUN = function(i){ i })
    #Create new names
    names(ls.seq) <- gsub(pattern = "\\.(\\d+\\-\\d+$)", replacement = ":\\1",
                          x = names(ls.seq))
  }
  return(ls.seq)
}


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
#' \itemize{
#'  \item{\href{https://cran.r-project.org/web/packages/hoardeR/index.html}{hoardeR: Collect and Retrieve Annotation Data for Various Genomic Data Using Different Webservices.}}
#'  \item{\href{https://academic.oup.com/nar/article/36/suppl_2/W5/2505810}{Johnson, M. et al. NCBI BLAST: a better web interface. Nucleic Acids Research 36, W5–W9 (2008).}}
#' }

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

#' Automatically submits all sequences and retrieves all results in a
#' data.table.
#' 
#' @param sequences   Either a \code{list} of DNA sequences, or a
#'                    \code{data.frame} in a BED-like format containing Genbank
#'                    accession IDs and coordinates to extract sequences from.
#'                    \cr For more information about the data.frame format see
#'                    'GBaccess.bed' argument in \link{prepare.gb.access}.
#'                    \cr For more information about the list format see the
#'                    example in \link{submit_NCBI_BLAST}.
#' @param db          A \code{character} string specifying an NCBI genome or
#'                    reference sequence set database on which the submitted
#'                    sequences will be BLASTed (Default using GRCh38 (hg38)
#'                    genome assembly database:
#'                    db = "genomic/9606/GCF_000001405.39").
#'                    For more supported databases see Details section in
#'                    \link{submit_NCBI_BLAST}.
#' @param res.dir     A \code{character} string specifying the directory where
#'                    NCBI BLAST submission results will be stored. Each BLAST
#'                    result is stored in a separate folder within the given
#'                    directory. Blast hits for each submission are available as
#'                    a XML file in each matching folder.
#' @param ncores      An \code{integer} specifying the number of cores or
#'                    threads to be used for parallel processing.
#' @param auto.rm.dir A \code{logical} specifying whether the results directory
#'                    should be immediately deleted after function execution end
#'                    (auto.rm.dir = TRUE) or not (auto.rm.dir = FALSE). If
#'                    FALSE, the folder will remain available after execution
#'                    with all raw results in it.
#' @param delay.req   Mandatory seconds between BLAST request submissions as an
#'                    \code{integer} (Default: delay.req = 10).
#' @param email       Mandatory E-mail adress for NCBI BLAST request, specified
#'                    as a \code{character} string.
#' @return A \code{data.table} aggregating all NCBI BLAST results.
#' @author Yoann Pageaud.
#' @export 
#' @examples
#' #Create 1 row data.frame for the example
#' df <- data.frame("GBaccessID" = "AC073318", "Start" = 71401, "End" = 72401)
#' #Submit the sequence extracted from AC073318 to NCBI BLAST API
#' example.dt <- get.NCBI.BLAST2DT(
#'   sequences = df, db = "genomic/9606/GCF_000001405.25",
#'   res.dir = "~/result_BLAST", ncores = 2, auto.rm.dir = FALSE,
#'   email = "myemailadress@dkfz.de")
#' #example.dt is a data.table containing all BLAST hits for the sequence submitted.
#' @references
#' \itemize{
#'  \item{\href{https://cran.r-project.org/web/packages/hoardeR/index.html}{hoardeR: Collect and Retrieve Annotation Data for Various Genomic Data Using Different Webservices.}}
#'  \item{\href{https://academic.oup.com/bioinformatics/article/35/3/526/5055127}{Paradis E. & Schliep K. 2019. ape 5.0: an environment for modern phylogenetics and evolutionary analyses in R. Bioinformatics 35: 526–528. doi:10.1093/bioinformatics/bty633. HAL: ird-01920132.}}
#'  \item{\href{https://www.ncbi.nlm.nih.gov/BLAST/tutorial/Altschul-1.html}{The Statistics of Sequence Similarity Scores.}}
#'  \item{\href{https://academic.oup.com/nar/article/36/suppl_2/W5/2505810}{Johnson, M. et al. NCBI BLAST: a better web interface. Nucleic Acids Research 36, W5–W9 (2008).}}
#' }

get.NCBI.BLAST2DT <- function(
  sequences, db = "genomic/9606/GCF_000001405.39", res.dir, ncores,
  auto.rm.dir = TRUE, delay.req = 10, email){
  #Convert data.frame into a list of sequences
  if(is.data.frame(sequences)){
    sequences <- NCBI.BLAST2DT::prepare.gb.access(
      GBaccess.bed = sequences, ncores = ncores)
  } else if(!is.list(sequences)){
    stop("Format not supported for 'sequences'.")
  }
  #Submit list of sequences to NCBI BLAST servers
  NCBI.BLAST2DT::submit_NCBI_BLAST(
    seq.list = sequences, res.dir = res.dir, delay.req = delay.req, email, db)
  #Get XML results into a data.table
  dt.res <- NCBI.BLAST2DT::aggregate_NCBI_BLAST_XMLs2DT(
    dir.to.xmls = res.dir, seq.names = names(sequences), ncores = ncores)
  #Remove result directory if requested
  if(auto.rm.dir){ unlink(x = res.dir) } else {
    message(paste("auto.rm.dir = FALSE: Raw results are stored in:", res.dir))
  }
  return(dt.res)
}