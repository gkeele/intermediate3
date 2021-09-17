#' @name Tmem68
#' @aliases Tmem68
#' 
#' @title Liver protein expression dataset
#'
#' @details
#' Data for mediation analysis of liver protein expression.
#' Dataset is a list of objects for `target`, `mediator`, `annotation`, `covar` (covariates),
#' `driver`.
#' There are 192 diversity outbred mice and 8050 genes with measured mRNA expression.
#' The 'target' is the level of expression of gene 'Tmem68'.
#' The mRNA transcript expression is measured for 8050 genes in 'mediator'.
#' The 'annotation' has positional information for each of the genes.
#' The 'covar' has Sex and Diet information for the 192 mice.
#' The 'driver' object has the driver as allele probabilities for the 8 founder genotypes for the 192 mice
#' at the location of the Tmem68 gene.
#' 
#' This version only contains data for chr 4 and 13, and has had the following changes
#' to elements:
#' `covar`: first two columns renamed to `Intercept` and `SexM`;
#' `driver`: renames from the original `qtl.geno`.
#' See example below for more details.
#'
#' @docType data
#'
#' @usage data(Tmem68)
#'
#' @format A list object.
#'
#' @keywords datasets
#'
#' @references Chick et al. (2016) 
#' (\href{https://dx.org/10.1038/nature18270}{Nature 534: 500-505}) 
#'
#' @source \href{https://github.com/churchill-lab/intermediate}{Churchill GitHub}
#'
#' @examples
#' data(Tmem68)
#' str(Tmem68)
#' 
#' \dontrun{
#' # Here are details for transforming data from original.
#' 
#' # Download file from Churchill Lab and load
#' tmpfile <- tempfile()
#' file <- paste0(file.path("https://raw.githubusercontent.com/churchill-lab",
#'                          "intermediate/master/data/Tmem68.rda"))
#' download.file(file, tmpfile)
#' load(tmpfile)
#' unlink(tmpfile)
#' 
#' # Fix column names of covariates.
#' colnames(Tmem68$covar)[1:2] <- c("Intercept", "SexM")
#'
#' # Reconstruct 8-allele genotype probabilities.
#' m <- match("qtl.geno", names(Tmem68))
#' if(!is.na(m)) names(Tmem68)[m] <- "driver"
#' if(ncol(Tmem68$driver) == 7) {
#'   Tmem68$driver <- cbind(A = 1 - apply(Tmem68$driver, 1, sum), Tmem68$driver)
#' }
#' 
#' # Subset to chr 4 and 13
#' id <- which(Tmem68$annotation$chr %in% c("4", "13"))
#' Tmem68$mediator <- object$mediator[, id, drop = FALSE]
#' Tmem68$annotation <- object$annotation[id,, drop = FALSE]
#' 
#' # Save file.
#' save(Tmem68, "Tmem68.rda")
#' 
#' # Resave file in more compressed format.
#' tools::resaveRdaFiles("Tmem68.rda")
#' }

NULL


