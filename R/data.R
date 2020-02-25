#' Partial results of a differential methylation test.
#'
#' Output from a limma analysis comparing DNA methylation status between 4 
#' healthy and 4 cancerous samples. Data are from an Illumina 450k chip and 
#' distributed with the \code{ChAMP} package. Output only includes CpGs from 
#' chromosome 13. See epigen.R for complete preprocessing pipeline.  
#'
#' @format A data frame with 10286 rows and 6 columns:
#' \describe{
#'   \item{cpg}{CpG island id.}
#'   \item{pos}{Genomic location.}
#'   \item{chr}{Chromosome.}
#'   \item{z}{z-statistic from differential methylation test.}
#'   \item{p.value}{Corresponding p-value.}
#'   \item{BH_q.value}{Corresponding Benjamini-Hochberg q-value.}
#' }
#' @source \url{https://bioconductor.org/packages/release/bioc/html/ChAMP.html}
"DNAm"