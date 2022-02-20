# -----------------------------------------------------------------------------
# E-size
# -----------------------------------------------------------------------------

#' Generate E-size metric from DNAStringSet
#' @param obj DNAStringSet
#' @importFrom magrittr %>%
.e_size <- function(obj) {
  G <- gsize(obj)

  tibble::as_tibble(Biostrings::seqinfo(obj),rownames = "seqnames")[-3:-4] %>%
    dplyr::mutate(LC2 = seqlengths^2) %>%
    dplyr::pull(LC2) %>%
    magrittr::divide_by(G) %>%
    sum()
}

#' Generate `E-size` metric.
#' For details see methods for https://dx.doi.org/10.1101%2Fgr.131383.111
#' @rdname e_size
#' @param obj character, DNAStringSet
#' @export
setGeneric("e_size", function(obj) standardGeneric("e_size"))

#' @rdname e_size
setMethod("e_size", signature(obj="DNAStringSet"),
  function(obj){
    .e_size(obj)
})

#' @rdname e_size
setMethod("e_size", signature(obj="character"),
          function(obj){
            stopifnot(file.exists(obj))
            obj <- Biostrings::readDNAStringSet(obj)
            .e_size(obj)
})

# -----------------------------------------------------------------------------
# N50 and related metrics
# -----------------------------------------------------------------------------

#' Generate N50 or related metrics.
#' @param obj DNAStringSet
#' @importFrom magrittr %>%
.NX0 <- function(obj, percentile = 50L) {
  G <- gsize(obj)

  contig_sizes <- Biostrings::width(obj) %>% sort(decreasing = T)

  cumulative_sizes <- cumsum(contig_sizes)

  GX0 <- 0.01 * percentile * G

  #message(paste("assembly size after filtering:",G))
  #message(paste0(percentile,"% of assembly size after filtering: ",GX0))

  max(contig_sizes[cumulative_sizes >= GX0])
}


#' Generate N50 or related metrics.
#' generate N50, N90, etc.
#' @rdname NX0
#' @param obj character, DNAStringSet
#' @export
setGeneric("NX0",
           function(obj, percentile = 50L) {standardGeneric("NX0")},
           signature=c("obj"))

#' @rdname NX0
setMethod("NX0", signature(obj="DNAStringSet"),
          function(obj, percentile=50L){
            .NX0(obj, percentile=percentile)
          })

#' @rdname NX0
setMethod("NX0", signature(obj="character"),
          function(obj, percentile=50L){
            stopifnot(file.exists(obj))
            obj <- Biostrings::readDNAStringSet(obj)
            .NX0(obj, percentile=percentile)
          })


# -----------------------------------------------------------------------------
# AUN curve generation
# -----------------------------------------------------------------------------

#' Generate auN curve
#' Inspired by http://lh3.github.io/2020/04/08/a-new-metric-on-assembly-contiguity
#'
#' Generates N10, N20 ... N100 rather than NG*.
#' @rdname nx_curve
#' @importFrom magrittr %>%
#' @param obj DNAStringSet, character
#' @export
nx_curve <- function(obj,percentiles = as.integer(seq(0L,100L,10L))) {

  magrittr::set_names(percentiles,percentiles) %>%
    purrr::map_int(~NX0(obj,percentile = .x)) %>%
    tibble::enframe("percentile","csize")

}

#' @rdname nx_curve
#' @param obj DNAStringSet, character
#' @export
multi_nx_curve <- function(objs, percentiles = seq(0L,100L,10L)) {
  if (is.character(objs)){
    stopifnot(all(file.exists(objs)))
  }
  purrr::map_df(objs,nx_curve,percentiles,.id="assembly")
}


