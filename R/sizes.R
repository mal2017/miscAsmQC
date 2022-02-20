gsize <- function(obj) {sum(Biostrings::width(obj))}

#' extract info about contig sizes
#' @rdname sizes
#' @export
n_largest_csizes <- function(obj, n=1) {sort(Biostrings::width(obj),decreasing = T)[1:n]}


