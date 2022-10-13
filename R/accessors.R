#' get the dataframe of ALL PFAM results from a `buscar` search object from `buscar()`
#'
#' @param b buscar search object returned from `buscar`
#' @return dataframe
#' @export
pfam_results <- function(b) {
  b$pfams
}

#' get the dataframe of passing deeptmhmm results from a `buscar` search object from `buscar()`
#'
#' @param b buscar search object returned from `buscar`
#' @return dataframe
#' @export
deeptmhmm_results <- function(b) {
  b$tm_signal_pep
}
#' get the dataframe of ALL blast results from a `buscar` search object from `buscar()`
#'
#' @param b buscar search object returned from `buscar`
#' @return dataframe
#' @export
blast_results <- function(b) {
  b$blasts
}


#' get the dataframe of classification results from a `buscar` search object from `buscar()`
#'
#' @param b buscar search object returned from `buscar`
#' @return dataframe
#' @export
groups <- function(b) {
  b$classes
}

#' get the computed classification matrix for the proteins
#'
#' @param b buscador object
#' @return dataframe
#' @export
classification_matrix <- function(b) {
  b$matrix
}
