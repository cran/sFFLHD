#' Split a matrix by rows, based on either the number of rows per group
#' or number of splits.
#'
#' @param mat A matrix to be split.
#' @param rowspergroup  Number of rows in a group.
#' @param nsplits Number of splits to make.
#' @param shuffle Should the splits be shuffled before returning?
#'
#' @return A list of the splits of the matrix.
#' @export
#'
#'
#' @examples
#' mat <- matrix(1:12, ncol=2)
#' split_matrix(mat, 4, shuffle=FALSE)
#' split_matrix(mat, 4, shuffle=TRUE)
#' split_matrix(mat, nsplits=3, shuffle=FALSE) # same as 4 rowspergroup
split_matrix <- function(mat,rowspergroup=NULL,nsplits=NULL,shuffle=TRUE) {
  if(is.null(rowspergroup)) {
    rowspergroup <- ceiling(nrow(mat) / nsplits)
  } else {
    nsplits <- ceiling(nrow(mat) / rowspergroup)
  }
  lapply(ifelse(shuffle,sample,identity)(1:nsplits),
         function(ii){mat[((ii-1)*rowspergroup+1):min(ii*rowspergroup, nrow(mat)), , drop=FALSE]})
}
