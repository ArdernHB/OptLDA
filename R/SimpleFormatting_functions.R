

#' Internal function for conversion of a matrix of landmarks to an array
#'
#' This function takes a matrix where rows are specimens and columns are
#' alternating landmark coordinate values and converts it to a 3 dimensional
#' array where rows are landmarks and columns correspond with the dimensions.
#' @return
#'
#' This function returns an array
#'
#' @param Mat is a matrix of landmark data
#' @param LMdim is the number of dimensions of the data, either 2 or 3
#'
#' @keywords internal
#' @author Ardern Hulme-Beaman
#' @export



Mat2Array <- function(Mat, LMdim){
  NewArray <- array(data = NA, dim = c(dim(Mat)[2]/LMdim, LMdim, dim(Mat)[1]))


  for (i in 1:dim(Mat)[1]){
    #i <- 1
    Mat4Array <- matrix(as.numeric(Mat[i,]), nrow = dim(Mat)[2]/LMdim, ncol = LMdim, byrow = TRUE)
    NewArray[,,i] <- Mat4Array
  }


  return(NewArray)
}



#' Internal function for conversion of an array to a matrix
#'
#' This function takes an array and converts it to a matrix where rows
#' specimens and columns are alternating landmark variables
#' (e.g. X1, X2, Y1, Y2...)
#' @return
#'
#' This function returns an array
#'
#' @param Array is an array of landmark data
#'
#' @keywords internal
#' @author Ardern Hulme-Beaman
#' @export


Array2Mat <- function(Array){
  Matrix <- matrix(NA, nrow = dim(Array)[3], ncol = length(c(t(Array[,,1]))))
  for (i in 1:dim(Array)[3]){
    #i <- 1
    Matrix[i,] <- c(t(Array[,,i]))
  }
  return(Matrix)
}

