
#' Rattini shape dataset
#'
#' A geometric morphometric derived dataset from the data of Hulme-Beaman et al.
#' 2019. The dataset contains an array of 395 specimens each with 105 landmarks.
#' Centroid size data is also included as is individual specimen info. A matrix
#' of full Procrustes distances calculated among the specimens is also provided.
#'
#' @format A list of 2 objects:
#' \describe{
#'   \item{LMs}{Specimen landmark data in array form, where each row is a landmark, each column a landmark dimension and each slice is an individual specimen. These landmarks have been scaled.}
#'   \item{Info}{A dataframe corresponding to the LM data which includes specimen genus, species and centroid size (CS) data.}
#' }
#' @source \url{https://doi.org/10.1007/s10914-017-9423-8}
"RatData"



#' Black Rat shape dataset
#'
#' A geometric morphometric derived dataset containing an data from 100 specimens
#' with 120 landmarks all collected from Mainland Southeast Asia.
#'
#' @format A list of 2 objects:
#' \describe{
#'   \item{LMArray}{Specimen landmark data in array form, where each row is a landmark, each column a landmark dimension and each slice is an individual specimen. These landmarks have been scaled.}
#' }
#' @source \url{https://doi.org/10.1007/s10914-017-9423-8}
"BlackRat"




