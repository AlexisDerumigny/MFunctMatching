

#' Standard version of the MFM algorithm
#'
#' This function uses the MFM algorithm with predefined functionals,
#' which are:
#' \itemize{
#'    \item the mean and standard deviation for each \code{numeric} or \code{integer} variable
#'    \item Kendall's tau for each pair of \code{numeric} / \code{integer} variable;
#'    \item the empirical proportion for each unique possibility of each variable of
#'    type \code{character} or \code{factor}.
#' }
#'
#' @param data the dataframe to use
#' @param nDraws the number of random samples to generate from the population
#' @param sizeSample the size of the sample
#' @param methodNormalization the method used for normalization
#'
#' @param coordX the vector of X-coordinates
#' @param coordY the vector of Y-coordinates
#' @param minDistance minimum distance between two observations
#' of the selected sample
#'
#' @examples
#' df = data.frame(X1 = rnorm(2000), X2 = rnorm(2000),
#'                 X3 = c(rep("A",500), rep("B", 500), rep("C", 1000)))
#' mfm_standard(data = df, nDraws = 500, sizeSample = 10,
#'              coordX = rnorm(2000), coordY = rnorm(2000))
#'
#' # agricultural example
#' df = agridat::gartner.corn
#' mfm_standard(data = df[,c("mass", "moist", "elev")],
#'              nDraws = 500, sizeSample = 10,
#'              coordX = df$long, coordY = df$lat)
#'
#' @export
#'
mfm_standard <- function(data, coordX = NULL, coordY = NULL,
                         minDistance = 0, nDraws, sizeSample,
                         methodNormalization ="ecdf")
{
  data = as.data.frame(data)
  listfunctionals = construct_fns(data)

  result = mfm(
    data = data,
    listfns = listfunctionals$listfns,
    listinputs = listfunctionals$listinputs,
    coordX = coordX, coordY = coordY, minDistance = minDistance,
    nDraws = nDraws, sizeSample = sizeSample,
    methodNormalization = methodNormalization,
    weightsfns = listfunctionals$weightsfns)

  return (result)
}


#' Construct the standard list of functionals
#'
#' These functionals are:
#' \itemize{
#'    \item the mean and standard deviation for each \code{numeric} or \code{integer} variable
#'    \item Kendall's tau for each pair of \code{numeric} / \code{integer} variable;
#'    \item the empirical proportion for each unique possibility of each variable of
#'    type \code{character} or \code{factor}.
#' }
#'
#' @param data the data frame of observations from which to construct the functionals.
#'
#' @return a list with
#' \itemize{
#'    \item listfns: the constructed functionals
#'    \item listinputs: their inputs
#'    \item weightsfns: the weights that they should have so that each variable
#'    (categorical and continuous) has the influence
#' }
#'
construct_fns <- function(data)
{

  d = ncol(data)
  vecTypes = sapply(data, class)
  d_num = length(which(vecTypes == "numeric" | vecTypes == "integer"))

  # 1 functional for mean + 1 functional for sd + 1/2 a functional for each other variable
  weightsnum = 2 + (d_num-1)/2

  listfns = list()
  listinputs = list()
  weightsfns = c()

  pos = 1
  for (ivar in 1:d) {
    if (vecTypes[ivar] == "numeric" || vecTypes[ivar] == "integer") {
      listfns[[pos]] <- base::mean
      listinputs[[pos]] <- ivar
      pos = pos + 1

      listfns[[pos]] <- stats::sd
      listinputs[[pos]] <- ivar
      pos = pos + 1

      weightsfns = c(weightsfns, rep( 1/weightsnum, 2))

      if (ivar < d){
        for (jvar in (ivar+1):d) {
          if (vecTypes[jvar] == "numeric" || vecTypes[jvar] == "integer"){
            listfns[[pos]] <- function(...) {return(pcaPP::cor.fk(...)[1,2])}
            listinputs[[pos]] <- c(ivar, jvar)
            pos = pos + 1
            weightsfns = c(weightsfns,  2/weightsnum)
          }
        }
      }

    } else if (vecTypes[ivar] == "factor" || vecTypes[ivar] == "character" ) {
      uniqueCategories = unique(data[, ivar])
      for (i in 1:length(uniqueCategories)){
        listfns[[pos]] <- function(x){return(
          base::mean(as.numeric(x == uniqueCategories[i])) )}
        listinputs[[pos]] <- ivar
        pos = pos + 1
      }

      weightsfns = c(weightsfns, rep( 1/length(uniqueCategories), length(uniqueCategories)))
    }

  }

  return( list(listfns=listfns, listinputs = listinputs, weightsfns = weightsfns))
}

