

#' Find a sample representative of a given population
#' by multifunctional matching
#'
#' @param data the dataframe to use
#' @param listfns the list of functionals
#' @param listinputs the list of inputs, where
#' each input is a vector of columns corresponding to each functional
#'
#' @param coordX the vector of X-coordinates
#' @param coordY the vector of Y-coordinates
#' @param minDistance minimum distance between two observations
#' of the selected sample
#'
#' @param nDraws the number of random samples to generate from the population
#' @param sizeSample the size of the sample
#' @param methodNormalization the method used for normalization.
#' Can be \code{ecdf} for normalization by the empirical cumulative distribution
#' function or \code{meansd} for the normalization \eqn{(X-E[X])/sd(X)}.
#'
#' @param weightsfns weights for each of the functionals. These weights can be used
#' to give more or less importance to some of the functionals.
#' It should be a vector of the same length as the \code{listfns}.
#'
#' @param progressBar \code{TRUE}: display a progressbar
#'
#' @return a list with two elements
#' \itemize{
#'    \item \code{sample}: the optimal set of (indices) of observations representing
#'    best the population.
#'    \item \code{crit}: the value of the MFM criteria.
#' }
#'
#' @examples
#' # matching the first mean of a two-dimensional dataset
#' df = data.frame(rnorm(2000), rnorm(2000))
#' mfm(data = df, listfns = list(mean), listinputs = list(1),
#'     nDraws = 500, sizeSample = 10, methodNormalization = "ecdf")
#'
#' @export
mfm <- function(data, listfns, listinputs,
                coordX = NULL, coordY = NULL, minDistance = 0,
                nDraws, sizeSample,
                methodNormalization, weightsfns = NULL, progressBar = TRUE)
{
  n = nrow(data)
  p = length(listfns)
  if (length(listinputs) != p){
    stop("the number of functionals should be the same as the number of inputs")
  }
  if (!is.null(weightsfns)){
    if ( length(weightsfns) != p){
      stop("the number of weights should be the sample as the number of functionals.")
    }
  }

  # First part: computing the statistics on each generated sample
  if (is.null(coordX) & is.null(coordY)){
    matrixResults = generateSample_unconstrained(
      data = data, listfns = listfns, listinputs = listinputs,
      nDraws = nDraws, sizeSample = sizeSample, progressBar = progressBar)
  } else {
    if (length(coordX) != n || length(coordY) != n) {
      stop("The vectors of coordinate should have the length as the number of observations.")
    }

    matrixResults = generateSample_spatial(
      data = data, listfns = listfns, listinputs = listinputs,
      coordX = coordX, coordY = coordY, minDistance = minDistance,
      nDraws = nDraws, sizeSample = sizeSample, progressBar = progressBar)
  }

  # Second part : computing the statistics on the population

  # We prepare the right method for normalization
  switch(
    methodNormalization,
    "ecdf" = { funNorm <- function (x) { return ( stats::ecdf(x)(x) )}},
    "meansd" = { funNorm <- function (x) {
      return ( (x - base::mean(x)) / stats::sd(x) )}},
    stop("Normalisation method not implemented yet.")
  )

  matrixNormsDiffs = matrix(nrow = nDraws, ncol = p)

  for (iFunctional in 1:p){
    populationFunct =
      do.call(what = listfns[[iFunctional]],
              args = list( data[ , listinputs[[iFunctional]] ]) )

    matrixNormsDiffs[ , iFunctional] =
      funNorm( abs(matrixResults[ , sizeSample + iFunctional] - populationFunct) )
  }

  if (is.null(weightsfns)) {
    vecCrit = apply(X = matrixNormsDiffs, MARGIN = 1, FUN = sum)
  } else {
    vecCrit = apply(X = matrixNormsDiffs, MARGIN = 1,
                    FUN = function(x){return(x %*% weightsfns)})
  }


  optimalChoice = which.min(vecCrit)
  optimalSetTrees = matrixResults[optimalChoice , 1:sizeSample]

  return( list(sample = optimalSetTrees, crit = vecCrit[optimalChoice]) )
}




#' Help function to construct the matrix of the generated samples
#' in the case where there is a spatial constraint.
#'
#' @keywords internal
#'
generateSample_spatial <- function(data, listfns, listinputs,
                           coordX, coordY, minDistance,
                           nDraws, sizeSample, progressBar)
{
  n = nrow(data)
  p = length(listfns)
  # Initalization of the data frame to store the results
  matrixResults = matrix(nrow = nDraws, ncol = sizeSample + p)

  if (progressBar) {pb <- pbapply::startpb(min = 0, max = nDraws)}
  for (k in 1:nDraws){

    # This is the vector that control if a given observation
    # has any chances to be selected or not.
    isPossibleObs = rep(TRUE, n)

    matrixResults[k, 1] = sample(n, size = 1)

    # We create the list of remaining trees, at the minimum distance
    distanceToObs = sqrt( (coordX - coordX[matrixResults[k, 1]])^2
                          + (coordY - coordY[matrixResults[k, 1]])^2 )

    isPossibleObs = isPossibleObs & (distanceToObs > minDistance)

    # We sample the other observations, starting with the observation n°2
    for (iObs in 2:sizeSample)
    {
      # We test if one more observation can be added given the constraints

      if (length(which(isPossibleObs)) == 0){
        stop(paste0("Too many constraints! ",
                    "There is no more observation available to complete our sets of ",
                    sizeSample,
                    ". Stopping at ", iObs - 1, "observations."))
      }

      # We sample the new observation among the elements of the list of possible observation
      matrixResults[k, iObs] = sample(x = which(isPossibleObs), size = 1)

      # Then we remove the trees that are too far away from it
      distanceToObs = sqrt( (coordX - coordX[matrixResults[k, iObs]])^2
                            + (coordY - coordY[matrixResults[k, iObs]])^2 )

      isPossibleObs = isPossibleObs & (distanceToObs > minDistance)
    }

    # Computation of the functionals
    currentSample = matrixResults[k, 1:sizeSample]
    for (iFunctional in 1:p){
      matrixResults[k, sizeSample + iFunctional] =
        do.call(what = listfns[[iFunctional]],
                args = list( data[currentSample ,
                                  listinputs[[iFunctional]] ]) )
    }
    if (progressBar) {pbapply::setpb(pb = pb, value = k)}
  }
  if (progressBar) {pbapply::closepb(pb)}

  return(matrixResults)
}



#' Help function to construct the matrix of the generated samples
#' in the unconstrained case.
#'
#' @keywords internal
generateSample_unconstrained <- function(
  data, listfns, listinputs, nDraws, sizeSample, progressBar)
{
  n = nrow(data)
  p = length(listfns)

  # Initalization of the data frame to store the results
  matrixResults = matrix(nrow = nDraws, ncol = sizeSample + p)

  if (progressBar) {pb <- pbapply::startpb(min = 0, max = nDraws)}
  for (k in 1:nDraws){

    # This is the vector that control if a given observation
    # has any chances to be selected or not.
    isPossibleObs = rep(TRUE, sizeSample)

    matrixResults[k, 1:sizeSample] = sample(n,
                                            size = sizeSample, replace = FALSE)

    # Computation of the functionals
    currentSample = matrixResults[k, 1:sizeSample]
    for (iFunctional in 1:p){
      matrixResults[k, sizeSample + iFunctional] =
        do.call(what = listfns[[iFunctional]],
                args = list( data[currentSample ,
                                  listinputs[[iFunctional]] ]) )
    }
    if (progressBar) {pbapply::setpb(pb = pb, value = k)}
  }
  if (progressBar) {pbapply::closepb(pb)}

  return(matrixResults)
}



