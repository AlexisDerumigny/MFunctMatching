
#' Choice of the optimal sample size for the MFM algorithm
#'
#' @param data the dataframe to use
#' @param rangeSizeSample the list of the sample sizes to tested
#' @param alpha the percentage of reduction of the criterion
#' @param listfns list of functionals to be used. If \code{NULL} (the default)
#' then the mean, the standard deviation and Kendall's tau between each variables
#' are used.
#' @param listinputs list of inputs to the functionals.
#' Only used when \code{listfns} is not NULL.
#' @param coordX X coordinates of the points. Used for spatial data.
#' @param coordY Y coordinates of the points. Used for spatial data.
#' @param minDistance minimum distance between two observations
#' of the selected sample
#'
#' @param nDraws number of draws to compute each MFM criterion.
#' @param nReplications number of replications of the MFM criterion computation.
#'
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
#' @examples
#' variables = c("mass", "moist", "elev")
#' df = agridat::gartner.corn
#' result = mfm_sampleSize(
#'   data = df[,variables], rangeSizeSample = c(15,25,30),
#'   nDraws = 500, nReplications = 20,
#'   coordX = df$long, coordY = df$lat)
#'
#' @export
#'
mfm_sampleSize <- function(
  data, rangeSizeSample,
  alpha = 0.9,
  listfns = NULL, listinputs = NULL,
  coordX = NULL, coordY = NULL, minDistance = 0,
  nDraws = 2000, nReplications = 1000,
  methodNormalization = "ecdf", weightsfns = NULL, progressBar = TRUE)
{
  # Converting the data
  data = as.data.frame(data)
  n = nrow(data)

  # If there is no functionals provided, use the standard ones
  if (is.null(listfns)){
    if (!is.null(listinputs)){
      stop("Inputs provided but without any functionals")
    }
    listFunctionals = construct_fns(data)
    listfns <- listFunctionals$listfns
    listinputs <- listFunctionals$listinputs
    weightsfns <- listFunctionals$weightsfns

    p = length(listfns)
  } else {
    # Checking the user-provided functionals
    p = length(listfns)
    if (length(listinputs) != p){
      stop("the number of functionals should be the same as the number of inputs")
    }
    if (!is.null(weightsfns)){
      if ( length(weightsfns) != p){
        stop("the number of weights should be the sample as the number of functionals.")
      }
    }
  }

  # Checking the spatial data and choosing the function for sample generation
  if (is.null(coordX) & is.null(coordY)){
    generation_sample <- function (sizeSample){
      return(generateSample_unconstrained(
        data = data, listfns = listfns, listinputs = listinputs,
        nDraws = nDraws, sizeSample = sizeSample, progressBar = FALSE))
    }
  } else {
    if (length(coordX) != n || length(coordY) != n) {
      stop("The vectors of coordinate should have the length as the number of observations.")
    }
    generation_sample <- function (sizeSample){
      return(generateSample_spatial(
        data = data, listfns = listfns, listinputs = listinputs,
        coordX = coordX, coordY = coordY, minDistance = minDistance,
        nDraws = nDraws, sizeSample = sizeSample, progressBar = FALSE))
    }
  }

  # We prepare the right method for normalization
  switch(
    methodNormalization,
    "ecdf" = { funNorm <- function (x) { return ( stats::ecdf(x)(x) )}},
    "meansd" = { funNorm <- function (x) {
      return ( (x - base::mean(x)) / stats::sd(x) )}},
    stop("Normalisation method not implemented yet.")
  )

  # Computation of the population functionals
  populationFunct = rep(NA, p)
  for (iFunctional in 1:p){
    populationFunct[iFunctional] =
      do.call(what = listfns[[iFunctional]],
              args = list( data[ , listinputs[[iFunctional]] ]) )
  }

  # Main loop
  matrix_Crit = matrix(nrow = length(rangeSizeSample), ncol = nReplications)

  if (progressBar) {pb <- pbapply::startpb(
    min = 0, max = length(rangeSizeSample) * nReplications)}
  for (isizeSample in 1:length(rangeSizeSample))
  {
    sizeSample = rangeSizeSample[isizeSample]

    for (iReplication in 1:nReplications){
      matrixResults = generation_sample(sizeSample)

      matrixNormsDiffs = matrix(nrow = nDraws, ncol = p)

      for (iFunctional in 1:p){
        matrixNormsDiffs[ , iFunctional] =
          funNorm( abs(matrixResults[ , sizeSample + iFunctional] - populationFunct[iFunctional]) )
      }

      if (is.null(weightsfns)) {
        vecCrit = apply(X = matrixNormsDiffs, MARGIN = 1, FUN = sum)
      } else {
        vecCrit = apply(X = matrixNormsDiffs, MARGIN = 1,
                        FUN = function(x){return(x %*% weightsfns)})
      }

      optimalChoice = which.min(vecCrit)
      matrix_Crit[isizeSample , iReplication] = vecCrit[optimalChoice]

      if (progressBar) {
        pbapply::setpb(pb = pb, value = (isizeSample-1) * nReplications + iReplication)
      }
    }
  }
  if (progressBar) {pbapply::closepb(pb)}

  rownames(matrix_Crit) <- rangeSizeSample
  vecMeanCrit = apply(X = matrix_Crit, MARGIN = 1, FUN = base::mean)

  optimalSampleSize = which( (vecMeanCrit[1] - vecMeanCrit) >=
                               alpha * (utils::tail(vecMeanCrit,1) - vecMeanCrit[1]) )[1]

  return(list( optimalSampleSize = optimalSampleSize, vecMeanCrit = vecMeanCrit,
               matrix_Crit = matrix_Crit))
}

