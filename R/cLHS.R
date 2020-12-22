

cLHS <- function(data, sizeSample, niter = 1000,
                 factorTemp = 0.9, prob_p = 0.8)
{
  data = data.frame(data)
  N = nrow(data)
  Temperature = 1

  vecTypes = sapply(data, class)
  num_variables = as.numeric(which(vecTypes == "numeric" | vecTypes == "integer"))

  # Step 1-
  matQuantiles = matrix(nrow = length(num_variables), ncol = sizeSample+1)
  for (k in 1:length(num_variables))
  {
    matQuantiles[k, ] = as.numeric(
      stats::quantile(x = data[, num_variables[k]],
                      probs = (0:sizeSample)/sizeSample))
  }

  matC = stats::cor(x = data[, num_variables])

  # Step 2-
  my_sample = sample(x = N, size = sizeSample)

  # Follow-up of the objective
  vec_obj = rep(NA, niter+1)
  vec_Metro = rep(NA, niter)
  vec_obj[1] = Inf

  for (i_iter in 1:niter)
  {
    # Step 3-

    mat_eta = matrix(nrow = sizeSample, ncol = length(num_variables))
    for (k in 1:length(num_variables))
    {
      data_sample_k = data[my_sample, num_variables[k]]
      for (i in 1:sizeSample)
      {
        mat_eta[i, k] = length(which(matQuantiles[k, i] <= data_sample_k &
                                       data_sample_k < matQuantiles[k, i+1]))
      }
    }
    Obj1 = sum(abs(mat_eta - 1))

    # Categorical data -> to be implemented later
    Obj2 = 0

    matT = stats::cor(x = data[my_sample, num_variables])
    Obj3 = sum(abs(matT - matC))

    Obj = Obj1 + Obj2 + Obj3

    # Step 4-
    vec_obj[i_iter+1] <- Obj
    Metro = exp( (vec_obj[i_iter] - vec_obj[i_iter+1]) / Temperature)
    vec_Metro[i_iter] = Metro

    # Step 5-
    if (stats::runif(1) < Metro){
      old_sample = my_sample
    } else {
      my_sample = old_sample
    }

    # Step 6
    if (stats::runif(1) < prob_p){
      my_sample[sample(x = sizeSample, size = 1)] <- sample(x = (1:N)[-my_sample], size = 1)
    } else {
      k = sample(1:length(num_variables), size = 1)
      bad_quantiles_i = which.max(mat_eta[, k])
      which_bad = which(matQuantiles[k, bad_quantiles_i] <= data_sample_k &
                          data_sample_k < matQuantiles[k, bad_quantiles_i+1])[1]
      my_sample[which_bad] <- sample(x = (1:N)[-my_sample], size = length(which_bad))

    }

    Temperature = Temperature * factorTemp
  }

  return(list(my_sample = my_sample, vec_obj = vec_obj, vec_Metro = vec_Metro))
}

