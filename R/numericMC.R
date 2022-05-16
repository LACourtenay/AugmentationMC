
#' Monte Carlo Simulations for Numeric Variables
#'
#' The present function can be used to simulate numerical (quantitative) data
#'
#' @param data_set An matrix or array of numeric data to simulate
#' @param n_simulations An integer defining the number of simulations
#' @param num_breaks An optional parameter that can be used to fine-tune the calculation
#' of the probability distribution from which to sample from.
#'
#' @section Details:
#' This function uses the statistical chain rule to calculate multivariate probability
#' vectors associated with each categorical variable from which to sample from
#'
#' @return simulated data
#'
#' @author Lloyd A. Courtenay
#'
#' @examples
#'
#' data(multiclass_dataset)
#' data_set <- as.matrix(multiclass_dataset[,c(1, 2)])
#'
#' # augment data
#'
#' augmented_data <- numericMC(data_set, 1000)
#'
#' #
#'
#' @export

numericMC <- function(data_set, n_simulations, num_breaks = 0) {

  if(!is.numeric(n_simulations) | n_simulations <= 0 | n_simulations %% 1 != 0) {
    stop("The number of simulations must be a non-negative integer, e.g. 2")
  }

  if ("matrix" %!in% class(data_set) & "array" %!in% class(data_set)) {
    stop("Input data must be a matrix or array")
  }

  if(missing_data_type(data_set, type = "Numeric")) {
    stop("The dataset does not contain any numeric variables")
  }

  if(!missing_data_type(data_set, type = "Character")) {
    stop("The dataset contains invalid input for this function")
  }

  if(!missing_data_type(data_set, type = "Factor")) {
    stop("The dataset contains invalid input for this function")
  }

  pb <- txtProgressBar(min = 0, max = n_simulations,
                       style = 3, width = 100, char = "=")

  if (dim(data_set)[2] == 1) {

    simulated_data <- c()

    while(length(simulated_data) < n_simulations) {

      simulated_data <- c(
        simulated_data,
        single_value_simulator(data_set,
                               num_breaks = num_breaks)$value
      )

      if(!is.null(length(simulated_data))) {
        setTxtProgressBar(pb, length(simulated_data))
      }

    }

  } else {

    simulated_data <- array(dim = c(0, dim(data_set)[2]))

    while(dim(simulated_data)[1] < n_simulations) {

      single_numeric_simulation <- multivariate_prob_sampling_v2(
        data_set, 1, num_breaks = num_breaks
      )

      if (NA %!in% single_numeric_simulation) {
        simulated_data <- abind::abind(simulated_data,
                                       single_numeric_simulation,
                                       along = 1)
      }

      if(!is.null(dim(simulated_data))) {
        setTxtProgressBar(pb, dim(simulated_data)[1])
      }

    }

  }

  cat("\n\n")

  return(simulated_data)

}
