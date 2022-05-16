
#' Monte Carlo Simulation of Bivariate Multimodal Data
#'
#' The present function can be used to simulate a set of 2 variables; 1 categorical
#' and 1 numerical.
#'
#' @param data_set An object of type data.frame containing a single categorical
#' variable and a single numerical variable that the user wishes to augment
#' @param n_simulations An integer defining the number of simulations
#' @param num_breaks An optional parameter that can be used to fine-tune the calculation
#' of the probability distribution from which to sample from.
#'
#' @section Details:
#' This function uses the statistical chain rule to calculate probability
#' vectors from which to sample from so as to sample from the multivariate
#' categorical-numerical distribution
#'
#' @return simulated data
#'
#' @author Lloyd A. Courtenay
#'
#' @examples
#'
#' bivariate_trial <- data.frame(x = c(rnorm(100, 0, 2), rnorm(100, -4, 3)),
#'                               sample = as.factor(c(rep("Group_1", 100),
#'                                                    rep("Group_2", 100))))
#'
#' augmented_trial <- bivariate_multimodalMC(bivariate_trial, 100)
#'
#' par(mfrow = c(2,2))
#' hist(bivariate_trial[bivariate_trial$sample == "Group_1",1],
#'      main = "Original Data Group 1",
#'      xlab = "Group 1")
#' hist(c(bivariate_trial[bivariate_trial$sample == "Group_1",1],
#'        augmented_trial[augmented_trial$sample == "Group_1",1]),
#'      main = "Augmented Data Group 1",
#'      xlab = "Group 1")
#' hist(bivariate_trial[bivariate_trial$sample == "Group_2",1],
#'      main = "Original Data Group 2",
#'      xlab = "Group 2")
#' hist(c(bivariate_trial[bivariate_trial$sample == "Group_2",1],
#'        augmented_trial[augmented_trial$sample == "Group_2",1]),
#'      main = "Augmented Data Group 2",
#'      xlab = "Group 2")
#' par(mfrow = c(1,1))
#'
#' #
#'
#' @export

bivariate_multimodalMC <- function(data_set, n_simulations, num_breaks = 0) {

  if(!is.numeric(n_simulations) | n_simulations <= 0 | n_simulations %% 1 != 0) {
    stop("The number of simulations must be a non-negative integer, e.g. 2")
  }

  if(num_breaks != 0) {
    if(num_breaks < 0 | !is.numeric(num_breaks) | num_breaks %% 1 != 0) {
      stop("The number of breaks must be a non-negative integer, e.g. 2, or left at the default value of 0")
    }
  }

  if ("data.frame" %!in% class(data_set)) {
    stop("Input data must be in the format of a data.frame")
  }

  if(!missing_data_type(data_set, type = "Character")) {
    stop("The dataset contains character variables. Consider converting using as.factor")
  }

  if (missing_data_type(data_set, type = "Factor")) {
    stop("The dataset does not contain any qualitative (factor) variables")
  }

  if (missing_data_type(data_set, type = "Numeric")) {
    stop("The dataset does not contain any quantitative (numeric) variables")
  }

  if (ncol(data_set) < 2) {
    stop("Datasets must contain at least 2 variables")
  }

  original_col_names <- colnames(data_set)
  original_column_orders <- c()

  factor_list <- find_factor_columns(data_set, original_column_orders)
  number_list <- find_numeric_columns(data_set, factor_list$original_column_orders)

  original_column_orders <- number_list$original_column_orders
  factor_variables <- factor_list$factor_variables
  numeric_variables <- number_list$numeric_variables

  n_factors <- factor_list$n_variables
  n_numbers <- number_list$n_variables

  if ((n_factors > 1) & (n_numbers > 1)) {
    stop(paste0(
      "This function only accepts 2 columns, 1 categorical and 1 numeric variable.",
      " If the user wishes to simulate more variables, try conditionalMC or multimodalMC",
      " instead."
    ))
  }

  simulated_data_set <- array(dim = c(0, ncol(data_set)))
  iteration = 0
  pb <- txtProgressBar(min = 0, max = n_simulations,
                       style = 3, width = 100, char = "=")

  while (dim(simulated_data_set)[1] < n_simulations) {

    factor_variable <- sample_qualitative_value(factor_variables)

    remaining_numeric <- conditioned_filter_quantitative_value(
      numeric_variables,
      factor_variables,
      factor_variable
    )

    remaining_numeric <- as.numeric(remaining_numeric)

    numeric_variable <- single_value_simulator(
      remaining_numeric, num_breaks = num_breaks
    )$value

    if (NA %!in% c(factor_variable, numeric_variable)) {
      iteration = iteration + 1
      single_iteration <- c(factor_variable, numeric_variable)
      simulated_data_set <- abind::abind(simulated_data_set, single_iteration,
                                         along = 1)
    }

    if(!is.null(dim(simulated_data_set))) {
      setTxtProgressBar(pb, dim(simulated_data_set)[1])
    }

  }

  augmented_data <- data.frame(simulated_data_set)
  augmented_data[,1] <- as.factor(augmented_data[,1])
  augmented_data[,2] <- as.numeric(augmented_data[,2])

  augmented_data_return <- augmented_data
  for (order in 1:ncol(augmented_data)) {
    augmented_data_return[, original_column_orders[order]] <- augmented_data[,order]
  }
  colnames(augmented_data_return) <- original_col_names

  cat("\n\n")

  return(augmented_data_return)

}
