
#' Conditional Monte Carlo Simulations for Categorical and Numeric Variables
#'
#' The present function can be used to simulate a set of numeric values related to
#' a qualitative variable, or a set of categorical values related to a quantitative variable
#'
#' @param data_set An object of type data.frame containing either a single categorical
#' variable associated with a series of numeric variables, or a single numeric variable associated
#' with a series of categorical variables.
#' @param n_simulations An integer defining the number of simulations
#' @param num_breaks An optional parameter that can be used to fine-tune the calculation
#' of the probability distribution from which to sample from.
#'
#' @section Details:
#' This function uses the statistical chain rule to calculate multivariate probability
#' vectors from which to sample from, firstly calculating the probability of
#' the conditioning variable, and then sampling values associated with these
#' conditioning variables.
#'
#' @return simulated data
#'
#' @author Lloyd A. Courtenay
#'
#' @examples
#'
#' # example 1 ---------------------
#'
#' data("multiclass_dataset")
#' trial <- conditionalMC(multiclass_dataset[,c(1,2,3)], 100)
#'
#' plot(multiclass_dataset[multiclass_dataset$Sample_1 == "A1",1],
#'      multiclass_dataset[multiclass_dataset$Sample_1 == "A1",2],
#'      asp = 1, pch = 1,
#'      xlim = c(min(multiclass_dataset$x), max(multiclass_dataset$x)),
#'      ylim = c(min(multiclass_dataset$y), max(multiclass_dataset$y)))
#' points(multiclass_dataset[multiclass_dataset$Sample_1 == "B1",1],
#'        multiclass_dataset[multiclass_dataset$Sample_1 == "B1",2],
#'        pch = 1, col = "red")
#' points(multiclass_dataset[multiclass_dataset$Sample_1 == "C1",1],
#'        multiclass_dataset[multiclass_dataset$Sample_1 == "C1",2],
#'        pch = 1, col = "blue")
#'
#' points(trial[trial$Sample_1 == "A1", 1],
#'        trial[trial$Sample_1 == "A1", 2],
#'        pch = 19, col = "black", cex = 2)
#' points(trial[trial$Sample_1 == "B1", 1],
#'        trial[trial$Sample_1 == "B1", 2],
#'        pch = 19, col = "red", cex = 2)
#' points(trial[trial$Sample_1 == "C1", 1],
#'        trial[trial$Sample_1 == "C1", 2],
#'        pch = 19, col = "blue", cex = 2)
#'
#' # example 2 ---------------------
#'
#' # generate toy data
#'
#' x1 <- rnorm(100, mean = 2, sd = 1)
#' y1 <- rnorm(100, mean = 2, sd = 1)
#' x2 <- rnorm(100, mean = -2, sd = 1)
#' y2 <- rnorm(100, mean = -2, sd = 1)
#' group_1 <- cbind(x1, y1, rep("Group_1", 100))
#' group_2 <- cbind(x2, y2, rep("Group_2", 100))
#' samples <- data.frame(rbind(group_1, group_2))
#' samples[,1] <- as.numeric(samples[,1])
#' samples[,2] <- as.numeric(samples[,2])
#' samples[,3] <- as.factor(samples[,3])
#'
#' # augment
#'
#' trial2 <- conditionalMC(samples, 100, num_breaks = 20)
#'
#' # plot
#'
#' plot(samples[samples$V3 == "Group_1", 1], samples[samples$V3 == "Group_1", 2],
#'      asp = 1, xlim = c(-5, 5), ylim = c(-5, 5), pch = 1)
#' points(samples[samples$V3 == "Group_2", 1], samples[samples$V3 == "Group_2", 2],
#'        col = "red", pch = 1)
#'
#' points(trial2[trial2$V3 == "Group_1", 1], trial2[trial2$V3 == "Group_1", 2],
#'        pch = 19)
#' points(trial2[trial2$V3 == "Group_2", 1], trial2[trial2$V3 == "Group_2", 2],
#'        pch = 19, col = "red")
#'
#' #
#'
#' @export

conditionalMC <- function(data_set, n_simulations, num_breaks = 0) {

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

  if (ncol(data_set) < 3) {
    stop("Datasets must contain at least 3 variables")
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

  if ((n_factors >= 2) & (n_numbers >= 2)) {
    stop(paste0(
      "This function is only for the augmentation of variables conditioned by either 1 categorcial",
      "or 1 numeric variable. Try multiodalMC instead."
    ))
  }

  if ((n_factors == 1) & (n_numbers == 1)) {
    stop(paste0(
      "This function is only for the augmentation of variables conditioned by either 1 categorcial",
      "or 1 numeric variable. Try bivariate_multimodalMC instead."
    ))
  }

  simulated_data_set <- array(dim = c(0, ncol(data_set)))
  iteration = 0
  pb <- txtProgressBar(min = 0, max = n_simulations,
                       style = 3, width = 100, char = "=")

  if (n_factors < n_numbers) {

    while (dim(simulated_data_set)[1] < n_simulations) {

      single_factor_simulation <- sample_qualitative_value(factor_variables)

      remaining_numeric <- conditioned_filter_quantitative_value(
        numeric_variables,
        factor_variables,
        single_factor_simulation
      )

      #if (length(remaining_numeric) == 0 | is.na(dim(remaining_numeric))) {
      #  single_numeric_simulation <- rep(NA, n_numbers)
      #} else {
      #  single_numeric_simulation <- multivariate_prob_sampling(
      #    remaining_numeric, num_breaks = num_breaks
      #  )
      #}

      if (length(remaining_numeric) == 0) {
        single_numeric_simulation <- rep(NA, n_numbers)
      } else {
        single_numeric_simulation <- multivariate_prob_sampling_v2(
          remaining_numeric, 1, num_breaks = num_breaks
        )
      }

      if (NA %!in% c(single_factor_simulation, single_numeric_simulation)) {
        iteration = iteration + 1
        single_iteration <- c(single_factor_simulation, single_numeric_simulation)
        simulated_data_set <- abind::abind(simulated_data_set, single_iteration,
                                           along = 1)
      }

      if(!is.null(dim(simulated_data_set))) {
        setTxtProgressBar(pb, dim(simulated_data_set)[1])
      }
    }

  } else { # (n_numbers < n_factors)

    while (dim(simulated_data_set)[1] < n_simulations) {

      single_numeric_simulation <- single_value_simulator(numeric_variables)$value

      remaining_categorical <- conditioned_filter_qualitative_value(
        numeric_variables,
        factor_variables,
        single_value_simulator(numeric_variables)
      )

      if (length(remaining_categorical) == 0) {
        single_factor_simulation <- rep(NA, n_factors)
      } else {
        single_factor_simulation <- sample_qualitative_value(remaining_categorical)
      }

      if (NA %!in% c(single_factor_simulation, single_numeric_simulation)) {
        iteration = iteration + 1
        single_iteration <- c(single_factor_simulation, single_numeric_simulation)
        simulated_data_set <- abind::abind(simulated_data_set, single_iteration,
                                           along = 1)
      }

      if(!is.null(dim(simulated_data_set))) {
        setTxtProgressBar(pb, dim(simulated_data_set)[1])
      }
    }

  }

  augmented_data <- data.frame(empty = rep(as.double(NA), n_simulations))

  if (n_factors < n_numbers) {
    column_factor <- as.factor(simulated_data_set[,1])
    augmented_data <- cbind(augmented_data, column_factor)
    for (nums in 2:dim(simulated_data_set)[2]) {
      column_factor <- as.numeric(simulated_data_set[,nums])
      augmented_data <- cbind(augmented_data, column_factor)
    }
  } else { # (n_numbers < n_factors)
    for (facs in 1:dim(factor_variables)[2]) {
      column_factor <- as.factor(simulated_data_set[,facs])
      augmented_data <- cbind(augmented_data, column_factor)
    }
    for (nums in (dim(factor_variables)[2] + 1):(dim(factor_variables)[2] + 1)) {
      column_number <- as.numeric(simulated_data_set[,nums])
      augmented_data <- cbind(augmented_data, column_number)
    }
  }

  augmented_data <- augmented_data[,-1]

  augmented_data_return <- as.data.frame(augmented_data)
  for (order in 1:ncol(augmented_data)) {
    augmented_data_return[, original_column_orders[order]] <- augmented_data[,order]
  }
  colnames(augmented_data_return) <- original_col_names

  cat("\n\n")

  return(augmented_data_return)

}
