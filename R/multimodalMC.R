
#' Monte Carlo Simulation of Multivariate Multimodal Data
#'
#' The present function can be used to simulate a mixture of both categorical (qualitative)
#' and numeric (quantitative) variables using a Monte Carlo approach
#'
#' @param data_set An object of type data.frame containing both qualitative
#' and quantitative variables that the user wishes to augment
#' @param n_simulations An integer defining the number of simulations
#' @param num_breaks An optional parameter that can be used to fine-tune the calculation
#' of the probability distribution from which to sample from.
#'
#' @section Details:
#' This function uses the statistical chain rule to calculate multivariate probability
#' vectors from which to sample from, firstly calculating the probability of
#' qualitative variables, and then sampling numeric values associated with these
#' qualitative variables.
#'
#' @return simulated data
#'
#' @author Lloyd A. Courtenay
#'
#' @examples
#'
#' data(multiclass_dataset)
#'
#' trial_simulation <- multimodalMC(multiclass_dataset, n_simulations = 100)
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
#' points(trial_simulation[trial_simulation$Sample_1 == "A1", 1],
#'        trial_simulation[trial_simulation$Sample_1 == "A1", 2],
#'        pch = 19, col = "black", cex = 2)
#' points(trial_simulation[trial_simulation$Sample_1 == "B1", 1],
#'        trial_simulation[trial_simulation$Sample_1 == "B1", 2],
#'        pch = 19, col = "red", cex = 2)
#' points(trial_simulation[trial_simulation$Sample_1 == "C1", 1],
#'        trial_simulation[trial_simulation$Sample_1 == "C1", 2],
#'        pch = 19, col = "blue", cex = 2)
#'
#' plot(multiclass_dataset[multiclass_dataset$Sample_2 == "A2",1],
#'      multiclass_dataset[multiclass_dataset$Sample_2 == "A2",2],
#'      asp = 1, pch = 1,
#'      xlim = c(min(multiclass_dataset$x), max(multiclass_dataset$x)),
#'      ylim = c(min(multiclass_dataset$y), max(multiclass_dataset$y)))
#' points(multiclass_dataset[multiclass_dataset$Sample_2 == "B2",1],
#'        multiclass_dataset[multiclass_dataset$Sample_2 == "B2",2],
#'        pch = 1, col = "red")
#'
#' points(trial_simulation[trial_simulation$Sample_2 == "A2", 1],
#'        trial_simulation[trial_simulation$Sample_2 == "A2", 2],
#'        pch = 19, col = "black", cex = 2)
#' points(trial_simulation[trial_simulation$Sample_2 == "B2", 1],
#'        trial_simulation[trial_simulation$Sample_2 == "B2", 2],
#'        pch = 19, col = "red", cex = 2)
#'
#' #
#'
#' @export

multimodalMC <- function (data_set, n_simulations, num_breaks = 0) {

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
    stop("Datasets must contain at least 2 qualitative and 2 quantitative variable")
  }

  original_col_names <- colnames(data_set)
  original_column_orders <- c()

  simulated_data_set <- array(dim = c(0, ncol(data_set)))

  factor_list <- find_factor_columns(data_set, original_column_orders)
  number_list <- find_numeric_columns(data_set, factor_list$original_column_orders)

  original_column_orders <- number_list$original_column_orders
  factor_variables <- factor_list$factor_variables
  numeric_variables <- number_list$numeric_variables

  if (factor_list$n_variables < 2 | number_list$n_variables < 2) {
    stop("dataset must contain at least 2 factor and 2 numeric columns")
  }

  rm(number_list, factor_list)

  iteration = 0

  pb <- txtProgressBar(min = 0, max = n_simulations,
                       style = 3, width = 100, char = "=")

  while (dim(simulated_data_set)[1] < n_simulations) {

    remaining_factor <- factor_variables
    sampling_dims <- ncol(factor_variables)
    single_simulation_fac <- array(dim = c(ncol(factor_variables)))

    for (i in 1:ncol(factor_variables)) {

      variable_to_sample_from <- remaining_factor[,i]

      prob_table <- as.data.frame(prop.table(table(variable_to_sample_from)) * 1000)

      vector_prob <- c()
      for (q_var in 1:nrow(prob_table)) {
        vector_prob <- c(
          vector_prob,
          rep(as.character(prob_table[q_var,1]), prob_table[q_var,2])
        )
      }

      vector_prob <- vector_prob[sample(length(vector_prob),
                                        size = length(vector_prob),
                                        replace = FALSE)]

      sim_var <- vector_prob[sample(1:length(vector_prob), size = 1)]

      single_simulation_fac[i] <- sim_var

      sampling_dims <- sampling_dims - 1

      if (sampling_dims != 0) {
        remaining_factor <- remaining_factor[
          remaining_factor[, i] == sim_var,
        ]
      } else {
        break
      }

    }

    sampling_dims <- ncol(factor_variables) # should it be numeric_variables?
    single_simulation_num <- array(dim = c(ncol(numeric_variables)))

    filter_data <- cbind(factor_variables, numeric_variables)
    for (filt_col in 1:length(single_simulation_fac)) {
      filter_data <- filter_data[
        filter_data[,filt_col] == single_simulation_fac[filt_col],
      ]
    }
    remaining_numeric <- filter_data[
      , (ncol(factor_variables) + 1):(ncol(factor_variables) + ncol(numeric_variables))
    ]

    for (i in 1:ncol(numeric_variables)) {

      prob_table <- dist_assoc_prob_table(remaining_numeric[,i], num_breaks = num_breaks)
      prob_table$frequencies <- (prob_table$frequencies * 1000)
      vector_prob <- c(); for (j in 1:nrow(prob_table)) {
        counts <- rep(j, prob_table[j, 3])
        vector_prob <- c(vector_prob, counts)
      }
      vector_prob <- vector_prob[sample(length(vector_prob),
                                        size = length(vector_prob),
                                        replace = FALSE)]
      target_bin <- vector_prob[sample(1:length(vector_prob), size = 1)]
      simulated_value <- runif(1,
                               min = prob_table[target_bin,]$lower,
                               max = prob_table[target_bin,]$upper)
      single_simulation_num[i] <- simulated_value

      sampling_dims <- sampling_dims - 1

      if (length(sampling_dims != 0)) {
        remaining_numeric <- remaining_numeric[
          (remaining_numeric[,i] > prob_table[target_bin,]$lower) &
            (remaining_numeric[,i] < prob_table[target_bin,]$upper),
        ]
      }

      if (dim(remaining_numeric)[1] == 0) {
        break
      }

    }

    if (NA %!in% c(single_simulation_fac, single_simulation_num)) {
      iteration = iteration + 1
      single_iteration <- c(single_simulation_fac, single_simulation_num)
      simulated_data_set <- abind::abind(simulated_data_set, single_iteration,
                                         along = 1)
    }

    if(!is.null(dim(simulated_data_set))) {
      setTxtProgressBar(pb, dim(simulated_data_set)[1])
    }

  }

  augmented_data <- data.frame(empty = rep(as.double(NA), n_simulations))
  for (facs in 1:ncol(factor_variables)) {
    column_factor <- as.factor(simulated_data_set[,facs])
    augmented_data <- cbind(augmented_data, column_factor)
  }
  for (nums in (ncol(factor_variables) + 1):(ncol(factor_variables) + ncol(numeric_variables))) {
    column_number <- as.numeric(simulated_data_set[,nums])
    augmented_data <- cbind(augmented_data, column_number)
  }
  augmented_data <- augmented_data[,-1]

  augmented_data_return <- augmented_data
  for (order in 1:ncol(augmented_data)) {
    augmented_data_return[, original_column_orders[order]] <- augmented_data[,order]
  }
  colnames(augmented_data_return) <- original_col_names

  cat("\n\n")

  return(augmented_data_return)

}

