
#' Monte Carlo Simulations for Categorical Variables
#'
#' The present function can be used to simulate categorical (qualitative) data
#'
#' @param data_set An object of type data.frame containing either uni- or multivariate
#' categorical variables
#' @param n_simulations An integer defining the number of simulations
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
#' data_set <- multiclass_dataset[,c(3, 4)]
#'
#' # augment data
#'
#' augmented_data <- categoricalMC(data_set, 1000)
#'
#' # visualise data
#'
#' par(mfrow = c(2,2))
#' barplot(table(data_set$Sample_1),
#'         main = "Original Sample 1")
#' barplot(table(augmented_data$Sample_1),
#'         main = "Augmented Sample 1")
#' barplot(table(data_set$Sample_2),
#'         main = "Original Sample 2")
#' barplot(table(augmented_data$Sample_2),
#'         main = "Augmented Sample 2")
#' par(mfrow = c(1,1))
#'
#' # statistical proportion test
#'
#' n_original <- nrow(data_set)
#' n_augmented <- nrow(augmented_data)
#'
#' # example A1
#'
#' sample_1_original_A1 <- table(data_set$Sample_1)[1]
#' sample_1_augmented_A1 <- table(augmented_data$Sample_1)[1]
#'
#' prop.test(x = c(sample_1_original_A1,
#'                 sample_1_augmented_A1),
#'           n = c(n_original, n_augmented))
#'
#' # example B1
#'
#' sample_1_original_B1 <- table(data_set$Sample_1)[2]
#' sample_1_augmented_B1 <- table(augmented_data$Sample_1)[2]
#'
#' prop.test(x = c(sample_1_original_B1,
#'                 sample_1_augmented_B1),
#'           n = c(n_original, n_augmented))
#'
#' # example C1
#'
#' sample_1_original_C1 <- table(data_set$Sample_1)[3]
#' sample_1_augmented_C1 <- table(augmented_data$Sample_1)[3]
#'
#' prop.test(x = c(sample_1_original_C1,
#'                 sample_1_augmented_C1),
#'           n = c(n_original, n_augmented))
#'
#' # example A2
#'
#' sample_2_original_A2 <- table(data_set$Sample_2)[1]
#' sample_2_augmented_A2 <- table(augmented_data$Sample_2)[1]
#'
#' prop.test(x = c(sample_2_original_A2,
#'                 sample_2_augmented_A2),
#'           n = c(n_original, n_augmented))
#'
#' # example B2
#'
#' sample_2_original_B2 <- table(data_set$Sample_2)[2]
#' sample_2_augmented_B2 <- table(augmented_data$Sample_2)[2]
#'
#' prop.test(x = c(sample_2_original_B2,
#'                 sample_2_augmented_B2),
#'           n = c(n_original, n_augmented))
#'
#' #
#'
#' @export

categoricalMC <- function(data_set, n_simulations) {

  if(!is.numeric(n_simulations) | n_simulations <= 0 | n_simulations %% 1 != 0) {
    stop("The number of simulations must be a non-negative integer, e.g. 2")
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

  if(!missing_data_type(data_set, type = "Numeric")) {
    stop("The dataset contains invalid input for this function")
  }

  pb <- txtProgressBar(min = 0, max = n_simulations,
                       style = 3, width = 100, char = "=")

  if (ncol(data_set) == 1) {

    simulated_variables <- c()

    for (sim in 1:n_simulations) {
      simulated_variables <- c(simulated_variables,
                               sample_qualitative_value(data_set))
      setTxtProgressBar(pb, sim)
    }

    simulated_variables <- as.factor(simulated_variables)

  } else {

    simulated_variables <- array(dim = c(0, ncol(data_set)))

    for (sim in 1:n_simulations) {

      simulated_variables <- abind::abind(
        simulated_variables,
        sample_qualitative_value(data_set),
        along = 1
      )
      setTxtProgressBar(pb, sim)

    }

    simulated_variables <- as.data.frame(simulated_variables,
                                         stringsAsFactors = T)

  }


  colnames(simulated_variables) <- colnames(data_set)

  cat("\n\n")

  return(simulated_variables)

}
