
# Miscellaneous Functions

`%!in%` = Negate(`%in%`)

dist_assoc_prob_table <- function(X, num_breaks = 0) {

  X <- as.numeric(X)

  if (num_breaks == 0) {
    histogram <- hist(X, plot = FALSE)
  } else {
    histogram <- hist(X, plot = FALSE, breaks = num_breaks)
  }
  hist_data <- data.frame(lower = histogram$breaks[-length(histogram$breaks)],
                          upper = histogram$breaks[-1])
  hist_data$frequencies <- histogram$counts / sum(histogram$counts)

  return(hist_data)

}

missing_data_type <- function(data_set, type) {

  if (type == "Factor") {

    n_variables <- ncol(data_set)
    factor_bool <- c()
    for (column in 1:n_variables) {
      factor_bool <- c(factor_bool, is.factor(data_set[,column]))
    }

    if (TRUE %!in% factor_bool) {
      return(TRUE)
    } else {
      return(FALSE)
    }

  } else if (type == "Numeric") {

    n_variables <- ncol(data_set)
    numeric_bool <- c()
    for (column in 1:n_variables) {
      numeric_bool <- c(numeric_bool, is.numeric(data_set[,column]))
    }

    if (TRUE %!in% numeric_bool) {
      return(TRUE)
    } else {
      return(FALSE)
    }

  } else if (type == "Character") {

    n_variables <- ncol(data_set)
    character_bool <- c()
    for (column in 1:n_variables) {
      character_bool <- c(character_bool, is.character(data_set[,column]))
    }

    if (TRUE %!in% character_bool) {
      return(TRUE)
    } else {
      return(FALSE)
    }

  }

}

sample_qualitative_value <- function(factor_variables) {

  factor_variables <- as.matrix(factor_variables)

  remaining_factor <- factor_variables
  sampling_dims <- dim(factor_variables)[2]
  single_simulation_fac <- array(dim = c(dim(factor_variables)[2]))

  for (variable_i in 1:dim(factor_variables)[2]) {

    if (is.null(dim(remaining_factor))) {
      variable_to_sample_from <- remaining_factor[variable_i]
    } else {
      variable_to_sample_from <- remaining_factor[,variable_i]
    }

    prob_table <- as.data.frame(prop.table(table(variable_to_sample_from)) * 1000)

    if (nrow(prob_table) != 0) {
      vector_prob <- c()
      for (q_var in 1:nrow(prob_table)) {
        vector_prob <- c(
          vector_prob,
          rep(as.character(prob_table[q_var,1]), prob_table[q_var,2])
        )
      }
    } else {
      return(NA)
    }

    vector_prob <- vector_prob[sample(length(vector_prob),
                                      size = length(vector_prob),
                                      replace = FALSE)]

    sim_var <- vector_prob[sample(1:length(vector_prob), size = 1)]

    single_simulation_fac[variable_i] <- sim_var

    sampling_dims <- sampling_dims - 1

    if (sampling_dims != 0) {
      remaining_factor <- remaining_factor[
        remaining_factor[, variable_i] == sim_var,
      ]
    } else {
      break
    }

  }

  return(single_simulation_fac)

}

conditioned_filter_quantitative_value2 <- function(numeric_variables,
                                                  factor_variables,
                                                  single_factor_simulation) {

  factor_variables <- as.matrix(factor_variables)
  numeric_variables <- as.matrix(numeric_variables)

  filter_data <- as.matrix(cbind(factor_variables, numeric_variables))
  for (filt_col in 1:length(single_factor_simulation)) {
    if(is.null(dim(filter_data))) {
      filter_data <- filter_data[
        filter_data[filt_col] == single_factor_simulation[filt_col]
      ]
    } else {
      filter_data <- filter_data[
        filter_data[,filt_col] == single_factor_simulation[filt_col],
      ]
    }
  }
  if (is.null(dim(filter_data))) {
    remaining_numeric <- filter_data[
      (dim(factor_variables)[2] + 1):(dim(factor_variables)[2] + dim(numeric_variables)[2])
    ]
  } else {
    remaining_numeric <- filter_data[
      , (dim(factor_variables)[2] + 1):(dim(factor_variables)[2] + dim(numeric_variables)[2])
    ]
  }


  return(remaining_numeric)

}

conditioned_filter_quantitative_value <- function(numeric_variables,
                                                 factor_variables,
                                                 single_factor_simulation) {
  single_factor_simulation <- single_factor_simulation[[1]]

  factor_variables <- as.matrix(factor_variables)
  numeric_variables <- as.matrix(numeric_variables)

  filter_data <- as.matrix(cbind(factor_variables, numeric_variables))

  filter_data <- filter_data[factor_variables == single_factor_simulation[[1]],]

  return(filter_data[,-1])

}

conditioned_filter_qualitative_value <- function(numeric_variables,
                                                 factor_variables,
                                                 single_value_simulator_class) {

  factor_variables <- as.matrix(factor_variables)
  numeric_variables <- as.matrix(numeric_variables)

  filter_data <- as.matrix(cbind(numeric_variables, factor_variables))

  filter_data <- filter_data[
    (filter_data[,1] > single_value_simulator_class$prob_table[
      single_value_simulator_class$target_bin,
    ]$lower) &
      (filter_data[,1] < single_value_simulator_class$prob_table[
        single_value_simulator_class$target_bin,
      ]$upper),
  ]

  if (length(filter_data) == 0) {
    return(NA)
  } else if (length(filter_data) == (ncol(factor_variables) + 1)) {
    return(NA)
  } else {
    return(filter_data[,-1])
  }

}

find_numeric_columns <- function(data_set, original_column_orders) {

  numeric_variables <- data.frame(empty = rep(as.double(NA), nrow(data_set)))

  n_variables <- 0

  for (i in 1:ncol(data_set)){
    if (is.numeric(data_set[,i])) {
      n_variables <- n_variables + 1
      numeric_variables <- cbind(numeric_variables, data_set[,i])
      original_column_orders <- c(original_column_orders, i)
    }
  }

  numeric_variables <- numeric_variables[,-1]

  return(list(numeric_variables = numeric_variables,
              n_variables = n_variables,
              original_column_orders = original_column_orders))

}

find_factor_columns <- function(data_set, original_column_orders) {

  factor_variables <- data.frame(empty = rep(as.double(NA), nrow(data_set)))

  n_variables <- 0

  for (i in 1:ncol(data_set)){
    if (is.factor(data_set[,i])) {
      n_variables <- n_variables + 1
      factor_variables <- cbind(factor_variables, data_set[,i])
      original_column_orders <- c(original_column_orders, i)
    }
  }

  factor_variables <- factor_variables[,-1]

  return(list(factor_variables = factor_variables,
              n_variables = n_variables,
              original_column_orders = original_column_orders))

}

multivariate_prob_sampling <- function(numeric_variables, num_breaks = 0) {

  single_simulation_num <- array(dim = c(dim(numeric_variables)[2]))
  sampling_dims <- dim(numeric_variables)[2]
  reduced_numbers <- numeric_variables
  reduced_numbers <- apply(reduced_numbers, 2, as.numeric)

  for (variable_i in 1:dim(numeric_variables)[2]) {

    simulation <- single_value_simulator(
      numeric_variables[,variable_i], num_breaks = num_breaks
    )

    single_simulation_num[variable_i] <- simulation$value

    sampling_dims <- sampling_dims - 1

    # something makes no sense at all here

    if (length(sampling_dims) != 0) {
      reduced_numbers <- reduced_numbers[
        (reduced_numbers[,variable_i] > simulation$prob_table[simulation$target_bin,]$lower) &
          (reduced_numbers[,variable_i] < simulation$prob_table[simulation$target_bin,]$upper),
      ]
    }

    if (NA %!in% single_simulation_num) {
      break
    } else {
      if (length(reduced_numbers) <= dim(
        numeric_variables
      )[2]) {
        break
      }
    }

  }

  return(single_simulation_num)

}

multivariate_prob_sampling_v2 <- function(distribution, values = 1, num_breaks = 0) {

  distribution <- apply(distribution, 2, as.numeric)

  simulated_distribution <- array(numeric(),
                                  dim = c(0, dim(distribution)[2]))
  available_dimensions <- seq(from = 1, to = dim(distribution)[2], by = 1)

  while (dim(simulated_distribution)[1] < values) {
    sampling_dimensions <- available_dimensions[sample(length(available_dimensions),
                                                       size = length(available_dimensions),
                                                       replace = FALSE)]
    matrix_values <- distribution
    single_simulation <- array(numeric(),
                               dim = c(1, dim(distribution)[2]))
    while (length(sampling_dimensions) != 0) {
      target_dimensions <- sample(sampling_dimensions, 1)
      sampling_dimensions <- sampling_dimensions[-which(target_dimensions == sampling_dimensions)]
      if (!is.null(dim(matrix_values))){
        vector_to_sample_from <- matrix_values[,target_dimensions]
      } else {
        vector_to_sample_from <- matrix_values
      }
      prob_table <- dist_assoc_prob_table(vector_to_sample_from, num_breaks = num_breaks)
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
      single_simulation[1, target_dimensions] <- simulated_value

      if (length(sampling_dimensions) != 0) {
        if (!is.null(dim(matrix_values))) {
          matrix_values <- matrix_values[
            (matrix_values[,target_dimensions] > prob_table[target_bin,]$lower) &
              (matrix_values[,target_dimensions] < prob_table[target_bin,]$upper),
          ]
        } else {
          break
        }
      }
    }
    if(complete.cases(single_simulation)) {
      simulated_distribution <- abind::abind(simulated_distribution, single_simulation, along = 1)
    }

  }

  return(simulated_distribution)

}

single_value_simulator <- function(X, num_breaks = 0) {

  prob_table <- dist_assoc_prob_table(X,
                                      num_breaks = num_breaks)
  prob_table$frequencies <- (prob_table$frequencies * 1000)

  vector_prob <- c(); for (q_var in 1:nrow(prob_table)) {
    counts <- rep(q_var, prob_table[q_var, 3])
    vector_prob <- c(vector_prob, counts)
  }

  vector_prob <- vector_prob[sample(length(vector_prob),
                                    size = length(vector_prob),
                                    replace = FALSE)]

  target_bin <- vector_prob[sample(1:length(vector_prob), size = 1)]

  simulated_value <- runif(1,
                           min = prob_table[target_bin,]$lower,
                           max = prob_table[target_bin,]$upper)

  return(list(value = simulated_value,
              prob_table = prob_table,
              target_bin = target_bin))

}

# doesn't work yet -----------

sample_quantitative_value <- function(numeric_variables, num_breaks = num_breaks) {

  if (is.null(dim(numeric_variables))) {
    return(rep(NA, length(numeric_variables)))
  }

  numeric_variables <- as.matrix(numeric_variables)

  sampling_dims <- dim(numeric_variables)[2] # before was set to factor_variabels but i'm not sure why
  single_simulation_num <- array(dim = c(dim(numeric_variables)[2]))

  for (variable_i in 1:dim(numeric_variables)[2]) {

    prob_table <- dist_assoc_prob_table(numeric_variables[,variable_i], num_breaks = num_breaks)
    prob_table$frequencies <- (prob_table$frequencies * 1000)

    vector_prob <- c(); for (q_var in 1:nrow(prob_table)) {
      counts <- rep(q_var, prob_table[q_var, 3])
      vector_prob <- c(vector_prob, counts)
    }

    vector_prob <- vector_prob[sample(length(vector_prob),
                                      size = length(vector_prob),
                                      replace = FALSE)]

    target_bin <- vector_prob[sample(1:length(vector_prob), size = 1)]

    simulated_value <- runif(1,
                             min = prob_table[target_bin,]$lower,
                             max = prob_table[target_bin,]$upper)

    single_simulation_num[variable_i] <- simulated_value

    sampling_dims <- sampling_dims - 1

    if (length(sampling_dims) != 0) {
      numeric_variables <- numeric_variables[
        (numeric_variables[,variable_i] > prob_table[target_bin,]$lower) &
          (numeric_variables[,variable_i] < prob_table[target_bin,]$upper),
      ]
    }

    if (is.null(dim(numeric_variables))) {
      break
    } else {
      if (dim(numeric_variables)[1] == 0) {
        break
      }
    }


  }

  return(single_simulation_num)

}

# doesn't work yet ------------

sample_quantitative_value2 <- function(numeric_variables, num_breaks = num_breaks) {

  numeric_variables <- as.matrix(numeric_variables)

  sampling_dims <- dim(numeric_variables)[2] # before was set to factor_variabels but i'm not sure why
  single_simulation_num <- array(dim = c(dim(numeric_variables)[2]))

  for (variable_i in 1:dim(numeric_variables)[2]) {

    prob_table <- dist_assoc_prob_table(numeric_variables[,variable_i], num_breaks = num_breaks)
    prob_table$frequencies <- (prob_table$frequencies * 1000)
    vector_prob <- c(); for (q_var in 1:nrow(prob_table)) {
      counts <- rep(q_var, prob_table[q_var, 3])
      vector_prob <- c(vector_prob, counts)
    }
    vector_prob <- vector_prob[sample(length(vector_prob),
                                      size = length(vector_prob),
                                      replace = FALSE)]
    target_bin <- vector_prob[sample(1:length(vector_prob), size = 1)]
    simulated_value <- runif(1,
                             min = prob_table[target_bin,]$lower,
                             max = prob_table[target_bin,]$upper)
    single_simulation_num[variable_i] <- simulated_value

    sampling_dims <- sampling_dims - 1

    if (length(sampling_dims) != 0) {
      numeric_variables <- numeric_variables[
        (numeric_variables[,variable_i] > prob_table[target_bin,]$lower) &
          (numeric_variables[,variable_i] < prob_table[target_bin,]$upper),
      ]
    }

    if (is.na(dim(numeric_variables))) {
      break
    }


  }

  return(single_simulation_num)

}

FPR <- function(p, priors = 0.5){

  bayes_factor_bound <- (1 / (-exp(1) * p * log(p)))

  pH = priors / (1 - priors)

  false_positive_risk <- 1 / (1 + (pH * bayes_factor_bound))

  return(false_positive_risk)

}

round_scientific_number <- function(scientific_number) {

  if (scientific_number < 0.0001) {
    scientific_number = format(scientific_number, scientific = TRUE, digits = 3)
  } else {
    scientific_number = round(scientific_number, 4)
  }

  return(scientific_number)

}
