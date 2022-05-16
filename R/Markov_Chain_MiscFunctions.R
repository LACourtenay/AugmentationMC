
bwmv <- function(X, sbwmv = TRUE) {

  U.i = (X-median(X))/(9*qnorm(.75)*mad(X))
  a.i = ifelse(U.i<=-1|U.i>=1,0,1)
  n = nrow(as.matrix(X))
  nx = sqrt(n)*sqrt(sum((a.i*((X-median(X))^2))*((1-U.i^2)^4)))
  dx = abs(sum(a.i*(1-U.i^2)*(1-5*U.i^2)))
  bwmv_value = (nx/dx)^2

  if (sbwmv == TRUE) {
    return(sqrt(bwmv_value))
  } else {
    return(bwmv_value)
  }
}

probability_association_table <- function(M_dimension, num_breaks = 0) {

  if (num_breaks == 0) {
    histogram <- hist(M_dimension, plot = FALSE)
  } else {
    histogram <- hist(M_dimension, plot = FALSE, breaks = num_breaks)
  }

  hist_data <- data.frame(lower = histogram$breaks[-length(histogram$breaks)],
                          upper = histogram$breaks[-1])
  hist_data$frequencies <- histogram$counts / sum(histogram$counts)

  return(hist_data)

}

grid_values <- function(new_value, hist_data) {

  prob <- 0
  lower_range <- 0
  upper_range <- 0

  if (new_value < min(hist_data$lower) | new_value > max(hist_data$upper)) {
    return(list(
      probability = prob, lower_limit = lower_range, upper_limit = upper_range
    ))
  } else {

    for(i in 1:nrow(hist_data)) {
      if (new_value >= hist_data[i,][["lower"]] & new_value < hist_data[i,][["upper"]]) {
        prob <- hist_data[i,][["frequencies"]]
        lower_range <- hist_data[i,][["lower"]]
        upper_range <- hist_data[i,][["upper"]]
      }
    }

    return(list(
      probability = prob, lower_limit = lower_range, upper_limit = upper_range
    ))

  }

}

multivariate_interpolated_probability <- function(X, matrix_distribution,
                                                   num_breaks = 0) {

  original_dimensions <- dim(matrix_distribution)[2]

  probabilities = c()

  for (i in 1:dim(matrix_distribution)[2]) {

    if (length(dim(matrix_distribution)[1] != 0) == 0) {
      p <- 0
      return(p)
    }

    dist_table <- dist_assoc_prob_table(matrix_distribution[,i], num_breaks = num_breaks)
    prob_params <- grid_values(X[i],
                               dist_table)

    density_function <- splinefun(
      x = seq(min(dist_table$lower), max(dist_table$upper),
              length.out = nrow(dist_table)),
      y = dist_table$frequencies
    )

    prob_point <- density_function(X[i])

    probabilities <- c(probabilities, prob_point)

    matrix_distribution <- matrix_distribution[
      (matrix_distribution[,i] > prob_params$lower_limit) &
        (matrix_distribution[,i] < prob_params$upper_limit),
    ]

    if(length(matrix_distribution) < original_dimensions) {
      p <- 0
      return(p)
    }

  }

  p <- prod(probabilities)

  if(is.na(p)) {
    p <- 0
    return(p)
  }

  return(p)

}
