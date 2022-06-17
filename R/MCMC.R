
#' @include MCMC_Trace.R
NULL

#' Markov Chain Monte Carlo for Multivariate Numeric Distributions
#'
#' The present function can be used to simulate numerical data using the
#' Metropolis-Hastings 1st Order Variant of the
#' Markov Chain Monte Carlo algorithm.
#'
#' @param data_set An matrix or array of numeric data
#' @param n_iterations The number of iterations to compute of the Markov Chain
#' @param step_size The step size of the Markov Chain. 0.5 is the default step size.
#' A vector of step sizes may also be provided, defining the step size for each dimension separately.
#' @param num_breaks An optional parameter that can be used to fine-tune the calculation
#' of the probability distribution from which to sample from.
#'
#' @section Details:
#' This function uses the the Metropolis-Hastings 1st Order Variant (Metropolis et al., 1953;
#' Hastings, 1970) of the Markov Chain Monte Carlo (MCMC) algorithm (Gamerman and Lopes, 2006),
#' to simulate numerical data from a multivariate distribution.
#'
#' MCMCs are stochastic algorithms that inherit from Monte Carlo based methods
#' by generating random values from a given probability distribution, while using
#' Markovian chains to explore the distribution via "random walks". The Metropolis-Hastings
#' algorithm is used to define an acceptance criterion for a proposed state of the Markov
#' Chain.
#'
#' State probabilities for this funciton are calculated using the chain rule in statistics.
#' Probability density functions for each dimension's distribution were additionally smoothed
#' using polynomial interpolation.
#'
#' @section Bibliography:
#'
#' Gamerman, D.; Lopes, H.R. (2006) Markov Chain Monte Carlo, Boca Raton: Chapman & Hall.
#'
#' Hastings, W.K. (1970) Monte-Carlo sampling methods using Markov Chains and their applications, Biometrika, 57:97-109
#'
#' Metropolis, N.; Rosenbluth, A.W.; Rosenbluth, M.N.; Teller, A.H.; Teller, E. (1953)
#' Equation of state calculations by fast computing machines, Journal of Chemical Physics.
#' 21:1087-1092
#'
#' @return Markov Chain Monte Carlo Trace S4 Object
#'
#' @seealso \code{\link{MCMC_Trace}}, \code{\link{burn_in}}, \code{\link{sample_from_trace}},
#' \code{\link{effective_sample_size}}, \code{\link{TOST}}, \code{\link{autocorrelation}},
#' \code{\link{trace_plot}}
#'
#' @author Lloyd A. Courtenay
#'
#' @examples
#'
#' data("single_moon")
#'
#' augmentation_trial <- MCMC(as.matrix(single_moon), 10000, step_size = 1)
#'
#' trace <- sample_from_trace(augmentation_trial, 100)
#'
#' plot(single_moon, asp = 1, pch = 19)
#' points(trace[,1], trace[,2], col = "red", pch  = 19)
#'
#' #
#'
#' @export

MCMC <- function(
  data_set,
  n_iterations,
  step_size = 0.5,
  num_breaks = 0
) {

  if ("matrix" %!in% class(data_set) & "array" %!in% class(data_set)) {
    stop("Input data must be a matrix or array")
  }

  if(missing_data_type(data_set, type = "Numeric")) {
    stop("This function is for numeric matrices only.")
  }

  if(!is.numeric(n_iterations) | n_iterations <= 0 | n_iterations %% 1 != 0) {
    stop("The number of simulations must be a non-negative integer, e.g. 2")
  }

  if(!is.numeric(step_size)) {
    stop("Step Size must be a non-negative number.")
  }

  if(length(step_size) != 1) {
    if (length(step_size) != dim(data_set)[2]) {
      stop(paste0(
        "If a vector of step sizes is provided, then the vector must be the same length",
        " as the number of dimensions to be sampled from."
      ))
    }
    if (TRUE %in% (step_size <= 0)) {
      stop("Step sizes must be non-negative numbers")
    }
  } else {
    if (step_size <= 0) {
      stop("Step size must be a non-negative number")
    }
  }

  if (num_breaks != 0) {
    if(!is.numeric(num_breaks) | num_breaks <= 0 | num_breaks %% 1 != 0) {
      stop("The number of breaks must be a non-negative integer, e.g. 2. If unsure, leave the default value of 0")
    }
  }

  #if (ncol(data_set) == 1) {
  #  stop("This function is for multivariate distributions only. Try univariateMC instead.")
  #}

  cat("\n")

  trace <- array(numeric(),
                 dim = c(0, dim(data_set)[2]))

  cat("Initialising " , dim(data_set)[2], " Chains", "\n")

  old_prob <- 0

  while (old_prob == 0) {

    old_x <- data_set[sample(nrow(data_set), 1),]
    old_prob <- multivariate_interpolated_probability(old_x, data_set,
                                                num_breaks = num_breaks)
  }

  delta = array(
    numeric(),
    dim = c(n_iterations, dim(data_set)[2])
  )

  if (length(step_size) == 1) {
    for (i in 1:dim(delta)[2]) {
      delta[,i] <- rnorm(dim(delta)[1], mean = 0, sd = step_size)
    }
  } else {
    for (i in 1:dim(delta)[2]) {
      delta[,i] <- rnorm(dim(delta)[1], mean = 0, sd = step_size[i])
    }
  }

  cat("Begin sampling...", "\n")

  pb <- txtProgressBar(min = 0, max = n_iterations,
                       style = 3, width = 100, char = "=")

  for (i in 1:n_iterations) {

    new_x = old_x + delta[i,]
    new_prob = multivariate_interpolated_probability(new_x, data_set,
                                                      num_breaks = num_breaks)
    acceptance = new_prob / old_prob
    if ( acceptance >= runif(1, 0, 1)) {
      trace <- abind::abind(trace, new_x, along = 1)
      old_x = new_x
      old_prob = new_prob
    } else {
      trace <- abind::abind(trace, old_x, along = 1)
    }
    setTxtProgressBar(pb, i)
  }
  cat("\n")
  cat("Algorithm Finished", "\n")

  trace_obj <- new("MCMC_Trace")
  trace_obj <- setTrace(trace_obj, trace)

  return(trace_obj)

}
