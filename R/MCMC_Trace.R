
# MCMC_Trace S4 Class ---------------------------

#'
#' An S4 Class for Markov Chain Monte Carlo traces
#'
#' @slot  Trace A numerical matrix containing the MCMC trace
#' @slot  burn_in A boolean value indicating whether the burn in period of the
#' trace has been removed
#'
#' @export

setClass(
  "MCMC_Trace",
  representation(
    Trace = "matrix",
    burn_in = "logical"
  ),
  prototype = list(
    Trace = NULL,
    burn_in = FALSE
  )
)
setValidity("MCMC_Trace", function(object){
  if(!is.numeric(object@Trace)) {
    "MCMC_Trace@Trace must be a numeric matrix"
  }
  if(!is.logical(object@burn_in)) {
    "MCMC_Trace@burn_in must be a boolean value"
  }
})

#' Define a MCMC_Trace object
#'
#' The is function is to define an MCMC_Trace object.
#'
#' @param traceObject an MCMC_Trace object
#' @param trace a numeric vector, matrix, or array
#'
#' @author Lloyd A. Courtenay
#'
#' @export

setGeneric(
  name = "setTrace",
  def = function(traceObject, trace) standardGeneric("setTrace")
)
setMethod(
  f = "setTrace",
  signature = "MCMC_Trace",
  definition = function(traceObject, trace) {
    traceObject@Trace <- trace
    return(traceObject)
  }
)

# Burn-In ---------------------------

#' Burn-In Sample
#'
#' The present function can be used to remove the burn in from a MCMC Trace Object
#'
#' @param traceObject an MCMC_Trace object
#' @param burn_in a numeric value between 0 and 1 defining the percentage of burn-in
#' to be removed. Default is 0.25
#'
#' @return An MCMC_Trace object excluding the burn-in period
#'
#' @author Lloyd A. Courtenay
#'
#' @examples
#'
#' data("single_moon")
#'
#' augmentation_trial <- MCMC(as.matrix(single_moon), 10000, step_size = 1)
#'
#' augmentation_trial <- burn_in(augmentation_trial, 0.25)
#'
#' #
#'
#' @export

setGeneric(
  name = "burn_in",
  def = function(traceObject,
                 burn_in = 0.25) standardGeneric("burn_in")
)
setMethod(
  f = "burn_in",
  signature = "MCMC_Trace",
  definition = function(traceObject,
                        burn_in = 0.25) {

    if ("MCMC_Trace" %!in% class(traceObject)) {
      stop("This function only accepts objects of class MCMC_Trace. See ??MCMC_Trace.")
    }

    if(!is.numeric(burn_in)) {
      stop("Burn in must be a numeric value between 0 and 1")
    }

    if(burn_in < 0 | burn_in > 1) {
      stop("Burn in must be a numeric value between 0 and 1")
    }

    trace <- as.matrix(traceObject@Trace)
    length_trace <- dim(trace)[1]


    trace <- trace[-c(1:(length_trace * burn_in)),]

    traceObject@Trace <- as.matrix(trace)
    traceObject@burn_in <- TRUE

    return(traceObject)

  }
)

# Sample from Trace ---------------------------

#' Sample from an MCMC Trace
#'
#' The present function can be used to sample values from an MCMC trace
#'
#' @param traceObject an MCMC_Trace object
#' @param n_simulations An integer defining the number of simulations to extract
#' from the MCMC_Trace
#'
#' @details This function first removes duplicate values from the MCMC's trace,
#' followed by a random sampling of the trace to produce a final matrix of new
#' simulated values.
#'
#' @return Matrix of simulated values sampled from the trace
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

setGeneric(
  name = "sample_from_trace",
  def = function(traceObject, n_simulations) standardGeneric("sample_from_trace")
)
setMethod(
  f = "sample_from_trace",
  signature = "MCMC_Trace",
  definition = function(traceObject, n_simulations) {

    if ("MCMC_Trace" %!in% class(traceObject)) {
      stop("This function only accepts objects of class MCMC_Trace. See ??MCMC_Trace.")
    }

    if(!is.numeric(n_simulations) | n_simulations <= 0 | n_simulations %% 1 != 0) {
      stop("The number of simulations must be a non-negative integer, e.g. 2")
    }

    if (traceObject@burn_in == FALSE) {
      warning(paste0(
        "This is a raw MCMC trace, optimal results are normally obtained",
        " taking into account a burn-in period. See help(burn_in)"
      ))
    }

    if (dim(traceObject@Trace)[2] == 1) {

      #res_sample <- traceObject@Trace[
      #  (dim(traceObject@Trace)[1] / 2):dim(traceObject@Trace)[1]
      #]

      res_sample <- traceObject@Trace

      res_sample <- res_sample[!duplicated(res_sample)]

      if (length(res_sample) < n_simulations) {
        warning(paste(
          "Trace not long enough to produce", n_simulations, "samples, only",
          length(res_sample), "have been obtained...",
          "\nThe user should consider;",
          "\n(A): Choosing a smaller number of samples",
          "\n(B): Increasing mcmc chain length",
          "\n(C): Improving mcmc parameters"
        ))

      } else {
        res_sample <- sample(res_sample, size = n_simulations, replace = FALSE)
      }

      return(res_sample)

    } else {
      #res_sample <- traceObject@Trace[
      #  (dim(traceObject@Trace)[1] / 2):dim(traceObject@Trace)[1],
      #]
      res_sample <- traceObject@Trace

      res_sample <- res_sample[
        !duplicated(res_sample[,1:dim(traceObject@Trace)[2]]),
      ]

      if (dim(res_sample)[1] < n_simulations) {
        warning(paste(
          "Trace not long enough to produce", n_simulations, "samples, only",
          dim(res_sample)[1], "have been obtained...",
          "\nThe user should consider;",
          "\n(A): Choosing a smaller number of samples",
          "\n(B): Increasing mcmc chain length",
          "\n(C): Improving mcmc parameters"
        ))

      } else {

        res_sample <- res_sample[sample(nrow(res_sample),
                                        size = n_simulations, replace = FALSE),]

      }

      return(res_sample)

    }

  }
)

# Effective Sample Size ---------------------------

#' Calculate Effective Sample Size of MCMC Trace
#'
#' Calculation of the effective sample size for mean estimation,
#' adjusted for autocorrelation
#'
#' @param traceObject an MCMC_Trace object
#'
#' @details This function calculates the effective sample size for mean estimation
#' from an MCMC trace by estimating the spectral density of the trace at frequency 0.
#' Said estimation is performed by fitting an autoregressive model to each dimension
#' of the trace, framed as a time-series. ESS thus quantifies the number of samples
#' with the same standard error as the original sample, so as to compute the reliability
#' of the posterior distribution drawn from the chain.
#'
#' @return Values of the Effective Sample Size (ESS) for each dimension of the trace
#'
#' @section Note:
#'
#' The smaller the ESS, the poorer the posterior distribution of the MCMC is considered
#' to be at estimating the original parameters of the distribution.
#'
#' A definition of a small ESS is hard to define, with not all authors coming to a
#' consensus. Likewise, ESS calculations for univariate chains tend to present
#' much higher values than individual ESS values for each dimension of multiple chains.
#' For a univariate chain, values above 1000 are recommendable,
#' while for multivariate chains, values above 100 (or preferably 200) are recommendable
#' per dimension.
#'
#' For an alternative ESS calculation, we strongly recommend the functions available
#' in the mcmcse R library
#'
#' @section Bibliography:
#'
#' Vats, D.; Flegal, J.M.; Jones, G.L. (2015) Multivariate output analysis for Markov
#' Chain Monte Carlo. arXiv: 1512.07713
#'
#' Elvira, V.; Martino, L.; Robert, C. (2018) Rethinking the Effective Sample Size.
#' arXiv: 1809.04129v1
#'
#' @author Lloyd A. Courtenay
#'
#' @examples
#'
#' data("single_moon")
#'
#' trial_1 <- MCMC(as.matrix(single_moon[,1]), 10000, step_size = 1)
#' trial_2 <- MCMC(as.matrix(single_moon[,1]), 10000, step_size = 0.5)
#' trial_3 <- MCMC(as.matrix(single_moon[,1]), 10000, step_size = 0.1)
#'
#'
#' # example of convergence
#' effective_sample_size(trial_1)
#'
#' # example of poor convergence
#' effective_sample_size(trial_2)
#'
#' # example of no convergence
#' effective_sample_size(trial_3)
#'
#' #
#'
#' @export

setGeneric(
  name = "effective_sample_size",
  def = function(traceObject) standardGeneric("effective_sample_size")
)
setMethod(
  f = "effective_sample_size",
  signature = "MCMC_Trace",
  definition = function(traceObject) {

    if ("MCMC_Trace" %!in% class(traceObject)) {
      stop("This function only accepts objects of class MCMC_Trace. See ??MCMC_Trace.")
    }

    if (traceObject@burn_in == FALSE) {
      warning(paste0(
        "This is a raw MCMC trace, optimal results are normally obtained",
        " taking into account a burn-in period. See help(burn_in)"
      ))
    }

    trace <- traceObject@Trace

    dim_variance <- numeric(dim(trace)[2])
    lag <- 1:dim(trace)[1]

    for (dimension in 1:dim(trace)[2]) {

      linear_model <- lm(trace[,dimension] ~ lag)
      if (identical(all.equal(sd(residuals(linear_model)), 0), TRUE)) {
        dim_variance[dimension] <- 0
      } else {
        autoregressive_model <- ar(trace[,dimension], aic = TRUE)
        dim_variance[dimension] <- autoregressive_model$var.pred / (
          1 - sum(autoregressive_model$ar)
        )^2
      }

    }

    ess <- ifelse(
      dim_variance == 0, 0, dim(trace)[1] * apply(trace, 2, var) / dim_variance
    )

    return(ess)

  }
)

# TOST ---------------------------

#' Two-One Sided Test (TOST) for Equivalence between Trace and Original Data
#'
#' Calculation of the statistical equivalence between the augmented data
#' sampled by the MCMC algorithm and the original data used to train the MCMC
#'
#' @param trace either an MCMC_Trace object, a matrix of numeric values,
#' or a vector of numeric values.
#' @param original_data the original matrix
#' @param dimension an integer defining which dimension of the distribution should
#' be used for TOST calculations (in the case of multivariate data)
#' @param robust a boolean value defining whether the robust calculation of TOST
#' should be performed. Default is set to TRUE.
#' @param plot a boolean value defining whether a density plot should be produced
#' comparing the original and the augmented data
#'
#' @details This function can be used to compare two distributions to check for
#' their equivalence, and the magnitude of said equivalence. As recommended by Courtenay
#' and Gonzalez-Aguilera (2020), a useful way of testing the quality of synthetically
#' produced data is by using a Two-One Sided Test for Equivalence. From this perspective
#' the null hypothesis (pH_{0}) is that the two samples are different. Low p-Values
#' thus represent a high level of equivalence between the synthetically produced data
#' and the original data.
#'
#' The way TOST is calculated depends on sample normality. If distributions are
#' Gaussian in nature, then a traditional TOST test using Welch's Test Statistic
#' can be used to calculate the magnitude of differences and similarities. Nevertheless,
#' if the distributions are not Gaussian, then a non-parametric test statistic should be
#' used. In this case we use the Yuen-Welch Trimmed test statistic, here referred to as
#' robust TOST (or rTOST)
#'
#' The calculated absolute value is a normalised value, where both distribution 1
#' and distribution 2 are scaled to contain values between 0 and 1.
#'
#' This test also returns (when p-Values are found to be under 0.3681) the probability
#' that the p-Value is capturing a false positive (Type I Statistical Error). The reported
#' probability uses a 0.5 prior, however is also accompanied by [0.7, 0.3] prior confidence
#' intervals
#'
#' @section Bibliography:
#'
#' Yuen, K.K.; Dixon, W.J. (1973) The Approximate Behaviour and Performance of the
#' Two-Sample Trimmed T. Biometrika. 60:369-374
#'
#' Yuen, K.K. (1974) The Two-Sample Trimmed T for Unequal Population Variances. Biometrika.
#' 61:165-170
#'
#' Hauk, D.W.W.; Anderson, S. (1984) A New Statistical Procedure for Testing Equivalence
#' in Two-Group Comparative Bioavailability Trials. Journal of Pharmacokinetics and
#' Biopharmaceutics. 12:83-91
#'
#' Schurimann, D.L. (1987) A comparison of the two one-sided test procedure and the
#' power approach for assessing the equivalence of average biovariability.
#' Journal of Pharmacokinetics and Biopharmaceutics. 15:657-680
#'
#' Cohen, J. (1988) Statistical Power Analysis for Behavioural Sciences. New York: Routledge
#'
#' Courtenay, L.A.; Gonzalez-Aguilera, D. (2020) Geometric Morphometric Data Augmentation using
#' Generative Computational Learning Algorithms, Applied Sciences.10:9133. DOI:
#' 10.3390/app10249133
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
#' TOST(trace, single_moon, dimension = 1, robust = TRUE)
#'
#' #
#'
#' @rdname TOST
#' @export
#'
setGeneric(
  name = "TOST",
  def = function(trace,
                 original_data,
                 dimension = NULL,
                 robust = TRUE,
                 plot = TRUE) standardGeneric("TOST")
)

#' @rdname TOST
setMethod(
  f = "TOST",
  signature = "MCMC_Trace",
  definition = function(trace,
                        original_data,
                        dimension = NULL,
                        robust = TRUE,
                        plot = TRUE) {

    if (missing(original_data) | missing(trace)) {
      stop("TOST can only be used to evaluate original_data with a simulated trace.")
    }

    trace <- as.matrix(trace@Trace)
    dim_1 <- dim(trace)[1]
    dim_2 <- dim(trace)[2]

    if (dim_2 != dim(as.matrix(original_data))[2]) {
      stop("Both the original data and the trace must have the same size")
    }

    if (missing_data_type(trace, "Numeric")) {
      stop("TOST can only be calculated with numeric variables")
    }

    if (is.null(dimension)) {
      if (dim_2 != 1) {
        stop("For non-vectorial data, the user must specify which variable they wish to compare")
      }
    } else {
      if (dimension > dim_2) {
        stop(paste0(
          "The user has specified a dimension or variable that does not exist.\n",
          "The data provided has ", dim_2, " dimensions/variables, but the user has",
          " specified ", dimension
        ))
      }
    }

    if (!is.logical(robust)) {
      stop("The robust argument only accepts boolean values")
    }

    if (!is.logical(plot)) {
      stop("The plot argument only accepts boolean values")
    }

    if (dim_2 >= 2) {
      trace <- trace[(dim_1 / 2):dim_1,]
      trace <- trace[!duplicated(trace),]
      trace <- trace[,dimension]
      original_data <- original_data[,dimension]
    } else {
      trace <- trace[(dim_1 / 2):dim_1]
      trace <- trace[!duplicated(trace)]
    }

    scaled_values <- scale_values(c(trace, original_data))
    scaled_trace <- scaled_values[1:length(trace)]
    scaled_original <- scaled_values[
      length(trace) + 1:length(scaled_values) - length(trace)
    ]

    if (robust == TRUE) {

      rTOST(scaled_trace, scaled_original)

    } else {

      tTOST(scaled_trace, scaled_original)

    }

    if (plot == TRUE) {

      trace_dense <- density(trace)
      original_dense <- density(original_data)

      ylim = c(
        0, max(c(trace_dense$y, original_dense$y))
      )
      xlim = c(
        min(c(trace_dense$x, original_dense$x)),
        max(c(trace_dense$x, original_dense$x))
      )

      par(mar = c(5.1, 5, 4.1, 2.))

      plot(original_dense$x, original_dense$y,
           xlim = xlim, ylim = ylim,
           type = "l", lwd = 2,
           xlab = "", ylab = "",
           xaxt = "none", yaxt = "none")
      lines(trace_dense$x, trace_dense$y,
            lwd = 2, col =" red", lty = 2)

      mtext(side = 1, line = 3, "Distribution", cex = 1.25, font = 2)
      mtext(side = 2, line = 3, "Density", cex = 1.25, font = 2)
      axis(1, round(seq(xlim[1], xlim[2], length.out = 5), 2), font = 1, cex.axis = 1)
      axis(2, round(seq(ylim[1], ylim[2], length.out = 5), 2), font = 1, cex.axis = 1)

      legend("topright", legend = c("Original", "Simulated"),
             inset = 0.015,
             col = c("Black", "Red"),
             lty = 1:2, lwd = 2)

      par(mar = c(5.1, 4.1, 4.1, 2.1))

    }

  }
)

#' @rdname TOST
setMethod(
  f = "TOST",
  signature = "matrix",
  definition = function(trace,
                        original_data,
                        dimension = NULL,
                        robust = TRUE,
                        plot = TRUE) {

    if (missing(original_data) | missing(trace)) {
      stop("TOST can only be used to evaluate original_data with a simulated trace.")
    }

    if (missing_data_type(trace, "Numeric")) {
      stop("TOST can only be calculated with numeric variables")
    }

    trace <- as.matrix(trace)
    dim_1 <- dim(trace)[1]
    dim_2 <- dim(trace)[2]

    if (dim_2 != dim(as.matrix(original_data))[2]) {
      stop("Both the original data and the trace must have the same size")
    }

    if (is.null(dimension)) {
      if (dim_2 != 1) {
        stop("For non-vectorial data, the user must specify which variable they wish to compare")
      }
    } else {
      if (dimension > dim_2) {
        stop(paste0(
          "The user has specified a dimension or variable that does not exist.\n",
          "The data provided has ", dim_2, " dimensions/variables, but the user has",
          " specified ", dimension
        ))
      }
    }

    if (!is.logical(robust)) {
      stop("The robust argument only accepts boolean values")
    }

    if (!is.logical(plot)) {
      stop("The plot argument only accepts boolean values")
    }

    if (dim_2 >= 2) {
      trace <- trace[!duplicated(trace),]
      trace <- trace[,dimension]
      original_data <- original_data[,dimension]
    } else {
      trace <- trace[(dim_1 / 2):dim_1]
      trace <- trace[!duplicated(trace)]
    }

    scaled_values <- scale_values(c(trace, original_data))
    scaled_trace <- scaled_values[1:length(trace)]
    scaled_original <- scaled_values[
      length(trace) + 1:length(scaled_values) - length(trace)
    ]

    if (robust == TRUE) {

      rTOST(scaled_trace, scaled_original)

    } else {

      tTOST(scaled_trace, scaled_original)

    }

    if (plot == TRUE) {

      trace_dense <- density(trace)
      original_dense <- density(original_data)

      ylim = c(
        0, max(c(trace_dense$y, original_dense$y))
      )
      xlim = c(
        min(c(trace_dense$x, original_dense$x)),
        max(c(trace_dense$x, original_dense$x))
      )

      par(mar = c(5.1, 5, 4.1, 2.))

      plot(original_dense$x, original_dense$y,
           xlim = xlim, ylim = ylim,
           type = "l", lwd = 2,
           xlab = "", ylab = "",
           xaxt = "none", yaxt = "none")
      lines(trace_dense$x, trace_dense$y,
            lwd = 2, col =" red", lty = 2)

      mtext(side = 1, line = 3, "Distribution", cex = 1.25, font = 2)
      mtext(side = 2, line = 3, "Density", cex = 1.25, font = 2)
      axis(1, round(seq(xlim[1], xlim[2], length.out = 5), 2), font = 1, cex.axis = 1)
      axis(2, round(seq(ylim[1], ylim[2], length.out = 5), 2), font = 1, cex.axis = 1)

      legend("topright", legend = c("Original", "Simulated"),
             inset = 0.015,
             col = c("Black", "Red"),
             lty = 1:2, lwd = 2)

      par(mar = c(5.1, 4.1, 4.1, 2.1))

    }

  }
)

#' @rdname TOST
setMethod(
  f = "TOST",
  signature = "numeric",
  definition = function(trace,
                        original_data,
                        robust = TRUE,
                        plot = TRUE) {

    if (missing(original_data) | missing(trace)) {
      stop("TOST can only be used to evaluate original_data with a simulated trace.")
    }

    if (!is.logical(robust)) {
      stop("The robust argument only accepts boolean values")
    }

    if (!is.logical(plot)) {
      stop("The plot argument only accepts boolean values")
    }

    trace <- as.matrix(trace)
    dim_1 <- dim(trace)[1]
    dim_2 <- dim(trace)[2]
    trace <- trace[!duplicated(trace)]

    if (dim_2 != dim(as.matrix(original_data))[2]) {
      stop("Both the original data and the trace must have the same size")
    }

    scaled_values <- scale_values(c(trace, original_data))
    scaled_trace <- scaled_values[1:length(trace)]
    scaled_original <- scaled_values[
      length(trace) + 1:length(scaled_values) - length(trace)
    ]

    if (robust == TRUE) {

      rTOST(scaled_trace, scaled_original)

    } else {

      tTOST(scaled_trace, scaled_original)

    }

    if (plot == TRUE) {

      trace_dense <- density(trace)
      original_dense <- density(original_data)

      ylim = c(
        0, max(c(trace_dense$y, original_dense$y))
      )
      xlim = c(
        min(c(trace_dense$x, original_dense$x)),
        max(c(trace_dense$x, original_dense$x))
      )

      par(mar = c(5.1, 5, 4.1, 2.))

      plot(original_dense$x, original_dense$y,
           xlim = xlim, ylim = ylim,
           type = "l", lwd = 2,
           xlab = "", ylab = "",
           xaxt = "none", yaxt = "none")
      lines(trace_dense$x, trace_dense$y,
            lwd = 2, col =" red", lty = 2)

      mtext(side = 1, line = 3, "Distribution", cex = 1.25, font = 2)
      mtext(side = 2, line = 3, "Density", cex = 1.25, font = 2)
      axis(1, round(seq(xlim[1], xlim[2], length.out = 5), 2), font = 1, cex.axis = 1)
      axis(2, round(seq(ylim[1], ylim[2], length.out = 5), 2), font = 1, cex.axis = 1)

      legend("topright", legend = c("Original", "Simulated"),
             inset = 0.015,
             col = c("Black", "Red"),
             lty = 1:2, lwd = 2)

      par(mar = c(5.1, 4.1, 4.1, 2.1))

    }

  }
)

# Autocorrelation ---------------------------

#' Autocorrelation Function (ACF) Plots
#'
#' Computation and plotting of a Autocorrelation Function plot (or Correlogram)
#'
#' @param trace an MCMC_Trace object
#' @param dimension an integer defining which dimension of the distribution should
#' be visualised (if the trace is multivariate)
#' @param main A string defining the title for the plot. If no main is defined then
#' "Correlogram" will be used as default.
#'
#' @details ACF Plots can be used in combination with trace plots to visualise
#' and measure the serial correlation of a chain. Effective use of MCMC
#' should see the curve approximate 0 quickly, while chains that take a while to
#' reach 0 indicate the step size of MCMC may not be ideal for this particular
#' case study. Likewise, another common issue in MCMCs may be that the size of the trace
#' is too small, and increasing the number of iterations will improve results.
#' From this perspective, the user can assess how well and how
#' quickly the MCMC chain has converged on the actual distribution.
#'
#' @author Lloyd A. Courtenay
#'
#' @examples
#'
#' data("single_moon")
#'
#' trial_1 <- MCMC(as.matrix(single_moon[,1]), 10000, step_size = 1)
#' trial_2 <- MCMC(as.matrix(single_moon[,1]), 10000, step_size = 0.5)
#' trial_3 <- MCMC(as.matrix(single_moon[,1]), 10000, step_size = 0.1)
#'
#'
#' # example of convergence
#' autocorrelation(trial_1, main = "Step Size = 1")
#'
#' # example of poor convergence
#' autocorrelation(trial_2, main = "Step Size = 0.5")
#'
#' # example of no convergence
#' autocorrelation(trial_3, main = "Step Size = 0.1")
#'
#' #
#'
#' @export
#'

setGeneric(
  name = "autocorrelation",
  def = function(trace,
                 dimension = NULL,
                 main = NULL) standardGeneric("autocorrelation")
)
setMethod(
  f = "autocorrelation",
  signature = "MCMC_Trace",
  definition = function(trace,
                        dimension = NULL,
                        main = NULL) {

    if ("MCMC_Trace" %!in% class(trace)) {
      stop("This function only accepts objects of class MCMC_Trace. See ??MCMC_Trace.")
    }

    trace <- trace@Trace

    n_variables <- dim(trace)[2]

    if (n_variables >= 2) {
      if (is.null(dimension)) {
        stop("The user must specify a dimension of the trace to visualise")
      } else {
        trace <- trace[,dimension]
      }
    }

    autocorrelation_features <- acf(
      trace, las = 1, plot = FALSE
    )$acf
    autocor_data <- data.frame(
      Lag = seq(1:length(autocorrelation_features)),
      Autocorrelation = autocorrelation_features
    )

    par(mar = c(5.1, 5, 4.1, 2))
    plot(
      autocor_data$Lag, autocor_data$Autocorrelation,
      col = NULL,
      xlab = "", ylab = "",
      ylim = c(-1, 1),
      xaxt = "none", yaxt = "none"
    )
    abline(h = 0, lty = 2, col = "grey")
    lines(autocor_data$Lag, autocor_data$Autocorrelation, lty = 1, col = "red",
          lwd = 2)
    points(autocor_data$Lag, autocor_data$Autocorrelation, pch = 19)
    mtext(side = 1, line = 3, "Lag", cex = 1.2, font = 2)
    mtext(side = 2, line = 3, "Autocorrelation Function", cex = 1.2, font = 2)
    axis(
      1, round(seq(
        0, max(autocor_data$Lag),
        by = 5
      )), font = 1, cex.axis = 1
    )
    axis(2, seq(-1, 1, 0.5), font = 1, cex.axis = 1)
    if (is.null(main)) {
      title("Correlogram", cex.main = 1.75)
    } else {
      title(main, cex.main = 1.75)
    }

    par(mar = c(5.1, 4.1, 4.1, 2.1))

  }
)

# Trace Plots ---------------------------

#' MCMC Trace Plots
#'
#' Plot the Trace of a MCMC
#'
#' @param trace an MCMC_Trace object
#' @param dimension an integer defining which dimension of the distribution should
#' be visualised (if the trace is multivariate)
#' @param main an optional plot title
#'
#' @details Trace plots are a visual and useful guide to ensuring that the MCMC
#' has converged on the original distribution. Traces that do not converge will
#' appear to slowly move around the distribution, while traces that do converge
#' will fixate on where the original distribution is.
#'
#' @author Lloyd A. Courtenay
#'
#' @examples
#'
#' data("single_moon")
#'
#' trial_1 <- MCMC(as.matrix(single_moon[,1]), 10000, step_size = 1)
#' trial_2 <- MCMC(as.matrix(single_moon[,1]), 10000, step_size = 0.5)
#' trial_3 <- MCMC(as.matrix(single_moon[,1]), 10000, step_size = 0.1)
#'
#'
#' # example of convergence
#' trace_plot(trial_1, main = "Step Size = 1")
#'
#' # example of poor convergence
#' trace_plot(trial_2, main = "Step Size = 0.5")
#'
#' # example of no convergence
#' trace_plot(trial_3, main = "Step Size = 0.1")
#'
#' #
#'
#' @export
#'

setGeneric(
  name = "trace_plot",
  def = function(trace,
                 dimension = NULL,
                 main = NULL) standardGeneric("trace_plot")
)
setMethod(
  f = "trace_plot",
  signature = "MCMC_Trace",
  definition = function(trace,
                        dimension = NULL,
                        main = NULL) {

    if ("MCMC_Trace" %!in% class(trace)) {
      stop("This function only accepts objects of class MCMC_Trace. See ??MCMC_Trace.")
    }

    trace <- trace@Trace

    n_variables <- dim(trace)[2]

    if (n_variables >= 2) {
      if (is.null(dimension)) {
        stop("The user must specify a dimension of the trace to visualise")
      } else {
        trace <- trace[,dimension]
      }
    }

    sd_trace <- sd(trace)
    trace_range <- range(trace)

    par(mar = c(5.1, 5, 4.1, 2))

    plot(
      1:length(trace), trace,
      col = NULL,
      xlab = "",
      ylab = "",
      ylim = c(
        trace_range[1] - sd_trace,
        trace_range[2] + sd_trace
      )
    )
    lines(1:length(trace), trace)
    mtext(side = 1, line = 3, "Iterations", cex = 1.2, font = 2)
    mtext(side = 2, line = 3, "Sampled Values", cex = 1.2, font = 2)
    if (is.null(main)) {
      title("Trace Plot", cex.main = 1.75)
    } else {
      title(main, cex.main = 1.75)
    }

    par(mar = c(5.1, 4.1, 4.1, 2.1))

  }
)

