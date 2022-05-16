
# Traditional TOST

tTOST <- function(dist_1, dist_2) {

  n1 <- length(dist_1)
  mu1 <- mean(dist_1)
  var1 <- var(dist_1)

  n2 <- length(dist_2)
  mu2 <- mean(dist_2)
  var2 <- var(dist_2)

  estimate <- c(mu1, mu2)

  std_error_dist_1 <- sqrt(var1 / n1)
  std_error_dist_2 <- sqrt(var2 / n2)
  std_error <- sqrt(std_error_dist_1^2 + std_error_dist_2^2)

  df <- std_error^4 / (std_error_dist_1^4 / (n1 - 1) + std_error_dist_2^4 / (n2 - 1))

  test_stat <- (mu1 - mu2) / std_error

  difference <- mu1 - mu2

  conf_interval <- c(test_stat - qt(0.95, df),
                     test_stat + qt(0.95, df)) * std_error

  p_value <- as.numeric(pt(
    (1 - abs(difference)) / std_error,
    df,
    lower.tail = FALSE
  ))

  tTOST_return <- list(
    robust = FALSE,
    p = p_value,
    abs_difference = abs(difference),
    test_statistic = test_stat,
    df = df,
    confidence_interval = conf_interval
  )

  class(tTOST_return) <- "TOST_obj"

  print_TOST(tTOST_return)

}

# Robust TOST

rTOST <- function(dist_1, dist_2) {

  alpha <- 1 - 0.9
  epsilon = 0.31

  m1 <- length(dist_1) - 2 * floor(0.2 * length(dist_1))
  m2 <- length(dist_2) - 2 * floor(0.2 * length(dist_2))
  d1 <- (length(dist_1) - 1) * gamma_winsorized_variance(dist_1) / (m1 * (m1 - 1))
  d2 <- (length(dist_2) - 1) * gamma_winsorized_variance(dist_2) / (m2 * (m2 - 1))
  yuen_df <- (d1 + d2)^2 / (d1^2 / (m1 - 1) + d2^2 / (m2 - 1))

  se <- sqrt(d1 + d2)

  difference <- mean(dist_1, trim = 0.2) - mean(dist_2, trim = 0.2)

  test_stat <- difference / se

  criteria <- qt(1 - alpha/2, yuen_df)
  lower <- difference - criteria * se
  upper <- difference + criteria * se
  yuen_estimate <- c(mean(dist_1, trim = 0.2), mean(dist_2, trim = 0.2))
  yuen_confidence_interval <- c(lower, upper)

  mean_difference <- ifelse(
    length(yuen_estimate) == 1,
    as.numeric(yuen_estimate),
    as.numeric(yuen_estimate[1] - yuen_estimate[2])
  )

  se_difference <- as.numeric(
    (yuen_confidence_interval[2] - yuen_confidence_interval[1]) /
      qt(0.95, df = yuen_df)
  ) / 2

  se_difference <- as.numeric((yuen_confidence_interval[2] - yuen_confidence_interval[1])/
                                qt(1 - 0.05, df = yuen_df))/2

  rtost_test_stat <- (epsilon - abs(mean_difference)) / se_difference

  p_value <- as.numeric(pt(
    rtost_test_stat, yuen_df,
    lower.tail = FALSE
  ))

  rTOST_return <- list(
    robust = TRUE,
    p = p_value,
    abs_difference = abs(mean_difference),
    test_statistic = rtost_test_stat,
    df = yuen_df,
    confidence_interval = yuen_confidence_interval
  )

  class(rTOST_return) <- "TOST_obj"

  print_TOST(rTOST_return)

}

# Misc Functions

gamma_winsorized_variance <- function(dist) {

  dist_ordered <- sort(dist)
  n <- length(dist)
  bottom_indx <- floor(0.2 * n) + 1
  top_indx <- length(dist) - bottom_indx + 1
  bottom_x <- dist_ordered[bottom_indx]
  top_x <- dist_ordered[top_indx]
  y <- ifelse(dist_ordered <= bottom_x, bottom_x, dist_ordered)
  y <- ifelse(y >= top_x, top_x, y)
  return(var(y))

}

# Print Function for Results

print_TOST <- function(TOST_obj) {

  robust <- TOST_obj$robust

  p <- TOST_obj$p
  d <- TOST_obj$abs_difference
  t <- TOST_obj$test_statistic
  df <- TOST_obj$df
  ci <- TOST_obj$confidence_interval
  false_positive_risk_mean <- FPR(p, priors = 0.5) * 100
  false_positive_risk_upper <- FPR(p, priors = 0.3) * 100
  false_positive_risk_lower <- FPR(p, priors = 0.7) * 100

  p <- round_scientific_number(p)

  if (robust == TRUE) {

    cat("\n Robust Two-One Sided Test of Equivalence (rTOST)\n")
    cat(" - Yuen-Welch Trimmed Test Statistic\n\n")
    null_hyp <- "Null Hypothesis (pH0): Trimmed Sample Means are Different"

  } else {

    cat("\n Two-One Sided Test of Equivalence (TOST)\n")
    cat(" - Welch Test Statistic\n\n")
    null_hyp <- "Null Hypothesis (pH0): Sample Means are Different"

  }

  cat(paste0(
    "\tt = ", round(t, 4), ", df = ", round(df, 2), ", p-value = ", p, "\n\n"
  ))
  cat(paste0(null_hyp, "\n\n"))
  cat(paste0(
    "95 Percent Confidence Intervals:\n\n",
    "\t", round(ci[1], 4), "\t", round(ci[2], 4), "\n\n"
  ))
  cat(paste0(
    "Normalised Absolute Difference:\n\n\t",
    round(d, 4), "\n\n"
  ))
  if (as.numeric(p) < 0.3681) {
    false_positive_risk_mean <- round_scientific_number(false_positive_risk_mean)
    false_positive_risk_upper <- round_scientific_number(false_positive_risk_upper)
    false_positive_risk_lower <- round_scientific_number(false_positive_risk_lower)
    cat(paste0(
      "Probability of Type I Statistical Error\n",
      "When accepting the Alternative Hypothesis:",
      "\n\n\t",
      false_positive_risk_mean, " +/- [ ", false_positive_risk_lower,
      " , ", false_positive_risk_upper, " ] %"
    ))
  }

  cat("\n\n")

}

#

scale_values <- function(x) {
  return((x - min(x)) / (max(x) - min(x)))
}

#
