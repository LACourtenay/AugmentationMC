#' @name AugmentationMC
#' @docType package
#' @aliases AugmentationMC
#' @title Monte Carlo based tools and algorithms for Data Augmentation
#' @author Lloyd A. Courtenay
#' @description This library provides a set of functions for the implementation of Monte Carlo
#' based algorithms for the augmentation of both qualitative and quantitative data.
#' This includes a series of Monte Carlo simulation functions, as well as an implementation
#' of a Metropolis-Hastings variant of the Markov Chain Monte Carlo algorithm, with
#' a series of functions for the evaluation of chains.
#'
#' @import stats
#' @import utils
#' @import graphics
#' @import grDevices
#' @import utils
#' @import methods
#'
NULL

#' Toy Circles Dataset
#'
#' A toy dataset containing x and y coordinates for a multivariate inverse-gaussian distribution
#'
#' @name circles
#' @docType data
#' @author Lloyd A. Courtenay
#' @keywords datasets
#' @format A data frame with 500 rows and 2 variables (x, y)
NULL

#' Toy Moons Dataset
#'
#' A toy dataset containing x and y coordinates of two multivariate moons
#'
#' @name moons
#' @docType data
#' @author Lloyd A. Courtenay
#' @keywords datasets
#' @format A data frame with 1000 rows and 2 variables (x, y)
NULL

#' Toy Single Moon Dataset
#'
#' A toy dataset containing x and y coordinates of a single multivariate moon
#'
#' @name single_moon
#' @docType data
#' @author Lloyd A. Courtenay
#' @keywords datasets
#' @format A data frame with 500 rows and 2 variables (x, y)
NULL

#' Toy Qualitative and Quantitative Moons Dataset
#'
#' A toy dataset containing x and y coordinates of two multivariate moons with an additional factor column
#'
#' @name moons_qual_quant
#' @docType data
#' @author Lloyd A. Courtenay
#' @keywords datasets
#' @format A data frame with 1000 rows and 3 variables; 2 quantitative (x, y), 1 qualitative
#' (category)
NULL

#' Toy Complex and Imbalanced Qualitative and Quantitative Dataset
#'
#' A toy dataset containing x and y coordinates of two multivariate distributions with
#' 2 seperate categorical classes
#'
#' @name multiclass_dataset
#' @docType data
#' @author Lloyd A. Courtenay
#' @keywords datasets
#' @format A data frame with 210 rows and 4 variables; 2 quantitative (x, y), 2 qualitative (Sample_1, Sample_2).
#' Sample_1 consistes of three classes, A1, B1, and C1, of which the latter is imbalanced (number of individuals = 10)
NULL

