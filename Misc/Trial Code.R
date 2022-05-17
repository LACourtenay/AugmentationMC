
# Trial -------

mc_qual_quant <- function (data_set, n_values, num_breaks = 0) {

  `%!in%` = Negate(`%in%`)

  original_col_names <- colnames(data_set)
  original_column_orders <- c()

  simulated_data_set <- array(dim = c(0, ncol(data_set)))

  factor_variables <- data.frame(empty = rep(as.double(NA), nrow(data_set)))
  numeric_variables <- data.frame(empty = rep(as.double(NA), nrow(data_set)))
  for (i in 1:ncol(data_set)){
    if (is.factor(data_set[,i])) {
      factor_variables <- cbind(factor_variables, data_set[,i])
      original_column_orders <- c(original_column_orders, i)
    }
  }
  for (i in 1:ncol(data_set)){
    if (is.numeric(data_set[,i])) {
      numeric_variables <- cbind(numeric_variables, data_set[,i])
      original_column_orders <- c(original_column_orders, i)
    }
  }

  factor_variables <- factor_variables[,-1]
  numeric_variables <- numeric_variables[,-1]

  iteration = 0

  pb <- txtProgressBar(min = 0, max = n_values,
                       style = 3, width = 100, char = "=")

  while (dim(simulated_data_set)[1] < n_values) {

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

    sampling_dims <- ncol(factor_variables)
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

    #if (complete.cases(c(single_simulation_fac, single_simulation_num))) {
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

  augmented_data <- data.frame(empty = rep(as.double(NA), n_values))
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

  return(augmented_data_return)
}

multivariate_prob_sampling2 <- function(distribution, n_values, num_breaks = 0) {

  if (dim(distribution)[1] == 0) {
    stop("No data to augment")
  }

  simulated_distribution <- array(numeric(),
                                  dim = c(0, dim(distribution)[2]))
  available_dimensions <- seq(from = 1, to = dim(distribution)[2], by = 1)

  pb <- txtProgressBar(min = 0, max = n_values,
                       style = 3, width = 100, char = "=")

  while (dim(simulated_distribution)[1] < n_values) {
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

    if(!is.null(dim(simulated_distribution))) {
      setTxtProgressBar(pb, dim(simulated_distribution)[1])
    }

  }

  cat("\n")

  return(simulated_distribution)

}

#

# work ---------

library(AugmentationMC)

data("moons_qual_quant")

left_moon <- moons[moons$category == "Left",]
right_moon <- moons[moons$category == "Right",]

# visualise original data

#trial <- mc_qualitative_quantitative2(moons_qual_quant, n_simulations = 100, num_breaks = 20)
trial <- conditionalMC(moons_qual_quant, 100)

left_moon_sim <- trial[trial$category == "Left",]
right_moon_sim <- trial[trial$category == "Right",]

#

plot(left_moon$x, left_moon$y,
     asp = 1, pch = 19,
     xlim = c(min(moons$x), max(moons$x)),
     ylim = c(min(moons$y), max(moons$y)))
points(right_moon$x, right_moon$y, pch = 19, col = "red")

points(left_moon_sim$x, left_moon_sim$y, pch = 19, col = "black")
points(right_moon_sim$x, right_moon_sim$y, pch = 19, col = "red")

#

a <- read.table(
  "C:\\Users\\Lloyd\\Desktop\\TIDOP\\Julia Tesis\\Augmentation and 3D Model Simulation\\Trial_QL_QN_Data.csv",
  sep = ";", head = T)[,2:7]
a$Taxon <- as.factor(a$Taxon)
a$Sex <- as.factor(a$Sex)
a$Region <- as.factor(a$Region)

trial1 <- mc_qualitative_quantitative2(a, n_simulations = 1000)
trial2 <- mc_qual_quant(a, n_values = 1000)

summary(a$Taxon)
summary(trial1$Taxon)
summary(trial2$Taxon)

mc_qualitative_quantitative2(a, n_simulations = 1000)

# toy example -----------------------

library(AugmentationMC)

data("multiclass_dataset")

sample_data_2 <- multiclass_dataset
sample_data_2 <- cbind(sample_data_2, sample_data_2[sample(nrow(sample_data_2),
                                                           size = nrow(sample_data_2),
                                                           replace = FALSE),1])

trial <- multimodalMC(multiclass_dataset, 100)

par(mfrow = c(1,2))

plot(
  multiclass_dataset[multiclass_dataset$Sample_1 == "A1",1],
  multiclass_dataset[multiclass_dataset$Sample_1 == "A1",2],
  asp = 1, pch = 1,
  xlim = c(min(multiclass_dataset$x), max(multiclass_dataset$x)),
  ylim = c(min(multiclass_dataset$y), max(multiclass_dataset$y)),
  xlab = "x", ylab = "y"
)
points(
  multiclass_dataset[multiclass_dataset$Sample_1 == "B1",1],
  multiclass_dataset[multiclass_dataset$Sample_1 == "B1",2],
  pch = 1, col = "red"
)
points(
  multiclass_dataset[multiclass_dataset$Sample_1 == "C1",1],
  multiclass_dataset[multiclass_dataset$Sample_1 == "C1",2],
  pch = 1, col = "blue"
)

plot(
  multiclass_dataset[multiclass_dataset$Sample_2 == "A2",1],
  multiclass_dataset[multiclass_dataset$Sample_2 == "A2",2], asp = 1, pch = 19,
  xlim = c(min(multiclass_dataset$x), max(multiclass_dataset$x)),
  ylim = c(min(multiclass_dataset$y), max(multiclass_dataset$y)),
  xlab = "x", ylab = "y"
)
points(
  multiclass_dataset[multiclass_dataset$Sample_2 == "B2",1],
  multiclass_dataset[multiclass_dataset$Sample_2 == "B2",2], pch = 19, col = "red"
)

par(mfrow = c(1,1))

points(trial[trial$Sample_1 == "A1", 1],
       trial[trial$Sample_1 == "A1", 2],
       pch = 19, col = "black", cex = 2)
points(trial[trial$Sample_1 == "B1", 1],
       trial[trial$Sample_1 == "B1", 2],
       pch = 19, col = "red", cex = 2)
points(trial[trial$Sample_1 == "C1", 1],
       trial[trial$Sample_1 == "C1", 2],
       pch = 19, col = "blue", cex = 2)

trial <- multimodalMC(multiclass_dataset, n_simulations = 100)
plot(
  multiclass_dataset[multiclass_dataset$Sample_2 == "A2",1],
  multiclass_dataset[multiclass_dataset$Sample_2 == "A2",2], asp = 1, pch = 19,
  xlim = c(min(multiclass_dataset$x), max(multiclass_dataset$x)),
  ylim = c(min(multiclass_dataset$y), max(multiclass_dataset$y))
)
points(
  multiclass_dataset[multiclass_dataset$Sample_2 == "B2",1],
  multiclass_dataset[multiclass_dataset$Sample_2 == "B2",2], pch = 19, col = "red"
)

points(trial[trial$Sample_2 == "A2", 1],
       trial[trial$Sample_2 == "A2", 2], pch = 19, col = "black", cex = 2)
points(trial[trial$Sample_2 == "B2", 1],
       trial[trial$Sample_2 == "B2", 2], pch = 19, col = "red", cex = 2)

#

# Conditional MC -------------------

#trial <- conditionalMC(multiclass_dataset[,c(1,2,3)], 100)
#trial <- multimodalMC(multiclass_dataset, n_simulations = 100)

data("multiclass_dataset")
trial <- conditionalMC(multiclass_dataset[,c(1,2,4)], 100)

plot(multiclass_dataset[multiclass_dataset$Sample_2 == "A2",1],
     multiclass_dataset[multiclass_dataset$Sample_2 == "A2",2], asp = 1, pch = 1,
     xlim = c(min(multiclass_dataset$x), max(multiclass_dataset$x)),
     ylim = c(min(multiclass_dataset$y), max(multiclass_dataset$y)))
points(multiclass_dataset[multiclass_dataset$Sample_2 == "B2",1],
       multiclass_dataset[multiclass_dataset$Sample_2 == "B2",2], pch = 1, col = "red")

points(trial[trial$Sample_2 == "A2", 1],
       trial[trial$Sample_2 == "A2", 2], pch = 19, col = "black", cex = 2)
points(trial[trial$Sample_2 == "B2", 1],
       trial[trial$Sample_2 == "B2", 2], pch = 19, col = "red", cex = 2)

# prop.test()

# univariate conditional MC ----------

data("multiclass_dataset")


library(sm)

trial <- conditionalMC(multiclass_dataset[,c(2,3,4)], 1000, num_breaks = 20)
sm.density.compare(multiclass_dataset[,2], multiclass_dataset[,3])
sm.density.compare(trial[,1], trial[,2])

#

trial <- multimodalMC(multiclass_dataset, n_simulations = 1000)
sm.density.compare(multiclass_dataset[,2], multiclass_dataset[,4])
sm.density.compare(trial[,2], trial[,4])

x1 <- rnorm(100, mean = 2, sd = 1)
y1 <- rnorm(100, mean = 2, sd = 1)
x2 <- rnorm(100, mean = -2, sd = 1)
y2 <- rnorm(100, mean = -2, sd = 1)

plot(x1, y1, xlim = c(-5, 5), ylim = c(-5, 5), asp = 1, pch = 19)
points(x2, y2, col = "red", pch = 19)

group_1 <- cbind(x1, y1, rep("Group_1", 100))
group_2 <- cbind(x2, y2, rep("Group_2", 100))
samples <- data.frame(rbind(group_1, group_2))
samples[,1] <- as.numeric(samples[,1])
samples[,2] <- as.numeric(samples[,2])
samples[,3] <- as.factor(samples[,3])

trial <- conditionalMC(samples, 100, num_breaks = 20)

plot(samples[samples$V3 == "Group_1", 1], samples[samples$V3 == "Group_1", 2],
     asp = 1, xlim = c(-5, 5), ylim = c(-5, 5), pch = 1)
points(samples[samples$V3 == "Group_2", 1], samples[samples$V3 == "Group_2", 2],
       col = "red", pch = 1)

points(trial[trial$V3 == "Group_1", 1], trial[trial$V3 == "Group_1", 2],
       pch = 19)
points(trial[trial$V3 == "Group_2", 1], trial[trial$V3 == "Group_2", 2],
       pch = 19, col = "red")

#

# bivariate MC -----------

bivariate_trial <- data.frame(x = c(rnorm(100, 0, 2), rnorm(100, -4, 3)),
                              sample = as.factor(c(rep("Group_1", 100),
                                                   rep("Group_2", 100))))

augmented_trial <- bivariate_multimodalMC(bivariate_trial, 100)

par(mfrow = c(2,2))
hist(bivariate_trial[bivariate_trial$sample == "Group_1",1],
     main = "Original Data Group 1",
     xlab = "Group 1")
hist(c(bivariate_trial[bivariate_trial$sample == "Group_1",1],
       augmented_trial[augmented_trial$sample == "Group_1",1]),
     main = "Augmented Data Group 1",
     xlab = "Group 1")
hist(bivariate_trial[bivariate_trial$sample == "Group_2",1],
     main = "Original Data Group 2",
     xlab = "Group 2")
hist(c(bivariate_trial[bivariate_trial$sample == "Group_2",1],
       augmented_trial[augmented_trial$sample == "Group_2",1]),
     main = "Augmented Data Group 2",
     xlab = "Group 2")
par(mfrow = c(1,1))

#

# Synthetic Data base ----------------

x <- rnorm(100, mean = 0, sd = 1)
y <- rnorm(100, mean = 0, sd = 2)
xy <- cbind(x, y)
alpha <- -atan((xy[1,2]-tail(xy,1)[,2])/(xy[1,1]-tail(xy,1)[,1]))
rotm <- matrix(c(cos(alpha),sin(alpha),-sin(alpha),cos(alpha)),ncol=2)
M2 <- t(rotm %*% (
  t(xy)-c(xy[1,1],xy[1,2])
)+c(xy[1,1],xy[1,2]))

plot(M2[,1], M2[,2], asp = 1, pch = 19)

x2 <- rnorm(100, mean = 2, sd = 1)
y2 <- rnorm(100, mean = 2, sd = 2)
xy2 <- cbind(x2, y2)
alpha <- -atan((xy2[1,2]-tail(xy2,1)[,2])/(xy2[1,1]-tail(xy2,1)[,1]))
rotm <- matrix(c(cos(alpha),sin(alpha),-sin(alpha),cos(alpha)),ncol=2)
M3 <- t(rotm %*% (
  t(xy2)-c(xy2[1,1],xy2[1,2])
)+c(xy2[1,1],xy2[1,2]))

plot(M2[,1], M2[,2], asp = 1, pch = 19, xlim = c(-5, 10), ylim = c(-5, 5))
points(M3[,1], M3[,2], asp = 1, pch = 19, col = "red")

coords <- rbind(M2, M3)
labels <- c(rep("Sample_1", 100),rep("Sample_2", 100))
labels <- as.factor(labels)
data_set <- data.frame(coords, Sample = labels)

x3 <- rnorm(10, mean = 6, sd = 0.5)
y3 <- rnorm(10, mean = 4, sd = 0.75)
xy3 <- cbind(x3, y3)

plot(M2[,1], M2[,2], asp = 1, pch = 19, xlim = c(-5, 10), ylim = c(-5, 5))
points(M3[,1], M3[,2], asp = 1, pch = 19, col = "red")
points(xy3[,1], xy3[,2], asp = 1, pch = 19, col = "blue")

coords <- rbind(M2, M3, xy3)
labels <- c(rep("A1", 100),rep("B1", 100), rep("C1", 10))
labels <- as.factor(labels)
data_set <- data.frame(coords, Sample = labels)

plot(M2[,1], M2[,2], asp = 1, pch = 19, xlim = c(-5, 10), ylim = c(-5, 5))
points(M3[,1], M3[,2], asp = 1, pch = 19, col = "red")
points(xy3[,1], xy3[,2], asp = 1, pch = 19, col = "blue")
plot(function(x) {(x * 0.5) - 1}, col = "red", from = -10, to = 15, add = TRUE)

data_set2 <- as_tibble(data_set) %>%
  mutate(class = ((x3 * 0.5) - 1),
         Sample_2 = ifelse(y3 < class, "A2", "B2")) %>% select(x3, y3, Sample, Sample_2) %>%

data_set2$Sample_2 <- as.factor(data_set2$Sample_2)
qualitative_quantitative_data <- data_set2
multiclass_dataset <- qualitative_quantitative_data
colnames(multiclass_dataset) <- c("x", "y", "Sample_1", "Sample_2")

usethis::use_data(multiclass_dataset, overwrite = TRUE)

#

# Categorical_MC ------------

data(multiclass_dataset)
data_set <- multiclass_dataset[,c(3, 4)]

# augment data

augmented_data <- categoricalMC(data_set, 1000)

# visualise data

par(mfrow = c(2,2))
barplot(table(data_set$Sample_1),
        main = "Original Sample 1")
barplot(table(augmented_data$Sample_1),
        main = "Augmented Sample 1")
barplot(table(data_set$Sample_2),
        main = "Original Sample 2")
barplot(table(augmented_data$Sample_2),
        main = "Augmented Sample 2")
par(mfrow = c(1,1))

# statistical proportion test

n_original <- nrow(data_set)
n_augmented <- nrow(augmented_data)

# example A1

sample_1_original_A1 <- table(data_set$Sample_1)[1]
sample_1_augmented_A1 <- table(augmented_data$Sample_1)[1]

prop.test(x = c(sample_1_original_A1,
                sample_1_augmented_A1),
          n = c(n_original, n_augmented))

# example B1

sample_1_original_B1 <- table(data_set$Sample_1)[2]
sample_1_augmented_B1 <- table(augmented_data$Sample_1)[2]

prop.test(x = c(sample_1_original_B1,
                sample_1_augmented_B1),
          n = c(n_original, n_augmented))

# example C1

sample_1_original_C1 <- table(data_set$Sample_1)[3]
sample_1_augmented_C1 <- table(augmented_data$Sample_1)[3]

prop.test(x = c(sample_1_original_C1,
                sample_1_augmented_C1),
          n = c(n_original, n_augmented))

# example A2

sample_2_original_A2 <- table(data_set$Sample_2)[1]
sample_2_augmented_A2 <- table(augmented_data$Sample_2)[1]

prop.test(x = c(sample_2_original_A2,
                sample_2_augmented_A2),
          n = c(n_original, n_augmented))

# example B2

sample_2_original_B2 <- table(data_set$Sample_2)[2]
sample_2_augmented_B2 <- table(augmented_data$Sample_2)[2]

prop.test(x = c(sample_2_original_B2,
                sample_2_augmented_B2),
          n = c(n_original, n_augmented))

#

# numeric_MC ----------------

data(multiclass_dataset)
data_set <- multiclass_dataset[,c(1, 2)]
data_set <- as.matrix(data_set)

# augment data

augmented_data <- numericMC(data_set, 100)

plot(data_set, asp = 1, pch = 19, xlim = c(-10, 15), ylim = c(-8, 8))
points(augmented_data[,1], augmented_data[,2], col = "red", pch = 19)

#

# MCMC ---------------------

library(AugmentationMC)
data("single_moon")

trial <- multivariateMCMC(as.matrix(single_moon), 20000, step_size = 1)
trial <- burn_in(trial, 0.1)

effective_sample_size(trial)

trace_plot(trial, dimension = 1)
trace_plot(trial, dimension = 2)

autocorrelation(trial, dimension = 1)
autocorrelation(trial, dimension = 2)

augmented_data <- sample_from_trace(trial, 100)
plot(single_moon, asp = 1, pch = 19)
points(augmented_data[,1], augmented_data[,2], col = "red", pch  = 19)

TOST(augmented_data[,1], single_moon[,1])
TOST(augmented_data[,2], single_moon[,2])


