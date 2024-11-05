# Simulation functions

# Load libraries
library(dplyr)
library(purrr)
library(compositions)
library(MASS)
library(tidyr)
library(Hmisc)

# Impute Missing or Zero RFV Values
impute_rfv_values <- function(row) {
  # This function imputes the rfv_l, rfv_r, and rfv_h columns with default or derived values if missing or zero.
  # If rfv_r is zero, it assigns default values to rfv_r, rfv_l, and rfv_h.
  # If rfv_r is between 0 and 2, it sets rfv_l to 0.01 and rfv_h to rfv_r + 2.
  # If rfv_r is greater than 2, it sets rfv_l to rfv_r - 2 and rfv_h to rfv_r + 2.
  #
  # @param row A dataframe row containing columns 'rfv_r', 'rfv_l', and 'rfv_h'.
  # @return The updated row with imputed 'rfv_r', 'rfv_l', and 'rfv_h' values.

  # Case 1: If 'rfv_r' is zero, assign default values
  if (row["rfv_r"] == 0) {
    row["rfv_r"] <- 0.02
    row["rfv_l"] <- 0.01
    row["rfv_h"] <- 0.03

    # Case 2: If 'rfv_r' is greater than 0 and less than or equal to 2
  } else if (row["rfv_r"] > 0 && row["rfv_r"] <= 2) {
    row["rfv_l"] <- 0.01  # Assign default low value if not already set
    row["rfv_h"] <- row["rfv_r"] + 2  # Set high value to 'rfv_r + 2'

    # Case 3: If 'rfv_r' is greater than 2
  } else if (row["rfv_r"] > 2) {
    row["rfv_l"] <- row["rfv_r"] - 2  # Set low value to 'rfv_r - 2'
    row["rfv_h"] <- row["rfv_r"] + 2  # Set high value to 'rfv_r + 2'
  }

  return(row)
}


remove_organic_layer <- function(df) {
  #' Remove Organic Layers and Adjust Depths (Including Depth Estimates)
  #'
  #' This function removes rows with organic horizons (where 'hzname' contains a capital 'O')
  #' and adjusts the depths ('hzdept_r', 'hzdepb_r', 'hzdept_l', 'hzdepb_l', 'hzdept_h', 'hzdepb_h')
  #' within each group (grouped by 'compname'). The depth adjustments ensure that the remaining horizons'
  #' depths are recalculated based on the removal of the organic layers, preserving the original horizon thicknesses.
  #' If the '_l' or '_h' values are missing or the columns do not exist, only the '_r' values will be adjusted.
  #'
  #' @param df A data frame containing soil horizon data, with columns 'hzname', 'hzdept_r', 'hzdepb_r',
  #'           and optionally 'hzdept_l', 'hzdepb_l', 'hzdept_h', 'hzdepb_h', and 'compname'.
  #' @return A data frame with organic layers removed and depth adjustments made within each group.
  #' @examples
  #' # Example dataframe
  #' df <- data.frame(
  #'   compname = c("P1", "P1", "P1", "P2", "P2"),
  #'   hzname = c("O", "A", "B", "O", "C"),
  #'   hzdept_r = c(0, 10, 20, 0, 15),
  #'   hzdepb_r = c(10, 20, 30, 15, 55)
  #' )
  #'
  #' # Apply the function to remove organic layers and adjust depths
  #' result <- remove_organic_layer(df)
  #' print(result)

  # Function to process each group
  process_group <- function(group) {
    # Filter out rows where 'hzname' contains 'O' (Organic layers)
    filtered_group <- group[!grepl("O", group$hzname), ]

    # If the group is empty after filtering, return it directly
    if (nrow(filtered_group) == 0) {
      return(filtered_group)
    }

    # After removing O horizons, take the first representative mineral horizon depth (baseline) and substract it from all horizon depths
    # Extract the first value of hzdept_r
    initial_hzdept_r <- filtered_group$hzdept_r[1]

    # Subtract the initial hzdept_r from all depth columns and set any negative values to 0
    filtered_group <- filtered_group %>%
      mutate(across(starts_with("hzdep"), ~ . - initial_hzdept_r))

    # Adjust 'hzdept_l' and 'hzdept_h' to 0 if 'hzdept_r' is 0
    for (i in seq_len(nrow(filtered_group))) {
      if (filtered_group$hzdept_r[i] == 0) {
        if ('hzdept_l' %in% names(filtered_group)) {
          if (!is.na(filtered_group$hzdept_l[i])) {
            filtered_group$hzdept_l[i] <- 0
          }
        }
        if ('hzdept_h' %in% names(filtered_group)) {
          if (!is.na(filtered_group$hzdept_h[i])) {
            filtered_group$hzdept_h[i] <- 0
          }
        }
      }
    }

    # Calculate the thickness of each horizon for 'hzdept_r' and 'hzdepb_r'
    thickness_r <- filtered_group$hzdepb_r - filtered_group$hzdept_r

    # Adjust depths for 'hzdept_r' and 'hzdepb_r'
    filtered_group$hzdept_r <- cumsum(c(0, thickness_r[-length(thickness_r)]))
    filtered_group$hzdepb_r <- cumsum(thickness_r)

    # Check if 'hzdept_l' and 'hzdepb_l' exist and adjust if they do
    if (all(c('hzdept_l', 'hzdepb_l') %in% names(filtered_group))) {
      if (!all(is.na(filtered_group$hzdept_l))) {
        thickness_l <- filtered_group$hzdepb_l - filtered_group$hzdept_l
        filtered_group$hzdept_l <- cumsum(c(0, thickness_l[-length(thickness_l)]))
        filtered_group$hzdepb_l <- cumsum(thickness_l)
      }
    }

    # Check if 'hzdept_h' and 'hzdepb_h' exist and adjust if they do
    if (all(c('hzdept_h', 'hzdepb_h') %in% names(filtered_group))) {
      if (!all(is.na(filtered_group$hzdept_h))) {
        thickness_h <- filtered_group$hzdepb_h - filtered_group$hzdept_h
        filtered_group$hzdept_h <- cumsum(c(0, thickness_h[-length(thickness_h)]))
        filtered_group$hzdepb_h <- cumsum(thickness_h)
      }
    }
    return(filtered_group)
  }
  # Group by 'cokey' and apply the processing function to each group
  result <- df %>%
    group_by(cokey) %>%
    group_modify(~ process_group(.x)) %>%
    ungroup()

  return(result)
}


infill_soil_data <- function(df) {
  #' Infill Soil Data for Missing Values
  #'
  #' This function imputes missing values for soil data, performing various replacements
  #' based on specific rules for particle size, bulk density, and water retention data.
  #' It handles grouped data by 'compname' and applies filters and imputations accordingly.
  #'
  #' @param df A dataframe containing soil data with columns such as 'sandtotal_r', 'claytotal_r',
  #'           'silttotal_r', 'dbovendry_r', 'wthirdbar_r', 'wfifteenbar_r', etc.
  #' @return A dataframe with missing values imputed according to the specified rules.

  # Step 1: Remove rows where all particle size 'r' values are missing
  df_filtered <- df %>%
    filter(
      !(is.na(sandtotal_r) & is.na(claytotal_r) & is.na(silttotal_r))
    )

  # Step 2: Group by 'compname'
  df_filtered <- df_filtered %>%
    group_by(compname) %>%
    # Step 3: Exclude groups where any 'hzdepb_r' <= 50 and any of the particle sizes are missing
    filter(
      !(
        any(hzdepb_r <= 50, na.rm = TRUE) &
          (any(is.na(sandtotal_r)) | any(is.na(claytotal_r)) | any(is.na(silttotal_r)))
      )
    ) %>%
    ungroup()

  # Step 4: Replace missing '_l' and '_h' for particle size values with '_r' values +/- 8 (row-wise)
  particle_cols <- c("sandtotal", "claytotal", "silttotal")
  for (col in particle_cols) {
    df_filtered <- df_filtered %>%
      mutate(
        # Replace missing lower bounds
        "{col}_l" := if_else(
          is.na(.data[[paste0(col, "_l")]]),
          pmax(.data[[paste0(col, "_r")]] - 8, 0),
          .data[[paste0(col, "_l")]]
        ),
        # Replace missing upper bounds
        "{col}_h" := if_else(
          is.na(.data[[paste0(col, "_h")]]),
          pmax(.data[[paste0(col, "_r")]] + 8, 0),
          .data[[paste0(col, "_h")]]
        )
      )
  }

  # Step 5: Replace missing 'dbovendry_l' and 'dbovendry_h' with 'dbovendry_r' +/- 0.01 (row-wise)
  df_filtered <- df_filtered %>%
    mutate(
      dbovendry_l = if_else(is.na(dbovendry_l), pmax(dbovendry_r - 0.01, 0), dbovendry_l),
      dbovendry_h = if_else(is.na(dbovendry_h), pmax(dbovendry_r + 0.01, 0), dbovendry_h)
    )

  # Step 6: Replace missing 'wthirdbar_l' and 'wthirdbar_h' with 'wthirdbar_r' +/- 1 (row-wise)
  df_filtered <- df_filtered %>%
    mutate(
      wthirdbar_l = if_else(is.na(wthirdbar_l), pmax(wthirdbar_r - 1, 0), wthirdbar_l),
      wthirdbar_h = if_else(is.na(wthirdbar_h), pmax(wthirdbar_r + 1, 0), wthirdbar_h)
    )

  # Step 7: Replace missing 'wfifteenbar_l' and 'wfifteenbar_h' with 'wfifteenbar_r' +/- 0.6 (row-wise)
  df_filtered <- df_filtered %>%
    mutate(
      wfifteenbar_l = if_else(is.na(wfifteenbar_l), pmax(wfifteenbar_r - 0.6, 0), wfifteenbar_l),
      wfifteenbar_h = if_else(is.na(wfifteenbar_h), pmax(wfifteenbar_r + 0.6, 0), wfifteenbar_h)
    )

  # Step 8: Impute 'rfv_l' and 'rfv_h' values based on 'rfv_r' (row-wise)
  df_filtered <- df_filtered %>%
    mutate(
      # Adjust 'rfv_r' if it is zero
      rfv_r = if_else(rfv_r == 0, 0.02, rfv_r),
      # Impute 'rfv_l' based on conditions
      rfv_l = case_when(
        rfv_r == 0 ~ 0.01,
        (rfv_r > 0 & rfv_r <= 2) & is.na(rfv_l) ~ 0.01,
        (rfv_r > 2) & is.na(rfv_l) ~ rfv_r - 2,
        TRUE ~ rfv_l  # Keep existing value if not missing
      ),
      # Impute 'rfv_h' based on conditions
      rfv_h = case_when(
        rfv_r == 0 ~ 0.03,
        (rfv_r > 0 & rfv_r <= 2) & is.na(rfv_h) ~ rfv_r + 2,
        (rfv_r > 2) & is.na(rfv_h) ~ rfv_r + 2,
        TRUE ~ rfv_h  # Keep existing value if not missing
      )
    )

  return(df_filtered)
}


slice_and_aggregate_soil_data <- function(df) {
  #' Slice and Aggregate Soil Data at set depth intervals (0-30, and 30-100)
  #'
  #' This function slices a data frame with soil data into 1 cm increments based on
  #' depth ranges provided in the 'hzdept_r' and 'hzdepb_r' columns, and calculates
  #' mean values for each depth increment across all other numeric data columns.
  #'
  #' @param df A data frame where each row represents a soil sample with 'hzdept_r'
  #'        and 'hzdepb_r' columns representing depth ranges.
  #' @return A data frame with depth ranges and mean values of soil properties for each range.

  # dplyr::select numeric columns for aggregation, excluding the depth range columns
  data_columns <- df %>%
    dplyr::select(where(is.numeric)) %>%
    dplyr::select(-hzdept_r, -hzdepb_r) %>%
    colnames()

  # Generate a data frame for each 1 cm increment within each row's depth range
  rows_list <- list()

  for (i in 1:nrow(df)) {
    row <- df[i, ]
    # Convert the dplyr::selected columns to a named vector
    row_data <- as.list(row[data_columns])
    for (depth in seq(row$hzdept_r, row$hzdepb_r - 1)) {
      # Add the 'Depth' value to the row data
      row_data$Depth <- depth
      # Append the modified row data to the list
      rows_list <- append(rows_list, list(as.data.frame(row_data, stringsAsFactors = FALSE)))
    }
  }

  # Combine the list of data frames into a single data frame
  aggregated_data <- do.call(rbind, rows_list)

  # Convert columns back to numeric where applicable
  aggregated_data$Depth <- as.numeric(aggregated_data$Depth)

  # # Calculate mean values for each depth increment
  # depth_increment_means <- aggregated_data %>%
  #   group_by(Depth) %>%
  #   summarise(across(everything(), mean, na.rm = TRUE)) %>%
  #   ungroup()

  # Define depth ranges (0-30 cm and 30-100 cm)
  depth_ranges <- list(c(0, 30), c(30, 100))
  results <- list()

  # Process each depth range
  for (range in depth_ranges) {
    top <- range[1]
    bottom <- range[2]

    # Subset data for the current depth range
    subset <- aggregated_data %>%
      dplyr::select(where(is.numeric)) %>%
      filter(Depth > top & Depth <= bottom)

    # Calculate the mean for each column in the subset
    layer_depth <- max(subset$Depth)
    mean_values <- colMeans(subset, na.rm = TRUE)
    mean_values["hzdept_r"] <- top
    if(layer_depth < bottom){
      mean_values["hzdepb_r"] <- layer_depth
    } else {
      mean_values["hzdepb_r"] <- bottom
    }

    results <- append(results, list(c(aggregated_data$compname[1], mean_values)))
  }

  # Convert results to a data frame
  result_df <- do.call(rbind, results) %>% as.data.frame()

  # Ensure missing ranges are added
  if (!(30 %in% result_df$hzdept_r)) {
    missing_row <- as.list(rep(NA, ncol(result_df)))
    names(missing_row) <- colnames(result_df)
    missing_row$hzdept_r <- 30
    missing_row$hzdepb_r <- 100
    result_df <- rbind(result_df, missing_row)
  }

  return(result_df)
}


# Triangular Distribution Sampling Function, code from the 'triangle' R package
tri_dist <- function(n = 1, a = 0, b = 1, c = (a + b)/2) {
  #' Triangular Distribution Sampling
  #'
  #' This function generates random samples from a triangular distribution
  #' defined by the minimum value (`a`), maximum value (`b`), and the mode (`c`).
  #'
  #' @param n Integer, number of random samples to generate (default is 1).
  #'          If `n` is a vector, the length of the vector is used.
  #' @param a Numeric, the minimum value of the distribution (default is 0).
  #' @param b Numeric, the maximum value of the distribution (default is 1).
  #' @param c Numeric, the mode or most likely value of the distribution.
  #'          It defaults to the midpoint `(a + b) / 2` if not provided.
  #' @return A vector of random samples from the triangular distribution.
  #'         If invalid inputs are provided, a vector of `NaN` is returned.
  #' @examples
  #' # Generate 10 random samples from a triangular distribution
  #' samples <- tri_dist(10, a = 0, b = 1, c = 0.5)
  #' print(samples)

  # If n is a vector, take its length as the number of samples
  if (length(n) > 1)
    n <- length(n)

  # Ensure n is a valid number and positive
  if (n < 1 | is.na(n))
    stop(paste("invalid argument: n =", n))

  # Ensure n is an integer (in case it isn't)
  n <- floor(n)

  # Check for any missing (NA) or invalid parameters in a, b, or c
  if (any(is.na(c(a, b, c))))
    return(rep(NaN, times = n))

  # Check if the mode (c) is within the range [a, b]
  if (a > c | b < c)
    return(rep(NaN, times = n))

  # Check for infinite values in a, b, or c
  if (any(is.infinite(c(a, b, c))))
    return(rep(NaN, times = n))

  # Generate n uniform random numbers between 0 and 1
  p <- runif(n)

  # Determine which values of p should be mapped to the left side of the triangle (a to c)
  # and which values should be mapped to the right side (c to b)
  if (a != c) {
    # For general case where mode (c) is not equal to the minimum (a)
    i <- which((a + sqrt(p * (b - a) * (c - a))) <= c)
    j <- which((b - sqrt((1 - p) * (b - a) * (b - c))) > c)
  } else {
    # Special case when a == c (the distribution starts from the mode)
    i <- which((a + sqrt(p * (b - a) * (c - a))) < c)
    j <- which((b - sqrt((1 - p) * (b - a) * (b - c))) >= c)
  }

  # Calculate the left side of the triangular distribution for indices in i
  if (length(i) != 0)
    p[i] <- a + sqrt(p[i] * (b - a) * (c - a))

  # Calculate the right side of the triangular distribution for indices in j
  if (length(j) != 0)
    p[j] <- b - sqrt((1 - p[j]) * (b - a) * (b - c))

  # Return the vector of random samples
  return(p)
}


# Function to simulate soil component composition from a dataframe
sim_component_comp <- function(data, n_simulations = 1000) {
  #' Simulate Soil Component Composition
  #'
  #' This function simulates soil component compositions from a given dataframe
  #' that includes the low, mode (representative), and high values for each component.
  #' It uses a triangular distribution to generate random values for each component
  #' and adjusts the results to ensure the total sum is 100%.
  #'
  #' @param data Dataframe containing the soil component information. It must contain
  #'        the columns: 'compname' (component names), 'comppct_r' (mode/representative value),
  #'        'comppct_l' (low value), and 'comppct_h' (high value).
  #' @param n_simulations Integer, number of simulations to run (default is 1000).
  #' @return A dataframe where each row corresponds to a simulated composition
  #'         and each column represents one component. The percentages in each row sum to 100.
  #'
  #' @examples
  #' # Create an example dataframe
  #' data_table <- data.frame(
  #'   compname = c("Component1", "Component2", "Component3", "Component4"),
  #'   comppct_r = c(30, 40, 20, 10),  # Mode (representative values)
  #'   comppct_l = c(25, 35, 15, 5),   # Low values
  #'   comppct_h = c(35, 45, 25, 15)   # High values
  #' )
  #'
  #' # Simulate 1000 compositions
  #' simulated_df <- simulate_soil_composition(data_table, n_simulations = 1000)
  #' head(simulated_df)

  data <- data %>%
    dplyr::select(mukey, cokey, compname, comppct_l, comppct_r, comppct_h) %>% dplyr::distinct()
  # Fill missing comppct_l and comppct_h based on comppct_r
  data$comppct_l[is.na(data$comppct_l)] <- data$comppct_r[is.na(data$comppct_l)] - 2
  data$comppct_h[is.na(data$comppct_h)] <- data$comppct_r[is.na(data$comppct_h)] + 2

  n_components <- nrow(data)  # Number of components
  compositions <- matrix(NA, nrow = n_simulations, ncol = n_components)  # Initialize matrix to store compositions for 1000 simulations

  # Generate 1000 random samples from triangular distributions for each component
  for (i in 1:n_components) {
    low <- data$comppct_l[i]  # Low value for the component
    mode <- data$comppct_r[i] # Mode (representative value) for the component
    high <- data$comppct_h[i] # High value for the component

    # Sample 1000 values from the triangular distribution for each component
    compositions[, i] <- tri_dist(n_simulations, a = low, b = high, c = mode)

  }

  # Normalize across each row to ensure the total percentage sums to 100
  sim_compositions <- apply(compositions, 2, function(col) {
    round( sum(col) / 100, 0)  # Normalize each row so it sums to 100
  })
  data$sim_comppct <- sim_compositions

  return(data)
}


# Function to simulate correlated samples from triangular distributions
simulate_correlated_triangular <- function(n, params, correlation_matrix, random_seed = NULL) {
  #' Simulate Correlated Samples from Triangular Distributions
  #'
  #' This function generates correlated random samples from multiple triangular distributions
  #' with specified lower limits, modes, and upper limits. The correlations between the variables
  #' are controlled by the input correlation matrix. The function uses a Cholesky decomposition
  #' to induce the specified correlations between the generated samples.
  #'
  #' @param n Integer, the number of samples to generate.
  #' @param params List of parameter sets for the triangular distributions, where each element
  #'        is a vector of three values:
  #'        - Lower limit (`a`), Mode (`b`), and Upper limit (`c`) for the triangular distribution.
  #'        For example: `params = list(c(a1, b1, c1), c(a2, b2, c2), ...)`.
  #' @param correlation_matrix A square matrix specifying the desired correlations between
  #'        the variables. It must be positive semi-definite and of size equal to the number of triangular distributions.
  #' @param random_seed Optional integer, to set a specific random seed for reproducibility.
  #' @return A matrix of correlated samples where each column corresponds to samples from one triangular distribution.
  #' @details
  #' The function first generates uncorrelated standard normal variables and then uses a Cholesky
  #' decomposition of the correlation matrix to introduce correlations. The correlated standard
  #' normal variables are then transformed into uniform variables using the cumulative distribution
  #' function (CDF) of the normal distribution. Finally, these uniform variables are transformed
  #' into samples from triangular distributions using the inverse CDF of the triangular distribution.
  #'
  #' @examples
  #' # Define parameters for three triangular distributions
  #' params <- list(c(0, 5, 10),  # Triangle distribution with limits 0, 10 and mode 5
  #'                c(1, 4, 6),   # Triangle distribution with limits 1, 6 and mode 4
  #'                c(2, 3, 5))   # Triangle distribution with limits 2, 5 and mode 3
  #'
  #' # Define a correlation matrix
  #' correlation_matrix <- matrix(c(1, 0.5, 0.3,
  #'                                0.5, 1, 0.2,
  #'                                0.3, 0.2, 1), nrow = 3)
  #'
  #' # Simulate 1000 correlated samples
  #' samples <- simulate_correlated_triangular(1000, params, correlation_matrix, random_seed = 42)
  #' head(samples)

  # Optionally set the seed for reproducibility
  if (!is.null(random_seed)) {
    set.seed(random_seed)
  }

  # Generate uncorrelated standard normal variables
  uncorrelated_normal <- MASS::mvrnorm(n, mu = rep(0, length(params)), Sigma = diag(length(params)))

  # Cholesky decomposition of the correlation matrix
  L <- chol(correlation_matrix)

  # Compute correlated normal variables by multiplying the uncorrelated normals with L
  correlated_normal <- uncorrelated_normal %*% L

  # Initialize a matrix to store the samples for each triangular distribution
  samples <- matrix(NA, nrow = n, ncol = length(params))

  # Loop through each triangular distribution
  for (i in seq_along(params)) {
    a <- params[[i]][1]  # Lower limit of the triangle distribution
    b <- params[[i]][2]  # Mode (peak) of the triangle distribution
    c <- params[[i]][3]  # Upper limit of the triangle distribution

    # Get the correlated normal variable for this triangular distribution
    normal_var <- correlated_normal[, i]

    # Transform the normal variable to a uniform distribution using the CDF of the normal distribution
    u <- pnorm(normal_var)

    # Handle the case of division by zero for degenerate triangular distributions
    if (c == a) {
      samples[, i] <- rep(a, n)
    } else {
      # Calculate samples for the triangular distribution
      condition <- u <= (b - a) / (c - a)

      # For values satisfying the condition (left side of the triangle)
      samples[condition, i] <- a + sqrt(u[condition] * (c - a) * (b - a))

      # For values not satisfying the condition (right side of the triangle)
      samples[!condition, i] <- c - sqrt((1 - u[!condition]) * (c - a) * (c - b))
    }
  }

  # Return the matrix of correlated samples
  return(samples)
}


#' Adjust Depthwise Property Values Using GP Predictions
#'
#' This function adjusts simulated property values across different depths to align with the relative changes predicted by a Gaussian Process (GP) model. The adjustment ensures that the simulated values reflect the trends indicated by the GP model while maintaining their statistical properties.
#'
#' @param simulated_values A numeric matrix of simulated property values, where each row corresponds to a depth and each column corresponds to a simulation.
#' @param gp_model A fitted Gaussian Process model used to predict property means at different depths.
#' @param depths A numeric vector of depth values corresponding to the rows of `simulated_values`.
#'
#' @return A numeric matrix of adjusted simulated property values with the same dimensions as `simulated_values`.
#'
#' @examples
#' # Assuming 'simulated_values' is your matrix of simulated data,
#' # 'gp_model' is your fitted GP model, and 'depths' is your depth vector:
#' adjusted_values <- adjust_depthwise_property_GP_Quant(simulated_values, gp_model, depths)
#'
#' @export
adjust_depthwise_property_GP_Quant <- function(simulated_values, gp_model, depths) {

  # Get GP-predicted means for each depth
  gp_pred_means <- predict.GP(gp_model, newdata = as.matrix(depths))$Y_hat

  # Initialize adjusted values matrix with the original simulated values
  adjusted_values <- simulated_values

  # Loop through each depth starting from the second one
  for (i in 2:length(depths)) {

    # Extract simulated values at the previous and current depths
    simulated_prev <- adjusted_values[i - 1, ]
    simulated_curr <- adjusted_values[i, ]

    # Compute the empirical cumulative distribution function (ECDF) of the previous depth's simulated values
    prev_ecdf <- ecdf(simulated_prev)
    prev_quantiles <- prev_ecdf(simulated_prev)

    # Calculate the ratio of GP-predicted means between current and previous depths
    gp_mean_ratio <- gp_pred_means[i] / gp_pred_means[i - 1]

    # Calculate the standard deviation adjustment factor
    sd_adjustment_factor <- sd(simulated_curr) / sd(simulated_prev)

    # Initialize a vector to store adjusted simulated values for the current depth
    adjusted_curr <- numeric(length(simulated_curr))

    # Adjust each simulated value at the current depth
    for (j in 1:length(simulated_curr)) {
      # Get the quantile for the simulated value at the previous depth
      q <- prev_quantiles[j]

      # Find the corresponding value at the same quantile in the current depth's distribution
      quantile_value <- quantile(simulated_curr, probs = q)

      # Adjust the current value based on the GP mean ratio and SD adjustment factor
      adjusted_curr[j] <- quantile_value + (simulated_prev[j] * gp_mean_ratio - quantile_value) * sd_adjustment_factor
    }

    # Store the adjusted simulated values back into the matrix
    adjusted_values[i, ] <- adjusted_curr
  }

  # Return the matrix of adjusted simulated values
  return(adjusted_values)
}


readjust_between_property_correlations <- function(adjusted_properties, correlation_matrix) {
  #'
  #' This function adjusts the between-property correlations for each depth/horizon
  #' after depth-wise adjustments using Cholesky decomposition. It assumes that
  #' the properties are already adjusted for depth-wise trends and now need to
  #' conform to a specific correlation matrix within each horizon.
  #'
  #' @param adjusted_properties A matrix of properties that have been depth-adjusted.
  #'        Each row is a horizon, and each column is a property (e.g., sand, silt, clay).
  #' @param correlation_matrix A square matrix specifying the desired correlations between
  #'        the variables (e.g., sand, silt, clay, bulk density, etc.). It must be positive
  #'        semi-definite and of size equal to the number of properties (columns).
  #'
  #' @return A matrix of adjusted properties with the desired between-property correlations.

  # Check if the correlation matrix is positive semi-definite
  if (!isSymmetric(correlation_matrix) || any(eigen(correlation_matrix)$values < 0)) {
    stop("Correlation matrix must be symmetric and positive semi-definite.")
  }

  # Cholesky decomposition of the correlation matrix
  L <- chol(correlation_matrix)

  # Generate correlated normal variables based on the adjusted properties
  n <- nrow(adjusted_properties)
  uncorrelated_normal <- MASS::mvrnorm(n, mu = rep(0, ncol(adjusted_properties)), Sigma = diag(ncol(adjusted_properties)))

  # Compute the correlated normal variables by multiplying with L
  correlated_normal <- uncorrelated_normal %*% L

  # Adjust each property's distribution using the cumulative distribution function (CDF)
  for (i in seq_len(ncol(adjusted_properties))) {
    # Use the CDF of the adjusted properties to match the new correlated normal variables
    u <- pnorm(correlated_normal[, i])

    # Re-map the correlated normal variable to the adjusted property range
    adjusted_properties[, i] <- quantile(adjusted_properties[, i], u)
  }

  # Return the re-adjusted properties matrix
  return(adjusted_properties)
}


simulate_soil_properties <- function(data) {
  #'
  #' This function simulates soil properties based on input data from a soil profile dataset.
  #' It calculates a local soil property correlation matrix and simulates soil properties.
  #'
  #' @param data A data frame containing input soil horizon data.
  #'
  #' @return A data frame with simulated soil properties.

  # Step 1: Calculate a local soil property correlation matrix
  # Subset data based on required soil inputs
  sim_columns <- c(
    "compname", "mukey", "cokey", "hzname", "sim_comppct", "hzdept_r", "hzdepb_r",
    "sandtotal_l", "sandtotal_r", "sandtotal_h",
    "silttotal_l", "silttotal_r", "silttotal_h",
    "claytotal_l", "claytotal_r", "claytotal_h",
    "dbovendry_l", "dbovendry_r", "dbovendry_h",
    "wthirdbar_l", "wthirdbar_r", "wthirdbar_h",
    "wfifteenbar_l", "wfifteenbar_r", "wfifteenbar_h",
    "rfv_l", "rfv_r", "rfv_h"
  )

  sim_data_columns <- c(
    "hzdept_r", "hzdepb_r", "sandtotal_l", "sandtotal_r", "sandtotal_h",
    "silttotal_l", "silttotal_r", "silttotal_h", "claytotal_l", "claytotal_r", "claytotal_h",
    "dbovendry_l", "dbovendry_r", "dbovendry_h", "wthirdbar_l", "wthirdbar_r", "wthirdbar_h",
    "wfifteenbar_l", "wfifteenbar_r", "wfifteenbar_h", "rfv_l", "rfv_r", "rfv_h"
  )

  sim_data <- data %>% dplyr::select(sim_columns)

  # Fill missing soil data
  sim_data <- infill_soil_data(sim_data)

  if (!"genhz" %in% colnames(sim_data)) {
    # Remove any digits from the beginning of hzname
    sim_data$genhz <- sub("^[0-9]+", "", sim_data$hzname)

    sim_data$genhz <- generalizeHz(
      sim_data$genhz,
      new = c('O', 'A', 'B', 'C', 'Cr', 'R'),
      pattern = c('O', '^A', '^B', '^C', '^Cr', '^R'),
      ordered = TRUE
    )
  }

  # Define global correlation matrix data (used if local matrix calculation fails)
  # parameter order: bulk_density_third_bar, water_retention_third_bar, water_retention_15_bar, ilr1, ilr2
  global_correlation_df <- data.frame(
    bulk_density_third_bar = c(1.00000000, -0.69590081, -0.52909188, -0.23033210, -0.10138811,  0.02442569),
    water_retention_third_bar = c(-0.69590081, 1.00000000, 0.65403992, 0.24722247, 0.19316597, -0.05426055),
    water_retention_15_bar = c(-0.52909188, 0.65403992, 1.00000000, 0.35081340, 0.52221620, -0.12133250),
    ilr1 = c(-0.23033210, 0.24722247, 0.35081340, 1.00000000, 0.60031030, -0.21691814),
    ilr2 = c(-0.10138811, 0.19316597, 0.52221620, 0.60031030, 1.00000000, -0.23230431),
    rfv = c(0.02442569, -0.05426055, -0.12133250, -0.21691814, -0.23230431, 1.00000000)
  )

  # Assign row names
  row.names(global_correlation_df) <- names(global_correlation_df)

  global_correlation_matrix <- as.matrix(global_correlation_df)

  # parameter order: sand_total, silt_total, clay_total
  texture_correlation_matrix <- matrix(
    c(1.0000000, -0.76231798, -0.67370589,
      -0.7623180, 1.00000000, 0.03617498,
      -0.6737059, 0.03617498, 1.00000000),
    nrow = 3, byrow = TRUE
  )

  # Define the is.positive.definite function
  is.positive.definite <- function(mat) {
    result <- tryCatch({
      chol(mat)
      TRUE
    }, error = function(e) {
      FALSE
    })
    return(result)
  }

  # Extract columns with names ending in '_r'
  sim_data_r <- sim_data %>% dplyr::select(dplyr::ends_with("_r"), genhz)

  # Compute the local correlation matrix (Spearman correlation)
  rep_columns <- sim_data_r %>% dplyr::select(-sandtotal_r, -silttotal_r, -claytotal_r)

  ilr_site_txt <- ilr(sim_data %>% dplyr::select(sandtotal_r, silttotal_r, claytotal_r)) %>% as.matrix()
  rep_columns$ilr1 <- ilr_site_txt[,1]
  rep_columns$ilr2 <- ilr_site_txt[,2]

  # Perform the replacement
  exclude_cols <- c("hzdept_r", "hzdepb_r", "genhz")
  include_cols <- setdiff(names(rep_columns), exclude_cols)

  rep_columns[include_cols] <- lapply(rep_columns[include_cols], function(column) {
    # Replace zeros with 0.01
    column[column == 0] <- 0.01
    return(column)
  })
  # Define columns for correlation calculation
  cor_cols <- c("dbovendry_r", "wthirdbar_r", "wfifteenbar_r", "ilr1", "ilr2", "rfv_r")

  # Reset factor levels for genhz
  rep_columns$genhz <- factor(rep_columns$genhz)  # Reset factor levels to the levels present in rep_columns$genhz

  # Initialize an empty list to store the correlation matrices
  correlation_matrices <- list()

  # Get the unique genetic horizons
  genhz_levels <- unique(rep_columns$genhz)

  # Loop through each genetic horizon and calculate the correlation matrix
  for (genhz in genhz_levels) {
    # Filter data for the current genetic horizon
    rep_columns_filtered <- rep_columns %>%
      filter(genhz == !!genhz, wthirdbar_r != 0.01)

    # Select columns needed for correlation calculation
    correlation_matrix_data <- rep_columns_filtered %>%
      dplyr::select(all_of(cor_cols))

    # Calculate the correlation matrix if there are more than 4 observations
    if (nrow(rep_columns_filtered) > 4) {
      local_correlation_matrix <- rcorr(as.matrix(correlation_matrix_data), type = "pearson")$r
    } else {
      # Use the global correlation matrix if there are insufficient observations
      local_correlation_matrix <- global_correlation_matrix
    }

    # Ensure symmetry in the correlation matrix
    local_correlation_matrix <- (local_correlation_matrix + t(local_correlation_matrix)) / 2

    # Check positive definiteness of the matrix
    if (is.positive.definite(local_correlation_matrix)) {
      # Keep the local matrix as is
    } else if (any(is.na(local_correlation_matrix)) || any(is.infinite(local_correlation_matrix))) {
      # Replace with the global matrix if NaNs or Infs are present
      local_correlation_matrix <- global_correlation_matrix
    } else {
      # Adjust to be near positive definite
      local_correlation_matrix <- as.matrix(nearPD(local_correlation_matrix)$mat)
    }

    # Save the resulting matrix in the list using the genetic horizon as the key
    correlation_matrices[[genhz]] <- local_correlation_matrix
  }

  # Initialize an empty list to store the results
  result_list <- list()

  # Get the unique 'cokey' values
  unique_cokeys <- unique(sim_data$cokey)

  # Loop over each unique 'cokey'
  for (i in seq_along(unique_cokeys)) {
    # Current 'cokey' value
    cokey_value <- unique_cokeys[i]

    # Subset the data for the current 'cokey'
    sim_subset <- sim_data[sim_data$cokey == cokey_value, ]

    # Apply the 'simulate_cokey' function to the subset
    sim_result <- simulate_cokey(sim_subset, correlation_matrices, texture_correlation_matrix)

    # Append the result to the list
    result_list[[i]] <- sim_result
  }

  # Combine all results into a single data frame
  sim_data_df <- do.call(rbind, result_list)

  return(sim_data_df)
}


simulate_cokey <- function(sim_cokey, correlation_matrices, texture_correlation_matrix) {

  # Step 2: Simulate data for each row
  sim_data_out <- list()

  for (index in seq_len(nrow(sim_cokey))) {
    row <- sim_cokey[index, ]
    local_correlation_matrix <- correlation_matrices[[as.character(row$genhz)]]

    # Step 2a: Simulate sand, silt, clay percentages
    params_txt <- list(
      c(row$sandtotal_l, row$sandtotal_r, row$sandtotal_h),
      c(row$silttotal_l, row$silttotal_r, row$silttotal_h),
      c(row$claytotal_l, row$claytotal_r, row$claytotal_h)
    )

    # Step 2b: Convert simulated data using acomp and ilr transformation
    simulated_txt <- acomp(simulate_correlated_triangular((as.integer(row$sim_comppct)), params_txt, texture_correlation_matrix))
    simulated_txt_ilr <- ilr(simulated_txt)

    # Step 2c: Extract l, r, h values for ilr1 and ilr2
    ilr1_values <- simulated_txt_ilr[, 1]
    ilr2_values <- simulated_txt_ilr[, 2]

    ilr1_l <- min(ilr1_values)
    ilr1_r <- median(ilr1_values)
    ilr1_h <- max(ilr1_values)
    ilr2_l <- min(ilr2_values)
    ilr2_r <- median(ilr2_values)
    ilr2_h <- max(ilr2_values)

    params <- list(
      c(row$dbovendry_l, row$dbovendry_r, row$dbovendry_h),
      c(row$wthirdbar_l, row$wthirdbar_r, row$wthirdbar_h),
      c(row$wfifteenbar_l, row$wfifteenbar_r, row$wfifteenbar_h),
      c(ilr1_l, ilr1_r, ilr1_h),
      c(ilr2_l, ilr2_r, ilr2_h),
      c(row$rfv_l, row$rfv_r, row$rfv_h)
    )

    # Step 2d: Simulate all properties and perform inverse ilr transformation
    tryCatch({
      n_sim <- as.integer(row$sim_comppct)
      sim_data <- simulate_correlated_triangular(
        n = n_sim,
        params = params,
        correlation_matrix = local_correlation_matrix
      )

      sim_data <- data.frame(sim_data)
      colnames(sim_data) <- c("bulk_density_third_bar", "water_retention_third_bar", "water_retention_15_bar", "ilr1", "ilr2", "rfv")

      sim_data$water_retention_third_bar <- sim_data$water_retention_third_bar / 100
      sim_data$water_retention_15_bar <- sim_data$water_retention_15_bar / 100
      sim_txt <- ilrInv(sim_data[, c("ilr1", "ilr2")])
      sim_txt <- data.frame(sim_txt)
      colnames(sim_txt) <- c("sand_total", "silt_total", "clay_total")
      sim_txt <- sim_txt * 100
      multi_sim <- cbind(sim_data[, !colnames(sim_data) %in% c("ilr1", "ilr2")], sim_txt)
      multi_sim$compname <- row$compname
      multi_sim$mukey <- row$mukey
      multi_sim$cokey <- row$cokey
      multi_sim$hzdept_r <- row$hzdept_r
      multi_sim$hzdepb_r <- row$hzdepb_r

      # Add simulation number
      multi_sim$simulation_number <- seq_len(nrow(multi_sim))

      # Create unique identifier: 'cokey' plus "_" and the simulation number
      multi_sim$unique_id <- paste0(row$cokey, "_", multi_sim$simulation_number)

      sim_data_out <- append(sim_data_out, list(multi_sim))
    }, error = function(e) {
      cat("Error: ", e$message, "\n")
    })
  }

  # Concatenate simulated values to a data frame
  sim_data_df <- bind_rows(sim_data_out) %>% replace(is.na(.), NA)

  return(sim_data_df)
}

# van Genuchten function with unit conversion
van_genuchten <- function(h, alpha, n, theta_r, theta_s) {
  #' van Genuchten Soil Water Retention Model
  #'
  #' This function computes the volumetric water content (`theta`) for a given
  #' matric potential (`h`) using the van Genuchten equation. The equation models
  #' the soil water retention curve, describing the relationship between soil water
  #' content and soil matric potential.
  #'
  #' @param h Numeric, the matric potential (pressure head) in cmH2O. Often negative
  #'          for unsaturated soils (e.g., -33 cmH2O for field capacity, -1500 cmH2O for wilting point).
  #' @param alpha Numeric, the van Genuchten parameter related to the inverse of the air-entry value.
  #'              Units are in 1/cm.
  #' @param n Numeric, the van Genuchten shape parameter related to pore size distribution.
  #'          Must be greater than 1.
  #' @param theta_r Numeric, the residual water content (minimum water content).
  #' @param theta_s Numeric, the saturated water content (maximum water content).
  #'
  #' @return Numeric, the predicted volumetric water content (`theta`) at the given matric potential (`h`).
  #'
  #' @details The van Genuchten model equation is given by:
  #'   \deqn{\theta(h) = \theta_r + \frac{(\theta_s - \theta_r)}{(1 + (|\alpha h|)^n)^{m}}}
  #' where \eqn{m = 1 - \frac{1}{n}}. This equation models the water retention curve
  #' in soils, describing how water content decreases with decreasing matric potential.
  #'
  #' @examples
  #' # Example: Calculate water content at field capacity (-33 kPa) using van Genuchten parameters
  #' h_fc <- -33 * 10.19716  # Field capacity matric potential in cmH2O (conversion from kPa)
  #' alpha <- 0.08  # Example alpha (1/cm)
  #' n <- 1.2  # Example n parameter (dimensionless)
  #' theta_r <- 0.05  # Residual water content
  #' theta_s <- 0.4   # Saturated water content
  #' theta_fc <- van_genuchten(h_fc, alpha, n, theta_r, theta_s)
  #' print(theta_fc)

  # Calculate m based on the n parameter
  m <- 1 - (1 / n)

  # van Genuchten equation for soil water retention
  return( as.numeric(theta_r + (theta_s - theta_r) / ((1 + (abs(alpha * h))^n)^m)))
}


# Function to simulate plant available water storage
simulate_vg_aws <- function(data, n_simulations = 100) {
  #' Simulate Available Water Holding Capacity (AWHC) using van Genuchten Parameters from the ROSETTA Model
  #'
  #' This function simulates the Available Water Holding Capacity (AWHC) for soil components
  #' based on the variability of van Genuchten parameters estimated from the ROSETTA model.
  #' The parameters for the van Genuchten model (alpha, n, theta_r, theta_s) are sampled using
  #' Monte Carlo simulations based on the provided means and standard deviations for each parameter.
  #' The function calculates water retention at specific matric potentials (Field Capacity and
  #' Permanent Wilting Point) and computes AWHC as the difference between these two points.
  #'
  #'
  #' @param data Dataframe containing soil parameters for each component. It must include the following columns:
  #'        - `alpha`: The van Genuchten alpha parameter (mean).
  #'        - `sd_alpha`: The standard deviation for alpha.
  #'        - `npar`: The van Genuchten n parameter (mean).
  #'        - `sd_npar`: The standard deviation for n.
  #'        - `theta_r`: The residual water content (mean).
  #'        - `sd_theta_r`: The standard deviation for theta_r.
  #'        - `theta_s`: The saturated water content (mean).
  #'        - `sd_theta_s`: The standard deviation for theta_s.
  #'        - `compname`: The component name.
  #'        - `hzname`: The horizon name.
  #' @param n_simulations Integer, number of Monte Carlo simulations to perform (default is 100).
  #' @return A list where each element contains simulated water retention parameters for a given soil component.
  #'         Each entry is a dataframe with the van Genuchten parameters and the corresponding
  #'         water content at Field Capacity (FC), Permanent Wilting Point (PWP), and Available Water Holding Capacity (AWHC).
  #'
  #' @examples
  #' # Example input data
  #' data <- data.frame(
  #'   compname = c("Component1", "Component2"),
  #'   hzname = c("A", "Bt"),
  #'   alpha = c(-1.84, -1.73),
  #'   sd_alpha = c(0.06, 0.06),
  #'   npar = c(0.13, 0.11),
  #'   sd_npar = c(0.01, 0.01),
  #'   theta_r = c(0.05, 0.07),
  #'   sd_theta_r = c(0.01, 0.015),
  #'   theta_s = c(0.4, 0.45),
  #'   sd_theta_s = c(0.02, 0.025)
  #' )
  #' # Run the water retention simulation
  #' results <- calculate_water_retention(data, n_simulations = 1000)
  #'
  #' # Access the results for a specific component
  #' head(results[["Component1_A"]])

  # Initialize result list to store the simulations for each component
  simulation_results <- list()

  # Iterate over each row in the data
  for (i in 1:nrow(data)) {

    # Extract van Genuchten parameters and their uncertainties for the current component
    alpha_mean <- data$alpha[i]
    alpha_sd <- data$sd_alpha[i]
    n_mean <- data$npar[i]
    n_sd <- data$sd_npar[i]
    theta_r_mean <- data$theta_r[i]
    theta_r_sd <- data$sd_theta_r[i]
    theta_s_mean <- data$theta_s[i]
    theta_s_sd <- data$sd_theta_s[i]

    # If any of the parameters are NA, skip this row
    if (is.na(alpha_mean) | is.na(n_mean) | is.na(theta_r_mean) | is.na(theta_s_mean)) {
      next
    }

    # Generate random samples for the van Genuchten parameters based on the provided means and standard deviations
    set.seed(123)  # Set seed for reproducibility
    alpha_samples <- rnorm(n_simulations, mean = alpha_mean, sd = alpha_sd)
    n_samples <- rnorm(n_simulations, mean = n_mean, sd = n_sd)
    theta_r_samples <- rnorm(n_simulations, mean = theta_r_mean, sd = theta_r_sd)
    theta_s_samples <- rnorm(n_simulations, mean = theta_s_mean, sd = theta_s_sd)

    # Convert matric potentials from kPa to cmH2O (1 kPa = 10.19716 cmH2O)
    h_fc <- -33 * 10.19716    # Field capacity (FC) matric potential in cmH2O
    h_pwp <- -1500 * 10.19716 # Permanent wilting point (PWP) matric potential in cmH2O

    # Create a dataframe to hold the sampled parameters
    results <- data.frame(
      alpha = 10^(alpha_samples),   # Convert alpha back from log scale
      n = 10^(n_samples),           # Convert npar back from log scale
      theta_r = theta_r_samples,
      theta_s = theta_s_samples,
      sim_num = seq(1:100)
    )

    # Calculate water retention at FC and PWP for each set of sampled parameters using the van Genuchten function
    results <- results %>%
      rowwise() %>%
      mutate(
        theta_fc = van_genuchten(h_fc, alpha, n, theta_r, theta_s),   # Water content at Field Capacity
        theta_pwp = van_genuchten(h_pwp, alpha, n, theta_r, theta_s)  # Water content at Permanent Wilting Point
      ) %>%
      ungroup()

    # Calculate Available Water Holding Capacity (AWHC) as the difference between FC and PWP
    results$AWHC <- results$theta_fc - results$theta_pwp

    # Store the results in the list with component name and horizon as the key
    simulation_results[[paste(data$layerID[i], i, sep = "_")]] <- results
  }

  # Return the list of simulation results for each component
  return(simulation_results)
}


calculate_aws <- function(sim_data_df) {
  #'
  #' This function calculates available water storage in the top 100 cm of the soil profile
  #' based on simulated soil properties.
  #'
  #' @param sim_data_df A data frame containing simulated soil properties.
  #'
  #' @return A list with available water storage prediction interval (aws_PIW90) and variable importance (var_imp).

  # Step 3: Run Rosetta and other Van Genuchten equations

  # Step 3a: Run Rosetta
  variables <- c(
    "sand_total", "silt_total", "clay_total",
    "bulk_density_third_bar", "water_retention_third_bar", "water_retention_15_bar"
  )
  rosetta_data <- soilDB::ROSETTA(sim_data_df, vars = variables, v="3", include.sd = TRUE)

  # Create layerID
  sim_data_df$layerID <- paste(sim_data_df$compname, sim_data_df$hzdept_r, sep = "_")
  rosetta_data$layerID <- sim_data_df$layerID

  aws <- simulate_vg_aws(rosetta_data)

  # Use Map to add new columns to each data frame in the list
  aws <- Map(function(df, top_value, bottom_value, compname_value) {
    df$top <- top_value
    df$bottom <- bottom_value
    df$compname <- compname_value
    return(df)
  }, aws, sim_data_df$hzdept_r, sim_data_df$hzdepb_r, sim_data_df$compname)

  # Combine the data frames into one
  sim_aws_df <- bind_rows(aws)

  aws_grouped <- sim_aws_df %>% group_by(top)
  data_len_depth <- summarise(aws_grouped, depth_len = n())
  sim_aws_df <- left_join(sim_aws_df, data_len_depth, by = "top")

  # Step 3b: Reshape aws data by ROI
  aws_grouped_bottom <- sim_aws_df %>% group_by(bottom)
  aws_quant_list <- summarise(aws_grouped_bottom, aws_quant = quantile(sim_aws_df, probs = c(0.05, 0.50, 0.95)))

  # Step 3c: Rename and group data
  aws05 <- calculate_aws(aws_quant_list, "0.05")
  aws95 <- calculate_aws(aws_quant_list, "0.95")

  # Width of 90th prediction interval for available water storage (aws_PIW90)
  aws_PIW90 <- round(aws95$aws0.95_100 - aws05$aws0.05_100, 2)

  return(list(aws_PIW90 = aws_PIW90, var_imp = "Data not available"))
}


# Calculate Available Water Storage (AWS) for Region of Interest (ROI)
calculate_aws <- function(df, quantile) {
  total <- sum(df[[quantile]] * df$depth * df$n)
  return(data.frame(!!paste0("aws", quantile, "_100") := total))
}


# Load necessary libraries
library(soilDB)

# Define the function
get_aws_data_by_mukey <- function(mukeys) {
  #' Query SSURGO database for soil data by mukey
  #'
  #' This function queries the NRCS SSURGO database for specific soil attributes based on the provided map unit key (mukey).
  #'
  #' @param mukeys A vector of string or numeric map unit keys (mukeys) representing the map units to retrieve data for.
  #'
  #' @return A data frame with soil horizon data for the given mukeys, including aggregated rock fragment data.

  # Step 1: Format mukey(s) for SQL
  formatted_mukey <- paste0("(", paste0("'", mukeys, "'", collapse = ","), ")")

  # Step 2: Construct the SQL query
  query <- sprintf("
    SELECT
    -- contextual data
    co.mukey, co.cokey, compname, comppct_l, comppct_r, comppct_h,
    -- horizon morphology
    ch.chkey, hzname, hzdept_l, hzdept_r, hzdept_h, hzdepb_l, hzdepb_r, hzdepb_h, hzthk_l, hzthk_r, hzthk_h,
    -- soil texture parameters
    sandtotal_l, sandtotal_r, sandtotal_h,
    silttotal_l, silttotal_r, silttotal_h,
    claytotal_l, claytotal_r, claytotal_h,
    -- bulk density and water retention parameters
    dbovendry_l, dbovendry_r, dbovendry_h,
    wthirdbar_l, wthirdbar_r, wthirdbar_h,
    wfifteenbar_l, wfifteenbar_r, wfifteenbar_h,
    -- rock fragment volume from chfrags table
    chf.fragvol_l AS rfv_l, chf.fragvol_r AS rfv_r, chf.fragvol_h AS rfv_h
    -- tables of interest
    FROM mapunit AS mu
    -- implicit filtering
    INNER JOIN component AS co ON mu.mukey = co.mukey
    INNER JOIN chorizon AS ch ON co.cokey = ch.cokey
    -- join with chfrags table using the correct chkey
    LEFT JOIN chfrags AS chf ON ch.chkey = chf.chkey
    -- filter by provided mukey
    WHERE mu.mukey IN %s
    -- order for visual inspection
    ORDER BY co.cokey, ch.hzdept_r ASC;", formatted_mukey)

  # Step 3: Submit the query and fetch results
  result <- SDA_query(query)

  # Step 4: Group by `chkey` and sum rock fragment values (rfv_l, rfv_r, rfv_h)
  rfv_sum <- result %>%
    dplyr::group_by(chkey) %>%
    summarise(
      rfv_l = sum(rfv_l, na.rm = TRUE),
      rfv_r = sum(rfv_r, na.rm = TRUE),
      rfv_h = sum(rfv_h, na.rm = TRUE)
    ) %>%
    ungroup()

  # Step 5: Merge the summed values back to the main dataset and remove duplicates
  result <- result %>%
    dplyr::select(-c(rfv_l, rfv_r, rfv_h)) %>%  # Remove original rfv columns
    left_join(rfv_sum, by = "chkey") %>%  # Merge summed values
    distinct()  # Remove any duplicates

  # Step 6: Return the final result
  return(result)
}

#--------------------------------------------------------------------------------------------------------
# Functions to run soil profile depth simulations

# Query OSD Data and Convert Horizon Distinctness to Offset
query_osd_distinctness <- function(horizon_data) {
  #'
  #' This function retrieves Official Series Description (OSD) data for a given list of soil series (component names),
  #' extracts horizon-related fields such as id, top, bottom, hzname, and distinctness, and converts the horizon
  #' distinctness codes into offset values using the `hzDistinctnessCodeToOffset` function.
  #'
  #' The function returns a data frame containing the extracted horizon data along with the distinctness code and the
  #' corresponding offset value.
  #'
  #' @param ssurgo_horizon_data Data frame containing SSURGO horizon data.
  #'
  #' @return A data frame containing the following columns:
  #'         - `id`: Unique identifier for the soil profile.
  #'         - `top`: The top depth of the horizon (in cm).
  #'         - `bottom`: The bottom depth of the horizon (in cm).
  #'         - `hzname`: The name of the soil horizon.
  #'         - `distinctness`: The distinctness code for the horizon (e.g., A, C, G, D).
  #'         - `bound_sd`: The calculated SD value from the distinctness code.
  #'
  #' @examples
  #' # Example soil series to query
  #' soils <- c('amador', 'pentz', 'pardee', 'auburn', 'loafercreek', 'millvilla')
  #'
  #' # Query OSD data and convert distinctness codes to offsets
  #' result <- query_osd_distinctness(soils)
  #'
  #' # View the result
  #' head(result)
  #'
  #' @export

  # subset ssurgo horizon data
  ssurgo_horizon_data <- horizon_data %>% dplyr::select(compname, hzname)

  # Step 1: Ensure `genhz` column exists in SSURGO horizon data (ssurgo_horizon_data)
  if (!"genhz" %in% colnames(ssurgo_horizon_data)) {
    ssurgo_horizon_data$genhz <- generalizeHz(
      ssurgo_horizon_data$hzname,
      new = c('O','A','B','C','Cr','R'),
      pattern = c('O', '^A','^B','^C','^Cr', '^R'),
      ordered = TRUE
    )
  }
  colnames(ssurgo_horizon_data)[colnames(ssurgo_horizon_data) == "compname"] <- "id"
  # Fetch OSD data for the list of soil component names
  s <- fetchOSD(unique(horizon_data$compname))

  # Check if the necessary fields are available in the horizon data
  if ("distinctness" %in% names(s@horizons) &&
      "top" %in% names(s@horizons) &&
      "bottom" %in% names(s@horizons) &&
      "hzname" %in% names(s@horizons) &&
      "id" %in% names(s@horizons)) {

    # Extract relevant fields: id, top, bottom, hzname, distinctness
    osd_horizon_data <- s@horizons[, c("id", "hzname", "distinctness")]

    osd_horizon_data$genhz <- generalizeHz(
      osd_horizon_data$hzname,
      new = c('O','A','B','C', 'Cr','R'),
      pattern = c('O', '^\\d*A','^\\d*B','^\\d*C', '^\\d*Cr', '^\\d*R'),
      ordered = TRUE
    )

    # Step 1: Identify missing genhz values in ssurgo_horizon_data that are not in osd_horizon_data$genhz
    missing_genhz <- setdiff(ssurgo_horizon_data$genhz, osd_horizon_data$genhz)

    # Step 2: If there are missing genhz values, create a data frame and append to osd_horizon_data
    if (length(missing_genhz) > 0) {
      # Create a new data frame with the missing genhz values
      missing_rows <- ssurgo_horizon_data[ssurgo_horizon_data$genhz %in% missing_genhz, c("id", "hzname", "genhz")]

      # Add placeholder values for the necessary fields (distinctness)
      missing_rows$distinctness <- NA

      # Append the missing rows to osd_horizon_data
      osd_horizon_data <- rbind(osd_horizon_data, missing_rows)

      # Print a message indicating the missing genhz values have been added
      message("Added missing genhz values: ", paste(missing_genhz, collapse = ", "))
    } else {
      message("No missing genhz values found.")
    }

    # Step 3: Infill missing distinctness values using the `infill_missing_distinctness` function
    osd_horizon_data <- infill_missing_distinctness(osd_horizon_data)

    # Convert the distinctness codes to offset values using hzDistinctnessCodeToOffset
    osd_horizon_data$bound_sd <- hzDistinctnessCodeToOffset(osd_horizon_data$distinctness)

    # Rename the columns for clarity
    names(osd_horizon_data) <- c("id", "hzname", "distinctness", "genhz", "bound_sd")

    return(osd_horizon_data)
  } else {
    stop("Necessary fields (distinctness, id, top, bottom, hzname) are not available in the OSD data.")
  }
}


infill_missing_distinctness <- function(horizon_data) {
  #' Infill Missing Distinctness Values for Horizons
  #'
  #' This function fills in missing `distinctness` values for common horizon names (O, A, B, C, R, Cr)
  #' based on typical boundary characteristics observed in the field. It assigns reasonable default
  #' values for common horizon types.
  #'
  #' ## Default Boundary Distinctness Assignments:
  #' - **O Horizons (Organic Layers):** Diffuse
  #' - **A Horizons (Topsoil):** Clear
  #' - **B Horizons (Subsoil):** Gradual
  #' - **C Horizons (Parent Material):** Gradual
  #' - **Cr Horizons (Weathered Bedrock):** Gradual
  #' - **R Horizons (Bedrock):** Abrupt
  #'
  #' @param horizon_data A data frame containing horizon data. It must include `hzname` and `distinctness` columns.
  #'
  #' @return A data frame with missing `distinctness` values infilled based on the horizon name or generalized horizon group.

  # List of default boundary distinctness values for common horizons
  default_distinctness <- list(
    "O" = "diffuse",      # Organic horizon
    "A" = "clear",        # Topsoil
    "B" = "gradual",      # Subsoil
    "C" = "gradual",      # Parent material
    "Cr" = "gradual",     # Weathered bedrock
    "R" = "abrupt"        # Bedrock
  )

  # Assign a generalized horizon group (genhz) if it's not already present
  if (!"genhz" %in% colnames(horizon_data)) {
    horizon_data$genhz <- generalizeHz(
      horizon_data$hzname,
      new = c('O','A','B','C','Cr','R'),
      pattern = c('O', '^A','^B','^C','^Cr', '^R'),
      ordered = TRUE
    )
  }

  # Function to assign default distinctness based on horizon name or generalized horizon group
  fill_default_distinctness <- function(hzname, genhz, distinctness) {
    if (!is.na(distinctness)) {
      return(distinctness)  # If distinctness exists, return it
    }

    # Use default distinctness based on hzname or genhz
    if (!is.na(genhz) && genhz %in% names(default_distinctness)) {
      return(default_distinctness[[genhz]])
    } else if (!is.na(hzname) && hzname %in% names(default_distinctness)) {
      return(default_distinctness[[hzname]])
    } else {
      return(NA)  # Return NA if no match is found
    }
  }

  # Apply the default distinctness infilling
  horizon_data$distinctness <- mapply(
    fill_default_distinctness,
    horizon_data$hzname,
    horizon_data$genhz,
    horizon_data$distinctness
  )

  return(horizon_data)
}


# Function to infill missing values with representative values ± 2, ensuring no negative values
infill_missing_depth_variability <- function(horizon_data) {
  # Infill missing top low values (hzdept_l) with hzdept_r - 2
  horizon_data$hzdept_l <- ifelse(is.na(horizon_data$hzdept_l),
                                  pmax(horizon_data$hzdept_r - 2, 0),  # Ensure no negative values
                                  horizon_data$hzdept_l)

  # Infill missing top high values (hzdept_h) with hzdept_r + 2
  horizon_data$hzdept_h <- ifelse(is.na(horizon_data$hzdept_h),
                                  pmax(horizon_data$hzdept_r + 2, 0),
                                  horizon_data$hzdept_h)

  # Infill missing bottom low values (hzdepb_l) with hzdepb_r - 2
  horizon_data$hzdepb_l <- ifelse(is.na(horizon_data$hzdepb_l),
                                  pmax(horizon_data$hzdepb_r - 2, 0),
                                  horizon_data$hzdepb_l)

  # Infill missing bottom high values (hzdepb_h) with hzdepb_r + 2
  horizon_data$hzdepb_h <- ifelse(is.na(horizon_data$hzdepb_h),
                                  pmax(horizon_data$hzdepb_r + 2, 0),
                                  horizon_data$hzdepb_h)

  return(horizon_data)
}


# Top-down simulation function
simulate_soil_profile_top_down <- function(horizon_data) {
  n <- nrow(horizon_data)

  # Initialize a result data frame with hzname, top, and bottom columns
  result <- data.frame(
    hzname = horizon_data$hzname,
    top = rep(NA, n),
    bottom = rep(NA, n)
  )

  # Initialize the top of the first horizon
  result$top[1] <- 0

  for (i in 1:n) {
    if (i > 1) {
      # Set the top of the current horizon equal to the bottom of the previous horizon
      result$top[i] <- result$bottom[i - 1]
    }

    # Simulate bottom depth using tri_dist instead of simulate_triangular
    bottom_low <- max(horizon_data$hzdepb_l[i], result$top[i])
    bottom_high <- max(bottom_low, horizon_data$hzdepb_h[i])
    bottom_mode <- min(max(horizon_data$hzdepb_r[i], bottom_low), bottom_high)

    # Generate a single random sample using tri_dist
    bottom_depth <- tri_dist(n = 1, a = bottom_low, b = bottom_high, c = bottom_mode)
    bottom_depth <- round(bottom_depth)

    # Ensure bottom depth is at least the top depth
    if (bottom_depth < result$top[i]) {
      bottom_depth <- result$top[i]
    }
    result$bottom[i] <- bottom_depth
  }

  return(result)
}



# Bottom-up simulation function
simulate_soil_profile_bottom_up <- function(horizon_data) {
  n <- nrow(horizon_data)

  # Initialize a result data frame with hzname, top, and bottom columns
  result <- data.frame(
    hzname = horizon_data$hzname,
    top = rep(NA, n),
    bottom = rep(NA, n)
  )

  # Process horizons from bottom to top
  for (i in n:1) {
    # Simulate bottom depth using tri_dist
    bottom_low <- horizon_data$hzdepb_l[i]
    bottom_high <- horizon_data$hzdepb_h[i]
    bottom_mode <- horizon_data$hzdepb_r[i]

    # Generate a single random sample using tri_dist
    bottom_depth <- tri_dist(n = 1, a = bottom_low, b = bottom_high, c = bottom_mode)
    bottom_depth <- round(bottom_depth)

    # For the bottom horizon, simulate top depth
    if (i == n) {
      top_low <- horizon_data$hzdept_l[i]
      top_high <- min(horizon_data$hzdept_h[i], bottom_depth)
      top_mode <- min(max(horizon_data$hzdept_r[i], top_low), top_high)

      # Generate a single random sample using tri_dist
      top_depth <- tri_dist(n = 1, a = top_low, b = top_high, c = top_mode)
      top_depth <- round(top_depth)

      # Ensure top depth is not greater than bottom depth
      if (top_depth > bottom_depth) {
        top_depth <- bottom_depth
      }
    } else {
      # For other horizons, set bottom depth to top depth of the horizon below
      bottom_depth <- result$top[i + 1]

      # Simulate top depth using tri_dist
      top_low <- horizon_data$hzdept_l[i]
      top_high <- min(horizon_data$hzdept_h[i], bottom_depth)
      top_mode <- min(max(horizon_data$hzdept_r[i], top_low), top_high)

      # Generate a single random sample using tri_dist
      top_depth <- tri_dist(n = 1, a = top_low, b = top_high, c = top_mode)
      top_depth <- round(top_depth)

      # Ensure top depth is not greater than bottom depth
      if (top_depth > bottom_depth) {
        top_depth <- bottom_depth
      }
    }

    # Assign the simulated depths to the result
    result$top[i] <- top_depth
    result$bottom[i] <- bottom_depth
  }

  # Ensure the top depth of the uppermost horizon is 0
  result$top[1] <- 0

  # Adjust the bottom depth of the first horizon if necessary
  if (result$bottom[1] < result$top[1]) {
    result$bottom[1] <- result$top[1]
  }

  # Ensure no gaps between horizons
  for (i in 1:(n - 1)) {
    if (result$bottom[i] != result$top[i + 1]) {
      result$bottom[i] <- result$top[i + 1]
    }
  }

  return(result)
}



#' Simulate Soil Profile Thickness
#'
#' This function simulates the thickness of soil horizons for a given soil profile, using top-down and bottom-up
#' methods, and calculates the standard deviation of horizon thickness across multiple simulations. The function
#' returns summarized results that include the horizon name, representative top and bottom depths, and the standard
#' deviation of thickness for each horizon.
#'
#' @param horizon_data A dataframe containing horizon data. It must include the following columns:
#'        - `hzname`: The name of each horizon (e.g., A, B, C).
#'        - `hzdept_r`: The representative top depth of each horizon.
#'        - `hzdepb_r`: The representative bottom depth of each horizon.
#'        - `hzdept_l`: The low estimate for the top depth of each horizon.
#'        - `hzdept_h`: The high estimate for the top depth of each horizon.
#'        - `hzdepb_l`: The low estimate for the bottom depth of each horizon.
#'        - `hzdepb_h`: The high estimate for the bottom depth of each horizon.
#' @param n_simulations Integer, number of simulations to perform (default is 500).
#'
#' @return A dataframe containing the following columns:
#'         - `hzname`: The horizon name.
#'         - `top`: The representative top depth of the horizon (from `hzdept_r`).
#'         - `bottom`: The representative bottom depth of the horizon (from `hzdepb_r`).
#'         - `sd_thickness`: The standard deviation of the horizon thickness across the simulations.
#'
#' @examples
#' # Example horizon data (with representative top and bottom depths)
#' horizon_data <- data.frame(
#'   hzname = c("A", "B", "C"),
#'   hzdept_r = c(0, 20, 35),  # Representative top depths
#'   hzdepb_r = c(20, 35, 50),  # Representative bottom depths
#'   hzdept_l = c(0, 15, 30),   # Low top depths
#'   hzdept_h = c(0, 25, 40),   # High top depths
#'   hzdepb_l = c(15, 30, 45),  # Low bottom depths
#'   hzdepb_h = c(25, 40, 60)   # High bottom depths
#' )
#'
#' # Simulate 500 soil profiles with thickness calculations
#' set.seed(123)
#' n_simulations <- 500
#' summarized_results <- simulate_soil_profile_thickness(horizon_data, n_simulations)
#'
#' # View the summarized results
#' print(summarized_results)
#'
#' @export
simulate_soil_profile_thickness <- function(horizon_data, n_simulations = 500) {
  # Step 1: Infill missing values
  horizon_data <- infill_missing_depth_variability(horizon_data)

  results_list <- list()

  # Loop over the number of simulations
  for (sim in 1:n_simulations) {
    # Randomly select the simulation direction
    if (runif(1) < 0.5) {
      # Simulate top-down direction
      result <- simulate_soil_profile_top_down(horizon_data)
      method <- "Top-Down"
    } else {
      # Simulate bottom-up direction
      result <- simulate_soil_profile_bottom_up(horizon_data)
      method <- "Bottom-Up"
    }

    # Calculate the thickness for each horizon as the difference between bottom and top depths
    result$Thickness <- result$bottom - result$top

    # Add the method used and the simulation number to the result
    result$Method <- method
    result$Simulation <- sim

    # Store the result of the current simulation in the results list
    results_list[[sim]] <- result
  }

  # Combine all simulation results into a single data frame
  all_results <- do.call(rbind, results_list)

  # Reorder columns for clarity: Simulation, Method, hzname, top, bottom, Thickness
  all_results <- all_results[, c("Simulation", "Method", "hzname", "top", "bottom", "Thickness")]

  # Group the results by hzname and calculate the standard deviation of thickness for each horizon
  horizon_sd <- all_results %>%
    group_by(hzname) %>%
    summarise(
      thickness_sd = sd(Thickness) %>% round(2)         # Calculate the standard deviation of thickness
    )

  # Merge the representative top and bottom depths from horizon_data
  summarized_results <- merge(horizon_data[, c("hzname", "hzdept_r", "hzdepb_r")], horizon_sd,
                              by.x = "hzname", by.y = "hzname")

  # Rename hzdept_r and hzdepb_r to top and bottom for consistency
  summarized_results <- summarized_results %>%
    rename(
      top = hzdept_r,
      bottom = hzdepb_r
    )

  # Return the summarized results
  return(summarized_results)
}


simulate_and_perturb_soil_profiles <- function(soil_profile) {

  # Step 1: Extract horizons
  soil_profile@horizons <- horizons(soil_profile) %>% dplyr::select(mukey, id, cokey, compname, hzname, sim_comppct, hzdept_l, hzdept_r, hzdept_h, hzdepb_l, hzdepb_r, hzdepb_h, hzthk_l, sim_comppct)
  horizon_data <- horizons(soil_profile)
  horizon_data$genhz <- generalizeHz(
    horizon_data$hzname,
    new = c('O','A','B','C', 'Cr','R'),
    pattern = c('O', '^\\d*A','^\\d*B','^\\d*C', '^\\d*Cr', '^\\d*R'),
    ordered = TRUE
  )
  n_simulations = unique(horizon_data$sim_comppct) #set number of simulations to simulated comppct metric
  # Check if there's only one horizon (e.g., R horizon only)
  if (nrow(horizon_data) == 1) {
    message("Only one horizon in the profile. Skipping perturbation for this profile.")

    # Return the original soil profile replicated n_simulations times (since no perturbation is needed)
    replicated_profiles <- list()
    for (i in 1:n_simulations) {
      # Create a copy of the soil profile and assign a unique profile ID
      new_profile <- soil_profile
      profile_id(new_profile) <- paste0(profile_id(new_profile), "_sim_", i)
      replicated_profiles[[i]] <- new_profile
    }

    # Combine replicated profiles into a single SoilProfileCollection
    simulated_profiles <- aqp::combine(replicated_profiles)

    return(simulated_profiles)
  }

  # Step 1: Simulate soil profile thickness
  simulated_thickness <- simulate_soil_profile_thickness(horizon_data, n_simulations)
  simulated_thickness$genhz <- generalizeHz(
    simulated_thickness$hzname,
    new = c('O','A','B','C', 'Cr','R'),
    pattern = c('O', '^\\d*A','^\\d*B','^\\d*C', '^\\d*Cr', '^\\d*R'),
    ordered = TRUE
  )

  # Step 2: Query OSD distinctness and get bound_sd values
  distinctness_data <- query_osd_distinctness(horizon_data)
  bound_lut <- distinctness_data %>%
    group_by(id, genhz) %>%
    summarise(bound_sd = mean(bound_sd, na.rm = TRUE)) %>%  # Calculate the mean of 'bound_sd'
    ungroup()  # Remove grouping

  # Step 3: Merge the simulated thickness data and distinctness data
  combined_data <- merge(simulated_thickness, bound_lut, by = "genhz") %>% dplyr::arrange(top)

  # Step 4: Add sd_thickness (thickness.attr) and bound_sd (boundary.attr) to the soil profile collection
  horizon_data <-  horizon_data %>% left_join(combined_data %>% dplyr::select(hzname, thickness_sd, bound_sd), by = "hzname")
  soil_profile@horizons <- horizon_data %>% dplyr::select( mukey,id,cokey,compname,hzname,sim_comppct, hzdept_r, hzdepb_r, thickness_sd, bound_sd)

  # Step 3: First perturbation (horizon thickness)
  perturbed_profiles_thickness <- aqp::perturb(
    soil_profile,
    n = n_simulations,
    thickness.attr = "thickness_sd"
  )

  # Step 4: Second perturbation (boundary distinctness)
  # Initialize an empty list to store profiles
  list_of_profiles <- list()

  for (i in 1:length(perturbed_profiles_thickness)) {
    # Extract individual profile
    single_profile <- perturbed_profiles_thickness[i]

    # Perturb boundary depths
    perturbed_profile <- aqp::perturb(
      single_profile,
      n = 1,
      boundary.attr = "bound_sd"
    )

    # Store the perturbed profile
    list_of_profiles[[i]] <- perturbed_profile
  }

  # Combine all perturbed profiles into a SoilProfileCollection
  simulated_profiles <- aqp::combine(list_of_profiles)

  # Step 5: Identify out-of-range horizons and adjust depths
  adjust_out_of_range_profiles <- function(simulated_profiles, horizon_data) {
    sim_horizons <- horizons(simulated_profiles) %>% dplyr::select(id,mukey,cokey,compname, hzname, top = hzdept_r, bottom = hzdepb_r)

    # Merge simulated horizons with the original horizon data for comparison
    merged_data <- merge(
      sim_horizons,
      horizon_data[, c("hzname", "hzdept_l", "hzdepb_h")],
      by = "hzname",
      all.x = TRUE
    )

    # Adjust the top and bottom depths rowwise
    merged_data <- merged_data %>%
      rowwise() %>%  # Apply the function rowwise
      mutate(
        top = if_else(is.na(hzdept_l), top, pmax(top, hzdept_l)),        # If hzdept_l is NA, keep top unchanged
        bottom = if_else(is.na(hzdepb_h), bottom, pmin(bottom, hzdepb_h))  # If hzdepb_h is NA, keep bottom unchanged
      ) %>%
      ungroup()  # Ungroup to return to the default structure

    # Update the adjusted data back to the profile
    simulated_profiles_hzdata <- horizons(simulated_profiles)
    simulated_profiles_hzdata <- simulated_profiles_hzdata %>% dplyr::left_join(merged_data %>% dplyr::select(hzname,id,top,bottom), by=c('id'='id', 'hzname'='hzname'))
    simulated_profiles@horizons <- simulated_profiles_hzdata

    return(simulated_profiles)
  }

  # Apply the adjustment for out-of-range horizons
  simulated_profiles <- adjust_out_of_range_profiles(simulated_profiles, horizon_data)

  # Step 6: Return the adjusted simulated profiles
  return(simulated_profiles)
}

# Run soil profile depth simulations for all profiles in soil profile collection
simulate_profile_depths_by_collection <- function(soil_collection, seed = 123) {

  # Check if input is a valid SoilProfileCollection
  if (!inherits(soil_collection, "SoilProfileCollection")) {
    stop("Input must be a SoilProfileCollection.")
  }

  # Step 1: Set up soil profile data for reproducibility
  set.seed(seed)

  # Initialize an empty list to store all simulated profiles
  all_simulated_profiles <- list()

  # Step 2: Loop over each profile in the SoilProfileCollection
  for (i in seq_along(soil_collection)) {
    # Select each profile
    soil_profile <- soil_collection[i, ]

    # Step 3: Simulate n_simulations for the current soil profile
    simulated_profiles <- simulate_and_perturb_soil_profiles(soil_profile)

    # Append the simulated profiles to the list
    all_simulated_profiles[[i]] <- simulated_profiles
  }

  # Step 4: Combine all simulated profiles into a single SoilProfileCollection
  combined_simulated_profiles <- aqp::combine(all_simulated_profiles)

  # Step 5: Return the combined simulated profiles
  return(combined_simulated_profiles)
}


# Run soil profile depth simulations for all profiles in soil profile collection using 'future.apply' packages
simulate_profile_depths_by_collection_parallel <- function(soil_collection,  seed = 123, n_cores = 6) {

  # Set up the parallel backend
  plan(multisession, workers = n_cores)  # Use multisession for parallel processing

  # Check if input is a valid SoilProfileCollection
  if (!inherits(soil_collection, "SoilProfileCollection")) {
    stop("Input must be a SoilProfileCollection.")
  }

  # Set seed for reproducibility
  set.seed(seed)

  # Function to simulate a single profile
  simulate_single_profile <- function(i) {
    soil_profile <- soil_collection[i, ]  # Select the ith profile
    cat("Processing profile ID:", soil_profile$id, "\n")
    simulate_and_perturb_soil_profiles(soil_profile)  # Run the simulation
  }
  tryCatch({
  # Apply future_lapply to run the simulation in parallel for each profile
  all_simulated_profiles <- future_lapply(seq_along(soil_collection), simulate_single_profile)
  # Add logging


  # Combine all simulated profiles into a single SoilProfileCollection
  combined_simulated_profiles <- aqp::combine(all_simulated_profiles)

  # Return the combined simulated profiles
  return(combined_simulated_profiles)
  }, error = function(e) {
    cat("Error in profile ID:", soil_profile$id, "\n")
    cat("Error message:", e$message, "\n")
    # Optionally return NA or NULL to indicate failure
    return(NULL)
  })
}
