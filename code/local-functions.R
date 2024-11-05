# Local Functions

# TODO:
#  * constant values?
#  * not enough values to fit a curve? dice() -> nls() | optim() -> coef()
#  * are exponential fits reasonable?

# Function to fit a decay function
fitDecayFunction <- function(z, p0, p) {

  # Exponential decay formula (solve for p)
  # res <- p0 * exp(-(z/p))

  return(res)
}

# Function to find dominant condition for a specific variable within components
dominantCondition <- function(i, v) {

  # Filter out miscellaneous areas
  i <- i[which(i$compkind != 'Miscellaneous area'), ]
  if(nrow(i) < 1) {
    return(NULL)
  }

  # Sum component percent by variable 'v'
  fm <- as.formula(sprintf("comppct_r ~ %s", v))
  a <- aggregate(fm, data = i, FUN = sum, na.rm = TRUE)

  # Get the most frequent class
  idx <- order(a[['comppct_r']], decreasing = TRUE)[1]

  # Create result data frame with the most frequent class
  res <- data.frame(
    mukey = i$mukey[1],
    cokey = i$cokey[1],
    compname = i$compname[1],
    source = i$source[1],
    v = a[[v]][idx],
    pct = a[['comppct_r']][idx]
  )

  # Rename columns
  names(res) <- c('mukey', 'cokey', 'compname', 'source', v, 'pct')

  return(res)
}

# Function to find the dominant value for a specific variable
dominantValue <- function(i, v) {

  # Filter out miscellaneous areas
  i <- i[which(i$compkind != 'Miscellaneous area'), ]
  if(nrow(i) < 1) {
    return(NULL)
  }

  # Find the component with the highest percentage
  idx <- order(i[['comppct_r']], decreasing = TRUE)[1]

  # Create result data frame
  res <- data.frame(
    mukey = i$mukey[1],
    v = i[[v]][idx],
    pct = i[['comppct_r']][idx]
  )

  # Rename columns
  names(res) <- c('mukey', v, 'pct')

  return(res)
}

# Function to build a soil parameter list from SSURGO/RSS component data
buildParameterList <- function(s, template = NULL) {

  # Create a parameter list or start from the provided template
  if(is.null(template)) {
    p <- list()
  } else {
    p <- template
  }

  # TODO: Decide how to handle organic horizons (possibly missing data).
  # Remove organic horizons
  s <- subsetHz(s, ! grepl('O', hzDesgn(s)))

  # Estimate soil depth
  .soildepth <- estimateSoilDepth(s)

  # Aggregate soil properties over soil depth
  a <- suppressMessages(
    slab(s, fm = ~ sandtotal_r + silttotal_r + claytotal_r, slab.structure = c(0, .soildepth), strict = FALSE, slab.fun = mean, na.rm = TRUE)
  )

  # Convert long format to wide format
  a.wide <- reshape2::dcast(a, top + bottom ~ variable, value.var = 'value')

  # Extract sand, silt, and clay percentages
  .clay <- a.wide$claytotal_r
  .sand <- a.wide$sandtotal_r
  .silt <- a.wide$silttotal_r

  # Ensure the sum of sand, silt, and clay is at most 100%
  if(.sand + .silt + .clay > 100) {
    .silt <- 100 - (.sand + .clay)
  }

  # Convert Ksat from um/s to m/d
  s$ksat_r <- s$ksat_r * 0.0864

  # Get Ksat from the first mineral horizon
  .ksat0 <- s[, , .FIRST]$ksat_r

  # TODO: Determine how to estimate the Ksat decay parameter.

  # Fill in the parameter list based on SSURGO/RSS data

  # Convert depth from cm to meters
  p$soil_depth <- .soildepth * 0.01

  # Set soil depth for the heat flux model (in meters)
  p$deltaZ <- .soildepth * 0.01

  # Set surface Ksat (in meters per day)
  p$Ksat_0 <- .ksat0

  # Set sand, silt, and clay as fractions
  p$sand <- .sand * 0.01
  p$silt <- .silt * 0.01
  p$clay <- .clay * 0.01

  # Return the parameter list
  return(p)
}

# Function to convert a soil parameter file into a list
soilParameterFileToList <- function(f) {

  # Read the file into a two-column data frame
  s <- read.table(f)

  # Assign column names and reorder columns
  names(s) <- c('value', 'parameter')
  s <- s[, c('parameter', 'value')]

  # Convert the data frame to a named list
  p <- as.list(setNames(s$value, s$parameter))

  return(p)
}

# Function to write the soil parameter list to a file
writeSoilParameterFile <- function(p, f = '') {

  # Get names of the parameters
  nm <- names(p)

  # Create a vector to hold the text lines
  textLines <- vector(mode = 'character', length = length(p))

  # Iterate through the parameters
  for(i in seq_along(p)) {
    .v <- p[[i]]  # value
    .n <- nm[i]   # name

    # Format each line as "[value] [parameter]"
    textLines[i] <- sprintf("%s %s", .v, .n)
  }

  # Write the text lines to the specified file
  cat(textLines, sep = '\n', file = f)
}
