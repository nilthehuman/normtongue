#' Filter ultrasound tongue imaging data by confidence
#'
#' A convenience function to throw away all data point with zero confidence and
#' keep the rest. Other functions may call this one as a necessary preprocessing
#' step before doing their actual job.
#'
#' @param data Input dataframe of a set of 2D positions recorded with ultrasound.
#' @param column_confi The name of the dataframe column containing confidence values.
#' @param confi_cutoff Data points up to this limit will be excluded (0 by default).
#' @returns The same dataframe with the zero confidence data filtered out.
#' @export

filterbyconfidence <- function(data, column_confi='confi', confi_cutoff=0) {
  filtered_data <- data[data[,column_confi] > confi_cutoff,]
  filtered_data
}


#' Automatically detect the straightedge in an ultrasound image
#'
#' Find and separate a range of contiguous points at the right end of a series of 2D positions
#' that probably represent the straightedge or other level object used for calibration.
#'
#' @param data Dataframe including several 2D points of an object on the occlusal plane.
#' @param column_x The name of the dataframe column containing a data point's position
#'   along the horizontal axis.
#' @param column_y The name of the dataframe column containing a data point's position
#'   along the vertical axis.
#' @returns A subset of the imaging data belonging to the straightedge.
#' @export

findstraightedge <- function(data, column_x='X', column_y='Y') {
  data <- filterbyconfidence(data)
  data <- data[order(data[,column_x]),]
  slope_deg_prev <- NULL
  stderror_prev <- NULL
  data_tail <- NULL
  # trust the rightmost two points to always be on the straightedge by default
  for (num_points in 2:nrow(data)) {
    data_tail <- tail(data, n=num_points)
    linear_regression <- summary(lm(data_tail[,column_y] ~ data_tail[,column_x], data=data_tail))
    slope <- linear_regression$coefficients[2,1]
    slope_deg <- atan(slope) / pi * 180
    stderror <- linear_regression$coefficients[2,2]
    # use reasonable heuristics to detect a sharp "knuckle" in the data
    MAX_SLOPE_DEVIATION <- 4.0  # degrees
    if (!is.null(slope_deg_prev) && abs(slope_deg - slope_deg_prev) > MAX_SLOPE_DEVIATION) {
      break
    }
    # N.B. stderror is NaN when you only have two points
    if (num_points > 3 && !is.null(stderror_prev) && stderror > stderror_prev) {
      break
    }
    slope_deg_prev <- slope_deg
    stderror_prev <- stderror
  }
  # backtrack one point to the state before the break statement
  tail(data_tail, n=-1)
}


#' Rotate a single tongue contour
#'
#' Fix the angle of one set of 2D tongue contour data based on previously identified
#' straightedge data points.
#'
#' @param tongue_data Dataframe containing the input positions of the tongue contour.
#' @param occlusal_data Dataframe containing positions on or near the occlusal plane.
#' @param column_x The name of the dataframe column containing a data point's position
#'   along the horizontal axis.
#' @param column_y The name of the dataframe column containing a data point's position
#'   along the vertical axis.
#' @returns The data frame with the rotated position values.
#' @export

rottongue <- function(tongue_data, occlusal_data, column_x='X', column_y='Y') {
  # extract the angle of the straightedge to be cancelled out
  linear_regression <- summary(lm(Y ~ X, data=occlusal_data))
  slope <- linear_regression$coefficients[2,1]
  angle <- atan(slope)
  # translate the whole dataset to be centered around (0, 0)
  center_x <- mean(tongue_data[[column_x]])
  center_y <- mean(tongue_data[[column_y]])
  tongue_data[,column_x] <- tongue_data[,column_x] - center_x
  tongue_data[,column_y] <- tongue_data[,column_y] - center_y
  # rotate around origin:
  # x' = cos(a) * x - sin(a) * y
  # y' = cos(a) * y + sin(a) * x
  old_x <- tongue_data[[column_x]]
  old_y <- tongue_data[[column_y]]
  tongue_data[,column_x] <- cos(-angle) * old_x - sin(-angle) * old_y
  tongue_data[,column_y] <- cos(-angle) * old_y + sin(-angle) * old_x
  # translate back
  tongue_data[,column_x] <- tongue_data[,column_x] + center_x
  tongue_data[,column_y] <- tongue_data[,column_y] + center_y
  tongue_data
}


#' Make all tongue position data level
#'
#' Rotate several instances of tongue imaging data based on a reference image
#' showing the occlusal plane.
#'
#' @param data Dataframe containing the input positions of several tongue contours
#'   and a straightedge contour for each "batch" of tongue contours.
#' @param keys Array of column names to uniquely identify a batch of data belonging
#'   to the same occlusal contour in 'data'.
#' @param column_word The name of the column in 'data' that contains the target
#'   word produced by the speaker when the tongue contour was captured.
#' @param occlusal_word The dummy word in the 'column_word' that identifies
#'   straightedge data as opposed to regular data.
#' @return The data frame with all position values rotated
#' @export

normtongue <- function(data, keys, column_word, occlusal_word) {
  data <- filterbyconfidence(data)
  occlusal_rows <- data[data[,column_word] == occlusal_word,]
  occlusal_keys <- occlusal_rows[!duplicated(occlusal_rows[, keys]),]
  occlusal_keys <- occlusal_keys[keys]
  data_rotated <- data.frame()
  for (row in 1:nrow(occlusal_keys)) {
    print(sprintf("processing batch #%d out of %d total", row, nrow(occlusal_keys)))
    keys_row <- occlusal_keys[row,]
    # extract straightedge contour
    occlusal_contour <- merge(data, keys_row)
    occlusal_contour <- occlusal_contour[occlusal_contour[,column_word] == occlusal_word,]
    occlusal_plane <- findstraightedge(occlusal_contour)
    # find all contours belonging to this straightedge contour
    contours <- merge(data, keys_row)
    contours <- rottongue(contours, occlusal_plane)
    # accumulate results in data_rotated
    data_rotated <- rbind(data_rotated, contours)
  }
  data_rotated
}


# ---- ---- ---- ---- ---- ---- ---- ----
# Usage examples:

data <- read.csv(file='ultrasound_splines.csv')

# Example 1: basic usage of the main function.

data_rotated <- normtongue(data, c("spk", "batch"), "word", "STRAIGHTEDGE")

# We can take a look at the results for a specific word and speaker, for example
contours <- data_rotated[data_rotated$spk == 'N01' & data_rotated$word == 'bricks',]
ggplot(contours, aes(x=X, y=Y, color=batch)) + geom_point()

# We can export the results if we are satisfied
write.csv(data_rotated, file='output.csv')

# Example 2: some of our data cannot be treated like the rest for some reason,
# so we decide to use a different straightedge from the dame data frame for it.

# The data we want to use a different straightedge for
which_spk <- "F01"
which_batch <- "1"
# The custom straightedge we want to use for the exceptional data
straightedge_to_use_spk <- "M01"
straightedge_to_use_batch <- "2"

# Find and process the data points of the straightedge we want to use
straightedge <- findstraightedge(data[data$spk == straightedge_to_use_spk & data$batch == straightedge_to_use_batch & data$word == "STRAIGHTEDGE",])
# Rotate selected exceptional data according to the custom straightedge
rotated_exceptional <- rottongue(data[data$spk == which_spk & data$batch == which_batch,], straightedge)

# Rotate all of the data frame using everyone's own straightedges
rotated_full <- normtongue(data, c("spk", "batch"), "word", "STRAIGHTEDGE")
# Swap out the subset that we rotated previously with the custom straightedge
rotated_full <- subset(rotated_full, spk != which_spk | batch != which_batch)
rotated_full <- rbind(rotated_full, rotated_exceptional)

# Example 3: it turns out the first straightedge contour is the most accurate
# for all speakers, so let's only use those straightedges for all rotations.

# For this we define a modified function based on the original normtongue
# function, adapted to our purpose. We process all data on a per-speaker basis,
# so instead of a vector of keys, we only use the speakers' identifiers here
normtongue_only_first <- function(data, column_word, occlusal_word) {
  data <- filterbyconfidence(data)
  data_rotated <- data.frame()
  for (speaker in unique(data$spk)) {
    print(sprintf("processing the contours of speaker %s...", speaker))
    # Extract first straightedge contour
    occlusal_contour <- data[data$spk == speaker & data$batch == 1 & data[,column_word] == occlusal_word,]
    occlusal_plane <- findstraightedge(occlusal_contour)
    # Find all contours from this speaker
    contours <- data[data$spk == speaker,]
    contours <- rottongue(contours, occlusal_plane)
    # accumulate results in data_rotated
    data_rotated <- rbind(data_rotated, contours)
  }
  data_rotated
}

data_rotated <- normtongue_only_first(data, "word", "STRAIGHTEDGE")
