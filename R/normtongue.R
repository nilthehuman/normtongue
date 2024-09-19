#' Filter ultrasound tongue imaging data by confidence
#'
#' A convenience function to throw away all data point with zero confidence and
#' keep the rest. Other functions may call this one as a necessary preprocessing
#' step before doing their actual job.
#'
#' @param data Input dataframe of a set of 2D positions recorded with ultrasound.
#' @param confi_column The name of the dataframe column containing confidence values.
#' @param confi_cutoff Data points up to this limit will be excluded (0 by default).
#' @returns The same dataframe with the zero confidence data filtered out.
#' @export

filterbyconfidence <- function(data, confi_column='confi', confi_cutoff=0) {
  filtered_data <- data[data[,confi_column] > confi_cutoff,]
  filtered_data
}


#' Automatically detect the straightedge in an ultrasound image
#'
#' Find and separate a range of contiguous points at the right end of a series of 2D positions
#' that probably represent the straightedge or other level object used for calibration.
#'
#' @param data Dataframe including several 2D points of an object on the occlusal plane.
#' @param x_column The name of the dataframe column containing a data point's position
#'   along the horizontal axis.
#' @param y_column The name of the dataframe column containing a data point's position
#'   along the vertical axis.
#' @returns A subset of the imaging data belonging to the straightedge.
#' @export

findstraightedge <- function(data, x_column='X', y_column='Y') {
  data <- filterbyconfidence(data)
  data <- data[order(data[,x_column]),]
  prev_slope_deg <- NULL
  prev_stderror <- NULL
  data_tail <- NULL
  # trust the rightmost two points to always be on the straightedge by default
  for (num_points in 2:nrow(data)) {
    data_tail <- tail(data, n=num_points)
    linear_regression <- summary(lm(data_tail[,y_column] ~ data_tail[,x_column], data=data_tail))
    slope <- linear_regression$coefficients[2,1]
    slope_deg <- atan(slope) / pi * 180
    stderror <- linear_regression$coefficients[2,2]
    # use reasonable heuristics to detect a sharp "knuckle" in the data
    MAX_SLOPE_DEVIATION <- 4.0  # degrees
    if (!is.null(prev_slope_deg) && abs(slope_deg - prev_slope_deg) > MAX_SLOPE_DEVIATION) {
      break
    }
    # N.B. stderror is NaN when you only have two points
    if (num_points > 3 && !is.null(prev_stderror) && stderror > prev_stderror) {
      break
    }
    prev_slope_deg <- slope_deg
    prev_stderror <- stderror
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
#' @param x_column The name of the dataframe column containing a data point's position
#'   along the horizontal axis.
#' @param y_column The name of the dataframe column containing a data point's position
#'   along the vertical axis.
#' @returns The data frame with the rotated position values.
#' @export

rottongue <- function(tongue_data, occlusal_data, x_column='X', y_column='Y') {
  # extract the angle of the straightedge to be cancelled out
  linear_regression <- summary(lm(Y ~ X, data=occlusal_data))
  slope <- linear_regression$coefficients[2,1]
  angle <- atan(slope)
  # translate the whole dataset to be centered around (0, 0)
  center_x <- mean(tongue_data[[x_column]])
  center_y <- mean(tongue_data[[y_column]])
  tongue_data[,x_column] <- tongue_data[,x_column] - center_x
  tongue_data[,y_column] <- tongue_data[,y_column] - center_y
  # rotate around origin:
  # x' = cos(a) * x - sin(a) * y
  # y' = cos(a) * y + sin(a) * x
  old_x <- tongue_data[[x_column]]
  old_y <- tongue_data[[y_column]]
  tongue_data[,x_column] <- cos(-angle) * old_x - sin(-angle) * old_y
  tongue_data[,y_column] <- cos(-angle) * old_y + sin(-angle) * old_x
  # translate back
  tongue_data[,x_column] <- tongue_data[,x_column] + center_x
  tongue_data[,y_column] <- tongue_data[,y_column] + center_y
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
#' @param word_column The name of the column in 'data' that contains the target
#'   word produced by the speaker when the tongue contour was captured.
#' @param occlusal_word The dummy word in the 'word_column' that identifies
#'   straightedge data as opposed to regular data.
#' @return The data frame with all position values rotated
#' @export

normtongue <- function(data, keys, word_column, occlusal_word) {
  data <- filterbyconfidence(data)
  occlusal_rows <- data[data[,word_column] == occlusal_word,]
  occlusal_keys <- occlusal_rows[!duplicated(occlusal_rows[, keys]),]
  occlusal_keys <- occlusal_keys[keys]
  data_rotated <- data.frame()
  for (row in 1:nrow(occlusal_keys)) {
    print(sprintf("processing batch #%d out of %d total", row, nrow(occlusal_keys)))
    keys_row <- occlusal_keys[row,]
    # extract straightedge contour
    occlusal_contour <- merge(data, keys_row)
    occlusal_contour <- occlusal_contour[occlusal_contour[,word_column] == occlusal_word,]
    occlusal_plane <- findstraightedge(occlusal_contour)
    # find all contours belonging to this straightedge contour
    contours <- merge(data, keys_row)
    contours <- rottongue(contours, occlusal_plane)
    # accumulate results in data_rotated
    data_rotated <- rbind(data_rotated, contours)
  }
  data_rotated
}
