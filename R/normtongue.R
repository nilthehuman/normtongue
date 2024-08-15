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
#' @return The data frame with all position values rotated
#' @export

normtongue <- function(tongue_data, occlusal_data=NULL) {
  # TODO
}


# ---- ---- ---- ---- ---- ---- ---- ----
# Usage example:
library(ggplot2)

spl_data <- read.csv(file='/home/nil/R/spl_vegl.csv', fileEncoding='latin1')

spl_sample <- spl_data[spl_data$spk == 'F03' & spl_data$szó == 'vonalzó' & spl_data$sorrend == 1& spl_data$confi > 0,]

ggplot(spl_sample, aes(x=X, y=Y)) + geom_point() + ylim(40, 70) + geom_abline(intercept=intercept, slope=slope)

straightedge <- findstraightedge(spl_sample)
rotated_sample <- rottongue(spl_sample, straightedge)

ggplot(rotated_sample, aes(x=X, y=Y)) + geom_point() + ylim(40, 70) + ggplot(spl_sample, aes(x=X, y=Y))

teklanak <- merge(spl_sample, rotated_sample, all=TRUE)
odd_rows <- seq_len(nrow(teklanak)) %% 2
teklanak$color <- "red"
for (i in 1:19) {
  if (odd_rows[i] == 0) { teklanak$color[i] <- 'red' }
  if (odd_rows[i] == 1) { teklanak$color[i] <- 'black' }
}
teklanak$color[19] <- 'red'

ggplot(teklanak, aes(x=X, y=Y)) + geom_point(colour=teklanak$color) + ylim(50, 70)

spl_sample_tail <- head(spl_sample, n=5)
linear_regression <- summary(lm(Y ~ X, data=spl_sample_tail))
intercept <- linear_regression$coefficients[1,1]
slope     <- linear_regression$coefficients[2,1]
ggplot(spl_sample, aes(x=X, y=Y)) + geom_point() + ylim(40, 70) + geom_abline(intercept=intercept, slope=slope)
