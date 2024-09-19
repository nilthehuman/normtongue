\name{normtongue}
\alias{normtongue}
\title{Level out ultrasound tongue contours automatically}
\description{
Ultrasound tongue imaging is a non-invasive technique used in phonetics
experiments to obtain data on the position and shape of a speaker's tongue at
certain points in time.  A thin, stiff object such as a straightedge may be
inserted into the speaker's mouth in order to provide an image of the occlusal
plane before the actual measurements, based on which the other data may be
oriented to the same plane.  This package provides a few basic, easy to use
functions to rotate all tongue image data into the plane indicated by the
straightedge.
}
\usage{
normtongue(dataframe, keys, word_column, straightedge_word)
}
\arguments{
  \item{dataframe}{the input data to be corrected.}
  \item{keys}{a vector of column names in the dataframe that identify sets of
    tongue contours belonging to the same occlusal plane contour.}
  \item{word_column}{the column name in the dataframe showing the target word
    produced by the speaker.}
  \item{straightedge_word}{special "dummy word" to identify rows that belong
    to occlusal plane contours.}
}
\value{
  A new dataframe containing all observations (with positive confidence values),
  including the occlusal plane observations, rotated so as to nullify the angle
  detected in the corresponding occlusal plane contours.
}
\examples{
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
}
