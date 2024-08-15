# normtongue

normtongue is a small self-contained package developed for a very specific purpose:
post-processing a certain kind of data.  Ultrasound tongue imaging, a non-invasive
technique used by phoneticians to find out more about the position of a speaker's
tongue during articulation, tends to produce tilted images.  However if you insert
a stiff flat object like a straightedge (i.e. a ruler) into the speaker's mouth
you can get reference images indicating the position and orientation of the
speaker's occlusal plane.  Once you have a good estimate of the tilt angle of the
occlusal plane itself you can rotate back the actual tongue images accordingly.

This package helps you fix the orientation of your tongue imaging data quickly
and painlessly.

## Example

Here's an example of how to use this package:

``` r
library(normtongue)
# TODO actual script
```
