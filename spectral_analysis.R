# Load necessary libraries
library(spatialwarnings)
library(jpeg)
library(png)

convert_rgb_into_gray_scale <- function(img_matrix){
  return(apply(img_matrix[,,1:3], c(1, 2), mean))
}

plot_radial_spectrum <- function(img_matrix){
  plot(rspectrum(img_matrix), type='l', xlab="Nombre d'onde",
       ylab="Puissance radiale")
}

# Load image
image_path <- "/figures/images/plants.png"
image <- readPNG(image_path)

plot_radial_spectrum(image)
