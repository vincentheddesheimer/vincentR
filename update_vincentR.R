library(devtools)
library(roxygen2)

rm(list = ls())

## Declare WDs

wd <- "~/Documents/GitHub/vincentR/"

## For some reasons some dependencies are not installed w/ the package
## This should fix this:

use_package("lfe")
use_package("stringdist")
use_package("stringr")
use_package("ggplot2")
use_package("dplyr")
use_package("fixest")
use_package("broom")
use_package("rdrobust")
use_package("clipr")
use_package("grDevices")
use_package("glue")
use_package("magick")
use_package("pdftools")
use_package("png")

document()

## PUSH BEFORE REINSTALLING

devtools::install_github("vincentheddesheimer/vincentR",
    upgrade = T, force = T, quiet = F
)
