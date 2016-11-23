# This script contains all the modular functions and a pipeline function to
# read, filter, and analyze DASL expression data per SCU pipeline.

# A custom S4 class is created to handle the object through the pipeline
# To access the data within the S4 classes, use object@data_type

library(beadarray)
library(limma)
library(org.Hs.eg.db)
library(scales)
library(toolkit)
source('code/def_classes.R')



