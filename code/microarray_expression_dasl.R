# This script contains all the modular functions and a pipeline function to
# read, filter, and analyze DASL expression data per SCU pipeline.

# A custom S4 class is created to handle the object through the pipeline
# To access the data within the S4 classes, use object@data_type

require(illuminaio)
require(limma)
require(org.Hs.eg.db)
require(scales)
require(DMRcate)
require(toolkit)
source('code/def_classes.R')
load('annot/il450k.rda')