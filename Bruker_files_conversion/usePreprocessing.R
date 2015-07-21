arg=(commandArgs(trailingOnly = TRUE));
writeLines( arg );
source("preprocessing.R");
MALDI_IMS_preprocessing( arg[1], arg[2], 0, 0 );