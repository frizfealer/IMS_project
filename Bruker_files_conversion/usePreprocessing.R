arg=(commandArgs(trailingOnly = TRUE));
writeLines( arg );
source("D:\\preprocessing.R");
MALDI_IMS_preprocessing( arg[1], arg[2] );