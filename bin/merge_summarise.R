#!/usr/bin/env Rscript
# Author: Etienne JEAN
# Version: 1.0.1
library(docopt)
"Usage:
    merge_summarise.R -o=FILE FILE...

Mandatory arguments:
    -o FILE, --output FILE      Output file
    
" -> doc
args = docopt(doc)
file_out=args$output
files_in=args$FILE
library(dplyr)


# Initialize all_data
all_data=NULL

# For all files
for (file in files_in) {
	# Read data
    data = readRDS(file)
    # Merge tables
    all_data = rbind(all_data, data)
}

all_data %>%
# Summarise
tibble::as_tibble() %>%
group_by(read_chrm, bin) %>% 
summarise(cov_nano = sum(cov_nano), cov_brdu=sum(cov_brdu), cov_track=sum(cov_track)) %>%
# Save data
saveRDS(file_out)
