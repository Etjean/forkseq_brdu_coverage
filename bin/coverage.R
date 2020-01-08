#!/usr/bin/env Rscript
# Version: 1.0.0
# Author: Etienne JEAN
library(docopt)
'Usage:
    coverage.R -f=FILE -r=FILE -o=FILE [options]

Mandatory arguments:
    -f FILE, --fasta FILE    Fasta input file
    -r FILE, --ratio FILE    Ratio input file
    -o FILE, --output FILE   Tibble output of coverage data in .RDS format

Parameters:
    --bin_size X             Bin size for coverage data (in bp) [default: 200]
    --brdu_thrs X            Threshold of mean proportion of BrdU to declare that BrdU/a track is present [default: 0.4]
    --min_read_length X      Minimum read length [default: 1000]
' -> doc
args = docopt(doc)
library(tidyverse)


######### Parameters #########
fa_path = args$fasta
fa_ratio_path = args$ratio
output_path = args$output
bin_size = as.numeric(args$bin_size)
brdu_thrs = as.numeric(args$brdu_thrs)
min_read_length = as.numeric(args$min_read_length)


######### Data #########
fa = read_lines(fa_path) %>% .[seq(1, length(.), by=2)]
fa_ratio = read_lines(fa_ratio_path) %>% .[seq(2, length(.), by=2)]


######### Coverage #########
# Get signals 
tibble(
    header = fa,
    ratio = fa_ratio
) %>% 
# extract info from header
mutate(read_id     = map_chr(header, ~ str_match(.x, "^.+\\/([^\\/\\s]+)")[2])) %>%
mutate(read_strand = map_chr(header, ~ str_match(.x, "'mapped_strand': '(.)'")[2])) %>%
mutate(read_start  = map_dbl(header, ~ str_match(.x, "'mapped_start': (\\d+)")[2] %>% as.numeric)) %>%
mutate(read_end    = map_dbl(header, ~ str_match(.x, "'mapped_end': (\\d+)")[2] %>% as.numeric)) %>%
mutate(read_chrm   = map_chr(header, ~ str_match(.x, "'mapped_chrom': '(\\S+)'")[2])) %>%
select(-header) %>%
# extract ratio & add genomic position
mutate(ratio = map(ratio, ~ str_split(.x, ' ')%>%unlist%>%as.numeric%>%na.omit)) %>%
filter(map_lgl(ratio, ~ length(.x)>min_read_length)) %>%   # Filter reads above read length threshold
mutate(ratio = case_when(
    read_strand == '+' ~ pmap(., function(read_start, ratio, ...) {tibble(x=read_start + 1:length(ratio), y=ratio)}),
    read_strand == '-' ~ pmap(., function(read_end, ratio, ...) {tibble(x=read_end - 1:length(ratio), y=ratio)})
)) %>%
select(-read_strand, -read_start, -read_end) %>%
# subsample the signal
mutate(ratio = map(
    ratio,
    ~ .x %>% slice(seq(1, nrow(.), by=10))
)) %>%
# bin the signal
unnest %>%
mutate(bin = round(x/bin_size)*bin_size) %>%
group_by(read_id, read_chrm, bin) %>%
arrange(read_id, bin) %>%
summarise(y = mean(y)) %>%
# detect brdu
group_by(read_id) %>%
mutate(is_brdu = y>brdu_thrs) %>%
# detect tracks (position of the tracks is defined as the max of the signal)
mutate(id_track = case_when(
    lag(y, default=0)<brdu_thrs & y>brdu_thrs ~ 1, 
    T ~ 0
)) %>%
mutate(id_track = case_when(
    y>brdu_thrs ~ cumsum(id_track), 
    T ~ 0
)) %>%
group_by(read_id, id_track) %>%
mutate(is_track = id_track>0 & y==max(y)) %>%
# summarise
group_by(read_chrm, bin) %>%
summarise(cov_nano = n_distinct(read_id), cov_brdu=sum(is_brdu), cov_track=sum(is_track)) -> coverage


######## Save ########
saveRDS(coverage, output_path)


