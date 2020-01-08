#!/usr/bin/env nextflow
// Version: 1.0.0
// Author: Etienne JEAN

//////////// PARAMETERS ////////////
// Help display
params.help = false
if (params.size()==1 || params.help) {
    println """
Mandatory arguments:
    --fasta_ratio_input      Fasta input file
    --output                 Tibble output of coverage data in .RDS format

Parameters:
    --nb_reads_process       Nb of reads per process
    --bin_size               Bin size for coverage data (in bp) [default: 200]
    --brdu_thrs              Threshold of mean proportion of BrdU to declare that BrdU/a track is present [default: 0.4]
    --min_read_length        Minimum read length [default: 1000]
"""
    exit 1
}

// Parameters default values
params.fasta_ratio_input = ''
params.output = ''
params.nb_reads_process = 5000
params.bin_size = 200
params.brdu_thrs = 0.4
params.min_read_length = 1000


// Parameters check
if (params.fasta_ratio_input==''){
    println('ERROR : Missing the input or output files.')
    exit 1
}
if (!(params.output =~ /\//)) {
    println('ERROR : The output file must contain a path (at least ./).')
    exit 1
}



// Parameters extraction
out_dir = (params.output =~ /(.+\/)(.+)/)[0][1]
out_file = (params.output =~ /(.+\/)(.+)/)[0][2]



//////////// PIPELINE ////////////
// Pair sequence and ratio files, and split them into small chunks
reads = Channel
    .fromFilePairs(params.fasta_ratio_input, flat:true)
    .splitFasta(by: params.nb_reads_process, file:true, elem:[1,2])

// Compute coverage
process Coverage {
    errorStrategy 'retry'
    maxRetries 3
    input:
        set val(id), file(reads_fa), file(reads_ra) from reads
    output:
        file "coverage*" into coverage
    script:
        name = reads_ra.getBaseName()
        """
        coverage.R \
        -f $reads_fa \
        -r $reads_ra \
        -o coverage_${name}.rds \
        --bin_size $params.bin_size \
        --brdu_thrs $params.brdu_thrs \
        --min_read_length $params.min_read_length
        """
}


// Merge coverage tables
process Merge {
    publishDir path: out_dir, mode: 'move'
    input:
        file all_coverage from coverage.collect()
    output:
        file "*"
    script:
        """
        merge_summarise.R \
        -o $out_file \
        $all_coverage
        """
}

// // Summarise them
// process Summarise {
//     publishDir path: params.output, mode: 'move'
//     input:
//         file coverage_binded
//     output:
//         file '*'
//     script:
//         """
//         #!/usr/bin/env R
//         library(dplyr)

//         readRDS($coverage_binded) %>%
//         tibble::as_tibble() %>%
//         group_by(read_chrm, bin) %>% 
//         summarise(cov_nano = sum(cov_nano), cov_brdu=sum(cov_brdu), cov_track=sum(cov_track)) %>%
//         saveRDS('cov_summarised.rds')
//         """
// }
