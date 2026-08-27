


process gp_register_stats {
    label 'single_cpu'
    label 'small_mem'
    input:
        path alignments
        path gencoll
        val name
    output:
        path "output/*.txt", emit: "stats"
    script:
    """
    mkdir -p output
    echo "${alignments.join('\n')}" > align.mft
    gp_register_stats -nogenbank -ifmt seq-align -gc-assembly $gencoll -input-manifest align.mft -stats-output output/align_${name}_stats.txt -by-run -collated-by-query -include-query-length -means-only -omit-empty-stats -tracking-server NONE -tracking-user Username -tracking-password Locator -tracking-database NONE
    """
    stub:
    """         
    mkdir -p output
    touch output/align_${name}_stats.txt
    """
}
