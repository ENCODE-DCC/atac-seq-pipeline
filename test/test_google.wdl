task pair_gen { # returns every pair of true replicate
        Int num_rep
        File genome_tsv
        String gensz = read_map(genome_tsv)['gensz']
        command {
                echo ${gensz} > out.txt
        }
        output {
                String pairs = num_rep
                File out = "out.txt"
        }
        runtime {
                disks: "local-disk 20 SSD"
        }
}

workflow test_wf {
        Int num_rep = 4
        String genome_tsv = "gs://atac-seq-pipeline-genome-data/mm10_google.tsv"
        call pair_gen { input : num_rep = num_rep, genome_tsv = genome_tsv }
}
