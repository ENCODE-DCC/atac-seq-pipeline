## How to download genome database

1. Choose `GENOME` from `hg19`, `hg38`, `mm9` and `mm10` and specify a destination directory.
    ```bash
    $ bash genome/download_genome_data.sh [GENOME] [DESTINATION_DIR]
    ```
2. Find a TSV file on the destination directory and use it for `"atac.genome_tsv"` in your input JSON.

# How to build genome database

1. Following [this instruction](install_conda.md) to install Conda and dependencies.

2. Choose `[GENOME]` from `hg19`, `hg38`, `mm9` and `mm10` and specify a destination directory. This will take several hours. We recommend not to run this installer on a login node of your cluster. It will take >8GB memory and >2h time.
    ```bash
    $ conda activate encode-atac-seq-pipeline
    $ bash conda/build_genome_data.sh [GENOME] [DESTINATION_DIR]
    ```

3. Find a TSV file on the destination directory and use it for `"atac.genome_tsv"` in your input JSON.


## How to build genome database for your own genome

1. You can build your own genome database if your reference genome has one of the following file types.
   * `.fasta.gz`
   * `.fa.gz`
   * `.fasta.bz2`
   * `.fa.gz2`
   * `.2bit`

2. Get a URL for your reference genome. You may need to upload it to somewhere on the internet.

3. Get a URL for a gzipped blacklist BED file for your genome. If you don't have one then skip this step. An example blacklist for hg38 is [here](http://mitra.stanford.edu/kundaje/genome_data/hg38/hg38.blacklist.bed.gz).

4. Find the following lines in `conda/build_genome_data.sh` and modify it. Give a good name `[YOUR_OWN_GENOME]` for your genome.
    ```bash
    ...

    elif [[ $GENOME == "YOUR_OWN_GENOME" ]]; then
      REF_FA="URL_FOR_YOUR_FASTA_OR_2BIT"
      BLACKLIST= # leave it empty if you don't have it

    ...
    ```

5. Specify a destination directory for your genome database and run the installer. This will take several hours.
    ```bash
    $ bash conda/build_genome_data.sh [YOUR_OWN_GENOME] [DESTINATION_DIR]
    ```

6. Find a TSV file in the destination directory and use it for `"atac.genome_tsv"` in your input JSON.
