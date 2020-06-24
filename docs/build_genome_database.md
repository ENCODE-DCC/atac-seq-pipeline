## How to download genome database

1. Choose `GENOME` from `hg19`, `hg38`, `mm9` and `mm10` and specify a destination directory.
    ```bash
    $ bash scripts/download_genome_data.sh [GENOME] [DESTINATION_DIR]
    ```
2. Find a TSV file on the destination directory and use it for `"atac.genome_tsv"` in your input JSON.

# How to build genome database

1. [Install Conda](https://conda.io/miniconda.html). Skip this if you already have equivalent Conda alternatives (Anaconda Python). Download and run the [installer](https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh). Agree to the license term by typing `yes`. It will ask you about the installation location. On Stanford clusters (Sherlock and SCG4), we recommend to install it outside of your `$HOME` directory since its filesystem is slow and has very limited space. At the end of the installation, choose `yes` to add Miniconda's binary to `$PATH` in your BASH startup script.
    ```bash
    $ wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    $ bash Miniconda3-latest-Linux-x86_64.sh
    ```

2. Install pipeline's Conda environment.
    ```bash
    $ bash scripts/uninstall_conda_env.sh  # to remove any existing pipeline env
    $ bash scripts/install_conda_env.sh
    ```

3. Choose `GENOME` from `hg19`, `hg38`, `mm9` and `mm10` and specify a destination directory. This will take several hours. We recommend not to run this installer on a login node of your cluster. It will take >8GB memory and >2h time.
    ```bash
    $ conda activate encode-atac-seq-pipeline
    $ bash scripts/build_genome_data.sh [GENOME] [DESTINATION_DIR]
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

3. Get a URL for a gzipped blacklist BED file for your genome. If you don't have one then skip this step. An example blacklist for hg38 is [here](https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz).

4. Find the following lines in `scripts/build_genome_data.sh` and modify them as follows. Give a good name `[YOUR_OWN_GENOME]` for your genome. For `MITO_CHR_NAME` use a correct mitochondrial chromosome name of your genome (e.g. `chrM` or `MT`). For `REGEX_BFILT_PEAK_CHR_NAME` Perl style regular expression must be used to keep regular chromosome names only in a blacklist filtered (`.bfilt.`) peaks files. This `.bfilt.` peak files are considered final peaks output of the pipeline and peaks BED files for genome browser tracks (`.bigBed` and `.hammock.gz`) are converted from these `.bfilt.` peaks files. Chromosome name filtering with `REGEX_BFILT_PEAK_CHR_NAME` will be done even without the blacklist itself.
    ```bash
    ...

    elif [[ $GENOME == "YOUR_OWN_GENOME" ]]; then
      # Perl style regular expression to keep regular chromosomes only.
      # this reg-ex will be applied to peaks after blacklist filtering (b-filt) with "grep -P".
      # so that b-filt peak file (.bfilt.*Peak.gz) will only have chromosomes matching with this pattern
      # this reg-ex will work even without a blacklist.
      # you will still be able to find a .bfilt. peak file
      REGEX_BFILT_PEAK_CHR_NAME="chr[\dXY]+"
      # mitochondrial chromosome name (e.g. chrM, MT)
      MITO_CHR_NAME="chrM"
      # URL for your reference FASTA (fasta, fasta.gz, fa, fa.gz, 2bit)
      REF_FA="https://some.where.com/your.genome.fa.gz"
      # 3-col blacklist BED file to filter out overlapping peaks from b-filt peak file (.bfilt.*Peak.gz file).
      # leave it empty if you don't have one
      BLACKLIST=
    ...
    ```

5. Specify a destination directory for your genome database and run the installer. This will take several hours.
    ```bash
    $ bash scripts/build_genome_data.sh [YOUR_OWN_GENOME] [DESTINATION_DIR]
    ```

6. Find a TSV file in the destination directory and use it for `"atac.genome_tsv"` in your input JSON.
