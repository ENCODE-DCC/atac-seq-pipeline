# How to install pipeline's Conda environment

> **WARNING**: DO NOT INSTALL CONDA 4.7 UNTIL WE FIX CONDA ENV INSTALLATION ISSUES. [4.6.14](https://repo.anaconda.com/miniconda/Miniconda3-4.6.14-Linux-x86_64.sh) IS RECOMMENDED.

1) Download [Conda](https://repo.anaconda.com/miniconda/Miniconda3-4.6.14-Linux-x86_64.sh).
  ```bash
  $ wget https://repo.anaconda.com/miniconda/Miniconda3-4.6.14-Linux-x86_64.sh
  $ bash Miniconda3-4.6.14-Linux-x86_64.sh
  ```

2) Install Conda environment for pipeline.

  ```bash
  $ conda/install_dependencies.sh
  ```

3) Initialize Conda and re-login.

  ```bash
  $ conda init bash
  $ exit
  ```

4) Configure pipeline's python2 and python3 environments.

  ```bash
  $ conda/config_conda_env.sh
  $ conda/config_conda_env_py3.sh
  ```

5) Update pipeline's Conda environment with pipeline's python source code. You need to run this step everytime you update (`git pull`) this pipeline.

  ```bash
  $ conda/update_conda_env.sh
  ```

## How to download genome database

1. Choose `GENOME` from `hg19`, `hg38`, `mm9` and `mm10` and specify a destination directory.
    ```bash
    $ bash genome/download_genome_data.sh [GENOME] [DESTINATION_DIR]
    ```
2. Find a TSV file on the destination directory and use it for `"chip.genome_tsv"` in your input JSON.
