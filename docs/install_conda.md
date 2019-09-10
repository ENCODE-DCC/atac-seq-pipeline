# How to install pipeline's Conda environment

1) Install [Conda](https://docs.conda.io/en/latest/miniconda.html).

2) Install Conda environment for pipeline.

  ```bash
  $ export CONDA_RESTORE_FREE_CHANNEL=1  # for Conda >= 4.7
  $ conda/install_dependencies.sh
  ```

3) Initialize Conda and re-login.

  ```bash
  $ # do nothing for Conda < 4.6
  $ conda init bash  # for Conda >= 4.6
  $ exit
  ```

4) Configure pipeline's python2 and python3 environments.

  ```bash
  $ # source activate encode-atac-seq-pipeline  # for Conda < 4.6
  $ conda activate encode-atac-seq-pipeline  # for >= 4.6
  $ conda/config_conda_env.sh
  ```

5) Update pipeline's Conda environment with pipeline's python source code. You need to run this step everytime you update (`git pull`) this pipeline.

  ```bash
  $ # source activate encode-atac-seq-pipeline  # for Conda < 4.6
  $ conda activate encode-atac-seq-pipeline  # for >= 4.6
  $ conda/update_conda_env.sh
  ```
