# How to install pipeline's Conda environment

> **IMPORTANT**: DO NOT USE A GLOBALLY INSTALLED CONDA. INSTALL YOUR OWN MINICONDA3 UNDER YOUR HOME DIRECTORY BY FOLLOWING THE INSTRUCTION HERE.

> **WARNING**: DO NOT INSTALL CONDA 4.7 OR LATER UNTIL WE FIX CONDA ENV INSTALLATION ISSUES. [4.6.14](https://repo.anaconda.com/miniconda/Miniconda3-4.6.14-Linux-x86_64.sh) IS RECOMMENDED.

> **WARNING**: DO NOT SKIP ANY OF THE FOLLOWING STEPS OR PIPELINE'S ENVIRONMENT WILL BE MESSED UP WITH YOUR LOCAL PYTHON/GLOBAL CONDA.

1) Download [Miniconda installer](https://repo.anaconda.com/miniconda/Miniconda3-4.6.14-Linux-x86_64.sh). Use default answers to all questions except for the first and last.
  ```bash
  $ wget https://repo.anaconda.com/miniconda/Miniconda3-4.6.14-Linux-x86_64.sh
  $ bash Miniconda3-4.6.14-Linux-x86_64.sh
  ```

  Type `yes` to the first question.
  ```bash
  Do you accept the license terms? [yes|no]
  [no] >>> yes
  ```

  Type `yes` to the last question.
  ```bash
  Do you wish the installer to initialize Miniconda3
  by running conda init? [yes|no]
  [no] >>> yes
  ```

2) **IMPORTANT**: Close your session and re-login. If you skip this step then pipeline's Conda environment will be messed up with base Conda environment.

3) **IMPORTANT**: Disable auto activation of base Conda enrivonment. 
  ```bash
  conda config --set auto_activate_base false
  ```

4) **IMPORTANT**: Close your session and re-login.

5) Install pipeline's Conda environment.

  ```bash
  $ scripts/uninstall_dependencies.sh  # uninstall it for clean-install
  $ scripts/install_dependencies.sh
  ```

6) Configure pipeline's python2 environment.

  ```bash
  $ conda activate encode-atac-seq-pipeline
  $ scripts/config_conda_env.sh
  ```

7) Update pipeline's Conda environment with pipeline's python wrappers. You need to repeat this step everytime you update (`git pull`) this pipeline.

  ```bash
  $ conda activate encode-atac-seq-pipeline
  $ scripts/update_conda_env.sh
  ```

8) To activate pipeline's Conda environment. Use `conda activate` instead of `source activate`.
  ```bash
  $ conda activate encode-atac-seq-pipeline
  ```
