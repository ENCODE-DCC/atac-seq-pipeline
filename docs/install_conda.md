# How to install pipeline's Conda environment

If you do not have miniconda (or anaconda) installed, follow the instructions below in steps 1 - 4 to install miniconda.

**IF YOU ALREADY HAVE ANACONDA OR MINICONDA INSTALLED, SKIP TO STEP 5** 

1) Download [Miniconda installer](https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh). Use default answers to all questions except for the first and last.
  ```bash
  $ wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
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
  $ conda config --set auto_activate_base false
  ```

4) **IMPORTANT**: Close your session and re-login.

5) Install pipeline's Conda environment.

  ```bash
  $ bash scripts/uninstall_conda_env.sh  # uninstall it for clean-install
  $ bash scripts/install_conda_env.sh
  ```

> **WARNING**: DO NOT PROCEED TO RUN PIPELINES UNTIL YOU SEE THE FOLLOWING SUCCESS MESSAGE OR PIPELINE WILL NOT WORK.
  ```bash
  === All done successfully ===
  ```

6) Activate pipeline's Conda environment before running a pipeline.
  ```bash
  $ conda activate encode-atac-seq-pipeline

  $ caper run ...
  $ caper server ...
  ```
