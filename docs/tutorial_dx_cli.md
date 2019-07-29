# Tutorial for DNAnexus Platform (CLI)

All test samples and genome data are shared on our public DNAnexus project. You don't have to download any data for testing our pipeline on DNAnexus platform.

There are two methods to run our pipeline on DNAnexus.

1) Building your own DX workflow from `atac.wdl` with dxWDL (CLI)
2) [Using a pre-built DX workflow on our public DX project (Web UI)](tutorial_dx_web.md)

This document describes instruction for the item 1).

1. Sign up for a [DNAnexus account](https://platform.DNAnexus.com/register).

2. Create a new [DX project](https://platform.DNAnexus.com/projects) with name `[YOUR_PROJECT_NAME]` by clicking on "+New Project" on the top left.

3. Download dxWDL.
    ```bash
    $ cd
    $ wget https://github.com/DNAnexus/dxWDL/releases/download/0.77/dxWDL-0.77.jar
    $ chmod +rx dxWDL-0.77.jar
    ```

4. Git clone this pipeline.
    ```bash
    $ cd
    $ git clone https://github.com/ENCODE-DCC/atac-seq-pipeline
    ```

5. Move to pipeline's directory.
    ```bash
    $ cd atac-seq-pipeline
    ```

6. Choose an appropriate input for your project (AWS or Azure):
    * AWS
      ```bash
      $ INPUT=example_input_json/dx/ENCSR356KRQ_subsampled_dx.json
      ```
    * Azure
      ```bash
      $ INPUT=example_input_json/dx_azure/ENCSR356KRQ_subsampled_dx_azure.json
      ```

7. Compile `atac.wdl` with an input JSON for the SUBSAMPLED (1/400) paired-end sample of [ENCSR356KRQ](https://www.encodeproject.org/experiments/ENCSR356KRQ/).
    ```bash
    $ PROJECT=[YOUR_PROJECT_NAME]
    $ OUT_FOLDER=/test_sample_atac_ENCSR356KRQ_subsampled

    $ java -jar dxWDL-0.77.jar compile atac.wdl -project ${PROJECT} -f -folder ${OUT_FOLDER} -defaults ${INPUT} -extras dev/workflow_opts/docker.json
    ```

8. Go to DNAnexus [project page](https://platform.DNAnexus.com/projects) and click on your project.

9. Move to the directory `/test_sample_atac_ENCSR356KRQ_subsampled`.

10. You will find a DX workflow `atac` with all parameters pre-defined. Click on it. 

11. Specify an output directory by clicking "Workflow Actions" on the top right. Click on "Set output folder" and choose an output folder.

12. Click on "Run as Analysis..." and you will be automatically redirected to the "Monitor" tab.

13. It will take about an hour. You will be able to find all outputs on your output folder. Final QC report (`qc.html`)/JSON (`qc.json`) will be found on it.

14. See full specification for [input JSON file](input.md).
