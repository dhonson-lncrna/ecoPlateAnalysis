# ecoPlateAnalysis
Scripts for analyzing and visualizing metabolism data from BioLog EcoPlates.

### Overview

This repo contains scripts for analyzing (`ecoplate_analysis.py`) and visualizing (`ecoplate_dataviz.py`) absorbance data from [Biolog EcoPlates](https://www.biolog.com/products/community-analysis-microplates/ecoplate/). The analysis script:

1. Combines data from all timepoints and samples into a tidy dataframe.
2. Blanks 590 values based on Water wells.
3. Averages replicates across the 96-well plate.
4. Performs a trapezoidal approximation of the area under the absorbance curve.
5. Performs principal component analysis on the trapezoidal approximations.

The data visualization script will generate some basic plots if executed as a script, but may be more useful as sample code for making more specific visualizations.

The importer assumes that the data were collected from a Cytation 5 and exported as Excel files with the data in a column format. Each Excel file should represent a timepoint and each sheet within the file should represent a sample. Names for each sheet in the Excel file are mapped to sample names using a JSON file. For example data files see `exampleAnalysis/data/` and for an example plate JSON see `exampleAnalysis/plateIDs.json`. The naming convention `NNh_anyOtherInfo.xlsx`, where `NN` represents the timepoint in hours, is required for the script to execute correctly.

All other parameters for the script are supplied by YAML files. The script can either be run with default settings, which requires only seven user-supplied parameter values, or with full user control which requires twenty one parameter values. If running without the default settings, consult the documentation in ecoplate_analysis.py for full descriptions of each parameter. An example of a default YAML file is located at `exampleAnalysis/defaultConfig.yaml`, and of the full YAML file at `exampleAnalysis/fullConfig.yaml.`


### Getting Started

This repo comes with a Conda environment that contains all packages necessary for executing the scripts at `resources/config.yaml`. I recommend using I recommend using [Mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html) for package management. After installing and activating Mamba, navigate to the home directory of this repo and run:

`mamba env create -f resources/config.yaml`

To activate the environment, run:

`conda activate ecoplate`

### Running the example

This repo comes with two datafiles for a sample analysis. Once the ecoplate conda environment is activated, simply execute:

`python ecoplate_analysis.py exampleAnalysis/defaultConfig.yaml`

If the pipeline executes correctly, the data will be deposited in a new directory called `defaultOutput`.

You can also run an example of the pipeline with all user-supplied parameters. In the example, the full output should be identical to the default pipeline. To run it, execute:

`python ecoplate_analysis.py exampleAnalysis/fullConfig.yaml`

If the pipeline executes correctly, the data will be deposited in a new directory called `fullOutput`.

To visualize the `defaultOutput` data, use the supplied config file by running:

`python ecoplate_dataviz.py exampleAnalysis/dataVizConfig.yaml`

If the script executes, plots should be saved in a new directory called `dataViz`. 

To see how to adapt these visualzations and rerun sections of the pipeline to explore the data more thoroughly, see the Example Exploratory Analysis. 

### Running your own data

The requirements for running your own data are:

1. A folder containing an Excel file for each timepoint. Each file should have the naming convention `NNh_anyOtherInfo.xlsx`, where `NN` represents the timepoint in hours.
2. A JSON mapping timepoint, to sheet names, to samples, as in `exampleAnalysis/plateIDs.json`.
3. A JSON mapping well position to metabolite, as in `resources/ecoPlate.json`. The JSON supplied with this repo is based on [this PDF](https://www.biolog.com/wp-content/uploads/2023/08/00A-012-Rev-F-EcoPlate-IFU.pdf).
4. A YAML file specifying parameters for the pipeline, as in `exampleAnalysis/defaultConfig.yaml` or `exampleAnalysis/fullConfig.yaml`.

