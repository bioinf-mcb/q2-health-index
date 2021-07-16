# q2-health-index

QIIME 2 plugin for calculating the Health Index from microbiome data.

## About

This plugin is based on the Gut Microbiome Health Index (GMHI) created by [Gupta et al. 2020](https://www.nature.com/articles/s41467-020-18476-8).

See also the paper's [GitHub repository](https://github.com/jaeyunsung/GMHI_2020).

## Installation

To install the most up to date version of the plugin:

- Install and activate conda environment with QIIME 2 (see [docs](https://docs.qiime2.org/2020.11/install/native/)), e.g. for Linux 64-bit:
    ```
    wget https://data.qiime2.org/distro/core/qiime2-2021.4-py38-linux-conda.yml
    conda env create -n qiime2-2021.4 --file qiime2-2021.4-py38-linux-conda.yml
    rm qiime2-2021.4-py36-linux-conda.yml
    source activate qiime2-2021.4
    qiime --help
    ```
Note that the plugin was tested with `qiime2-2021.4` but please use the latest version if available.  

- Fetch the repository and go to main folder:
    ```
    git clone https://github.com/bioinf-mcb/q2-health-index
    cd q2-health-index
    ```
- Install plugin: `make install`
- Test plugin e.g.: `qiime health-index --help`
## Tutorials

This is a QIIME 2 plugin. For details on QIIME 2 and tutorials demonstrating how to use this plugin, see the [QIIME 2 documentation](https://qiime2.org/).

**Note:** in the examples below all paths are related to the main repository directory.

### Calculate GMHI

In order to compute the GMHI (as a `qza` artifact) you need to provide the abundance table (`qza` artifact of the type 
`FeatureTable[Frequency] or FeatureTable[RelativeFrequency]`) and output file name. 

- Example:
  ```
  qiime health-index calculate-gmhi \
  --i-table q2_health_index/tests/data/input/abundances/4347_final_relative_abundances.qza \
  --o-gmhi-results q2_health_index/tests/data/gmhi_output
  ```
### Calculate and visualize GMHI

In order to compute and visualize the GMHI (in the form of  `qza` and `qzv` artifacts) you need to provide the 
abundance table (`qza` artifact of the type `FeatureTable[Frequency] or FeatureTable[RelativeFrequency]`),
the metadata file (e.g. `tsv` file) which contains description of **all** samples in the abundance table
and output file names. 

- Example:
  ```
  qiime health-index calculate-gmhi-viz \
  --i-table q2_health_index/tests/data/input/abundances/4347_final_relative_abundances.qza \
  --m-metadata-file q2_health_index/tests/data/input/metadata/4347_final_metadata.tsv \
  --o-gmhi-results q2_health_index/tests/data/gmhi_output \
  --o-gmhi-plot q2_health_index/tests/data/gmhi_plot
  ```
  
## Contributing

QIIME 2 is an open-source project, and we are very interested in contributions from the community.  
Please see the [contributing guidelines](https://github.com/qiime2/q2-sample-classifier/blob/master/.github/CONTRIBUTING.md) if you would like to get involved.
