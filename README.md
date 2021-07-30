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

## CLI parameters description

### Calculate GMHI
**Usage:** `qiime health-index calculate-gmhi [OPTIONS]`  
GMHI calculates the gut microbiome health index for each sample in the abundance table. 

**Inputs:**  

`--i-table	ARTIFACT	FeatureTable[Frequency] or FeatureTable[RelativeFrequency]`  
Abundance table artifact on which GMHI will be computed.

**Parameters:**  

| Parameter   |  Type  |  Optional / required / default      |  Description |
|:-----|:-----:|:-------------:|:------|
| `--p-healthy-species-fp` | TEXT |  optional | Path to file with healthy species (taxonomy is based on MetaPhlAn 2). |
| `--p-non-healthy-species-fp` | TEXT |    optional   |   Path to file with non-healthy species (taxonomy is based on MetaPhlAn 2). |
| `--p-mh-prime`  | INTEGER | default: 7 |  Median from the top 1% healthy samples in training dataset (see Gupta et al. 2020 Methods section). |
| `--p-rel-thresh` | NUMBER  | default: 1e-05 | Median from the top 1% non-healthy samples in training dataset (see Gupta et al. 2020 Methods section).  |
| `--p-rel-thresh` | NUMBER | default: 1e-05 | Relative frequency based threshold for discarding insignificant OTU. |
| `--p-log-thresh` | NUMBER | default: 1e-05 | Normalization value for `log10` in the last step of GMHI calculation.  |

**Outputs:**

`--o-gmhi-results	ARTIFACT SampleData[AlphaDiversity]` Calculated GMHI in tabular form.

### Calculate and visualize GMHI
**Usage:** `qiime health-index calculate-gmhi-viz [OPTIONS]`  
GMHI calculates and visualize the gut microbiome health index for each sample in the abundance table. 

**Inputs:**  

`--i-table	ARTIFACT	FeatureTable[Frequency] or FeatureTable[RelativeFrequency]`  
Abundance table artifact on which GMHI will be computed.

**Parameters:**  

| Parameter   |  Type  |  Optional / required / default      |  Description |
|:-----|:-----:|:-------------:|:------|
| `--m-metadata-file METADATA` | METADATA |  required | Metadata used for visualization. |
| `--p-healthy-species-fp` | TEXT |  optional | Path to file with healthy species (taxonomy is based on MetaPhlAn 2). |
| `--p-non-healthy-species-fp` | TEXT |    optional   |   Path to file with non-healthy species (taxonomy is based on MetaPhlAn 2). |
| `--p-mh-prime`  | INTEGER | default: 7 |  Median from the top 1% healthy samples in training dataset (see Gupta et al. 2020 Methods section). |
| `--p-rel-thresh` | NUMBER  | default: 1e-05 | Median from the top 1% non-healthy samples in training dataset (see Gupta et al. 2020 Methods section).  |
| `--p-rel-thresh` | NUMBER | default: 1e-05 | Relative frequency based threshold for discarding insignificant OTU. |
| `--p-log-thresh` | NUMBER | default: 1e-05 | Normalization value for `log10` in the last step of GMHI calculation.  |

**Outputs:**

`--o-gmhi-results	ARTIFACT SampleData[AlphaDiversity]` Calculated GMHI in tabular form.  
`--o-gmhi-plot VISUALIZATION` Bar plot showing calculated GMHI distribution.

## Tutorials

This is a QIIME 2 plugin. For details on QIIME 2 see [documentation](https://docs.qiime2.org/2021.4/).

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
  
**Important:** feature table must contain at least one healthy and non-healthy species.

### Calculate and visualize GMHI

In order to compute and visualize the GMHI (in the form of  `qza` and `qzv` artifacts) you need to provide the 
abundance table (`qza` artifact of the type `FeatureTable[Frequency] or FeatureTable[RelativeFrequency]`),
the metadata file (e.g. `tsv` file) and output file names.

- Example:
  ```
  qiime health-index calculate-gmhi-viz \
  --i-table q2_health_index/tests/data/input/abundances/4347_final_relative_abundances.qza \
  --m-metadata-file q2_health_index/tests/data/input/metadata/4347_final_metadata.tsv \
  --o-gmhi-results q2_health_index/tests/data/gmhi_output \
  --o-gmhi-plot q2_health_index/tests/data/gmhi_plot
  ```

**Important:** must contain description of **all** samples in the abundance table.

The visualization is generated using the `alpha-group-significance` 
function from the `q2-diversity` plugin (i.e. nonparametric Kruskal–Wallis test for healthy/non-healthy group comparison). 
      
## Contributing

QIIME 2 is an open-source project, and we are very interested in contributions from the community.  
Please see the [contributing guidelines](https://github.com/qiime2/q2-sample-classifier/blob/master/.github/CONTRIBUTING.md) if you would like to get involved.
