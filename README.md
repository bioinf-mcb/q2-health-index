# q2-health-index

QIIME 2 plugin for calculating the Health Index from microbiome data.

## About

This plugin is based on the Gut Microbiome Health Index (GMHI) created by Gupta et al. 2020.

*A Predictive Index for Health Status using Species-level Gut Microbiome Profiling*
Nature Communications (2020) https://www.nature.com/articles/s41467-020-18476-8

Vinod K. Gupta, Minsuk Kim, Utpal Baksh, Kevin Y. Cunningham, John M. Davis III, Konstantinos N. Lazaridis, Heidi Nelson, Nicholas Chia, & Jaeyun Sung

Paper's GitHub repository: https://github.com/jaeyunsung/GMHI_2020

## Installation

To install the most up to date version of the plugin:

- Install and activate conda environment with QIIME 2 (see [docs](https://docs.qiime2.org/2020.11/install/native/)), e.g. for Linux 64-bit:
```
wget https://data.qiime2.org/distro/core/qiime2-2020.11-py36-linux-conda.yml
conda env create -n qiime2-2020.11 --file qiime2-2020.11-py36-linux-conda.yml
rm qiime2-2020.11-py36-linux-conda.yml
source activate qiime2-2020.8
qiime --help
```
Note that the plugin was tested with `qiime2-2020.8` but please use the latest version.  

- Fetch the repository and go to main folder:
```
git clone https://github.com/bioinf-mcb/q2-health-index
cd q2-health-index
```
- Install plugin: `make install`
- Run `calculate-gmhi` action e.g. (note: merge all lines below into one!):
```
qiime health-index calculate-gmhi
--i-table q2_health_index/tests/data/input/abundances/4347_final_relative_abundances.qza
--m-metadata-file q2_health_index/tests/data/input/metadata/4347_final_metadata.tsv
--o-gmhi-results q2_health_index/tests/data/gmhi_output
--o-gmhi-plot q2_health_index/tests/data/gmhi_plot
--p-healthy-states Healthy
--p-non-healthy-states rest
--p-healthy-column phenotype
```

## Tutorials

This is a QIIME 2 plugin. For details on QIIME 2 and tutorials demonstrating how to use this plugin, see the [QIIME 2 documentation](https://qiime2.org/).

### How to use this plugin


## Contributing

QIIME 2 is an open-source project, and we are very interested in contributions from the community. Please see the [contributing guidelines](https://github.com/qiime2/q2-sample-classifier/blob/master/.github/CONTRIBUTING.md) if you would like to get involved.

## NOTES from Mehrbod

The plugin would super useful if it has a primary action which uses their paper's initial process to calculate the health and non-health prevalent taxa, which then  can be output as a metadata file. And then this metadata file can be fed into another action which then calculates the actual index values using their formula and outputs another metadata.

The benefit of this is then that it's resolution/pipeline agnostic, meaning you can use ASV, gOTUs, species, genus, proteins, whatever..., and you're not bound by those taxa that they discovered which may or may not be pipeline specific.
When we re-run their methods with our pipelines it's likely we'll have different taxa identified and so this plugin would be super convenient to just plug and go with new taxa.

## NOTES from Tomasz

Example usage:

* calculate Health Index from (any) data:
    `qiime health-index create-index --i-table FeatureTable.qza --m-metadata-file HealthySick.txt --o-index TSVTaxonomyFormat.qza`
* run Health Index:
    `qiime health-index classify-samples --i-table FeatureTable.qza --i-index TSVTaxonomyFormat.qza --o-visualization SampleHealthIndex.qzv`

