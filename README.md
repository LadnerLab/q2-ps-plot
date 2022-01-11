# q2-ps-plot
Qiime2 Plug-in for the creation of visualizations from PepSIRF outputs.

https://ladnerlab.github.io/pepsirf-q2-plugin-docs/

### Table of Contents
- [Installation](https://github.com/LadnerLab/q2-ps-plot#installation)

- [Tutorial](https://github.com/LadnerLab/q2-ps-plot#tutorial)

- [Usage](https://github.com/LadnerLab/q2-ps-plot#usage)

## Installation

### Qiime2 Installation:
Visit: https://docs.qiime2.org/2021.8/install/ for intallation documentation on Qiime2

### PepSIRF Installation:
Visit: https://github.com/LadnerLab/PepSIRF for installation documentation on PepSIRF

### q2-pepsirf1 Installation:
Visit: https://github.com/LadnerLab/q2-pepsirf1 for installation documentation on q2-pepsirf1

### q2-ps-plot Installation:
#### Dependencies:
- `qiime2-2021.11 +`
- `PepSIRF`
- `q2-pepsirf1`
- `altair`
- `altair_saver`
- `jsonschema 3.2`

#### Directions:

Visit: https://ladnerlab.github.io/pepsirf-q2-plugin-docs/Pepsirf_Plugins/q2-ps-plot/ for installation documentation on q2-ps-plot.

#### Update ps-plot:

Visit: https://ladnerlab.github.io/pepsirf-q2-plugin-docs/Pepsirf_Plugins/q2-ps-plot/#updating for updating documentation on q2-ps-plot.

## Tutorial
File inputs for tutorial located in q2-ps-plot/example

### Running zenrich
Download the following files:

- IM0032-pA_PV1_subset_CS.qza
- IM0032-pA_PV1_subset_Z-HDI95.qza
- pairs_source.tsv

Run the following command:

```
qiime ps-plot zenrich --i-data IM0032-pA_PV1_subset_CS.qza \
--i-zscores IM0032-pA_PV1_subset_Z-HDI95.qza \
--m-source-file pairs_source.tsv \
--m-source-column Source \
--p-negative-controls SB_pA_A SB_pA_B SB_pA_D \
--o-visualization testing_zenrich
```

*If you cannot run PepSIRF by simply calling 'pepsirf' in the command line, you need to add `--p-pepsirf-binary` to the above command, followed by how you can call PepSIRF on your system*

Once ps-plot has finished running you should see: `Saved Visualization to: testing_zenrich.qzv`

You can view this visualization by dropping the `testing_zenrich.qzv` file into https://view.qiime2.org/.

### Running zenrich-tsv
Download the following files:

- IM0032-pA_PV1_subset_CS.tsv
- IM0032-pA_PV1_subset_Z-HDI95.tsv
- pairs_source.tsv

Run the following command:

```
qiime ps-plot zenrich-tsv --p-data-filepath IM0032-pA_PV1_subset_CS.tsv \
--p-zscores-filepath IM0032-pA_PV1_subset_Z-HDI95.tsv \
--m-source-file pairs_source.tsv \
--m-source-column Source \
--p-negative-controls SB_pA_A SB_pA_B SB_pA_D \
--o-zenrich-vis testing_zenrich_tsv
```

*If you cannot run PepSIRF by simply calling 'pepsirf' in the command line, you need to add `--p-pepsirf-binary` to the above command, followed by how you can call PepSIRF on your system*

Once ps-plot has finished running you should see: `Saved Visualization to: testing_zenrich_tsv.qzv`

You can view this visualization by dropping the `testing_zenrich_tsv.qzv` file into https://view.qiime2.org/.

### `More visualizers coming soon...`
