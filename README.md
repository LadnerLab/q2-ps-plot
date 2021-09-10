# q2-ps-plot
Qiime2 Plug-in for the creation of visualizations from PepSIRF outputs.

## Installation

### Qiime2 Installation:
Visit: https://docs.qiime2.org/2021.8/install/ for intallation documentation on Qiime2

### PepSIRF Installation:
Visit: https://github.com/LadnerLab/PepSIRF for installation documentation on PepSIRF

### q2-ps-plot Installation:
#### Dependencies:
- `qiime2`
- `PepSIRF`
- `altair`

## Usage

### zenrich
| Optional/Required | Argument(s) | Description | Example | Default |
| :--: | :------: | --- | --- | -- |
| **Required** | `--i-data` | Featuretable[Normed] - FeatureTable containing normalized read counts of samples and peptides. First column header must be 'Sequence Name' as produced by pepsirf. | **--i-data** some_file.qza | N/A |
| **Required** | `--i-zscores` | Featuretable[Zscore] - FeatureTable containing z scores of the normalized read counts. Fist column header must be 'Sequence Name' as produced by pepsirf. | **--i-zscores** some_file.qza | N/A |
| *Optional* | `--i-negative-data` | Featuretable[Normed] - FeatureTable containing normalized read counts of controls and peptides. First column header must be 'Sequence Name' as produced by pepsirf. | **-i-negative-data** some_file.qza | None |
| **Required** | `--m-source-file` | Metadata file containing all sample names and their source groups. Used to create pairs tsv to run pepsirf enrich module. | **--m-source-file** some_file.tsv | N/A |
| **Required** | `--m-source-column` | Name of column to collect source info from. | **--m-source-column** Source | N/A |
| *Optional* | `--m-peptide-metadata-file` | Filename of file that contains peptide metadata related to data to be plotted. | **--m-peptide-metadata-file** some_file.tsv | None |
| **Required** | `--p-negative-controls` | Sample names of the negative controls to be used. | **--p-negative-controls** sample1 sample2 sample3 | N/A |
| *Optional* | `--p-tooltip` | List of title names found in the peptide metadata file to be added to the hover tooltip (Parameter is case sensitive).  'Peptide' and 'Zscores' will always be added to the list of titles, if peptide metadata is not provided just 'Peptide' and 'Zscores' will be shown. | **--p-tooltip** Species SpeciesID | **['Species', 'SpeciesID']** or None |
| *Optional* | `--p-step-z-thresh` | Integar to increment z-score thresholds. | **--p-step-z-thresh** 1 | **5** |
| *Optional* | `--p-upper-z-thresh` | Upper limit of z-score thresholds (non-inclusive). | **--p-upper-z-thresh** 10 | **30** |
| *Optional* | `--p-lower-z-thresh` | Lower limit of z-score thresholds (inclusive). | **--p-lower-z-thresh** 5 | **5** |
| *Optional* | `--p-pepsirf-binary` | The binary to call pepsirf on your system. Used to call pepsirf enrich module. | **--p-pepsirf-binary** pepsirf_executable | **'pepsirf'** |
| **Required** | `--o-visualization` | Visualization output name | **--o-visualization** file_name.qzv | N/A |
| *Optional* | `--output-dir` | Output unspecified results to a directory | **--output-dir** directory_name | **Current Working Directory** |

### zenrich_tsv
| Optional/Required | Argument(s) | Description | Example | Default |
| :--: | :------: | --- | --- | -- |
| **Required** | `--p-data-filepath` | Filepath of .tsv file containing normalized read counts of samples and peptides. First column header must be 'Sequence Name' as produced by pepsirf. | **--p-data-filepath** some_file.tsv | N/A |
| **Required** | `--p-zscores-filepath` | Filepath of .tsv file containing z scores of the normalized read counts. Fist column header must be 'Sequence Name' as produced by pepsirf. | **--p-zscores-filepath** some_file.tsv | N/A |
| *Optional* | `--p-negative-data-filepath` | Filepath of .tsv file containing normalized read counts of controls and peptides. First column header must be 'Sequence Name' as produced by pepsirf. | **-p-negative-data-filepath** some_file.tsv | None |
| **Required** | `--m-source-file` | Metadata file containing all sample names and their source groups. Used to create pairs tsv to run pepsirf enrich module. | **--m-source-file** some_file.tsv | N/A |
| **Required** | `--m-source-column` | Name of column to collect source info from. | **--m-source-column** Source | N/A |
| *Optional* | `--m-peptide-metadata-file` | Filename of file that contains peptide metadata related to data to be plotted. | **--m-peptide-metadata-file** some_file.tsv | None |
| **Required** | `--p-negative-controls` | Sample names of the negative controls to be used. | **--p-negative-controls** sample1 sample2 sample3 | N/A |
| *Optional* | `--p-tooltip` | List of title names found in the peptide metadata file to be added to the hover tooltip (Parameter is case sensitive).  'Peptide' and 'Zscores' will always be added to the list of titles, if peptide metadata is not provided just 'Peptide' and 'Zscores' will be shown. | **--p-tooltip** Species SpeciesID | **['Species', 'SpeciesID']** or None |
| *Optional* | `--p-step-z-thresh` | Integar to increment z-score thresholds. | **--p-step-z-thresh** 1 | **5** |
| *Optional* | `--p-upper-z-thresh` | Upper limit of z-score thresholds (non-inclusive). | **--p-upper-z-thresh** 10 | **30** |
| *Optional* | `--p-lower-z-thresh` | Lower limit of z-score thresholds (inclusive). | **--p-lower-z-thresh** 5 | **5** |
| *Optional* | `--p-pepsirf-binary` | The binary to call pepsirf on your system. Used to call pepsirf enrich module. | **--p-pepsirf-binary** pepsirf_executable | **'pepsirf'** |
| **Required** | `--o-zenrich-vis` | Visualization output name | **--o-zenrich-vis** file_name.qzv | N/A |
| *Optional* | `--output-dir` | Output unspecified results to a directory | **--output-dir** directory_name | **Current Working Directory** |

### `More visualizers coming soon...`
