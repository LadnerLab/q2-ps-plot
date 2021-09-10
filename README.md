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
| **Required** | `--i-data` | Featuretable[Normed] | **--i-data** some_file.qza | N/A |
| **Required** | `--i-zscores` | Featuretable[Zscore] | **--i-zscores** some_file.qza | N/A |
| *Optional* | `--i-negative-data` | Featuretable[Normed] | **-i-negative-data** some_file.qza | None |
| *Optional* | `--p-step-z-thresh` | Integer | **--p-step-z-thresh** 1 | **5** |
| *Optional* | `--p-upper-z-thresh` | non-inclusive | **--p-upper-z-thresh** 10 | **30** |

### zenrich_tsv

### `More visualizers coming soon...`
