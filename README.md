# artMS

___Analytical R Tools for Mass Spectrometry___

---

[![Travis build status](https://travis-ci.org/biodavidjm/artMS.svg?branch=master)](https://travis-ci.org/biodavidjm/artMS)
[![codecov](https://codecov.io/gh/biodavidjm/artMS/branch/master/graph/badge.svg)](https://codecov.io/gh/biodavidjm/artMS)


## Overview

`artMS` is an R package that provides a set of tools for the analysis and 
integration of large-scale proteomics datasets (mass-spectrometry-based) 
obtained using the popular proteomics software package 
[MaxQuant](http://www.biochem.mpg.de/5111795/maxquant).

`artMS` facilitates and simplifies the analysis on R by creating a 
configuration file where multiple options can be selected, including:

- Multiple quality control (QC) functions.
- Relative quantification using [MSstats](http://msstats.org/).
- Post-MSstats processing and integration.
- A set of miscellaneous functions to generate input files for other popular 
tools, including:
    - APMS: MIST and [CompPASS](http://besra.hms.harvard.edu/ipmsmsdbs/cgi-bin/tutorial.cgi)
    - APMS: [SAINTq](http://saint-apms.sourceforge.net/Main.html) and SAINT express
    - [Photon](https://github.com/jdrudolph/photon)
    - [Phosfate](http://phosfate.com/)

***Required input files***

`artMS` performs the full pack of analyses by taking the following files:

- `evidence.txt` file: MaxQuant's output
- `keys.txt` (tab-delimited) file describing the whole experimental designed 
(templates available)
- `contrast.txt` (tab-delimited) file with the comparisons between conditions 
that the user wished to be quantified (templates available).

Finally, among the most relevant feature of `artMS`, a configuration file 
(in `yaml` format) is used to specify all the set of available 
(generated by the user, templates available).


## How to install

Until availability in BioConductor, install as:

```
install.packages("devtools")
library(devtools)
install_github('biodavidjm/artMS', build_vignettes=TRUE)
```

## The `artMS` input files

Three basic tab-delimited files are required to performed the full pack of 
operations:

### `evidence.txt`

The output of the quantitative proteomics software package `MaxQuant`. 
It combines all the information about the identified peptides.

### `keys.txt`

It contains the experimental design. 
This file will be merged with the `evidence.txt`
file (see above) through the "Raw.file" column. 
Each raw file corresponds to a unique individual technical replicate / 
biological replicate / Condition / Run.

Example: (see `test/example-keys.txt`)

**RawFile**|**IsotopeLabelType**|**Condition**|**BioReplicate**|**Run**
:-----:|:-----:|:-----:|:-----:|:-----:
FU20170922-17|L|H1N1\_03H|H1N1\_03H-1|9
FU20170922-19|L|H1N1\_03H|H1N1\_03H-2|10
FU20170922-21|L|H1N1\_06H|H1N1\_06H-1|11
FU20170922-23|L|H1N1\_06H|H1N1\_06H-2|12
FU20170922-35|L|H1N1\_12H|H1N1\_12H-1|13
FU20170922-37|L|H1N1\_12H|H1N1\_12H-2|14
FU20170922-39|L|H1N1\_18H|H1N1\_18H-1|15
FU20170922-41|L|H1N1\_18H|H1N1\_18H-2|16
FU20170922-01|L|MOCK\_03H|MOCK\_03H-1|1
FU20170922-03|L|MOCK\_03H|MOCK\_03H-2|2
FU20170922-05|L|MOCK\_06H|MOCK\_06H-1|3
FU20170922-07|L|MOCK\_06H|MOCK\_06H-2|4
FU20170922-09|L|MOCK\_12H|MOCK\_12H-1|5
FU20170922-11|L|MOCK\_12H|MOCK\_12H-2|6
FU20170922-13|L|MOCK\_18H|MOCK\_18H-1|7
FU20170922-15|L|MOCK\_18H|MOCK\_18H-2|8

### `contrast.txt`

The comparisons between conditions to be quantified. The written
comparisons must follow the following consensus:

```
Condition_A-Condition_B_mutant
```

- The two conditions to be compared must be separated by a dash symbol (`-`)
- The condition on the left will take the positive log2FC sign 
(if it is more abundant) and the one on the right the negative log2FC 
(if it is more abundant)
- The only special character allowed for the condition names is the
underscore (`_`)

Example (see also `test/example-contrast.txt`):

```
H1N1_03H-MOCK_03H
H1N1_06H-MOCK_06H
H1N1_12H-MOCK_12H
H1N1_18H-MOCK_18H
```


## The configuration file (`.yaml`)

The configuration file in `yaml` format contains the details of most of analyses 
performed by `artMS`. Check the folder [`data-raw`](./data-raw/artms_config.yaml) for a sample configuration file depending on the experiment. It currently covers the quantification of protein abundance, phosphorylation (ph), ubiquitination (ub), and acetylation (ac).

*Although it seems complex, the default options work very well*. Anyway, a detailed explanation of every section is explained next:

The configuration (`yaml`) file contains the following sections:

### `files`

```
files :
  evidence : /path/to/the/evidence.txt
  keys : /path/to/the/keys.txt
  contrasts : /path/to/the/contrast.txt
  output : /path/to/the/output/results/results.txt
```

The file `path/name` of the required files. 

---

### `qc`

```
qc:
  enabled: 1 # 1 = yes; 0 = no
```
- `1` to run QC analysis
- `0` otherwise

---

### `data`

```
  enabled : 1 # 1 = yes; 0 = no
  fractions: 
    enabled : 0 # 1 for protein fractionation
  silac: 
    enabled : 0 # 1 for SILAC experiments
  filters: 
    enabled : 1
    contaminants : 1
    protein_groups : remove #remove, keep
    modifications : ab # PH, UB, AB, APMS
  sample_plots : 1 # correlation plots
```

Let's break it down `data`:

- `enabled`:

    - `1`: to pre-process the data provided in the *files* section.
    - `0`: won't process the data (and a pre-generated MSstats file will be expected)

- `fractions`:
Multiple fractionation or separation methods are often combined in proteomics 
to improve signal-to-noise and proteome coverage and to reduce interference
between peptides in quantitative proteomics.
    - `enabled : 1` is a fractionation dataset
    - `enabled : 0` no fractions

- `silac`:
    - `enabled : 1`: if the files belong to a SILAC experiment. See **Special case: SILAC** below for details
    - `enabled : 0`: it does not

- `filters`: the following filtering options are available
    -  `enabled : 1` Enables filtering
    -  `contaminants : 1` Removes contaminants (`CON__` and `REV__` labeled by MaxQuant)
    - `protein_groups : ` choose whether `remove` or `keep` protein groups
    -  `modifications : ` any of the proteomics experiments, `PH`, `UB`, or `AC` for posttranslational modifications, `AB` or `APMS` otherwise.

- `sample_plots` 
    - `1` Generate correlation plots
    - `0` otherwise

---

### `msstats`

```
msstats :
  enabled : 1
  msstats_input : 
  profilePlots : none # before, after, before-after, none
  normalization_method : equalizeMedians # globalStandards (include a reference protein(s) ), equalizeMedians, quantile, 0
  normalization_reference :  #should be a value in the Protein column
  summaryMethod : TMP # "TMP"(default) means Tukey's median polish, which is robust estimation method. "linear" uses linear mixed model. "logOfSum" conducts log2 (sum of intensities) per run.
  censoredInt : NA  # Missing values are censored or at random. 'NA' (default) assumes that all 'NA's in 'Intensity' column are censored. '0' uses zero intensities as censored intensity. In this case, NA intensities are missing at random. The output from Skyline should use '0'. Null assumes that all NA intensites are randomly missing.
  cutoffCensored : minFeature  # Cutoff value for censoring. only with censoredInt='NA' or '0'. Default is 'minFeature', which uses minimum value for each feature.'minFeatureNRun' uses the smallest between minimum value of corresponding feature and minimum value of corresponding run. 'minRun' uses minumum value for each run.
  MBimpute : 1 # only for summaryMethod="TMP" and censoredInt='NA' or '0'. TRUE (default) imputes 'NA' or '0' (depending on censoredInt option) by Accelated failure model. FALSE uses the values assigned by cutoffCensored.
  feature_subset: all # all|highQuality  : highQuality seems to be buggy right now
```

Let's break it down:

- `enabled : ` Choose `1` to run MSstats, `0` otherwise. To find out more about the meaning of all the options, please, check the MSstats documentation (`?MSstats`)
- `msstats_input :` blank if MSstats is going to be run. Provide the path to the previously generated `evidence-mss.txt` if already available.
- `profilePlots :` Several profile plots available. 
    * `before` plots only before normalization
    * `after` plots only after normalization
    * `before-after`: recommended, although computational expensive (time consuming)
    * `none` no normalization plots (convenient if time limitations)

- `normalization_method :` available options:
    - `equalizeMedians`
    - `globalStandards` if selected, specified the reference protein in `normalization_reference` (next)
    - `quantile`
    - `0`: no normalization (not recommended)
- `normalization_reference :` an UniProt id if `globalStandards` is chosen as the `normalization_method` (above)
- `summaryMethod :` TMP # "TMP"(default) means Tukey's median polish, which is robust estimation method. "linear" uses linear mixed model. "logOfSum" conducts log2 (sum of intensities) per run.
- `censoredInt :` 
    - `NA`  Missing values are censored or at random. 'NA' (default) assumes that all 'NA's in 'Intensity' column are censored. 
    - `0` uses zero intensities as censored intensity. In this case, NA intensities are missing at random. The output from Skyline should use `0`. Null assumes that all `NA` intensities are randomly missing.
- `cutoffCensored :` 
    - `minFeature` Cutoff value for censoring. Only with `censoredInt='NA'` or `0`. Default is 'minFeature', which uses minimum value for each feature.
    - `minFeatureNRun` uses the smallest between minimum value of corresponding feature and minimum value of corresponding run. 
    - `minRun` uses minimum value for each run.
- `MBimpute :` 
    - `TRUE` only for `summaryMethod="TMP"` and `censoredInt='NA'` or `0`. TRUE (default) imputes 'NA' or '0' (depending on censoredInt option) by Accelerated failure model. 
    - `FALSE` uses the values assigned by cutoffCensored.
- `feature_subset :` 
    - `all` : default
    - `highQuality`  : this option seems to be buggy right now

---

### output

```
output_extras :
  enabled : 1
  msstats_output : # if the -evidence-mss.txt is already available.
  annotate : 0 # 1|0 whether to annotate the proteins in the results or not
  species : HUMAN # if annotate = 1, provide specie name. It can use multiple species, but separate with a "-" eg. HUMAN-MOUSE-HIV-...
  annotation_dir : /path/to/the/files # Required if "annotate = 1"
  comparisons : all # or any grep expression that returns a subset of the contrasts file
  LFC : -1 1
  FDR : 0.05
  heatmap : 1 
  heatmap_cluster_cols : 0
  heatmap_display : log2FC #or pvalue
  volcano : 1
```

Extra actions to perform based on the MSstats results.

#### Special case: Protein fractionation

To handle protein fractionation experiments, two options need to be activated

1. The keys' file must contain and additional column named "`FractionKey`" with 
the information of fractions. For example:

**Raw.file**|**IsotopeLabelType**|**Condition**|**BioReplicate**|**Run**|**FractionKey**
:-----:|:-----:|:-----:|:-----:|:-----:|:-----:
S9524\_Fx1|L|AB|AB-1|1|1
S9524\_Fx2|L|AB|AB-1|1|2
S9524\_Fx3|L|AB|AB-1|1|3
S9524\_Fx4|L|AB|AB-1|1|4
S9524\_Fx5|L|AB|AB-1|1|5
S9524\_Fx6|L|AB|AB-1|1|6
S9524\_Fx7|L|AB|AB-1|1|7
S9524\_Fx8|L|AB|AB-1|1|8
S9524\_Fx9|L|AB|AB-1|1|9
S9524\_Fx10|L|AB|AB-1|1|10
S9525\_Fx1|L|AB|AB-2|2|1
S9525\_Fx2|L|AB|AB-2|2|2
S9525\_Fx3|L|AB|AB-2|2|3
S9525\_Fx4|L|AB|AB-2|2|4
S9525\_Fx5|L|AB|AB-2|2|5
S9525\_Fx6|L|AB|AB-2|2|6
S9525\_Fx7|L|AB|AB-2|2|7
S9525\_Fx8|L|AB|AB-2|2|8
S9525\_Fx9|L|AB|AB-2|2|9
S9525\_Fx10|L|AB|AB-2|2|10
S9526\_Fx1|L|AB|AB-3|3|1
S9526\_Fx2|L|AB|AB-3|3|2
S9526\_Fx3|L|AB|AB-3|3|3
S9526\_Fx4|L|AB|AB-3|3|4
S9526\_Fx5|L|AB|AB-3|3|5
S9526\_Fx6|L|AB|AB-3|3|6
S9526\_Fx7|L|AB|AB-3|3|7
S9526\_Fx8|L|AB|AB-3|3|8
S9526\_Fx9|L|AB|AB-3|3|9
S9526\_Fx10|L|AB|AB-3|3|10

Internally, the function `getMSstatsFormat` handles the key step 
(just a simple `sum` aggregation)

2. Enable *fractions* in the configuration file as follow:

```
fractions: 
  enabled : 1 # 1 for protein fractions, 0 otherwise
```

#### Special case: SILAC

One of the most widely used techniques that enable relative protein 
quantification is stable isotope labeling by amino acids in cell culture (SILAC). 
The keys file will capture the typical SILAC experiment. 
For example, let's show a SILAC experiment with two conditions, 
two biological replicates and two technical replicates:

**RawFile**|**IsotopeLabelType**|**Condition**|**BioReplicate**|**Run**
:-----:|:-----:|:-----:|:-----:|:-----:
QE20140321-01|H|iso|iso-1|1
QE20140321-02|H|iso|iso-1|2
QE20140321-04|L|iso|iso-2|3
QE20140321-05|L|iso|iso-2|4
QE20140321-01|L|iso\_M|iso\_M-1|1
QE20140321-02|L|iso\_M|iso\_M-1|2
QE20140321-04|H|iso\_M|iso\_M-2|3
QE20140321-05|H|iso\_M|iso\_M-2|4

It is also required to activate the *silac* option in the yaml file to be 
activated as follow:

```
silac: 
  enabled : 1 # 1 for SILAC experiments
```


## Tips

Do you need to remember the basics of markdown? [Check out this fantastic link](https://commonmark.org/help/tutorial/index.html).
