
# NetCoMi <img src="man/figures/logo.png" align="right" height="200" />

[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![install with
bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](https://anaconda.org/bioconda/r-netcomi)
[![DOI](https://zenodo.org/badge/259906607.svg)](https://zenodo.org/badge/latestdoi/259906607)
[![Citation
Badge](https://api.juleskreuer.eu/citation-badge.php?doi=10.1093/bib/bbaa290)](https://juleskreuer.eu/citation-badge/)
[![DOI
paper](https://img.shields.io/badge/doi-10.1093/bib/bbaa290-yellow.svg)](https://doi.org/10.1093/bib/bbaa290)

NetCoMi (**Net**work **Co**nstruction and Comparison for **Mi**crobiome
Data) is an R package designed to facilitate the construction, analysis,
and comparison of networks tailored to microbial compositional data. It
implements a comprehensive workflow introduced in [Peschel et
al. (2020)](https://academic.oup.com/bib/advance-article/doi/10.1093/bib/bbaa290/6017455),
which guides users through each step of network generation and analysis
with a strong emphasis on reproducibility and computational efficiency.

With NetCoMi, users can construct microbial association or dissimilarity
networks directly from sequencing data, typically provided as a read
count matrix. The package includes a broad selection of methods for
handling zeros, normalizing data, computing associations between
microbial taxa, and sparsifying the resulting matrices. By offering
these components in a modular format, NetCoMi allows users to tailor the
workflow to their specific research needs, creating highly customizable
microbial networks.

The package supports both the construction, analysis, and visualization
of a **single network** and the **comparison of two networks** through
graphical and quantitative approaches, including statistical testing.
Additionally, NetCoMi offers the capability of constructing
**differential networks**, where only differentially associated taxa are
connected.

<img src="man/figures/soilrep_networks.png" width=100% />

> Exemplary network comparison using soil microbiome data ([‘soilrep’
> data from phyloseq
> package](https://github.com/joey711/phyloseq/blob/master/data/soilrep.RData)).
> Microbial associations are compared between the two experimantal
> settings ‘warming’ and ‘non-warming’ using the same layout in both
> groups.

## Website

Please visit [netcomi.de](https://netcomi.de/) for a complete reference.

## Installation

``` r
# Required packages
install.packages("devtools")
install.packages("BiocManager")

# Since two of NetCoMi's dependencies are only available on GitHub, 
# it is recommended to install them first:
devtools::install_github("zdk123/SpiecEasi")
devtools::install_github("GraceYoon/SPRING")

# Install NetCoMi
devtools::install_github("stefpeschel/NetCoMi", 
                         repos = c("https://cloud.r-project.org/",
                                   BiocManager::repositories()))
```

If there are any errors during installation, please install the missing
dependencies manually.

Packages that are optionally required in certain settings are not
installed together with NetCoMi. These can be installed automatically
using:

``` r
installNetCoMiPacks()
```

If not installed via `installNetCoMiPacks()`, the required package is
installed by the respective NetCoMi function when needed.

## Bioconda

Thanks to [daydream-boost](https://github.com/daydream-boost), NetCoMi
can also be installed from conda bioconda channel with

``` bash
# You can install an individual environment firstly with
# conda create -n NetCoMi
# conda activate NetCoMi
conda install -c bioconda -c conda-forge r-netcomi
```

## Development version

Everyone who wants to use new features not included in any releases is
invited to install NetCoMi’s development version:

``` r
devtools::install_github("stefpeschel/NetCoMi", 
                         ref = "develop",
                         repos = c("https://cloud.r-project.org/",
                                   BiocManager::repositories()))
```

Please check the
[NEWS](https://github.com/stefpeschel/NetCoMi/blob/develop/NEWS.md)
document for features implemented on develop branch.

## References

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-peschel2020netcomi" class="csl-entry">

Peschel, Stefanie, Christian L Müller, Erika von Mutius, Anne-Laure
Boulesteix, and Martin Depner. 2020. “<span class="nocase">NetCoMi:
network construction and comparison for microbiome data in R</span>.”
*Briefings in Bioinformatics* 22 (4): bbaa290.
<https://doi.org/10.1093/bib/bbaa290>.

</div>

</div>
