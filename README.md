# ![wombat-p pipelines](docs/images/wombatp-logo.png) 

## Introduction

**wombat-p pipelines** is a bioinformatics analysis pipeline that bundles different workflow for the analysis of label-free proteomics data with the purpose of comparison and benchmarking. It allows using files from the [proteomics metadata standard SDRF](https://github.com/bigbio/proteomics-metadata-standard).

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. We used one of the [nf-core](https://nf-co.re/) templates. 

<!-- TODO add continuous integration, preferably with statistics -->

## Pipeline summary

This work contains four major different workflows for the analysis or label-free proteomics data, originating from LC-MS experiments.
1. [MaxQuant](https://www.maxquant.org/) + [Normalyzer](http://normalyzer.immunoprot.lth.se/)
2. [SearchGui](http://compomics.github.io/projects/searchgui) + [Proline](https://www.profiproteomics.fr/proline/) + [PolySTest](https://bitbucket.org/veitveit/polystest)
3. [Compomics tools](http://compomics.github.io/) + [FlashLFQ](https://github.com/smith-chem-wisc/FlashLFQ) + [MSqRob](https://github.com/statOmics/MSqRob)
4. Tools from the [Transproteomic Pipeline](http://tools.proteomecenter.org/TPP.php) + [ROTS](https://bioconductor.org/packages/release/bioc/html/ROTS.html)

Initialization and parameterization of the workflows is based on tools from the [SDRF pipelines](https://github.com/bigbio/sdrf-pipelines), the [ThermoRawFileParser](http://compomics.github.io/projects/ThermoRawFileParser) with our own contributions and additional programs from the wombat-p organizaion [https://github.com/wombat-p/Utilities] as well as our [fork](https://github.com/elixir-proteomics-community/sdrf-pipelines). This includes setting a generalized set of data analysis parameters and the calculation of a multiple benchmarks.

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.10.3`)

2. Install [`Docker`](https://docs.docker.com/engine/installation/) or [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) (you can follow [this tutorial](https://singularity-tutorial.github.io/01-installation/)) _(you can use [`Conda`](https://conda.io/miniconda.html) both to install Nextflow itself and also to manage software within pipelines. Please only use it within pipelines as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_.

3. Download the pipeline and test it on a minimal dataset with a single command:

   ```console
   wget https://github.com/wombat-p/WOMBAT-Pipelines
   nextflow run main.nf3 -profile test,YOURPROFILE 
   ```
  Substitute `wget` with `curl`or alike.

4. Setup of system for running the analysis

Note that some form of configuration will be needed so that Nextflow knows how to fetch the required software. This is usually done in the form of a config profile (`YOURPROFILE` in the example command above). You can chain multiple config profiles in a comma-separated string.

> - The pipeline comes with config profiles called `docker`, `singularity` and `conda` which instruct the pipeline to use the named tool for software management. For example, `-profile test,docker`.
> - If you are using `docker`, your host system might need to set a parameters that should stop the mono-based programs from failing when running large data sets on multiple threads. For that please set `sudo sysctl -w vm.max_map_count=262144`
> - If you are using `singularity`, setting the [`NXF_SINGULARITY_CACHEDIR` or `singularity.cacheDir`](https://www.nextflow.io/docs/latest/singularity.html?#singularity-docker-hub) Nextflow options enables you to store and re-use the images from a central location for future pipeline runs.
> - If you are using `conda`, it is highly recommended to use the [`NXF_CONDA_CACHEDIR` or `conda.cacheDir`](https://www.nextflow.io/docs/latest/conda.html) settings to store the environments in a central location for future pipeline runs.

5. Start running your own analysis!

For a detailed explanation of the parameters, see below. Not all parameters are needed. 

```console
nextflow run nf-core/wombat --sdrf experimental_metadata.sdrf --fasta your_fasta_file.fasta --parameters your_parameters_yaml --raws thermo_raw_files --exp_design simple_experimental_design --workflow [other more specific parameters] -profile <docker/singularity/conda>
```

<!-- TODO TODO: ADD DIAGRAM -->

WOMBAT-P can run workflows using different (minimal) input, such as
_1) with SDRF file (raw files can be given as parameter or are download from the location specified in the sdrf file):_
a) SDRF file + fasta file
b) SDRF file + fasta file + experimental design file (will overwrite experimental design in sdrf)
c) SDRF file + fasta file + experimental design file + yaml parameter file (will overwrite default and sdrf parameters)

_2) without SDRF file:_
a) Raw files + fasta file + yaml parameter file
b) Raw file + fasta file + yaml parameter file + experimental design file

__-profile__ Set the profile and environment as described above

__--sdrf__ This is a tab-delimited file containng details about experimental design and can also include all paramters given in the `--parameters` yaml file. Several data sets on the PRIDE repository come with an sdrf file which is can then be found toghether with the other deposited files. For the PXD001819, this would be https://ftp.pride.ebi.ac.uk/pride/data/archive/2015/12/PXD001819/sdrf.tsv
See also the URL for SDRF files: https://github.com/bigbio/proteomics-metadata-standard/tree/master/annotated-projects and the description of the extended SDRF including data analysis parameters: https://github.com/bigbio/proteomics-metadata-standard/blob/master/sdrf-proteomics/Data-analysis-metadata.adoc

__--fasta__ You also need a fasta database to run the database search in the workflows. Standard databases can be downloaded from [UniProt](http:///uniprot.org)

__--parameters__ When deviating from the standard settings, use a yaml file containing new parameters settings. For more details about the different parameters and an example file, see https://github.com/bigbio/proteomics-metadata-standard/blob/master/sdrf-proteomics/Data-analysis-metadata.adoc
As not all of these parameters are avaailable for all workflows, see <!-- TODO Provide link -->  this table for an overview

__--raws__ Without given sdrf file containing the paths to the raw data files (Thermo raw format) or if you have the files already downloaded, specify the wildcard (e.g. "*" or "?") to access the files on your system. We recommend putting this parameters in 'single quotes' as you might run into an error when using wildcards.

__--exp_design__ An experimental design is automatically calculated from differences in the samples in the SDRF file. Alternatively, provide a tab-separated file with the two columns _raw_file_ and _exp_condition_. _raw_file_: raw file names without path. Incorrect or incomplete names will lead to errors. _exp_condition_: arbitrary names for the sample groups. Files with the same sample group name will be considered replicates. See [example](https://github.com/wombat-p/WOMBAT-Pipelines/blob/dev/docs/examples/pxd001819.txt).

__--workflow__ Instead of running 'all' workflows (default), run only one of 'maxquant', 'proline', 'compomics' or 'tpp'

__other parameters__:  

--comps (only maxquant workflow): Provide contrasts (specific comparisons) for the statistical tests. This is a list of comma-separated group names, e.g. "B-A,C-A" when having the three sample groups A, B and C 

--proline_engine (only proline workflow): Define the search engine for the database search. Can be one or multiples of "xtandem", "msgf", "ms-amanda", "tide", "comet", "myrimatch", "meta_morpheus" and "andromeda". Note that not all engines necessarily work well with each data set.

You can add other NextFlow parameters as described extensively [here](https://www.nextflow.io/docs/edge/cli.html?highlight=timeline#)


## Valid data analysis parameters per workflow

TODOOO

## Credits

nf-core/wombat was originally written by Veit Schwaemmle.

We thank the following people for their extensive assistance in the development of this pipeline:

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#wombat` channel](https://nfcore.slack.com/channels/wombat) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use  nf-core/wombat for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
