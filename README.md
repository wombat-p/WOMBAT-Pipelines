# WOMBAT-P Pipelines
<img src="docs/images/wombatp-logo.png" width="200"></img>

## Summary

**wombat-p pipelines** is a bioinformatics analysis pipeline that bundles different workflow for the analysis of label-free proteomics data with the purpose of comparison and benchmarking. It allows using files from the [proteomics metadata standard SDRF](https://github.com/bigbio/proteomics-metadata-standard).

It aims both for experienced end-users that want to test different workflows and configurations and developers that want to e.g. test a new software in a workflow setting. Contributions to this project are most welcome.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. We used one of the [nf-core](https://nf-co.re/) templates.

<!-- TODO add continuous integration, preferably with statistics -->

This work contains four major different workflows for the analysis or label-free proteomics data, originating from LC-MS experiments.

1. [MaxQuant](https://www.maxquant.org/) + [NormalyzerDE](https://normalyzerde.immunoprot.lth.se/)
2. [SearchGui](http://compomics.github.io/projects/searchgui) + [Proline](https://www.profiproteomics.fr/proline/) + [PolySTest](https://bitbucket.org/veitveit/polystest)
3. [Compomics tools](http://compomics.github.io/) + [FlashLFQ](https://github.com/smith-chem-wisc/FlashLFQ) + [MSqRob](https://github.com/statOmics/MSqRob)
4. Tools from the [Trans-Proteomic Pipeline](http://tools.proteomecenter.org/TPP.php) + [ROTS](https://bioconductor.org/packages/release/bioc/html/ROTS.html)

Initialization and parameterization of the workflows is based on tools from the [SDRF pipelines](https://github.com/bigbio/sdrf-pipelines), the [ThermoRawFileParser](http://compomics.github.io/projects/ThermoRawFileParser) with our own contributions and additional programs from the wombat-p organizaion [https://github.com/wombat-p/Utilities] as well as our [fork](https://github.com/elixir-proteomics-community/sdrf-pipelines). This includes setting a generalized set of data analysis parameters and the calculation of a multiple benchmarks.

## Usage

### Installation and testing

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=23.0.4`)

2. Install [`Docker`](https://docs.docker.com/engine/installation/) or [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) (you can follow [this tutorial](https://singularity-tutorial.github.io/01-installation/)) _(you can use [`Conda`](https://conda.io/miniconda.html) both to install Nextflow itself and also to manage software within pipelines. Please only use it within pipelines as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_.

3. Download the pipeline and test it on a minimal dataset with a single command:

   ```console
   wget https://github.com/wombat-p/WOMBAT-Pipelines
   nextflow run main.nf -profile test,YOURPROFILE
   ```

   Substitute `wget` with `curl`or alike.

### Configuration and execution

4. Setup of system for running the analysis

Note that some form of configuration will be needed so that Nextflow knows how to fetch the required software. This is usually done in the form of a config profile (`YOURPROFILE` in the example command above). You can chain multiple config profiles in a comma-separated string.

> - The pipeline comes with config profiles called `docker`, `singularity` and `conda` which instruct the pipeline to use the named tool for software management. For example, `-profile test,docker`.
> - If you are using `docker`, your host system might need to set a parameters that should stop the mono-based programs from failing when running large data sets on multiple threads. For that please set `sudo sysctl -w vm.max_map_count=262144`
> - If you are using `singularity`, setting the [`NXF_SINGULARITY_CACHEDIR` or `singularity.cacheDir`](https://www.nextflow.io/docs/latest/singularity.html?#singularity-docker-hub) Nextflow options enables you to store and re-use the images from a central location for future pipeline runs.
> - If you are using `conda`, it is highly recommended to use the [`NXF_CONDA_CACHEDIR` or `conda.cacheDir`](https://www.nextflow.io/docs/latest/conda.html) settings to store the environments in a central location for future pipeline runs.

5. Start running your own analysis!

For a detailed explanation of the parameters, see below. Not all parameters are needed.

```console
nextflow run main.nf --sdrf experimental_metadata.sdrf --fasta your_fasta_file.fasta --parameters your_parameters_yaml --raws thermo_raw_files --exp_design simple_experimental_design --workflow [other more specific parameters] -profile <docker/singularity/conda>
```

### Input options and parameters

WOMBAT-P can run workflows using different (minimal) input, such as
_1) with SDRF file (raw files can be given as parameter or are download from the location specified in the sdrf file):_
a) SDRF file + fasta file
b) SDRF file + fasta file + experimental design file (will overwrite experimental design in sdrf)
c) SDRF file + fasta file + experimental design file + yaml parameter file (will overwrite default and sdrf parameters)

_2) without SDRF file:_
a) Raw files + fasta file + yaml parameter file
b) Raw file + fasta file + yaml parameter file + experimental design file

**-profile** Set the profile and environment as described above

**--sdrf** This is a tab-delimited file containng details about experimental design and can also include all paramters given in the `--parameters` yaml file. Several data sets on the PRIDE repository come with an sdrf file which is can then be found toghether with the other deposited files. For the PXD001819, this would be https://ftp.pride.ebi.ac.uk/pride/data/archive/2015/12/PXD001819/sdrf.tsv
See also the URL for SDRF files: https://github.com/bigbio/proteomics-metadata-standard/tree/master/annotated-projects and the description of the extended SDRF including data analysis parameters: https://github.com/bigbio/proteomics-metadata-standard/blob/master/sdrf-proteomics/Data-analysis-metadata.adoc

**--fasta** You also need a fasta database to run the database search in the workflows. Standard databases can be downloaded from [UniProt](http:///uniprot.org)

**--parameters** When deviating from the standard settings, use a yaml file containing new parameters settings. For more details about the different parameters and an example file, see https://github.com/bigbio/proteomics-metadata-standard/blob/master/sdrf-proteomics/Data-analysis-metadata.adoc
As not all of these parameters are available for all workflows, see <!-- TODO Provide link --> this table for an overview

**--raws** Without given sdrf file containing the paths to the raw data files (Thermo raw format) or if you have the files already downloaded, specify the wildcard (e.g. "\*" or "?") to access the files on your system. We recommend putting this parameters in 'single quotes' as you might run into an error when using wildcards.

**--exp_design** An experimental design is automatically calculated from differences in the samples in the SDRF file. Alternatively, provide a tab-separated file with the five columns _raw_file_, _exp_condition_, _biorep_, _fraction_ and _techrep_.<br/>
_raw_file_: raw file names without path. Incorrect or incomplete names will lead to errors.<br/>
_exp_condition_: arbitrary names for the sample groups. Files with the same sample group name will be considered replicates.<br/>
_biorep_: biological replicate with numbering starting with 1.<br/>
_fraction_: fraction number with numbering starting with 1.<br/>
_techrep_: technical replicate with numbering starting with 1.<br/>
__Note:__ The numbering needs to be consistent _and_ each line needs to be unique in the combination of _exp_condition_, _biorep_, _fraction_ and _techrep_.
See [example](https://github.com/wombat-p/WOMBAT-Pipelines/blob/dev/docs/examples/pxd001819.txt).

**--workflow** Instead of running 'all' workflows (default), run only one of 'maxquant', 'proline', 'compomics' or 'tpp'

**other parameters**:

--comps (only maxquant workflow): Provide contrasts (specific comparisons) for the statistical tests. This is a list of comma-separated group names, e.g. "B-A,C-A" when having the three sample groups A, B and C

--proline_engine (only proline workflow): Define the search engine for the database search. Can be one or multiples of "xtandem", "msgf", "ms-amanda", "tide", "comet", "myrimatch", "meta_morpheus" and "andromeda". Note that not all engines necessarily work well with each data set.

You can add other NextFlow parameters as described extensively [here](https://www.nextflow.io/docs/edge/cli.html?highlight=timeline#)

### Valid data analysis parameters per workflow
The following parameters can be provided via the parameter yaml file (`--parameters` flag)

As not all parameters are available for each workflow, the last columns describe their applicability. Here, _TRUE_ means that the parameter is available and can be modified.

| parameter                | type     | sdrf name                                 | default                                      | maxquant                 | proline                  | compomics                | tpp                      |
| ------------------------ | -------- | ----------------------------------------- | -------------------------------------------- | ------------------------ | ------------------------ | ------------------------ | ------------------------ |
| fixed_mods               | ontology | modification parameters                   | NT=Carbamidomethyl;TA=C;MT=fixed;AC=UNIMOD:4 | TRUE                     | TRUE                     | TRUE                     | TRUE                     |
| variable mods            | ontology | modification parameters                   | NT=oxidation;MT=variable;TA=M;AC=UNIMOD:35   | TRUE                     | TRUE                     | TRUE                     | TRUE                     |
| precursor_mass_tolerance | string   | precursor mass tolerance                  | 30 ppm                                       | TRUE                     | TRUE                     | TRUE                     | TRUE                     |
| fragment_mass_tolerance  | string   | fragment mass tolerance                   | 0.05 Da                                      | TRUE                     | TRUE                     | TRUE                     | TRUE                     |
| enzyme                   | ontology | cleavage agent details                    | Trypsin                                      | TRUE                     | TRUE                     | TRUE                     | TRUE                     |
| fions                    | class    | forward ions                              | b                                            | FALSE                    | TRUE                     | TRUE                     | TRUE                     |
| rions                    | class    | reverse ions                              | y                                            | FALSE                    | TRUE                     | TRUE                     | TRUE                     |
| isotope_error_range      | integer  | isotope error range                       | 0                                            | FALSE                    | TRUE                     | TRUE                     | TRUE                     |
| add_decoys               | boolean  | add decoys                                | true                                         | FALSE                    | TRUE                     | TRUE                     | TRUE                     |
| num_hits                 | integer  | num peptide hits                          | 1                                            | FALSE                    | FALSE                    | FALSE                    | FALSE                    |
| allowed_miscleavages     | integer  | allowed miscleavages                      | 1                                            | TRUE                     | TRUE                     | TRUE                     | TRUE                     |
| min_precursor_charge     | integer  | minimum precursor charge                  | 2                                            | FALSE                    | TRUE                     | TRUE                     | TRUE                     |
| max_precursor_charge     | integer  | maximum precursor charge                  | 3                                            | TRUE                     | TRUE                     | TRUE                     | TRUE                     |
| min_peptide_length       | integer  | minimum peptide length                    | 8                                            | TRUE                     | TRUE                     | TRUE                     | TRUE                     |
| max_peptide_length       | integer  | maximum peptide length                    | 12                                           | FALSE                    | TRUE                     | TRUE                     | TRUE                     |
| max_mods                 | integer  | maximum allowed modification              | 4                                            | TRUE                     | TRUE                     | TRUE                     | TRUE                     |
| ident_fdr_psm            | float    | fdr on psm level                          | 0.01                                         | TRUE                     | TRUE                     | TRUE                     | TRUE                     |
| ident_fdr_peptide        | float    | fdr on peptide level                      | 0.01                                         | TRUE                     | TRUE                     | TRUE                     | TRUE                     |
| ident_fdr_protein        | float    | fdr on protein level                      | 0.01                                         | TRUE                     | TRUE                     | Not clear                | Not clear                |
| match_between_runs       | boolean  | run match between runs                    | true                                         | TRUE                     | FALSE                    | TRUE                     | Not available            |
| protein_inference        | class    | protein inference method                  | unique                                       | TRUE                     | FALSE                    | TRUE                     | TRUE                     |
| quantification_method    | class    | quantification method                     | intensity                                    | FALSE                    | FALSE                    | FALSE                    | FALSE                    |
| summarization_proteins   | class    | summarization of proteins method          | sum_abs                                      | FALSE                    | FALSE                    | FALSE                    | FALSE                    |
| min_num_peptides         | integer  | minimum number of peptides per protein    | 2                                            | TRUE                     | TRUE                     | TRUE                     | TRUE                     |
| summarization_psms       | class    | summarization of psms method              | sum_abs                                      | FALSE                    | FALSE                    | FALSE                    | FALSE                    |
| quant_transformation     | class    | transformation of quantitative values     | log                                          | FALSE                    | FALSE                    | FALSE                    | FALSE                    |
| normalization_method     | class    | normalization method                      | median                                       | TRUE                     | FALSE                    | FALSE                    | FALSE                    |
| run_statistics           | boolean  | run statistical tests                     | true                                         | TRUE                     | TRUE                     | TRUE                     | TRUE                     |
| fdr_method               | class    | method for correction of multiple testing | benjamini-hochberg                           | FALSE                    | FALSE                    | FALSE                    | FALSE                    |
| fdr_threshold            | float    | threshold for statistical test fdr        | 0.01                                         | By filtering the results | By filtering the results | By filtering the results | By filtering the results |

## Workflow output

Intermediate and final files are provided in the _results_ folder or the folder specified via the `outdir` parameter.

On top of the workflow-specific output, a standardized tabular format on both peptide (`stand_pep_quant_merged.csv`) and protein (`stand_prot_quant_merged.csv`) level is given.

For each of the workflows, WOMBAT-Pipelines calculated the same set of benchmarks for more systematic and thorough comparison between workflows and/or between different values of the data analysis parameters. For details about the benchmarks, see the following table:

| Category      | Aspect       | Subgroup       | Name                                 | Name in JSON file                  | Definition                                                                         | Value          |
| ------------- | ------------ | -------------- | ------------------------------------ | ---------------------------------- | ---------------------------------------------------------------------------------- | -------------- |
| Functionality | Traceability | Spectra        | Tracable spectra                     | TraceableSpectra                   | Results tracable to original spectra                                               | Y/N            |
| Functionality | Traceability | Spectra        | Universal spectrum identifiers       | UniversalSpectumIdentifiers        | Workflow generates USIs (Universal Spectrum Identifier)                            | Y/N            |
| Functionality | Traceability | Spectra        | Peptide to spectra                   | PeptideToSpectra                   | Corresponding spectrum numbers/ids available from peptide level                    | Y/N            |
| Functionality | Traceability | Spectra        | Protein to spectra                   | ProteinToSpectra                   | Corresponding spectrum numbers/ids available from protein level                    | Y/N            |
| Functionality | Traceability | File names     | Results to raw files                 | ResultsToRawFiles                  | Raw input file names preserved in tables on PSM/peptide/protein level              | Y/N            |
| Functionality | Traceability | File names     | Public raw files                     | PublicRawFiles                     | Raw files publicly available                                                       | Y/N            |
| Functionality | Traceability | Parameters     | Experimental design                  | ExperimentalDesign                 | Biological and technical replicates can be identified in results                   | Y/N            |
| Functionality | Performance  | Identification | PSM number                           | PSMNumber                          | Number of identified PSMs passing preset FDR                                       | Integer        |
| Functionality | Performance  | Identification | Peptide number                       | PeptideNumber                      | Number of uniquely identified peptide identifications passing preset FDR           | Integer        |
| Functionality | Performance  | Identification | Protein number                       | ProteinNumber                      | Number of uniquely identified protein identifications passing preset FDR           | Integer        |
| Functionality | Performance  | Identification | Protein group number                 | ProteinGroupNumber                 | Number of different protein groups passing preset FDR                              | Integer        |
| Functionality | Performance  | Identification | Peptide coverage                     | PeptideCoverage                    | Percentage of peptides identified in all samples                                   | Double         |
| Functionality | Performance  | Identification | Protein coverage                     | ProteinCoverage                    | Percentage of proteins identified in all samples                                   | Double         |
| Functionality | Performance  | Identification | Peptides per protein                 | PeptidesPerProtein                 | Distribution of peptides per protein group                                         | Set of Integer |
| Functionality | Performance  | Quantification | Correlation peptides                 | CorrelationPeptides                | Mean of Pearsson correlation of protein abundances between replicates (log2-scale) | Double         |
| Functionality | Performance  | Quantification | Correlation proteins                 | CorrelationProteins                | Mean of Pearsson correlation of peptide abundances between replicates (log2-scale) | Double         |
| Functionality | Performance  | Quantification | Number peptides                      | NumberOfPeptides                   | Number of quantified peptides with at least 50% coverage                           | Integer        |
| Functionality | Performance  | Quantification | Number protein groups                | NumberOfProteinGroups              | Number of quantified proteins groups with at least 50% coverage                    | Integer        |
| Functionality | Performance  | Quantification | Dynamic peptide range                | DynamicPeptideRange                | Difference of peptide abundance (top 5% versus bottom 5% quantile)                 | Double         |
| Functionality | Performance  | Quantification | Dynamic protein range                | DynamicProteinRange                | Difference of protein abundance (top 5% versus bottom 5% quantile)                 | Double         |
| Functionality | Performance  | Statistics     | Differentially regulated peptides 5% | DifferentialRegulatedPeptides5Perc | Number of differentially regulated peptides with FDR below 5%                      | Set of Double  |
| Functionality | Performance  | Statistics     | Differentially regulated proteins 5% | DifferentialRegulatedProteins5Perc | Number of differentially regulated proteins with FDR below 5%                      | Set of Double  |
| Functionality | Performance  | Statistics     | Differentially regulated peptides 1% | DifferentialRegulatedPeptides1Perc | Number of differentially regulated peptides with FDR below 1%                      | Set of Double  |
| Functionality | Performance  | Statistics     | Differentially regulated proteins 1% | DifferentialRegulatedProteins1Perc | Number of differentially regulated proteins with FDR below 1%                      | Set of Double  |
| Functionality | Performance  | Statistics     | Missing peptide values               | MissingPeptideValues               | Percentage of missing values in entire peptide set                                 | Double         |
| Functionality | Performance  | Statistics     | Missing protein values               | MissingProteinValues               | Percentage of missing values in entire protein set                                 | Double         |
| Functionality | Performance  | Digestion      | Digestion efficiency                 | Efficiency                         | Distribution of number of miscleavages                                             | Set of Double  |
| Functionality | Performance  | PTMs           | PTM Distribution                     | PTMDistribution                    | Percentage of peptides with PTM xyz                                                | Set of Double  |
| Functionality | Performance  | PTMs           | PTM Occupancy                        | PTMOccupancy                       | Distribution of peptides with 1,2,... PTMs                                         | Set of Double  |
| Functionality | Parameter    | Identification | Database size                        | DatabaseSize                       | Number of entries in fasta file                                                    | Integer        |
| Functionality | Parameter    | Identification | Canonical sequences                  | CanonicalSequences                 | Database includes canonical sequences                                              | Y/N            |
| Functionality | Parameter    | identification | PTM localization                     | PTMLocalization                    | Is PTM localization scoring software included in the workflow                      | Y/N            |

## Contribute

Contributions to change and modify the workflows are most welcome. For this, please create a fork and add your changes. We strongly recommend reaching out to us either via email or by creating an issue in this repository, as we then help adding the new implementation.


## Credits

nf-core/wombat was originally written by the members of the ELIXIR Implementation study [Comparison, benchmarking and dissemination of proteomics data analysis pipelines](https://elixir-europe.org/internal-projects/commissioned-services/proteomics-pipelines) under the lead of Veit Schwämmle and major participation of David Bouyssié and Fredrik Levander.

## Citations

Preprint available: https://www.biorxiv.org/content/10.1101/2023.10.02.560412v1

As the workflows are using an nf-core template, we refer to the publication:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
