[![Snakemake](https://img.shields.io/badge/snakemake-≥5.14.0-brightgreen.svg)](https://snakemake.bitbucket.io)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
![Tests](https://github.com/jounikuj/snakemake-qiaseq-mirna/workflows/Tests/badge.svg?branch=master)

# snakemake-qiaseq-mirna

`snakemake-qiaseq-mirna` is a Snakemake pipeline for the quantification of known miRNAs and piRNAs from sequencing data that is obtained from sequencing libraries prepared with Qiagen's [QIAseq miRNA Library Kit](https://www.qiagen.com/fi/products/discovery-and-translational-research/next-generation-sequencing/metagenomics/qiaseq-mirna-ngs/?clear=true#orderinginformation).

## Authors
* Jouni Kujala, University of Eastern Finland, Institute of Clinical Medicine

## Usage

#### Step 1: Install dependencies

This pipeline has two dependencies that must be installed to your system before running this pipeline. These dependencies are 

* [Snakemake ≥5.22](https://snakemake.readthedocs.io/en/stable/)
* [Conda (either Miniconda or Anaconda version)](https://docs.conda.io/en/latest/)

We recommend you to install Conda package management system and then install Snakemake with it. After the succesfull installation of Conda, you can [create virtual environment](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html) and run Snakemake within them. Use of virtual environments is a recommended way to run pipeline because it isolates pipeline from your operating system and helps to avoid problems with conflicting packages. This can be achieved by running following commands in terminal;

```
# Install Mamba package manager.
conda install -c conda-forge mamba

# Create virtual environment with Snakemake.
mamba create -c conda-forge -c bioconda -n [name of the environment] snakemake

# Activate virtual environment.
conda activate [name of the environment]
```

Pipeline utilizes a number of bioinformatic tools that are installed by Snakemake itself. If interested, you can find complete list of used tools from configuration file `workflow/envs/env.yaml`.

#### Step 2: Obtain a copy of this pipeline

First step is obviously to obtain a copy of this pipeline to your system. This can be achieved in a number of ways;

1. Clone this repository by running following commands in your terminal;
   ```
   # Define directory to which the repository should be cloned.
   cd [path to directory]  

   # Clone repository
   git clone https://github.com/jounikuj/snakemake-qiaseq-mirna.git
   ```
2. Download this repository as [zipped file](https://github.com/jounikuj/snakemake-qiaseq-mirna/archive/master.zip), place it to desired location and unzip it.

#### Step 3: Set up configuration

Although this pipeline is designed to be as portable as possible, there are some things that must be provided by the user before running this pipeline.
1. Place your input files to `fastq/` directory. Input files are expected to be in .fastq.gz format.
2. Modify file `samples.tsv` according to your needs. This file specifies which .fastq.gz file belongs to which sample. File should be provided in following format;
    ```
    > head -n 5 config/samples.tsv

    sample  fastq
    sampleA fastq/sampleA_S1_L001_R1_001.fastq.gz
    sampleA fastq/sampleA_S1_L002_R1_001.fastq.gz
    sampleB fastq/sampleA_S2_L001_R1_001.fastq.gz
    sampleB fastq/sampleA_S2_L002_R1_001.fastq.gz
    ...
    ```
   Please note that pipeline renames files according to sample list. For example, .fastq.gz files `sampleA_S1_L001_R1_001.fastq.gz` and `sampleA_S1_L002_R1_001.fastq.gz` will produce filenames like `sampleA.bam` because indicated .fastq.gz files are associated with `sampleA` in the sample list.
3. Check file `config.yaml`. This file contains various parameters such as used adapter sequences and source URLs for used miRNA and piRNA sequences. You are free  to try different settings for your data but it is always a good option to start with the defaults.

#### Step 4: Execute pipeline

After completing the configuration you are ready to proceed to pipeline execution.

1. Test your configuration by running Snakemake in dryrun mode by executing a following command in your terminal;
   ```
   snakemake --dryrun
   ```
2. If no error messages are shown in `--dryrun` mode, you are ready to proceed into pipeline execution. Execute pipeline by running following command in your   terminal;
   ```
   snakemake --use-conda --cores [available cores]
   ```
   
   Please note that `--use-conda` parameter is required to run pipeline as it allows Snakemake to install all required dependencies. Snakemake should now print status messages to console and start running. You can follow the execution of pipeline from these messages.

## Pipeline description

This section aims to describe the workflow of this pipeline. 

#### Step 1: Read trimming and UMI detection

Raw reads contain a number of different sequences that must be removed before the actual analysis. For this purpose we will use [`cutadapt`](https://cutadapt.readthedocs.io/en/stable/index.html#) and [`umi_tools`](https://umi-tools.readthedocs.io/en/latest/). Due to structure of reads prepared with the QIAseq miRNA Library Kit the read trimming and UMI detection is done in 
sequential order.

1. First we remove the possible remnants of Illumina Universal Adapter sequence `AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC` from the 3' end of the read. At this step, reads shorter than 20 bp are discarded from further analysis.
2. A unique molecular identifier (UMI) of 12 random basepairs is removed from from the 3' end of the read and added to read header. Each RNA molecule is tagged with an unique UMI before sequencing library is amplified with a PCR. Each PCR product originating from the same RNA molecule is expected to have same UMI and this information can be later used to reduce the sequencing errors and quantitative bias caused by PCR amplification.
3. Adapter sequence `GTTCAGAGTTCTACAGTCCGACGATC` is trimmed from the 5' end of the reads.
4. Adapter sequence `AACTGTAGGCACCATCAAT` is trimmed from the 3' end of the read. By default, pipeline performs some quality-trimming for reads and discards reads shorter than 20 bp.

#### Step 2: Read alignment

Pipeline downloads FASTA sequences of known miRNAs and piRNAs from [miRBase](http://www.mirbase.org/) and [piRBase](http://www.regulatoryrna.org/database/piRNA/). 
Pipeline uses [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) with `--very-sensitive-local` parameter to map reads to FASTA sequences. Alignment 
is done in sequential order:

1. Reads are first mapped to known miRNA sequences. Unaligned reads are collected to new .fastq.gz file and forwarded to next step.
2. Unaligned reads are mapped to known piRNA sequences. Again, unaligned reads are collected to new .fastq.gz file `sampleA_unaligned.fastq.gz` etc.

After each step, aligned reads are converted to BAM files to recude the pipeline's memory usage. 

#### Step 3: Deduplication

Aligned reads are merged to one .bam file to reduce the number of processed files. We use `umi_tools` to deduplicate raw reads on the basis of their UMI sequences 
stored in the read header. Reads with the same UMI are collapsed into one read to reflect the state of the sample before library amplification.

#### Step 4: Read counting

Aligned reads are calculated with [`samtools idxstats`](http://www.htslib.org/doc/samtools.html) and merged into one file. Read counts are calculated both for the raw reads and deduplicated reads - you can freely choose which one you use in your downstream analyses. Resulting count table should look like this;

```
> head -n 5 results/counts.tsv

Molecule       sampleA sampleB sampleC
hsa-let-7a-1   0       0       4 
hsa-let-7a-2   0       4       0
hsa-let-7a-3   4       4       0
hsa-let-7b     12      4       4
...
```

Please note that read counts are not normalized.

## Citing this pipeline

At the moment, this pipeline has not been published anywhere. If you use this pipeline in your research, just cite the URL of this repository until we can provide a better option.