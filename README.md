# *FXN* Methylation Pipeline

This is a pipeline to analyze DNA methylation of *FXN* CPG island / shore which is found at the 5' region of *FXN* 
which includes the *FXN* promoter, exon 1 and intron 1. To learn more about the origin of this process and other 
related information, you can read the following:

*Lists papers*

Using this pipeline, you will be able to generate figures like the ones in [Paper name](). Also included are stats on
the methylation of each amplicon (1-5) included in these studies.

## Installation

### Prerequisites
* Docker: https://docs.docker.com/engine/install/
* Git: https://git-scm.com/book/en/v2/Getting-Started-Installing-Git

### Build Docker Image

Navigate to the directory where you would like to keep this project and run the following commands to build the project:

```bash
git clone <repo link>
cd fxn-methylation
docker build . -t fxn-methylation
```

This created a docker image with the tag *(-t)* `fxn-methylation`. 

### Run Docker Container

In order to run the script, you must first create 
a docker container with that same image that we just made. Here is the command that will open up an interactive terminal 
in our new container:

```bash
docker run -it -v <input directory>:/app/in -v <output directory>:/app/out fxn-methylation
```

Replace the following in the command above before running:

`<input directory>`: The directory on your computer that you have the FASQs that you would like to be processed.

`<output directory>`: The directory on your computer where you want your results.

## Usage

Now you should be in a new shell in the /app directory. From here, you are able to run the script:

```bash
./run_pipeline.sh in/*.fastq.gz
```

The command above will run the pipeline on amplicon 3 of all files which are in `<input directory>`. There are
other options as well for running the script. You can specify other amps by using `-a[1-5]`. for instance,
you could run the above command using amplicons 2 to 5 as well:

```bash
./run_pipeline.sh -a2 -a3 -a4 -a5 in/*.fastq.gz
```

> 💡 In [Research paper name]() we found that only methylation in amplicon 3 had a significant difference between those 
> with and without FRDA. However, if your sample contains amplicon 1, you can still run the pipeline to verify those 
> results by adding *-a1* to the command. Running with *-a1* when the sample doesn't contain amplicon 1 will result 
> in an error. Amplicon 1 has been excluded when using Next-generation sequencing on modern samples to reduce costs.

To specify the number of DNA strands to be observed for the analysis, you can use the `-n` option. By default 300 
strands are used. After 300 there is no statistical difference.

For example, you can analyze 500 strands of amplicon3 using the command below:

```bash
./run_pipeline.sh -n 500 in/*.fastq.gz
```

If you want to run the command on one sample, you can also do that:
```bash
./run_pipeline.sh in/SAMPLE1_R1.fastq.gz in/SAMPLE1_R2.fastq.gz
```

For a full list of options run with the help flag:
```bash
./run_pipeline.sh -h
```

## Analysis

All artifacts and significant files produced while running the script will be found in the `<output directory>`
that you gave in above steps. At the top level you will find bowtie indexes that are generated and used during
the process. There will also be a directory for each sample that you have run *(SAMPLE1)*. Inside you will find some 
files of interest:

| File path                                             | Usage                                                                                                            |
|-------------------------------------------------------|------------------------------------------------------------------------------------------------------------------|
| \<output directory>/SAMPLE1/amplicon[1..5].aligned    | Files which are used for the generation of figures and stats which contain methylation information per amplicon. |
| \<output directory>/SAMPLE1/analysis/stats[1..5].json | JSON file of methylation statistics for each amplicon.                                                           |
| \<output directory>/SAMPLE1/analysis/figure[1..5].svg | SVG image showing a visual like the one at the top of the README.                                                |

