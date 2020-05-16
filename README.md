# snakeRNAseq
RNAseq read processing and alignment Snakemake pipeline



# Before running
In the `config.yaml` file you should include details of samples to be analyzed and required index paths as per specifications. You can also specify custom options for the trimming and alignment steps. 

Samples can be run in single or paired ended modes, and the corresponding option can be specified in the `config.yaml` file as `se` or `pe` respectively. All input files must have the format `.fq.gz` and be placed a folder called `input` . Paired ended samples must be specified in the following format: `sample_1.fq.gz, sample_2.fq.gz`. Multiple samples can be specified under the keyword `samples` in the `config.yaml` file in the following format (not the space between `-` and `sample`):

```
samples:
  - sample1
  - sample2
  ...
```

Also, you can specify any one of the following alignment methods in the `config.yaml` file:
1) bowtie2
2) hisat2
3) kallisto

If you select `bowtie2` or `hisat2` the final output would be the sorted and indexed bam files. If you select `kallisto`, the final output would be the abundance files.


### Dependiencies
* python 3
* [snakemake](https://snakemake.readthedocs.io/en/stable/)
* [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [trimgalore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
* [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
* [hisat2](http://daehwankimlab.github.io/hisat2/)
* [samtools/1.9](http://www.htslib.org/)


### How to run?

```bash
mkdir fastQC_output
Dry run: snakemake -n
Actual run: snakemake --cores [number of cores]

Run snakemake --unlock if directory is locked
```
Dry run displays the rules that will be executed
