# ChIP-seq
Basic processing for collections of ChIP-seq data. Currently only configured for single-end reads but can be updated to accomodate paired-end reads.

## Prerequisites
The prerequisites for this pipeline are minimal due to Singularity. All you need is to install Singularity. This package installs all prerequisites for software. [Singularity Installation](http://singularity.lbl.gov/docs-installation)

Other than singularity, there are a few requirements. The first is your directory structure. This pipeline takes advantage of a regularized directory structure. Conveniently, most of the directories are generated by Snakemake automatically so your input is minimal. You need to define a parent directory with a folder in it labeled data. Data will have one subdirectory called raw_data. As the name implies, that is where you put your fastq.gz files. NOTE: this pipeline does no QC so be sure that your files are trimmed and gzipped beforehand. Here is how the setup of your directories can look:
```
mkdir -p PARENT_DIRECTORY/data/raw_data 
cp RAW_FILES raw_data
```
After depositing raw files into the raw_data file, it's time to load other prequisites. For now, this script does not automatically generate the 2 prerequisite files necessary (perhaps something for future modifications). These files are the chromosome sizes file and the Bowtie2 reference genome. You should maintain these independently of the pipeline and copy them into your data directory after. Naming is important for this pipeline so be sure to name anything that is lowercase as specified. Here is how you would create these files and then copy them into your data folder:
```
samtools faidx GENOME.fa
cut -f1,2 GENOME.fa.fai > GENOME.chrom.sizes
cp GENOME.chrom.sizes PARENT_DIRECTORY/data/chrom.sizes

bowtie2-build --threads NUM_THREADS GENOME.fa GENOME
mkdir PARENT_DIRECTORY/data/bowtie_index/
cp GENOME.* PARENT_DIRECTORY/data/bowtie_index/
```
After creating the Bowtie2 index and copying it, you have to rename each of the files by changing their prefix to idx. This may be updated for flexibility in the future (perhaps a config command with a Bowtie2 index). Here's an example for an hg38 prefix:
```
mv hg38.1.bt2 idx.1.bt2
```
Do this for the files in the reference folder.

## Getting Started
There are three fields in the config file (config.yaml) to alter: tag_dir_cmds, peak_cmds and paired_end. 
For tag_dir_cmds, this will allow you to customize your HOMER tag directories. Check out the [HOMER Tag Directory documentation](http://homer.ucsd.edu/homer/ngs/tagDir.html) Note: if you want to look at genome characteristics such as GC content add something like the following:
```
tag_dir_cmds: -genome PATH_TO_FASTA -checkGC
```
It is important to specify your path as this version of the pipeline does not load the HOMER version of the genome. This is due to an onerous amount of space required for constructing the Singularity container.

For peak_cmds use customize peak finding as done by HOMER. For more info, check out the [HOMER Peak finding documentation](http://homer.ucsd.edu/homer/ngs/peaks.html). One of the most common things to change is the peak style. For example, for histone data, the following style is recommend:
```
peak_cmds: -style histone
```
One tip that I have for peak finding is to include an input HOMER tag directory (not generated here). You would move it into the data directory like other prerequisites with the specified name like so:
```
cp -r INPUT_TAG_DIRECTORY PARENT_DIRECTORY/data/INPUT_TAG_DIRECTORY
```
To include the input into peak finding, you would include the following flag in the config filE:
```
peak_cmds: -style histone -i /scif/data/INPUT_TAG_DIRECTORY
```
The location "/scif/data/" is essential. Otherwise, this will NOT work.

The other option is to specify whether the raw files are paired end or not. Set this as True or False. For example, for paired-end data, do the following:
```
paired_end: True
```

The second part of this is generating the Singularity simg file. Use the Singularity manual for installation. It is required that you do this on a computer where you have sudo privileges and after you alter your config file. The command is simple:
```
sudo singularity build chipseq.simg Singularity
```
Do this inside the directory and it will generate the chipseq.simg file. This takes some time. Then, scp or cp the file to your directory of interest.
```
scp chipseq.simg PARENT_DIRECTORY
```
Once you have the simg file, do a dry run in order to test that everything is in the right place.
```
singularity run --bind data/:/scif/data chipseq.simg run snakemake '-n'
```
This will deposit the config file in your data folder so that you can edit it with the shell.

## Running the pipeline
It is recommended that you run this pipeline outside the Singularity environment. You also need to specify the number of cpus to use on this job. With the current setting, the maximum is 50. Here is how you run the pipeline:
```
singularity run --bind data/:/scif/data chipseq.simg run snakemake '-j NUM_CPU'
```
