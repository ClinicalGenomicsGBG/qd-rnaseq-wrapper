# qd-rnaseq-wrapper
This is a wrapper for starting nf-core/rnaseq and nf-core/rnafusion on a set of fastq files

## Dependencies
The wrapper itself is written in python and has a few python modules as dependencies. The main dependency is nextflow.
All these are handled by a conda environment which can be activated before running the wrapper.

## Pipeline installation
The wrapper use two nf-core pipelines, [nf-core/rnaseq](https://nf-co.re/rnaseq) and [nf-core/rnafusion](https://nf-co.re/rnafusion/2.0.0).
Both pipelines were downloaded using the [nf-core tools suite](https://nf-co.re/tools/).

### nf-core/rnaseq

The pipeline is installed in `/apps/bio/repos/nf-core/nf-core-rnaseq-X.X/`.

Paths to where singularity images the pipe is dependent on can be found is set
using the env variable `$NXF_SINGULARITY_CACHEDIR`. This is best to have set in a .bashrc

Code for installation

```
export NXF_SINGULARITY_CACHEDIR="/apps/bio/dependencies/nf-core/singularities"
cd /apps/bio/repos/nf-core/nf-core-rnaseq-X.Y/
nf-core download rnaseq --container singularity
```

**Pre-downloaded references**

Some dependencies, such as genome references can be downloaded on the fly using the 
--genome flag of the rnaseq pipeline. However, this download sometimes takes quite a bit of time, and can cause
the pipeline to crash in the case of network dropouts.

In order to circumvent this, the references can be downloaded from AWS iGenomes. For this to work you need to
use the aws s3 tool-suite. The following code snippet shows how to download the fasta, gtf, bed and STARIndex 
for GRCh38 from NCBI. To construct other download options, use the website [AWS iGenomes](https://ewels.github.io/AWS-iGenomes/).

To use them, specify paths in the config.ini under '[rnaseq-references]'

```
module load aws-cli/2.7.7
aws s3 --no-sign-request --region eu-west-1 sync s3://ngi-igenomes/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/ /apps/bio/dependencies/nf-core/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/
aws s3 --no-sign-request --region eu-west-1 sync s3://ngi-igenomes/igenomes/Homo_sapiens/NCBI/GRCh38/Annotation/Genes/ /apps/bio/dependencies/nf-core/igenomes/Homo_sapiens/NCBI/GRCh38/Annotation/Genes/ --exclude "*" --include "genes.gtf"
aws s3 --no-sign-request --region eu-west-1 sync s3://ngi-igenomes/igenomes/Homo_sapiens/NCBI/GRCh38/Annotation/Genes/ /apps/bio/dependencies/nf-core/igenomes/Homo_sapiens/NCBI/GRCh38/Annotation/Genes/ --exclude "*" --include "genes.bed"
aws s3 --no-sign-request --region eu-west-1 sync s3://ngi-igenomes/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/STARIndex/ /apps/bio/dependencies/nf-core/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/STARIndex/
```

### nf-core/rnafusion
The pipeline is installed in `/apps/bio/repos/nf-core/nf-core-rnafusion-X.Y.Z/`.

Paths to where singularity images the pipe is dependent on can be found is set
using the env variable `$NXF_SINGULARITY_CACHEDIR`. This is best to have set in a .bashrc

Dependencies in the form of reference data and genome index was built using the --build_references flag. 
In order to be able to build all dependencies, a [COSMIC](https://cancer.sanger.ac.uk/cosmic/) account is required.

**Code for installation**

```
export NXF_SINGULARITY_CACHEDIR="/apps/bio/dependencies/nf-core/singularities"
cd /apps/bio/repos/nf-core/nf-core-rnafusion-X.Y.Z/
nf-core download rnafusion --container singularity
nextflow /apps/bio/repos/nf-core/nf-core-rnafusion-2.0.0/workflow/main.nf 
    -c /apps/bio/repos/nf-core-configs/conf/medair.config 
    -profile singularity,byss 
    --build_references 
    --all
    --genomes_base /apps/bio/dependencies/nf-core/nf-core-rnafusion-X.Y.Z
    --outdir /apps/bio/dependencies/nf-core/nf-core-rnafusion-X.Y.Z
    ---cosmic_username <username>
    --cosmic_passwd <password>  
```

## Usage
### Input data
* The runner takes a directory with gzipped fastq files as input. There should only be two files per sample. If there are more, please concatenate them before running the wrapper

### The Wrapper
The pipeline is meant to be initiated by the wrapper which queries SLIMS for any samples marked for QD-RNA as its secondary analysis. 
Analysed samples are stored in a flat-text file for now in order to keep track of which samples have been analysed

### Output
The output from both pipelines is set in the config file. From this a subset 
of files is copied to a reporting directory (also set in config). Logs for the pipelines end up in the 
output directory. One exception to this is the wrapper log, which ends up in a log dir set in the config.

## TODO
* Allow for more than 2 input files per sample
* Move mapping files to IGV
* Report back to SLIMS if a samples has been analysed already or not