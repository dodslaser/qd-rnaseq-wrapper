# qd-rnaseq-wrapper
This is a wrapper for starting nf-core/rnaseq and nf-core/rnafusion on a set of fastq files

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

### nf-core/rnafusion
The pipeline is installed in `/apps/bio/repos/nf-core/nf-core-rnafusion-X.Y.Z/`.

Paths to where singularity images the pipe is dependent on can be found is set
using the env variable `$NXF_SINGULARITY_CACHEDIR`. This is best to have set in a .bashrc

Dependencies in the form of reference data and genome index was built using the --build_references flag. 
In order to be able to build all dependencies, a [COSMIC](https://cancer.sanger.ac.uk/cosmic/) account is required.

Code for installation

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
