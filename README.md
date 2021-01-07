# hierarchical_autophagy_regulation
Research compendium for the preprint "Hierarchical Regulation of Autophagy During Adipocyte Differentiation"

## Setting up the docker environment

The analysis was run on a [docker](https://hub.docker.com/r/bcmslab/adiporeg/)
image based on the the latest **rocker/verse**.
Other R packages were added to the image and were made available as an image 
that can be obtained and launched on any local machine running
[docker](https://hub.docker.com/r/bcmslab/adiporeg/).

```bash
$ docker pull bcmslab/adiporeg:latest
$ docker run -it bcmslab/adiporeg:latest bash
```

## Obtaining the source code

The source code is hosted publicly on this repository in the form of a research
compendium. This includes the scripts to reproduce the figures and tables in 
this manuscript.

```bash
$ git clone https://github.com/BCMSLab/hierarchical_autophagy_regulation
```

## Runing the analysis

In the directory `hierarchical_autophagy_regulation`, run `make`

```bash
$ cd hierarchical_autophagy_regulation
$ make all
```

## Details of the R environment
The version of **R** that was used to perform this analysis is the 4.02
(2020-06-22) on `x86\_64-pc-linux-gnu`.

## More

This manuscript was released as a preprint under the title [Hierarchical Regulation of Autophagy During Adipocyte Differentiation](https://doi.org/10.1101/2021.01.06.425505)
