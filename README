# SARIFA Analysis scripts

This repository houses all scripts need to reproduce the RNAseq analysis of *insert link here*.

## How to run

### Docker installation

The scripts here are packaged as a Docker container which Downloads the RNAseq 
data from TCGA GDC and executes the analysis on the data. 
For Information on how to install Docker on your system refer to https://docs.docker.com/get-docker/.
You will also need to have docker-compose available on your system. 

### Clone repo and execute Container

```bash
git clone https://github.com/schlesnerlab/sarifa-analysis
cd sarifa-analysis
docker-compose up -d --build
```

Once the container has finished execution results can be found in `./data/plots`.

### Hardware requirements

This container will download RNAseq counts for ~280 samples from TCGA which can be around 2GB.  We recommend around 4 GB of Hard drive space available. 

Further to ensure that the container has sufficient RAM, ensure that the container has access to at least 8 GB of RAM if not more. This can be set depending on your Docker installation. 

