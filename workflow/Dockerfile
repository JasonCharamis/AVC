## Dockerfile to build Docker image for this Snakemake workflow.
FROM ubuntu:20.04

LABEL maintainer = "Jason Charamis" 
LABEL contact    = "jason.charamis@gmail.com"
LABEL build_date = "2024-01-02"
LABEL version    = "v.0.0.1-dev"

# Set environment variables for non-interactive installation of R
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Europe

# Step 1: Install essential dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        libssl-dev \
        libxml2-dev \
        git \
        python3-pip \
	gzip \
	emacs \
	wget

## Load mambaforge image
FROM condaforge/mambaforge:latest

RUN mkdir -p /conda-envs/dnaseq
COPY envs/dnaseq.yaml /conda-envs/dnaseq/environment.yaml

## Install required packages through conda
RUN mamba env create --prefix /conda-envs/dnaseq --file /conda-envs/dnaseq/environment.yaml && \
    mamba clean --all -y

RUN apt-get update && \
    apt-get install --no-install-recommends -y build-essential curl

## Install Python packages, including Snakemake, through pip
RUN pip3 install wheel && \
    pip3 install datrie && \
    pip3 install snakemake

RUN git clone https://github.com/JasonCharamis/BiS.git
  
ENV PATH="~/.local:$PATH"
ENV PATH="BiS:$PATH"
ENV PATH="BiS/workflow/scripts:$PATH"
ENV PATH="/conda-envs/varscan/bin:$PATH"
ENV PATH="/usr/lib/R/bin:$PATH"
