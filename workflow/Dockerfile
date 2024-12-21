# Dockerfile to build Docker image for this workflow
FROM ubuntu:20.04

LABEL maintainer = "Jason Charamis" 
LABEL contact    = "jason.charamis@gmail.com"
LABEL build_date = "2024-12-21"
LABEL version    = "v.0.0.1-dev"

# Install essential Ubuntu dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    libssl-dev \
    libxml2-dev \
    git \
    gzip \
    tar \
    wget
    
# Load mambaforge image
FROM condaforge/mambaforge:latest

RUN mkdir -p /conda-envs/dnaseq
COPY envs/dnaseq.yaml /conda-envs/dnaseq/environment.yaml

# Install required packages through conda
RUN mamba env create --prefix /conda-envs/dnaseq --file /conda-envs/dnaseq/environment.yaml && \
    mamba clean --all -y

# Set the default environment path
ENV PATH /conda-envs/dnaseq/bin:$PATH

RUN git clone https://github.com/JasonCharamis/AVC.git

ENV PATH="/root/.local/bin:$PATH"
ENV PATH="/AVC:$PATH"
ENV PATH="/AVC/workflow/scripts:$PATH"

# Set the default command to run when starting the container
SHELL ["/bin/bash", "-c"]
CMD ["/bin/bash"]