FROM ubuntu:22.04

USER root

# Install base utilities
RUN apt-get update \
	&& apt-get install -y build-essential \
	&& apt-get install -y wget \
	&& apt-get clean \
	&& rm -rf /var/lib/apt/lists/*

# Install Miniconda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py311_23.5.2-0-Linux-x86_64.sh -O ~/miniconda.sh && \
	bash ~/miniconda.sh -b -p /opt/conda && \
	rm ~/miniconda.sh
ENV PATH /opt/conda/bin:$PATH

# Set the default environment for the container
ENV CONDA_DEFAULT_ENV=SComatic

# Create conda environment and install packages
RUN /opt/conda/bin/conda init bash && \
    /opt/conda/bin/conda create -n ${CONDA_DEFAULT_ENV} -c bioconda python=3.7 r-base=3.6.1 samtools datamash bedtools -y && \
    /bin/bash -c "source /opt/conda/bin/activate ${CONDA_DEFAULT_ENV}" && \
    /opt/conda/bin/conda clean -a -y
ENV PATH /opt/conda/envs/$CONDA_DEFAULT_ENV/bin:$PATH

# Ensure that the environment is activated by default
RUN echo "source activate ${CONDA_DEFAULT_ENV}" >> ~/.bashrc

# Install Python dependencies using pip in the Conda environment
COPY requirements.txt /opt/
RUN /bin/bash -c "source activate ${CONDA_DEFAULT_ENV} && pip install -r /opt/requirements.txt"

# Install R dependencies using Rscript in the Conda environment
COPY r_requirements_install.v3_6.R /opt/
RUN /bin/bash -c "source activate ${CONDA_DEFAULT_ENV} && Rscript /opt/r_requirements_install.v3_6.R"

# Copy panel of normals (PoNs)
COPY PoNs /opt/SComatic/PoNs
RUN gunzip /opt/SComatic/PoNs/*.gz

# Copy RNAediting
COPY RNAediting /opt/SComatic/RNAediting
RUN gunzip /opt/SComatic/RNAediting/*.gz

# Copy BED files of interest
COPY bed_files_of_interest /opt/SComatic/bed_files_of_interest

# Copy example data
COPY example_data /opt/SComatic/example_data

# Generate samtools index
RUN /bin/bash -c "gunzip /opt/SComatic/example_data/chr10.fa.gz && \
    source activate ${CONDA_DEFAULT_ENV} && \
    which samtools && \
    samtools faidx /opt/SComatic/example_data/chr10.fa"

# Copy all scripts to /usr/local/bin to make them easily accessible
COPY scripts/*/* /usr/local/bin
RUN find /usr/local/bin -type f -name "*.py" -exec sed -i '1i#!/usr/bin/env python' {} +
RUN chmod +x /usr/local/bin/*.py

## USER CONFIGURATION
RUN adduser --disabled-password --gecos '' ubuntu && chsh -s /bin/bash && mkdir -p /home/ubuntu

# Set user and working directory
USER    ubuntu
WORKDIR /home/ubuntu

# Load what we need 
CMD ["/bin/bash", "-c", "source ~/.bashrc"]