FROM ubuntu:latest
LABEL authors="akhorgad"

ENV SRC=/usr/local/src
ENV BIN=/usr/local/bin
ENV CONDA_DIR=/opt/conda

WORKDIR $SRC

RUN apt-get update --allow-releaseinfo-change
RUN apt-get update -m
RUN apt-get install -y curl git-lfs time datamash wget bzip2 zlib1g-dev libssl-dev lbzip2 make gcc libbz2-dev libncurses5-dev libncursesw5-dev liblzma-dev


# Install miniconda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    bash ~/miniconda.sh -b -p /opt/conda

# Adding conda in path so we can use conda activate
ENV PATH=$CONDA_DIR/bin:$PATH

# install gsutil
RUN curl https://sdk.cloud.google.com | bash
ENTRYPOINT ["top", "-b"]

ENV PATH=/root/google-cloud-sdk/bin:$PATH

#install skera
RUN conda install -c bioconda pbskera

# install lima
RUN conda install -c bioconda lima

# install ISOSEQ
RUN conda install -c bioconda isoseq

#adding packages to support pb plots
RUN conda install -y matplotlib -c conda-forge
RUN pip install pysam
RUN conda install -y tqdm -c conda-forge
RUN conda install -y numpy -c conda-forge

RUN wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && \
	tar jxf samtools-1.9.tar.bz2 && \
	rm samtools-1.9.tar.bz2 && \
	cd samtools-1.9 && \
	./configure --prefix $(pwd) && \
	make

ENV PATH=${PATH}:/usr/local/src/samtools-1.9

SHELL ["bash", "-l" ,"-c"]