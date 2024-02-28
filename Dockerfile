FROM mambaorg/micromamba:0.14.0
LABEL version=0.1.0
LABEL maintainer="ccl40@sfu.ca"
LABEL source="https://github.com/jimmyliu1326/SamnSorter"

ARG CONDA_DIR="/opt/conda"
ARG GIT_REPO="https://github.com/jimmyliu1326/SamnSorter"
ARG BINARY_PATH="/usr/local/bin/SamnSorter"
ARG SCHEMA_URL="https://object-arbutus.cloud.computecanada.ca/cidgohshare/eagle/jimmyliu/enterobase_senterica_cgmlst_3.2.2.tar.gz"
ARG REF_DATA="https://object-arbutus.cloud.computecanada.ca/cidgohshare/eagle/jimmyliu/SamnSorter/data/samnsorter_v1.0.tar.gz"

# install linux dependencies
RUN apt-get update && \
    apt-get install -y procps wget git build-essential && \
    rm -rf /var/lib/apt/lists/*

# compile cgmlst-dists
RUN git clone https://github.com/jimmyliu1326/cgmlst-dists && \
    cd /cgmlst-dists && \
    make && \
    make check && \
    make PREFIX=/usr/local install

# clone repo and add main script to PATH
RUN git clone ${GIT_REPO} && \
    chmod +x /SamnSorter/SamnSorter.R

# install env dependencies
RUN micromamba install -n base -y -c bioconda -c conda-forge -f /SamnSorter/conda.yml && \
    pip install -r /SamnSorter/requirements.txt && \
    micromamba clean --all --yes && \
    rm -rf $CONDA_DIR/conda-meta && \
    rm -rf $CONDA_DIR/include && \
    rm -rf $CONDA_DIR/lib/python3.*/site-packages/pip && \
    find $CONDA_DIR -name '__pycache__' -type d -exec rm -rf '{}' '+'

# download reference data
RUN wget ${SCHEMA_URL} -O /enterobase_senterica_cgmlst_3.2.2.tar.gz && \
    tar -xzvf /enterobase_senterica_cgmlst_3.2.2.tar.gz -C / && \
    wget ${REF_DATA} -O /ref_data.tar.gz && \
    mkdir /SamnSorter/ref/ && \
    tar -xzvf /ref_data.tar.gz -C /SamnSorter/ref/

# set env variables for container
ENV PATH="${PATH}:/SamnSorter"
ENV MODEL_DIR="/SamnSorter/models"
ENV REF_CLUSTERS="/SamnSorter/ref/clusters.tsv"
ENV REF_NWK="/SamnSorter/ref/reference.rooted.anon.nwk"
ENV REF_TAXONOMY="/SamnSorter/ref/taxonomy.tsv"
ENV CGMLST_DISTS="/cgmlst-dists/cgmlst-dists-query"
ENV REF_ALLELES="/SamnSorter/ref/allele_hashed_full_outgroup_anon.tsv"
ENV REF_PATH="/enterobase_senterica_cgmlst_3.2.2/"
ENV TRAINING_PATH="/enterobase_senterica_cgmlst_3.2.2/Salmonella_enterica.trn"