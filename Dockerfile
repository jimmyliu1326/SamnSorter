FROM mambaorg/micromamba:0.14.0
LABEL version=0.1.0
LABEL maintainer="ccl40@sfu.ca"
LABEL source="https://github.com/jimmyliu1326/SamnSorter"

ARG CONDA_DIR="/opt/conda"
ARG GIT_REPO="https://github.com/jimmyliu1326/SamnSorter"
ARG BINARY_PATH="/usr/local/bin/SamnSorter"
ARG SCHEMA_URL="https://object-arbutus.cloud.computecanada.ca/cidgohshare/eagle/jimmyliu/enterobase_senterica_cgmlst_3.2.2.tar.gz"

# include env dependencies
ADD ./conda.yml /

# install env dependencies
RUN apt-get update && \
    apt-get install -y procps wget git build-essential && \
    rm -rf /var/lib/apt/lists/* && \
    micromamba install -n base -y -c bioconda -c conda-forge -f /conda.yml && \
    micromamba clean --all --yes && \
    rm -rf $CONDA_DIR/conda-meta && \
    rm -rf $CONDA_DIR/include && \
    rm -rf $CONDA_DIR/lib/python3.*/site-packages/pip && \
    find $CONDA_DIR -name '__pycache__' -type d -exec rm -rf '{}' '+'

# install cgmlst-dists
RUN git clone https://github.com/jimmyliu1326/cgmlst-dists && \
    cd /cgmlst-dists && \
    make && \
    make check && \
    make PREFIX=/usr/local install

# clone repo and add main script to PATH
RUN cd / && \
    git clone ${GIT_REPO} && \
    ln -s /SamnSorter/src/SamnSorter.R ${BINARY_PATH} && \
    chmod +x $BINARY_PATH

# download schema
RUN wget ${SCHEMA_URL} -O /enterobase_senterica_cgmlst_3.2.2.tar.gz && \
    tar -xzvf /enterobase_senterica_cgmlst_3.2.2.tar.gz -C /

# set env variables for container
ENV MODEL_DIR="/SamnSorter/models"
ENV REF_CLUSTERS="/SamnSorter/ref/clusters.tsv"
ENV REF_NWK="/SamnSorter/ref/reference.rooted.nwk"
ENV REF_TAXONOMY="/SamnSorter/ref/taxonomy.tsv"
ENV CGMLST_DISTS="/cgmlst-dists/cgmlst-dists-query"
ENV REF_ALLELES="/SamnSorter/ref/results_allele_hashed_full_outgroup.tsv"
ENV REF_PATH="/enterobase_senterica_cgmlst_3.2.2/"
ENV TRAINING_PATH="/enterobase_senterica_cgmlst_3.2.2/Salmonella_enterica.trn"