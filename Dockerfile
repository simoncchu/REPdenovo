FROM ubuntu:16.04


# general updates & installing necessary Linux components
RUN apt-get update -y && apt-get install -y \
    build-essential \
    bzip2 \
    gcc \
    make \
    cmake \
    libbz2-dev \
    liblzma-dev \
    libjsoncpp-dev \
    libncurses-dev \
    zlib1g-dev \
    python-software-properties \
    time \
    git \
    wget


# download tools
WORKDIR /usr/local/bin
RUN mkdir jellyfish \
    && cd jellyfish \
    && wget https://github.com/gmarcais/Jellyfish/releases/download/v2.3.0/jellyfish-linux \
    && mv jellyfish-linux jellyfish \
    && chmod +x jellyfish

RUN wget https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2 \
    && tar -xjf bwa-0.7.17.tar.bz2 \
    && cd bwa-0.7.17 \
    && make \
    && cd .. \
    && ln -s bwa-0.7.17 bwa

RUN wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 \
    && tar -xjf samtools-1.9.tar.bz2 \
    && cd samtools-1.9 \
    && make \
    && cd .. \
    && ln -s samtools-1.9 samtools

RUN wget https://www.ebi.ac.uk/~zerbino/velvet/velvet_1.2.10.tgz \
    && tar -xzf velvet_1.2.10.tgz \
    && cd velvet_1.2.10 \
    && make ’MAXKMERLENGTH=60’ \
    && cd .. \
    && ln -s velvet_1.2.10 velvet 

RUN wget https://github.com/pezmaster31/bamtools/archive/v2.5.1.tar.gz \
    && tar -xzf v2.5.1.tar.gz \
    && cd bamtools-2.5.1 \
    && mkdir build \
    && cd build \
    && cmake .. \
    && make \
    && make install


# set path
ENV PATH=/usr/local/bin/bwa/:$PATH
ENV PATH=/usr/local/bin/samtools/:$PATH
ENV PATH=/usr/local/bin/jellyfish/:$PATH
ENV PATH=/usr/local/bin/velvet/:$PATH

# supporting UTF-8
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8

# wrapper
WORKDIR /home/REPdenovo/
COPY ./ /home/REPdenovo/
RUN cd TERefiner && make
RUN cd ContigsCompactor-v0.2.0/ContigsMerger/ && make
RUN cp ./TERefiner/TERefiner_1 ./ && cp ./ContigsCompactor-v0.2.0/ContigsMerger/ContigsMerger ./
RUN chmod +x *

RUN cd /home/REPdenovo/



