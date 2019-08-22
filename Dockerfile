FROM ubuntu:16.04

# general updates & installing necessary Linux components
RUN apt-get update -y && apt-get install -y \
    bzip2 \
    gcc \
    libncurses-dev \
    make \
    libbz2-dev \
    liblzma-dev \
    time \
    zlib1g-dev \
    python-software-properties \
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
    tar -xjf samtools-1.9.tar.bz2 \
    && cd samtools-1.9 \
    && make \
    && cd .. \
    && ln -s samtools-1.9 samtools

RUN wget https://www.ebi.ac.uk/~zerbino/velvet/velvet_1.2.10.tgz \
    tar -xzf velvet_1.2.10.tgz \
    && cd velvet_1.2.10 \
    && make \
    && ./velveth \
    && cd .. \
    && ln -s velvet_1.2.10 velvet 

RUN cd /home/ \
    && mkdir REPdenovo \
    && cd REPdenovo


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
COPY * /home/REPdenovo/
RUN chmod +x *.py



