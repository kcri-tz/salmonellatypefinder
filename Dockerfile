FROM debian:stretch

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update -qq; \
    apt-get install -y -qq git \
    apt-utils \
    wget \
    python-pip \
    ncbi-blast+ \
    libz-dev \
    ; \
    mkdir workdir; \
    mkdir database

ENV DEBIAN_FRONTEND Teletype

# Install KMA
RUN cd /usr/src; \
    git clone --branch 1.2.1 https://bitbucket.org/genomicepidemiology/kma.git; \
    cd kma && make; \
    mv kma* /bin/; \
    cd /usr/src; \
    rm -r kma

#
# Install Python 3.7 (apt-get only gets 3.5)
#
RUN apt-get install -y -qq libffi-dev libssl-dev;

RUN cd /usr/src; \
    wget https://www.python.org/ftp/python/3.7.2/Python-3.7.2.tgz; \
    tar xzf Python-3.7.2.tgz; \
    rm Python-3.7.2.tgz; \
    cd Python-3.7.2; \
    ./configure --enable-optimizations; \
    make altinstall; \
    cd /usr/src; \
    rm -r Python-3.7.2; \
    cd /usr/bin/; \
    ln -s /usr/local/bin/python3.7 python3; \
    cd /usr/bin; \
    ln -s /usr/local/bin/pip3.7 pip3

# Install python 3 dependencies
RUN pip3 install -U biopython==1.73 tabulate cgecore==1.4.4

#
# Install dependencies for SeqSero
#

# BLAST
RUN apt-get install cpio liblmdb-dev; \
    cd /usr/src; \
    wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.7.1/ncbi-blast-2.7.1+-src.tar.gz; \
    tar -zxf ncbi-blast-2.7.1+-src.tar.gz; \
    rm ncbi-blast-2.7.1+-src.tar.gz; \
    cd ncbi-blast-2.7.1+-src/c++; \
    ./configure && make && make install; \
    cd /usr/src; \
    rm -r ncbi-blast-2.7.1+-src

# Samtools (Compile without the ncurses lib)
RUN cd /usr/src; \
    wget https://sourceforge.net/projects/samtools/files/samtools/0.1.18/samtools-0.1.18.tar.bz2; \
    bzip2 -d samtools-0.1.18.tar.bz2; \
    tar -xf samtools-0.1.18.tar; \
    rm samtools-0.1.18.tar; \
    cd samtools-0.1.18; \
    sed -ri 's/^DFLAGS.+/DFLAGS=    -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -D_USE_KNETFILE -D_CURSES_LIB=0/g' Makefile; \
    sed -ri 's/^LIBCURSES/#  LIBCURSES/g' Makefile; \
    make; \
    mv samtools /usr/bin/; \
    cd /usr/src; \
    rm -r samtools-0.1.18

# BWA
RUN cd /usr/src; \
    wget https://sourceforge.net/projects/bio-bwa/files/bwa-0.7.17.tar.bz2; \
    bzip2 -d bwa-0.7.17.tar.bz2; \
    tar -xf bwa-0.7.17.tar; \
    rm bwa-0.7.17.tar; \
    cd bwa-0.7.17; \
    make; \
    mv bwa /usr/bin/; \
    cd /usr/src; \
    rm -r bwa-0.7.17

# SeqSero, SeqSero2 and biopython for python 2.7
RUN cd /usr/src; \
    pip install -U biopython; \
    git clone https://github.com/denglab/SeqSero2.git; \
    chmod 755 SeqSero2/SeqSero2_package.py; \
    git clone https://github.com/denglab/SeqSero.git; \
    chmod 755 SeqSero/SeqSero.py; \
    cd /usr/bin; \
    ln -s /usr/src/SeqSero2/SeqSero2_package.py SeqSero2_package.py; \
    ln -s /usr/src/SeqSero/SeqSero.py SeqSero.py

#
# SalmonellaTypeFinder
#
RUN cd /usr/src; \
    git clone https://bitbucket.org/genomicepidemiology/salmonellatypefinder.git; \
    chmod 755 salmonellatypefinder/SalmonellaTypeFinder.py; \
    cd /usr/bin; \
    ln -s /usr/src/salmonellatypefinder/SalmonellaTypeFinder.py SalmonellaTypeFinder.py

#
# Install MLST tool
#
# This section is left out of the light build from Dockerfile_light
RUN cd /usr/src; \
    git clone -b 2.0.4 https://bitbucket.org/genomicepidemiology/mlst.git; \
    chmod 755 mlst/mlst.py; \
    cd /usr/bin; \
    ln -s /usr/src/mlst/mlst.py mlst.py

# Setup .bashrc file for convenience during debugging
ENV PATH $PATH:/usr/src
RUN rm -rf /var/cache/apt/* /var/lib/apt/lists/*;
RUN echo "alias ls='ls -h --color=tty'\n" \
    "alias ll='ls -lrt'\n" \
    "alias l='less'\n" \
    "alias du='du -hP --max-depth=1'\n" \
    "alias cwd='readlink -f .'\n" \
    "PATH=$PATH\n" >> ~/.bashrc

WORKDIR /workdir

ENTRYPOINT ["/usr/src/salmonellatypefinder/SalmonellaTypeFinder.py"]
