SalmonellaTypeFinder
===================

This project documents SalmonellaTypeFinder service

## What is it?

Predicts Salmonella serotypes from WGS data, using SeqSero (https://github.com/denglab/SeqSero) or SeqSero2 (https://github.com/denglab/SeqSero2) and MLST type. The MLST type is found using the CGE MLST tool (https://bitbucket.org/genomicepidemiology/mlst). The MLST type is used to predict a serotype based on the known relation between serotypes found in the lab and their corresponding MLST type.

## Content of the repository
1. SalmonellaTypeFinder.py - the programme.
2. salmonellatypefinder    - folder containing the code for the programe.
3. data/db.json            - JSON formatted dictionary used to infer serotypes based on MLST.
4. data/output_headers.txt - Headers attached to output.
5. scripts                 - folder containing some in-house scripts.
6. Dockerfile              - dockerfile for running the service in docker (Note: Databases need to be mounted)
7. Dockerfile_light        - dockerfile for running the service in docker, but without the capability of finding MLST type.

## Installation using Docker

**Build time:** Docker image takes approximately 1 hour to build.

Setting up SalmonellaTypeFinder
```bash
# Go to wanted location for SalmonellaTypeFinder
cd /path/to/some/dir
# Clone and enter the salmonellatypefinder directory
git clone https://bitbucket.org/genomicepidemiology/salmonellatypefinder.git
cd salmonellatypefinder
```

The installation can either be with MLST typing capability or without. Choose without, if you intend to provide SalmonellaTypeFinder with an MLST type each time you run it, or if you have the CGE MLST tool installed, and intend to provide the path to it when invoking SalmonellaTypeFinder. Choose with, if you sometimes or always need SalmonellaTypeFinder to find the MLST type.

Build Docker container without MLST typing capability
```bash
# Rename Docker files
mv Dockerfile Dockerfile_complete
mv Dockerfile_light Dockerfile
# Build container
docker build -t salmonellatypefinder .
```

Build Docker container with MLST typing capability
Dependencies: KMA. If kma and kma_index has not bin installed please install
kma_index from the kma repository: https://bitbucket.org/genomicepidemiology/kma
```bash
# Build container
docker build -t salmonellatypefinder .
# Go to wanted location for MLST database
# Note: The location needs to be able to be attached to your docker container.
cd /path/to/db_dir
# Clone database from git repository
git clone https://bitbucket.org/genomicepidemiology/mlst_db.git
# Install MLST database with executable kma_index program
cd mlst_db
python3 INSTALL.py kma_index non_interactive
```

## Installation without using Docker

**Dependencies:**

SalmonellaTypeFinder per default assumes all dependencies are found in your path. If they are not, you need to specify the path to the dependency using the appropriate flag(s).

1. Python 3.7 with modules biopython, tabulate, and cgecore

2. Python 2.7 with module biopython (\*SeqSero)

3. CGE MLST tool version 2.0.4\*\*. https://bitbucket.org/genomicepidemiology/mlst

4. KMA version 1.2. https://bitbucket.org/genomicepidemiology/kma

5. BWA (\*SeqSero). (https://sourceforge.net/projects/bio-bwa/files/)

6. SAMtools version 0.1.18 (\*SeqSero). https://sourceforge.net/projects/samtools/files/samtools/0.1.18/samtools-0.1.18.tar.bz2

\*SeqSero: If you only want to use SeqSero2 these dependencies can be ignored.

\*\* CGE MLST tool Can be omitted if you intend to provide the service with the MLST type using the -st flag.


Install SalmonellaTypeFinder
```bash
# Go to wanted location for SalmonellaTypeFinder
cd /path/to/some/dir
# Clone the salmonellatypefinder directory
git clone https://bitbucket.org/genomicepidemiology/salmonellatypefinder.git
```

## Usage

The program can be invoked with the -h option to get help and more information
of the service.

```bash
usage: SalmonellaTypeFinder.py [-h] [-s {paired,single,assembled}]
                               [-o OUTPUT_TXT] [-t TMP_DIR] [-d JSON_MLST_DB]
                               [-m INT] [-f FRAC] [-st ST]
                               [--seromethod {seqsero,seqsero2}] [-p1 CGEMLST]
                               [-d1 CGEMLSTDB] [--python3 PYTHON3]
                               [--python2 PYTHON2] [--python2_env TXT]
                               [--seqsero SEQSERO] [--blastn BLASTN]
                               [--makeblastdb MAKEBLASTDB]
                               [--samtools SAMTOOLS] [--bwa BWA]
                               [--seqsero2 SEQSERO2]
                               FAST(Q|A) [FAST(Q|A) ...]

Given raw data or an assembly files outputs the Serotype and MLST.

positional arguments:
  FAST(Q|A)             Raw data in FASTQ format. Takes 1 (single-end) or 2
                        (paired-end) arguments.

optional arguments:
  -h, --help            show this help message and exit
  -s {paired,single,assembled}, --seq_type {paired,single,assembled}
                        Type of sequence: paired, single or assembled
  -o OUTPUT_TXT, --output OUTPUT_TXT
                        Path to file in which the results will be stored.
  -t TMP_DIR, --tmp_dir TMP_DIR
                        Temporary directory for storage of the results from
                        the external software.
  -d JSON_MLST_DB, --mlst_db JSON_MLST_DB
                        JSON formatted database used to predict serotypes from
                        MLST type. This option defaults to a database named
                        'db.json' located in the 'data' directory.
  -m INT, --mask_low_count_mlst INT
                        In the mlst<-->serovar database, ignore entries with
                        this number of isolates or fewer. This influences the
                        detailed MLST serovar output, but also the threshold
                        calculation, because the ignored entries will not be
                        included in the total sum of isolates with the given
                        MLST type and serovar. Default: 2
  -f FRAC, --fraction FRAC
                        Fraction of entries in mlst<-->serovar database that
                        needs to agree in order to call a serovar based on a
                        MLST type. Default: 0.75
  -st ST, --mlst ST     Optional. MLST type written as an integer. If given,
                        the programme will not find an MLST type but use the
                        one provided.
  --seromethod {seqsero,seqsero2}
                        Determines which version of SeqSero to use. Options
                        are 'seqsero' and 'seqsero2'. Note SeqSero2 is not yet
                        published and is currently still in the development
                        stage. Default: seqsero
  -p1 CGEMLST, --cgemlst_path CGEMLST
                        Path to cge mlst tool. Default: mlst.py
  -d1 CGEMLSTDB, --cgemlstdb_path CGEMLSTDB
                        Path to mlst database used for cge mlst tool.
  --python3 PYTHON3     Path to python3. Default: path to calling interpreter.
  --python2 PYTHON2     Path to python2.7. Default: python2.7
  --python2_env TXT     Path to a list of commands that will be executed just
                        before executing SeqSero. On most systems this will not
                        be necessary. The commands to be executed can set the
                        environment needed for SeqSero to run (e.g., the
                        python 2.7 environment).
  --seqsero SEQSERO     Path to SeqSero.py. Default: SeqSero.py
  --blastn BLASTN       Path to blastn. Default: blastn
  --makeblastdb MAKEBLASTDB
                        Path to makeblastdb. Default: makeblastdb
  --samtools SAMTOOLS   Path to samtools v. 0.18. Default: samtools
  --bwa BWA             Path to bwa. Default: bwa
  --seqsero2 SEQSERO2   Path to SeqSero2_package.py. Default:
                        SeqSero2_package.py
```

#### Example of use without Docker
```bash
SalmonellaTypeFinder.py -s paired -o results.txt \
    -d1 /path/to/mlst_db/ /path/to/data/sample*.gz
```

#### Example of use with Docker

```bash
# Running with a single fastq file
docker run --rm -v /path/to/files/:/path/to/files/ \
    -v $(pwd):/workdir \
    salmonellatypefinder -s paired \
    -o results.txt \
    -d1 /path/to/files/mlst_db/ /path/to/files/data/sample*.gz
```

## Web-server

A webserver implementing the methods is available at the [CGE
website](http://www.genomicepidemiology.org/) and can be found here:
https://cge.cbs.dtu.dk/services/SalmonellaTypeFinder/


## The Latest Version


The latest version can be found at
https://bitbucket.org/genomicepidemiology/salmonellatypefinder

## Documentation


The documentation available as of the date of this release can be found at
https://bitbucket.org/genomicepidemiology/salmonellatypefinder


Citation
=======

When using the method please cite:

Publication in final stages of writing


License
=======

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
