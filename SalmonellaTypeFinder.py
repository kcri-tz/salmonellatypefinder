#!/usr/bin/env python3

import argparse
import os.path
import re
import sys

from salmonellatypefinder.kauffmanwhite import KauffmanWhite
from salmonellatypefinder.mlst import MLST
from salmonellatypefinder.mlst2serotype import MLST2Serotype
from salmonellatypefinder.typingprofile import TypingProfile
from salmonellatypefinder.outputparser import Parser


def expand_and_check_path(path, force_dir=False):
    """ Returns path unaltered/unchecked if it doesn't contain a dir part.
        If path contains a dir part it finds the absolute path, checks the path
        and returns the absolute path if it exists.
    """
    dir_path = os.path.dirname(path)
    if(dir_path or force_dir):
        path = os.path.abspath(path)
        if(not os.path.exists(path)):
            sys.exit("Path not found: {}".format(path))
    return path


#
# Handling arguments
#
parser = argparse.ArgumentParser(description="Given raw data or an assembly\
    files outputs the Serotype and MLST.")
# Posotional arguments
parser.add_argument("input_files",
                    help="Raw data in FASTQ format. Takes 1 (single-end) or 2\
                          (paired-end) arguments.",
                    nargs='+',
                    metavar='FAST(Q|A)')
parser.add_argument("-s", "--seq_type",
                    help="Type of sequence: paired, single or assembled",
                    choices=["paired", "single", "assembled"],
                    default="paired")
parser.add_argument("-o", "--output",
                    help="Path to file in which the results will be stored.",
                    default=None,
                    metavar='OUTPUT_TXT')
parser.add_argument("-t", "--tmp_dir",
                    help="Temporary directory for storage of the results\
                          from the external software.",
                    default="SalmonellaTypeFinder_tmp_dir")
parser.add_argument("-d", "--mlst_db",
                    help="JSON formatted database used to predict serotypes\
                          from MLST type. This option defaults to a database\
                          named 'db.json' located in the 'data' directory.",
                    metavar='JSON_MLST_DB',
                    default=None)
parser.add_argument("-m", "--mask_low_count_mlst",
                    help="In the mlst<-->serovar database, ignore entries with\
                          this number of isolates or fewer. This influences\
                          the detailed MLST serovar output, but also the\
                          threshold calculation, because the ignored entries\
                          will not be included in the total sum of isolates\
                          with the given MLST type and serovar.\
                          Default: 2",
                    metavar="INT",
                    type=int,
                    default=2)
parser.add_argument("-f", "--fraction",
                    help="Fraction of entries in mlst<-->serovar database that\
                          needs to agree in order to call a serovar based on\
                          a MLST type.\
                          Default: 0.75",
                    default=0.75,
                    metavar="FRAC",
                    type=float)
parser.add_argument("-st", "--mlst",
                    help="Optional. MLST type written as an integer. If\
                          given, the programme will not find an MLST type\
                          but use the one provided.",
                    default=None,
                    metavar="ST",
                    type=int)
parser.add_argument("--seromethod",
                    help="Determines which version of SeqSero to use. Options\
                          are 'seqsero' and 'seqsero2'. Note SeqSero2 is not\
                          yet published and is currently still in the\
                          development stage.\
                          Default: seqsero",
                    choices=["seqsero", "seqsero2"],
                    default="seqsero")
parser.add_argument("-p1", "--cgemlst_path",
                    help="Path to cge mlst tool. Default: mlst.py",
                    metavar='CGEMLST',
                    default="mlst.py")
parser.add_argument("-d1", "--cgemlstdb_path",
                    help="Path to mlst database used for cge mlst tool.",
                    metavar='CGEMLSTDB',
                    default=None)
parser.add_argument("--python3",
                    help="Path to python3.\
                          Default: path to calling interpreter.",
                    default=None)
parser.add_argument("--python2",
                    help="Path to python2.7. Default: python2.7",
                    default="python2.7")
parser.add_argument("--python2_env",
                    help="Path to a list of commands that will be executed\
                          just before executing SeqSero. On most systems this\
                          won't be necessary. The commands to be executed can\
                          set the environment needed for SeqSero to run (e.g.,\
                          the python 2.7 environment).",
                    metavar='TXT',
                    default=None)
parser.add_argument("--seqsero",
                    help="Path to SeqSero.py. Default: SeqSero.py",
                    default="SeqSero.py")
parser.add_argument("--blastn",
                    help="Path to blastn. Default: blastn",
                    default="blastn")
parser.add_argument("--makeblastdb",
                    help="Path to makeblastdb. Default: makeblastdb",
                    default="makeblastdb")
parser.add_argument("--samtools",
                    help="Path to samtools v. 0.18. Default: samtools",
                    default="samtools")
parser.add_argument("--bwa",
                    help="Path to bwa. Default: bwa",
                    default="bwa")
parser.add_argument("--seqsero2",
                    help="Path to SeqSero2_package.py.\
                          Default: SeqSero2_package.py",
                    default="SeqSero2_package.py")

args = parser.parse_args()

# Check input files
input_files = []
if(args.input_files):
    if(len(args.input_files) > 2):
        sys.exit("! ERROR: Too many input arguments.")

    for filepath in args.input_files:
        filepath = os.path.abspath(filepath)
        if(not os.path.isfile(filepath)):
            sys.exit("! ERROR: Unable to locate input file: {}"
                     .format(filepath))
        input_files.append(filepath)
else:
    sys.exit("! ERROR: Too few input arguments.")

# Check tmp dir
args.tmp_dir = os.path.abspath(args.tmp_dir)
# Create tmp dir unless it already exists
os.makedirs(args.tmp_dir, exist_ok=True)

# Check JSON database.
if(not args.mlst_db):
    args.mlst_db = ("{}/data/db.json"
                    .format(os.path.dirname(os.path.realpath(__file__))))
if(not os.path.isfile(args.mlst_db)):
    print("JSON MLST database file not found:", args.mlst_db)
    sys.exit(1)

# Check programme paths

if(args.python3 is None):
    args.python3 = sys.executable

args.cgemlst_path = expand_and_check_path(args.cgemlst_path)
args.cgemlstdb_path = expand_and_check_path(args.cgemlstdb_path,
                                            force_dir=True)
args.python3 = expand_and_check_path(args.python3)
args.python2 = expand_and_check_path(args.python2)
args.seqsero = expand_and_check_path(args.seqsero)
args.seqsero2 = expand_and_check_path(args.seqsero2)
args.blastn = expand_and_check_path(args.blastn)
args.makeblastdb = expand_and_check_path(args.makeblastdb)
args.samtools = expand_and_check_path(args.samtools)
args.bwa = expand_and_check_path(args.bwa)

# Load database and create mlst2serotype object.
serotyper = MLST2Serotype(json_file=args.mlst_db,
                          min_sero_count=3,
                          min_frac=args.fraction,
                          mask_low_count=args.mask_low_count_mlst)

# SeqSero dependencies
seqsero_dependencies = {
    "seqsero": args.seqsero,
    "blastn": args.blastn,
    "makeblastdb": args.makeblastdb,
    "samtools": args.samtools,
    "bwa": args.bwa,
    "python2": args.python2
}

profile = TypingProfile(files=input_files,
                        mlst2serotype=serotyper,
                        seqtype=args.seq_type,
                        mlst=args.mlst,
                        tmp_dir=args.tmp_dir,
                        python2_env=args.python2_env,
                        cgemlst_path=args.cgemlst_path,
                        cgemlstdb_path=args.cgemlstdb_path,
                        python3=args.python3,
                        seqsero2=args.seqsero2,
                        seromethod=args.seromethod,
                        **seqsero_dependencies)

txt_output = Parser.output_txt(typing_profiles=[profile])

if(args.output):
    with open(args.output, "w", encoding="utf-8") as out_fh:
        out_fh.write(txt_output)
else:
    print(txt_output)

quit(0)
