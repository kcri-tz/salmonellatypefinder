#!/usr/bin/env python3

import argparse
import os.path
import gzip
import subprocess
import sys

from .kauffmanwhite import KauffmanWhite
from .mlst import MLST
from .mlst2serotype import MLST2Serotype, PredictedSerotype
from .outputparser import Parser


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


class TypingProfile():
    '''
    '''

    def __init__(self, files, mlst2serotype=None, seqtype="paired", mlst=None,
                 tmp_dir="tmp_dir", python2_env=None, cgemlst_path="mlst.py",
                 cgemlstdb_path=None, python3="python3", seqsero="SeqSero.py",
                 blastn="blastn", makeblastdb="makeblastdb",
                 samtools="samtools", bwa="bwa", python2="python2.7",
                 seqsero2="SeqSero2_package.py", seromethod="seqsero"):
        ''' Constructor.
        '''
        # SeqSero dependencies
        seqsero_dependencies = {
            "seqsero": seqsero,
            "blastn": blastn,
            "makeblastdb": makeblastdb,
            "samtools": samtools,
            "bwa": bwa,
            "python2": python2
        }

        self.mlst = ""
        self.kauffmanwhite = ""
        self.mlst_serotype = ""
        self.files = files
        self.serotype = ""
        self.uncertain_sero = False

        self.cgemlst_path = cgemlst_path
        self.cgemlstdb_path = cgemlstdb_path
        self.python3 = python3

        os.makedirs(tmp_dir, exist_ok=True)

        self.mlst = MLST((files[0], files[1]), seqtype=seqtype, mlst=mlst,
                         tmp_dir=tmp_dir, cgemlst_path=cgemlst_path,
                         cgemlstdb_path=cgemlstdb_path, python3_path=python3)

        self.kauffmanwhite = KauffmanWhite((files[0], files[1]),
                                           seqtype=seqtype, tmp_dir=tmp_dir,
                                           method=seromethod,
                                           python2_env=python2_env,
                                           seqsero2=seqsero2,
                                           **seqsero_dependencies)

        # Get serotype from MLST.
        if(mlst2serotype):
            self.mlst_serotype = mlst2serotype.mlst2serotype(self.mlst.st)

        # Get serotype from in silico KauffmanWhite.
        kauffwhite_sero = self.kauffmanwhite.serotype2string()

        # Both kauffmanwhite and MLST serotype are found.
        if(self.mlst_serotype.result and self.kauffmanwhite.serotypes):
            # MLST and kauffmanwhite agrees.
            if(self.mlst_serotype.result in self.kauffmanwhite.serotypes):
                self.serotype = self.mlst_serotype.result
            # The methods disagree, report kauffmanwhite
            else:
                self.uncertain_sero = True
                self.serotype = kauffwhite_sero
        # Only kauffmanwhite serotype is predicted.
        elif(self.kauffmanwhite.serotypes):
            self.uncertain_sero = True
            self.serotype = kauffwhite_sero
        # Only MLST serotype is predicted
        elif(self.mlst_serotype.result):
            self.uncertain_sero = True
            self.serotype = self.mlst_serotype.result
        # No serotype could be predicted
        else:
            self.uncertain_sero = True
            self.serotype = "n/a"


if __name__ == '__main__':

    #
    # Handling arguments
    #
    parser = argparse.ArgumentParser(description="Given raw data or an\
        assembly outputs the Serotype and MLST.")
    # Posotional arguments
    parser.add_argument("input_files",
                        help="Raw data in FASTQ format or an assembly in FASTA\
                              format.",
                        nargs='+',
                        metavar='FAST(Q|A)')
    parser.add_argument("-o", "--output",
                        help="Results will be written to this file. If no file\
                              is specified the results will be printed in a\
                              file called <input>.profile.txt",
                        default=None,
                        metavar="OUTPUT_TXT")
    parser.add_argument("-s", "--seq_type",
                        help="Type of sequence: paired, single or assembled",
                        choices=["paired", "single", "assembled"],
                        default="paired")
    parser.add_argument("-t", "--tmp_dir",
                        help="Temporary directory for storage of the results\
                              from the external software.",
                        default="TypingProfile_tmp_dir")
    parser.add_argument("-d", "--json_db",
                        help="",
                        metavar='JSON_DB')
    parser.add_argument("-m", "--mask_low_count_mlst",
                        help="In the mlst<-->serovar database, ignore entries\
                              with this number of isolates or fewer. This\
                              influences the detailed MLST serovar output, but\
                              also the threshold calculation, because the\
                              ignored entries will not be included in the\
                              total sum of isolates witht the given MLST type\
                              and serovar.",
                        metavar='INT',
                        type=int,
                        default=2)
    parser.add_argument("-st", "--mlst",
                        help="Optional. MLST type written as an integer. If\
                              given, the programme will not find an MLST type\
                              but use the one provided.",
                        default=None,
                        metavar="ST",
                        type=int)
    parser.add_argument("-p1", "--cgemlst_path",
                        help="Absolute path to cge mlst tool.\
                              Default: mlst.py.",
                        default="mlst.py")
    parser.add_argument("-d1", "--cgemlstdb_path",
                        help="Absolute path to cge mlst database.",
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
                              just before executing SeqSero. On most systems\
                              this won't be necessary. The commands to be\
                              executed can set the environment needed for\
                              SeqSero to run (e.g., the python 2.7\
                              environment).",
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

    args.tmp_dir = os.path.abspath(args.tmp_dir)

    if(args.python3 is None):
        args.python3 = sys.executable

    serotyper = MLST2Serotype(json_file=args.json_db,
                              min_sero_count=3,
                              min_frac=0.75,
                              mask_low_count=args.mask_low_count_mlst)

    eprint("Typing Profile:")
    subprocess.call("echo $0", shell=True)
    subprocess.call("printenv", shell=True)

    # SeqSero dependencies
    seqsero_dependencies = {
        "seqsero": args.seqsero,
        "blastn": args.blastn,
        "makeblastdb": args.makeblastdb,
        "samtools": args.samtools,
        "bwa": args.bwa,
        "python2": args.python2
    }

    profile = TypingProfile(files=(args.input_files[0], args.input_files[1]),
                            mlst2serotype=serotyper,
                            seqtype=args.seq_type,
                            mlst=args.mlst,
                            tmp_dir=args.tmp_dir,
                            python2_env=args.python2_env,
                            cgemlst_path=args.cgemlst_path,
                            cgemlstdb_path=args.cgemlstdb_path,
                            python3=args.python3,
                            **seqsero_dependencies)

    output_txt = Parser.output_txt(typing_profiles=[profile], headers=True)

    # Write output to file
    if(args.output):
        with open(args.output, "w", encoding="utf-8") as out_fh:
            out_fh.write(output_txt)
    else:
        args.output = os.path.basename(profile.files[0]) + ".profile.txt"
        with open(args.output, "w", encoding="utf-8") as out_fh:
            out_fh.write(output_txt)

    # Give "all" reading permissions to output file.
    try:
        subprocess.call(["chmod", "a+r", args.output])
    except subprocess.CalledProcessError:
        eprint("Warning: Couldn't change permissions of output file:\n"
               "\t{}\n".format(args.output))

    sys.exit()
