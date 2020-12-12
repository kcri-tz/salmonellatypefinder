#!/usr/bin/env python3

import subprocess
import re
import argparse
import os
import json
import gzip
import sys


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


class MLST():
    '''
    '''

    def __init__(self, files, method="default", seqtype="paired", mlst=None,
                 tmp_dir="tmp_dir", cgemlst_path="mlst.py",
                 cgemlstdb_path=None, python3_path="python3"):
        ''' Constructor.
            method: specifies what software to use in order to find the MLST
                    type. "default" is to employ SRST2 to reads and CGEMLST to
                    assembled genomes.
                    Other options are: cgemlst or srst2.
            seqtype: of data can be either: paired, single, or assembled.
            files: Path to file(s) are given as a list.
        '''
        self.cgemlst_path = cgemlst_path
        self.cgemlstdb = cgemlstdb_path
        self.python3 = python3_path
        self.alleles = {}
        self.st = None
        self.files = files
        self.method = None
        self.score = None  # Score depends on the method.
        self.cmd = None  # The exact cmd executed to run external software.

        os.makedirs(tmp_dir, exist_ok=True)

        if(mlst):
            self.st = mlst
        elif(method == "default"):
            if(seqtype == "paired" or seqtype == "single"):
                self.method = "CGE MLST"
                self.cgemlst(tmp_dir)
        elif(method == "cgemlst"):
            self.method = "CGE MLST"
            self.cgemlst(tmp_dir)

    def cgemlst(self, output):
        '''
        '''
        files = " ".join(self.files)
        if(os.path.dirname(self.cgemlst_path)):
            cmd = ("{python3} {mlst} -i {files} -o {out} -s senterica -p {db}"
                   .format(python3=self.python3, mlst=self.cgemlst_path,
                           files=files, out=output, db=self.cgemlstdb))
        else:
            cmd = ("{mlst} -i {files} -o {out} -s senterica -p {db}"
                   .format(mlst=self.cgemlst_path, files=files, out=output,
                           db=self.cgemlstdb))
        try:
            result_json = subprocess.run(cmd, capture_output=True, text=True,
                                         check=True, shell=True)
        except subprocess.CalledProcessError:
            eprint("ERROR: CGE MLST call failed")
            eprint("CMD that failed: " + cmd)
            quit(1)

        result_dict = json.loads(result_json.stdout)
        st = result_dict["mlst"]["results"]["sequence_type"]

        try:
            st = int(st)
        except ValueError:
            st = "unknown"

        self.st = st
        self.cmd = cmd


if __name__ == '__main__':

    #
    # Handling arguments
    #
    parser = argparse.ArgumentParser(description="Given raw data or an\
        assembly outputs the MLST type.")
    # Posotional arguments
    parser.add_argument("input_files",
                        help="Raw data in FASTQ format or an assembly in FASTA\
                              format.",
                        nargs='+',
                        metavar='FAST(Q|A)')
    parser.add_argument("-s", "--seq_type",
                        help="Type of sequence: paired, single or assembled",
                        choices=["paired", "single", "assembled"],
                        default="paired")
    parser.add_argument("-t", "--tmp_dir",
                        help="Temporary directory for storage of the results\
                              from the external software.",
                        default="MLST_tmp_dir")
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

    args = parser.parse_args()

    if(args.python3 is None):
        args.python3 = sys.executable

    profile = MLST(files=args.input_files,
                   seqtype=args.seq_type,
                   tmp_dir=args.tmp_dir,
                   cgemlst_path=args.cgemlst_path,
                   cgemlstdb_path=args.cgemlstdb_path,
                   python3_path=args.python3_path)

    print("ST " + str(profile.st))

    quit(0)
