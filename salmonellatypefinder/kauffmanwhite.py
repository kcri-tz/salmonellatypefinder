#!/usr/bin/env python3

import subprocess
import re
import argparse
import os.path
import json
import sys
import tempfile


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


class KauffmanWhite():
    '''
    '''

    def __init__(self, files, method="seqsero", seqtype="paired",
                 tmp_dir="tmp_dir", python2_env=None, seqsero="SeqSero.py",
                 blastn="blastn", makeblastdb="makeblastdb",
                 samtools="samtools", bwa="bwa", python2="python2.7",
                 seqsero2="SeqSero2_package.py", python3="python3"):
        ''' Constructor.
            method: specifies what software to use in order to find the
                    Kauffman-White serotype profile. Only seqsero is currently
                    implemented.
            seqtype: of data can be either: paired, single, or assembled.
            files: Path to file(s) are given as a list.
        '''
        # SeqSero dependencies
        self.seqsero_path = seqsero
        self.blastn = blastn
        self.makeblastdb = makeblastdb
        self.samtools = samtools
        self.bwa = bwa
        self.python2 = python2

        # SeqSero2 dependencies
        self.seqsero2_path = seqsero2
        self.python3 = python3

        self.o_type = ""
        self.h1_type = ""
        self.h2_type = ""
        self.sdf = ""
        self.profile = ""
        self.serotypes = {}
        self.files = files
        self.method = None
        self.cmd = None  # The exact cmd executed to run external software.

        os.makedirs(tmp_dir, exist_ok=True)

        self.pre_cmd = ""
        if(python2_env):
            with open(python2_env, "r") as py_env_fh:
                for line in py_env_fh:
                    line = line.strip()
                    # Skip empty lines.
                    if(not line):
                        continue
                    # Skip lines starting with "#"
                    if(re.search(r"^#", line)):
                        continue
                    self.pre_cmd += line + "\n"

        if(method == "seqsero"):
            self.method = "seqsero"
            self.seqsero(tmp_dir, seqtype)
        elif(method == "seqsero2"):
            self.method = "seqsero2"
            self.seqsero(tmp_dir, seqtype)

    def serotype2string(self):
        ''' returns a string with the serotype result.
        '''
        serotypes = list(self.serotypes.keys())
        # More than one serotype.
        if(len(serotypes) > 1):
            return ", ".join(serotypes)
        # One serotype
        elif(len(serotypes) == 1):
            return serotypes[0]
        # No serotype
        else:
            return ""

    @staticmethod
    def add_prgdir_to_envpath(env_path, prg_path):
        """
        Adds a prg directory to the environment PATH, if a directory is found.
        """
        dir = os.path.dirname(prg_path)
        if(dir):
            dir = os.path.abspath(dir)
            return "{}:{}".format(env_path, dir)
        else:
            return env_path

    def seqsero2(self, working_dir, seqtype):
        """
        """
        # A temp directory is created, in which SeqSero will run.
        tmp_dir = tempfile.mkdtemp(prefix='seqsero2_tmp', dir=working_dir)

        if(os.path.dirname(self.seqsero2_path)):
            seqsero2_pre_cmd = ("{python3} {seqsero2}"
                                .format(python3=self.python3,
                                        seqsero2=self.seqsero2_path))
        else:
            seqsero2_pre_cmd = ("{seqsero2}"
                                .format(python3=self.python3,
                                        seqsero2=self.seqsero2_path))

        # Adding the paired-end specific options for SeqSero2
        if(seqtype == "paired"):
            seqsero2_post_cmd = ("-t 2 -i {} {}"
                                 .format(self.files[0], self.files[1]))
        # Adding the single-end specific options for SeqSero2
        elif(seqtype == "single"):
            seqsero2_post_cmd = ("-t 3 -i {}".format(self.files[0]))
        # Adding the assembled specific options for SeqSero2
        elif(seqtype == "assembled"):
            seqsero2_post_cmd = ("-t 4 -i {}".format(self.files[0]))

        seqsero2_cmd = ("{pre} {post}"
                        .format(pre=seqsero2_pre_cmd,
                                post=seqsero2_post_cmd))

        self.cmd = seqsero2_cmd

        eprint(seqsero2_cmd)

        # SeqSero creates files in the current working directory, with no
        # option to change output dir it is necessary to change the cwd.
        working_dir = os.getcwd()
        os.chdir(tmp_dir)
        try:
            result_raw = subprocess.run(seqsero2_cmd, capture_output=True,
                                        text=True, check=True, shell=True)
        except subprocess.CalledProcessError as e:
            eprint("ERROR: SeqSero2 call failed")
            eprint("CMD that failed: " + seqsero2_cmd)
            eprint("ERROR MSG: " + str(e))
            quit(1)
        # Change working dir back
        os.chdir(working_dir)

        print("ERROR: " + result_raw.stderr)
        print("OUT: " + result_raw.stdout)

        self.load_seqsero_result(result_raw)

    def load_seqsero_result(self, result_raw):
        """
        """
        # Parse SeqSero results
        re_o_type = re.compile(r"^O antigen prediction:	(.+)")
        re_h1_type = re.compile(r"^H1 antigen prediction\(fliC\):\s+(.+)")
        re_h2_type = re.compile(r"^H2 antigen prediction\(fljB\):\s+(.+)")
        re_sdf_type = re.compile(r"^Sdf prediction:(.+)")
        re_profile = re.compile(r"^Predicted antigenic profile:\s+(.+)")
        re_serotype = re.compile(r"^Predicted serotype\(s\):\s+([^*]+)")
        re_serotype_NA = re.compile(r"See comments below")
        re_serotype_NA2 = re.compile(r"N\/A")
        # It seems easiest to parse the screen output
        for line in result_raw.stdout.splitlines():
            match_o = re_o_type.search(line)
            if(match_o):
                self.o_type = match_o.group(1)
                continue
            match_h1 = re_h1_type.search(line)
            if(match_h1):
                self.h1_type = match_h1.group(1)
                continue
            match_h2 = re_h2_type.search(line)
            if(match_h2):
                self.h2_type = match_h2.group(1)
                continue
            match_profile = re_profile.search(line)
            if(match_profile):
                self.profile = match_profile.group(1)
                continue
            match_sdf = re_sdf_type.search(line)
            if(match_sdf):
                self.sdf = match_sdf.group(1)
                continue
            match_serotype = re_serotype.search(line)
            if(match_serotype):
                match_serotype_NA = re_serotype_NA.search(line,
                                                          re.IGNORECASE)
                # No serotype found
                if(match_serotype_NA):
                    serotype = "NF*"
                    self.serotypes[serotype] = 1
                    return
                match_serotype_NA2 = re_serotype_NA2.search(line,
                                                            re.IGNORECASE)
                # No serotype found
                if(match_serotype_NA):
                    serotype = self.profile
                    self.serotypes[serotype] = 1
                    return
                # Serotype(s) found
                serotypes = match_serotype.group(1).split(" ")
                for serotype in serotypes:
                    serotype = serotype.lower()
                    serotype = serotype.strip()
                    # serotype = serotype.casefold()
                    if(serotype != "or"):
                        self.serotypes[serotype] = 1
                return

    def seqsero(self, working_dir, seqtype):
        """
        """
        # A temp directory is created, in which SeqSero will run.
        tmp_dir = tempfile.mkdtemp(prefix='seqsero_tmp', dir=working_dir)

        # Create environment for SeqSero
        new_path_env = os.environ["PATH"]
        new_path_env = self.add_prgdir_to_envpath(new_path_env,
                                                  self.blastn)
        new_path_env = self.add_prgdir_to_envpath(new_path_env,
                                                  self.makeblastdb)
        new_path_env = self.add_prgdir_to_envpath(new_path_env,
                                                  self.samtools)
        new_path_env = self.add_prgdir_to_envpath(new_path_env,
                                                  self.bwa)
        new_path_env = self.add_prgdir_to_envpath(new_path_env,
                                                  self.python2)
        os.environ["PATH"] = new_path_env

        # Create SeqSero command.
        if(os.path.dirname(self.seqsero_path)):
            seqsero_pre_cmd = ("{python2} {seqsero} -b sam"
                               .format(python2=self.python2,
                                       seqsero=self.seqsero_path))
        else:
            seqsero_pre_cmd = ("{seqsero} -b sam"
                               .format(python2=self.python2,
                                       seqsero=self.seqsero_path))

        # Adding the paired-end specific options for SeqSero
        if(seqtype == "paired"):
            seqsero_post_cmd = ("-m 2 -i {} {}"
                                .format(self.files[0], self.files[1]))
        # Adding the single-end specific options for SeqSero
        elif(seqtype == "single"):
            seqsero_post_cmd = ("-m 3 -i {}".format(self.files[0]))
        # Adding the assembled specific options for SeqSero
        elif(seqtype == "assembled"):
            seqsero_post_cmd = ("-m 4 -i {}".format(self.files[0]))

        seqsero_cmd = ("{pre} {post}"
                       .format(pre=seqsero_pre_cmd,
                               post=seqsero_post_cmd))

        # Add pre-SeqSero command if necessary
        if(self.pre_cmd):
            seqsero_cmd = self.pre_cmd + seqsero_cmd

        self.cmd = seqsero_cmd

        eprint(seqsero_cmd)

        # SeqSero creates files in the current working directory, with no
        # option to change output dir it is necessary to change the cwd.
        working_dir = os.getcwd()
        os.chdir(tmp_dir)

        try:
            result_raw = subprocess.run(seqsero_cmd, capture_output=True,
                                        text=True, check=True, shell=True)
        except subprocess.CalledProcessError as e:
            eprint("ERROR: SeqSero call failed")
            eprint("CMD that failed: " + seqsero_cmd)
            eprint("ERROR MSG: " + str(e))
            eprint("ERROR STDOUT: " + e.stdout)
            eprint("ERROR STDERR: " + e.stderr)
            quit(1)
        # Change working dir back
        os.chdir(working_dir)

        self.load_seqsero_result(result_raw)


if __name__ == '__main__':

    #
    # Handling arguments
    #
    parser = argparse.ArgumentParser(description="Given raw data or an\
        assembly outputs the Serotype type.")
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
                        default="KauffmanWhite_tmp_dir")
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
    parser.add_argument("--python2",
                        help="Path to python2.7. Default: python2.7",
                        default="python2.7")
    parser.add_argument("--python3",
                        help="Path to python3.\
                              Default: path to calling interpreter.",
                        default=None)
    parser.add_argument("--seqsero2",
                        help="Path to SeqSero2_package.py.\
                              Default: SeqSero2_package.py",
                        default="SeqSero2_package.py")

    args = parser.parse_args()

    args.tmp_dir = os.path.abspath(args.tmp_dir)

    seqsero_dependencies = {
        "seqsero": args.seqsero,
        "blastn": args.blastn,
        "makeblastdb": args.makeblastdb,
        "samtools": args.samtools,
        "bwa": args.bwa,
        "python2": args.python2
    }

    if(args.python3 is None):
        args.python3 = sys.executable

    seqsero = KauffmanWhite((args.input_files[0], args.input_files[1]),
                            seqtype=args.seq_type,
                            tmp_dir=args.tmp_dir,
                            method="seqsero", seqsero2=args.seqsero2,
                            python3=args.python3,
                            **seqsero_dependencies)
    print("Serotype: " + seqsero.serotype2string())
    quit(0)
