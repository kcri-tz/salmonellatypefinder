#!/usr/bin/env python3

import argparse
import os.path
import gzip
import sys

from .kauffmanwhite import KauffmanWhite
from .mlst import MLST
from .mlst2serotype import MLST2Serotype, PredictedSerotype


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


class Parser():
    '''
    '''

    @staticmethod
    def output_txt(typing_profiles, headers=True):
        if(headers):
            output_txt = "Sample\tPredicted Serotype\tST\tST mismatches\t"
            output_txt += "ST sero prediction\tSeqSero prediction\tO-type\t"
            output_txt += "H1-type\tH2-type\tMLST serotype details\tFlagged\n"
        else:
            output_txt = ""

        for profile in typing_profiles:
            # Sample
            filename = os.path.basename(profile.files[0])
            output_txt += filename
            # Predicted Serotype
            if(profile.serotype):
                output_txt += "\t" + profile.serotype
            else:
                output_txt += "\t" + "Unable to predict"
            # ST
            if(isinstance(profile.mlst.st, str)):  # Is it a string?
                output_txt += "\t" + profile.mlst.st
            elif(profile.mlst.st is None):
                output_txt += "\tNone"
            elif(profile.mlst.st):
                output_txt += "\t" + str(profile.mlst.st)
            else:
                output_txt += "\tNone"
            # ST mismatches
            output_txt += "\t"  # + profile.mlst.score
            # output_txt += str(profile.mlst.score, encoding='utf-8')
            # ST sero prediction
            if(profile.mlst_serotype.result):
                output_txt += "\t" + profile.mlst_serotype.result
            else:
                output_txt += "\tUnable to predict"
            # SeqSero prediction
            output_txt += "\t" + profile.kauffmanwhite.serotype2string()
            # O-type
            output_txt += "\t" + profile.kauffmanwhite.o_type
            # H1-type
            output_txt += "\t" + profile.kauffmanwhite.h1_type
            # H2-type
            output_txt += "\t" + profile.kauffmanwhite.h2_type
            # MLST serotype details
            mlst_serotype_details = profile.mlst_serotype.serotype2string()
            eprint("DETAILS: \"" + mlst_serotype_details + "\"")
            if(mlst_serotype_details):
                output_txt += "\t" + mlst_serotype_details
            else:
                output_txt += "\tNo serotypes"
            # Flagged
            if(profile.uncertain_sero):
                output_txt += "\t*"

            output_txt += "\n"

        return output_txt

    @staticmethod
    def output_html(output_txt):
        start_output = \
            '''
            <table style="width:100%">
            <tr>
                <th>Sample</th>
                <th>Predicted Serotype</th>
                <th>ST</th>
                <th>ST mismatches</th>
                <th>ST sero prediction</th>
                <th>SeqSero prediction</th>
                <th>O-type</th>
                <th>H1-type</th>
                <th>H2-type</th>
            </tr>
            '''

        lines = output_txt.split("\n")
        table_rows = []

        for line in lines:
            entries = line.split()

            row_str = "<tr>\n"

            # Sample name
            row_str += "\t<td>"
            row_str += entries[0]
            row_str += "</td>\n"
            # Predicted Serotype
            row_str += "\t<td>"
            row_str += entries[1]
            row_str += "</td>\n"
            # ST
            row_str += "\t<td>"
            row_str += entries[2]
            row_str += "</td>\n"
            # ST mismatches
            row_str += "\t<td>"
            row_str += entries[3]
            row_str += "</td>\n"
            # ST sero prediction
            row_str += "\t<td>"
            row_str += entries[4]
            row_str += "</td>\n"
            # SeqSero prediction
            row_str += "\t<td>"
            row_str += entries[5]
            row_str += "</td>\n"
            # O-type
            row_str += "\t<td>"
            row_str += entries[6]
            row_str += "</td>\n"
            # H1-type
            row_str += "\t<td>"
            row_str += entries[7]
            row_str += "</td>\n"
            # H2-type
            row_str += "\t<td>"
            row_str += entries[8]
            row_str += "</td>\n"

            row_str += "</tr>\n"

            table_rows.append(row_str)

        end_output = \
            '''
            </table>
            '''

        return start_output + "\n".join(table_rows) + end_output


if __name__ == '__main__':

    #
    # Handling arguments
    #
    parser = argparse.ArgumentParser(description="Given raw data or an\
        assembly outputs the Serotype and MLST.")
    # Posotional arguments
    parser.add_argument("input_files",
                        help="Input files are the txt output from the Parser\
                              class or salmonella_typefinder without headers.",
                        nargs='+',
                        metavar='TXT')
    parser.add_argument("-of", "--out_format",
                        help="Output format. Only HTML is implemented.",
                        choices=["html"],
                        default="html")
    parser.add_argument("-o", "--output",
                        help="Path to output file.")

    args = parser.parse_args()
    args.output = os.path.abspath(args.output)

    conc_txt = ""

    for file_txt in args.input_files:
        if(not os.path.isfile(file_txt)):
            eprint("Input file not found:", file_txt)
            quit(1)

        with open(file_txt, "r", encoding="utf-8") as fh:
            for line in fh:
                conc_txt += line

    if(args.out_format == "html"):
        print(Parser.output_html(conc_txt))

    quit(0)
