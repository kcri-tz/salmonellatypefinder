#! /tools/bin/python3

import subprocess
import re
import argparse
import os.path
import hashlib
import json
import textwrap
import gzip
from itertools import groupby


if __name__ == '__main__':

    #
    # Handling arguments
    #
    parser = argparse.ArgumentParser(description="Create JSON db for use with\
                                                  the MLST serotyper.")
    # Posotional arguments
    parser.add_argument("tab_db_file",
                        help="File containing the database in a tab seperated\
                              text file. Texfile is assumed to be encoded in\
                              utf-8.",
                        metavar='TAB_FILE')
    parser.add_argument("json_out",
                        help="Output path to write the JSON database/hash to.",
                        metavar='JSON_OUT')

    args = parser.parse_args()

    #
    # Load DB file.
    #
    output_hash = {}
    output_hash["ebg"] = {}  # stores eBG types.
    try:
        with open(args.tab_db_file, "r", encoding="utf-8") as tab_db_fh:
            header_line = tab_db_fh.readline()
            header_line = header_line.strip()
            headers = header_line.split("\t")

            index_st = ""
            index_serotype = ""
            index_ebg = ""
            index_subspecies = ""

            for i, header in enumerate(headers):
                header = header.lower()

                if(header == "st"):
                    index_st = i
                if(header == "ebg"):
                    index_ebg = i
                elif(header == "serovar"):
                    index_serotype = i
                elif(header == "subspecies"):
                    index_subspecies = i

            # Regexp to detect serotypes written as ex. "Salmonella Typhi"
            re_serotype_w_salm = re.compile(r"^salmonella\s*(.+)",
                                            re.IGNORECASE)

            for line in tab_db_fh:
                entries = line.split("\t")
                st = entries[index_st]
                ebg = entries[index_ebg]
                serotype = entries[index_serotype].lower()
                serotype = serotype.strip()

                # If serotype is blank, skip.
                if(not serotype):
                    continue

                # Filtering out an error in the enterodatabase.
                if(serotype == "shigella flexneri"):
                    continue

                # If required removes "Salmonella" from serotype.
                salm_match = re_serotype_w_salm.search(serotype)
                if(salm_match):
                    serotype = salm_match.group(1)

                # Ex.: changes saint-paul to saintpaul.
                serotype = serotype.replace("-", "")

                # Store the ST --> Serotype data
                if(st in output_hash):
                    if(serotype in output_hash[st]):
                        count = output_hash[st][serotype]
                        count = count+1
                        output_hash[st][serotype] = count
                    else:
                        output_hash[st][serotype] = 1
                else:
                    output_hash[st] = {}
                    output_hash[st][serotype] = 1

                # Store the eBG --> Serotype data
                if(ebg in output_hash["ebg"]):
                    if(serotype in output_hash["ebg"][ebg]):
                        count = output_hash["ebg"][ebg][serotype]
                        count = count+1
                        output_hash["ebg"][ebg][serotype] = count
                    else:
                        output_hash["ebg"][ebg][serotype] = 1
                else:
                    output_hash["ebg"][ebg] = {}
                    output_hash["ebg"][ebg][serotype] = 1

    except FileNotFoundError:
        print("The input file "+args.tab_db_file+" was not found\n")
        quit(1)

    #
    # Save JSON file
    #

    with open(args.json_out, 'w', encoding="utf-8") as json_fh:
        json.dump(output_hash, json_fh)

    print("# Wrote JSON hash to: "+args.json_out)

    quit(0)
