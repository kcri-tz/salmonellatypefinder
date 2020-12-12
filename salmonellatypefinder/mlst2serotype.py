#!/usr/bin/env python3

import subprocess
import re
import argparse
import os.path
import hashlib
import json
import textwrap
import gzip
from itertools import groupby


class PredictedSerotype(dict):
    ''' Key: Serovar Val: (isolate_count, total_isolate_count, isolate_frac)
    '''

    def __init__(self):

        '''Constructor
        '''
        self.result = None

    def serotype2string(self):
        '''
        '''
        out_list = []

        for key in self:
            isolate_data = self[key]
            partial = key + " ({0:d}, {1:.2f})"
            partial = partial.format(isolate_data[0], isolate_data[2] * 100)
            out_list.append(partial)

        output_txt = " | ".join(out_list)
        print("JOIN VAR: " + output_txt)

        return output_txt


class MLST2Serotype():
    ''' An object that will predict Salmonella serotypes from MLST ST types.
    '''

    def __init__(self, json_file, min_sero_count=3, min_frac=0.75,
                 mask_low_count=0):

        ''' Constructor
        '''
        # Checking validity of options.
        if(min_sero_count <= mask_low_count):
            print(("ERROR: Min. serovar count must be greater than mask low"
                   " count."), file=sys.stderr)
            print("ERROR: Min. serovar count: " + str(min_sero_count),
                  file=sys.stderr)
            print("ERROR: Mask low count: " + str(mask_low_count),
                  file=sys.stderr)
            quit(1)

        self.min_sero_count = min_sero_count
        self.min_frac = min_frac
        self.mask_low_count = mask_low_count
        try:
            with open(json_file, "r", encoding="utf-8") as json_fh:
                self.data = json.load(json_fh)
        except FileNotFoundError:
            print("The JSON file {} was not found\n".format(json_file))
            quit(1)

    def mlst2serotype(self, st):
        ''' Given a ST type and thresholds predicts a serotype.
            RETURN: Dictionary. If a result above the thresholds was found a
                    key="result" is found with value=<predicted serotype>. The
                    Dictionary will always contain the serotypes found as keys:
                        key:<serotype> val:<List>
                        The List contains:
                            [0] isolates found with the <serotype>
                            [1] total isolates with the ST type
                            [2] fraction: [0]/[1]

            NOTE: If min_frac is set to 0.5 or less, then it is possible to
                  encounter an ST type that should yield 2 serovars. The
                  function will only output the first encountered as a result.
        '''
        out_result = PredictedSerotype()

        # TODO: Convert json datase to store integers and not strings.
        st = str(st)

        if(st in self.data):
            st_hash = self.data[st]
        else:
            return out_result

        # Get total number of isolates with given ST
        total_isolate_count = 0
        for key in st_hash.keys():
            key_count = st_hash[key]
            # Only include serovar count above the mask threshold.
            if(key_count <= self.mask_low_count):
                continue
            total_isolate_count += st_hash[key]

        if(total_isolate_count <= self.mask_low_count):
                return out_result

        # If more than one serotype is found
        if(len(st_hash.keys()) > 1):
            predicted_serotype = None
            max_count = 0
            max_frac = 0
            for key in st_hash.keys():
                isolate_count = st_hash[key]
                # Ignore serovars below the mask threshold.
                if(isolate_count <= self.mask_low_count):
                    continue
                isolate_frac = isolate_count / total_isolate_count
                out_result[key] = \
                    (isolate_count, total_isolate_count, isolate_frac)
                # Max counts can be equal if fraction is 0.5 or lower.
                if(isolate_count > max_count):
                    max_count = isolate_count
                    max_frac = isolate_frac
                    predicted_serotype = key
            # Check if results are above thresholds
            if(max_count >= self.min_sero_count and max_frac >= self.min_frac):
                out_result.result = predicted_serotype

        # At least "min_sero_count" isolates with one serotype
        elif(total_isolate_count >= self.min_sero_count):
            predicted_serotype = list(st_hash.keys())[0]
            out_result.result = predicted_serotype
            out_result[predicted_serotype] = \
                (total_isolate_count, total_isolate_count, 1.0)
        # No prediction due to too few isolates found. If count is higher than
        # mask_low_count the serotype is still recorded in the
        # PredictedSerotype object.
        else:
            low_count_serotype = list(st_hash.keys())[0]
            out_result[low_count_serotype] = \
                (total_isolate_count, total_isolate_count, 1.0)

        return out_result


if __name__ == '__main__':

    #
    # Handling arguments
    #
    parser = argparse.ArgumentParser(description="")
    # Posotional arguments
    parser.add_argument("-st", "--mlst",
                        help="ST type",
                        metavar='ST',
                        type=int)
    parser.add_argument("-d", "--json_db",
                        help="",
                        metavar='JSON_DB')

    args = parser.parse_args()

    serotyper = MLST2Serotype(json_file=args.json_db)
    results = serotyper.mlst2serotype(args.mlst)
    if(results.result):
        print("Predicted serotype: " + results.result)

    print("Details:")
    print("\tSerotype\tCount\tTotal\tFrac")
    for serotype in results:
        (count, total, frac) = results[serotype]
        print("\t{sero}\t{count:d}\t{tot:d}\t{frac:.2f}"
              .format(sero=serotype, count=count, tot=total, frac=frac))

    quit(0)
