#! /home/data1/tools/bin/Anaconda3-2.5.0/bin/python3

import argparse
import sys


if __name__ == '__main__':

    #
    # Handling arguments
    #
    parser = argparse.ArgumentParser(description="Create JSON db for use with\
                                                  the MLST serotyper.")
    # Posotional arguments
    parser.add_argument("input",
                        help="Text output from salmonella_typefinder.",
                        metavar='TAB_FILE')
    parser.add_argument("--web",
                        help="Web ID from the wrapper.",
                        metavar='WEB_ID')
    parser.add_argument("--no_qc",
                        help="No QC data.",
                        action="store_true",
                        default=False)
    parser.add_argument("--button",
                        help="HTML txt to create download button.",
                        metavar='HTML')
    parser.add_argument("--failed",
                        help="File containing failed analysis.",
                        metavar='FILE',
                        default=None)

    args = parser.parse_args()

    output_html = ("<center><h2>SalmonellaTypeFinder Results</h2></center>\n"
                   "<p>Notice: 'Flagged' samples contain an uncertain or no "
                   "prediction</p>\n"
                   "<hr>\n"
                   "<br>\n")

    with open(args.input, "r", encoding="utf-8") as input_fh:
        header_line = input_fh.readline()
        headers = header_line.split("\t")

        output_html += ("<table>\n" +
                        "    <thead>\n" +
                        "        <tr>\n")

        for header in headers:
            output_html += ("            <th><strong>" + header + "</strong>"
                            "</th>\n")

        output_html += ("        </tr>\n" +
                        "   </thead>\n" +
                        "   <tbody>\n")

        # Parses each Sample.
        counter_results = 0
        for line in input_fh:
            counter_results += 1
            entries = line.split("\t")
            output_html += "        <tr>"
            # Parses each result of the sample.
            for entry in entries:
                entry = entry.rstrip()
                output_html += "<td>" + entry + "</td>"
            output_html += "</tr>\n"

    output_html += ("    </tbody>\n" +
                    "</table>\n")
    output_html += "<br><br>\n"

    downloadScript = "https://cge.cbs.dtu.dk/cge/download_files2.php"

    if(args.button and counter_results > 0):
        output_html += "<td>" + args.button + "</td>\n"
    elif(counter_results > 0):
        output_html += "<td><form action='" + downloadScript + "' method='post'><input type='hidden' name='service' value='SalmonellaTypeFinder'><input type='hidden' name='version' value='1.4'><input type='hidden' name='filename' value='typeFinderResults.txt'><input type='hidden' name='pathid' value='" + args.web + "'><input type='submit' value='Text'></form></td>\n"

        if(not args.no_qc):
            output_html += "<br><br>\n"
            output_html += "<td><form action='" + downloadScript + "' method='post'><input type='hidden' name='service' value='SalmonellaTypeFinder'><input type='hidden' name='version' value='1.4'><input type='hidden' name='filename' value='qc_summary.txt'><input type='hidden' name='pathid' value='" + args.web + "'><input type='submit' value='QC table'></form></td>\n"

    if(args.failed):
        output_html += "<br><br>\n"
        output_html += "<h3>Failed samples:</h3>\n"
        output_html += "<table>\n"
        output_html += "  <thead>\n"
        output_html += "    <tr>"
        output_html += "      <th>Filename</th>"
        output_html += "      <th>Failed Analysis</th>"
        output_html += "    </tr>\n"
        output_html += "  </thead>\n"
        with open(args.failed, "r") as fh:
            fh.readline()  # Skip header
            for line in fh:
                line = line.rstrip()
                filename, failed_analysis, result = line.split("\t")
                output_html += "  <tr>\n"
                output_html += "    <td>" + filename + "</td>\n"
                output_html += "    <td>" + failed_analysis + "</td>\n"
                output_html += "  </tr>\n"
        output_html += "</table>\n"
        output_html += "<br><br>\n"
        output_html += "<td><form action='" + downloadScript + "' method='post'><input type='hidden' name='service' value='SalmonellaTypeFinder'><input type='hidden' name='version' value='1.4'><input type='hidden' name='filename' value='failed.tar.gz'><input type='hidden' name='pathid' value='" + args.web + "'><input type='submit' value='Partial results'></form></td>\n"

    print(output_html)
    print("Done", file=sys.stderr)

    quit(0)
