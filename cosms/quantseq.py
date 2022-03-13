#!/usr/bin/env python3
# coding: utf-8
"""combine quantitative MS1 results and sequencing results
"""
import argparse
from cosms.core import inout

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = """
        combine quantitative MS1 results and sequencing results 
        """,
        formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('filevscosresults',
                        type = str,
                        help = "path to ms1 results file (filevscos)"
                        )
    parser.add_argument('seqresults',
                        type = str,
                        help = "path to sequencing results file"
                        )
    parser.add_argument('fileoverview',
                        type = str,
                        help = "path to csv file with MS1 and MS2 file names"
                        )
    parser.add_argument('-m', '--monomers',
                        nargs = '+',
                        type = str,
                        choices = ['A', 'D', 'R', 'G', 'M'],
                        default = ['A', 'R'],
                        help = "monomeric units, e.g. 'A D' or 'A R'"
                        )
#     parser.add_argument('--plot',
#                         action = 'store_true',
#                         default = False,
#                         help = "create additional output file with plots (chromatograms, spectra); requires rpy2"
#                         )
    args = parser.parse_args()
    
    quantseq_file = inout.QuantSeqFile(args.filevscosresults,
                                       args.seqresults,
                                       args.fileoverview)
    for pa_length in (1, 2, 3):
        quantseq_file.pa(monomers = args.monomers,
                         pa_length = pa_length)

