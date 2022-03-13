#!/usr/bin/env python3
# coding: utf-8
"""post-processing of MS1 results
"""
import argparse
from cosms.core import inout


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = """
        Process MS1 results from msanalysis.py: This script writes a table with
        the file names in rows and the COS in columns 
        """,
        formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('ms1results',
                        type = str,
                        help = "path to ms1 results file (output of msanalysis.py)"
                        )
    parser.add_argument('--snr',
                        nargs = '?',
                        type = float,
                        default = 3.0,
                        help = "minimum snr value"
                        )
    parser.add_argument('-a', '--adducts',
                        nargs = '+',
                        type = str,
                        choices = ['4H', '4Na', '4K', '3H', '3Na', '3K', '2H', 
                                   '2Na', '2K', 'H+NH4', 'H+Na', 'H+K', 'H', 
                                   'Na', 'K', 'NH4', '2Na-H', '2K-H', 'ACN+H', 
                                   'ACN+Na', '3H-', '2H-', 'H-', 'Na-2H', 
                                   'FA-H', 'Cl', 'HAc-H', 'TFA-H'],
                        default = False,
                        help = "considered adduct ions; if argument is not given all adduct ions are considered"
                        )
#     parser.add_argument('-c', '--charges',
#                         nargs = '+',
#                         type = int,
#                         default = False,
#                         help = "considered charge states; if argument is not given all charge states are considered"
#                         )
    parser.add_argument('--dp',
                        nargs = 2,
                        type = int,
                        default = False,
                        help = "minimum and maximum degree of polymerization"
                        )
    parser.add_argument('--mpl',
                        nargs = '?',
                        type = float,
                        default = False,
                        help = "maximum peak length [min]; everything eluting for a longer period of time is neglected"
                        )
    parser.add_argument('--minint',
                        nargs = '?',
                        type = float,
                        default = False,
                        help = "minimum relative intensity"
                        )
    parser.add_argument('--candidates',
                        nargs = '?',
                        type = str,
                        default = False,
                        help = "file with candidate oligomers (output from candidates.py); if given, all COS will be written to output file header"
                        )
#     parser.add_argument('--plot',
#                         action = 'store_true',
#                         default = False,
#                         help = "create additional output file with plots (chromatograms, spectra); requires rpy2"
#                         )
    args = parser.parse_args()
    
    ms1_out_file = inout.MS1File(args.ms1results)
    ms1_out_file.file_vs_cos(min_snr=args.snr, 
                             adducts=args.adducts, 
#                              charges=args.charges,
                             max_peak_length=args.mpl,
                             minint=args.minint,
                             minmaxdp=args.dp,
                             candidates_file=args.candidates)

