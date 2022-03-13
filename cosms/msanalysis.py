#!/usr/bin/env python3
# coding: utf-8
"""Analyze MS data in mzML file format
Targeted MS1 and MS2 analysis
"""
import argparse
import os
from cosms.core import mzml
from cosms.core import inout

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = """
        Analyze MS data in mzML file format; MS1 and MS2 analysis; 
        requires pymzml
        """,
        formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('mzml',
                        type = str,
                        help = "path to single mzML file or directory containing multiple mzML files"
                        )
    parser.add_argument('candidates',
                        type = str,
                        help = "file with candidate oligomers (output from candidates.py)"
                        )
    parser.add_argument('mstype',
                        nargs = '+',
                        type = str,
                        choices = ['ms1', 'ms2', 'ms2screen'],
                        help = "type(s) of MS analysis; usually 'ms1' OR 'ms2'"
                        )
    parser.add_argument('-p', '--precision',
                        nargs = '?',
                        type = float,
                        default = 0.12,
                        help = "m/z precision [Da]"
                        )
    parser.add_argument('--precursor',
                        nargs = '?',
                        type = float,
                        default = 0.3,
                        help = "m/z precision [Da] of precursor ion (only for MS2 measurements; fragment ion precision is controlled by --precision argument)"
                        )
    parser.add_argument('--mz-correction',
                        nargs = '?',
                        type = str,
                        default = False,
                        help = "csv file with correction values for each file in case of bad calibration of MS"
                        )
    parser.add_argument('--snr',
                        nargs = '?',
                        type = float,
                        default = 3.0,
                        help = "minimum snr value"
                        )
    parser.add_argument('-o', '--output',
                        nargs = '?',
                        type = str,
                        default = False,
                        help = "output directory"
                        )
    parser.add_argument('--plot',
                        action = 'store_true',
                        default = False,
                        help = "create additional output file with plots (chromatograms, spectra); requires rpy2"
                        )
    parser.add_argument('--deconv',
                        action = 'store_true',
                        default = False,
                        help = "isotopologue-(not charge-)deconvolution and correction for overlapping isotopologue peaks; use only for testing!"
                        )
    parser.add_argument('-f', '--force-rt',
                        action = 'store_true',
                        default = False,
                        help = "Use retention times exactly as given in candidates file"
                        )
    args = parser.parse_args()
    
    if not args.output:
        head, tail = os.path.split(args.mzml)
        output_dir = os.path.join(head, 
                                  '{0}_results'.format(tail.replace('.mzML', '')
                                                       ))
    else:
        output_dir = args.output
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    ms1out = os.path.split(args.candidates)[1].replace('.csv', '_MS1-results.csv')
    ms2out = os.path.split(args.candidates)[1].replace('.csv', '_MS2-results.csv')

    if 'ms1' in args.mstype:    
        mzml.ms1_analysis(mzml_files = args.mzml,
                          cos_candidates_file_name = args.candidates,
                          out_file_name = os.path.join(output_dir, ms1out),
                          precision_da = args.precision,
                          plotting = args.plot,
                          mz_correction = args.mz_correction,
                          min_snr = args.snr,
                          force_rt = args.force_rt,
                          deconv = args.deconv
                          )
        if not args.mz_correction:
            ms1_out_file = inout.MS1File(os.path.join(output_dir, ms1out))
            ms1_out_file.ms_precision()
    if 'ms2' in args.mstype:    
        # TODO sequencing options
        mzml.ms2_analysis(mzml_files = args.mzml,
                          cos_candidates_file_name = args.candidates,
                          out_file_name = os.path.join(output_dir, ms2out),
                          precision_da = args.precision,
                          precursor_precision_da = args.precursor,
                          plotting = args.plot,
                          mz_correction = args.mz_correction
                          )
    if 'ms2screen' in args.mstype:    
        mzml.ms2_screen(mzml_files = args.mzml,
                        cos_candidates_file_name = args.candidates,
                        out_file_name = os.path.join(output_dir, ms2out),
                        precision_da = args.precision,
                        precursor_precision_da = args.precursor,
                        plotting = args.plot,
                        mz_correction = args.mz_correction
                        )
        if not args.mz_correction:
            ms2_out_file = inout.FragmentFile(os.path.join(output_dir, ms2out))
            ms2_out_file.ms_precision()
    
    
    print("""
    --- Important Note ---
    Do not blindly trust the results of this script!!!
    It is crucial to check the raw MS data in Bruker DataAnalysis.
    Always check the results file of this script before proceeding with further 
    analyses. Results can be improved by setting time windows for each COS. 
    """)
