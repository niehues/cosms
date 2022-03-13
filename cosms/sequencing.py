#!/usr/bin/env python3
# coding: utf-8
"""post-processing of MS1 results
"""
import argparse
from cosms.core import inout


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = """
        Sequencing of COS with given FragmentIonFile 
        """,
        formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('ms2results',
                        type = str,
                        help = "path to ms2 results file (fragment ions; output of msanalysis.py)"
                        )
    args = parser.parse_args()
    
    ms2_out_file = inout.FragmentFile(args.ms2results)
    ms2_out_file.sequencing()


