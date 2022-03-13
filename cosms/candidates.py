#!/usr/bin/env python3
# coding: utf-8
"""Create a COS candidates file
"""
import argparse
import os
from cosms.core import inout

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = """
        Create csv file where candidate oligomers for MS analysis including 
        corresponding standards are defined.
        """,
        formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('filename',
                        type = str,
                        help = "file name or directory"
                        )
    parser.add_argument('--dp',
                        nargs = 2,
                        type = int,
                        default = [1, 6],
                        help = "minimum and maximum degree of polymerization"
                        )
    parser.add_argument('-m', '--monomers',
                        nargs = '+',
                        type = str,
                        choices = ['GlcNAc', 'A',
                                   'GlcN', 'D',
                                   'GlcNSuc', 'S',
                                   'GlcNGal', 'DGal',
                                   'GalA', 'G',
                                   'R',
                                   'Rstar',
                                   'GalAac',
                                   'GalAm', 'M',
                                   'GalAmd',
                                   'GlcNf', 'F',
                                   'Rha', 'Rhamnose',
                                   'RhaS', 'RhamnoseS',
                                   'Xyl', 'Xylose',
                                   'XylS', 'XyloseS',
                                   'Fuc','Fucose',
                                   'GalN', 'Galactosamin',
                                   'GalA', 'Galacturonic acid',
                                   'GlcA', 'Glucuronic acid',
                                   'ManN', 'Mannosamine',
                                   'Rib', 'Ribose',
                                   'Ara', 'Arabinose',
                                   'Man', 'Mannose',
                                   'Glc', 'Glucose',
                                   'Gal', 'Galactose',        
                                   'Ara', 'Arabinose',	
                                   'Fru', 'Fructose',
                                   'All', 'Allulose',
                                   'Suc', 'Sucrose',
                                   'Koji', 'Kojibiose'
				   'Hex'],
                        default = ['A', 'D'],
                        help = "monomeric units (separated by space), e.g. 'A D' or 'A R'"
                        )
    parser.add_argument('-a', '--adducts',
                        nargs = '+',
                        type = str,
                        choices = ['4H', '4Na', '4K', '3H', '3Na', '3K', '2H', 
                                   '2Na', '2K', 'H+NH4', 'H+Na', 'H+K', 'H', 
                                   'Na', 'K', 'NH4', '2Na-H', '2K-H', 'ACN+H', 
                                   'ACN+Na', '3H-', '2H-', 'H-', 'Na-2H', 
                                   'FA-H', 'Cl', 'HAc-H', 'TFA-H'],
#                         choices = ['H+', 'Na+', 'K+', 'NH4+'],
                        default = ['H'],
                        help = "Adduct ions"
                        )
#     parser.add_argument('-c', '--charges',
#                         nargs = '+',
#                         type = int,
#                         default = [1],
#                         help = "charge states"
#                         )
    parser.add_argument('-n', '--number-of-molecules',
                        nargs = '+',
                        type = int,
                        default = [1],
                        help = "can be used to check for dimeric ions -> '1 2'"
                        )
#     parser.add_argument('-l', '--label',
#                         nargs = '?',
#                         type = str,
#                         choices = ['18O', 'Ox', 'UDP', 'Meth', 'Methd'],
#                         default = None,
#                         help = "label at reducing end"
#                         )
    parser.add_argument('--redend_label',
                        nargs = '+',
                        type = str,
                        choices = ['18O', 'Ox', 'UDP', 'Meth', 'Methd', 'None',
                                   'Ox-2', 'Ox+16', 'PMP'],
                        default = [None],
                        help = "label at reducing end"
                        )
    parser.add_argument('--nonredend_label',
                        nargs = '+',
                        type = str,
                        choices = ['H2Oloss', 'None'],
                        default = [None],
                        help = "label at non-reducing end"
                        )
    parser.add_argument('-s', '--standard-mass',
                        nargs = '?',
                        type = float,
                        default = None,
                        help = "Added amount of standard in ng"
                        )
    parser.add_argument('-q', '--quant-standard',
                        nargs = '?',
                        type = str,
                        choices = ['Rstar', 'DR', 'A', 'R'],
                        default = 'Rstar',
                        help = "type of standard added for quantification"
                        )
    parser.add_argument('--mod',
                        nargs = '+',
                        choices = ['H2Oloss', 'HSO4', 'None'],
                        default = [None],
                        help = "Adduct ions"
                        )
    parser.add_argument('--redend_ftype',
                        nargs = '?',
                        type = str,
                        default = 'y',
                        choices = ['y', 'z'],
                        help = "Reducing end fragment type"
                        )
    parser.add_argument('--nonredend_ftype',
                        nargs = '?',
                        type = str,
                        default = 'b',
                        choices = ['b', 'c'],
                        help = "Non-reducing end fragment type"
                        )
    parser.add_argument('--min-rt',
                        nargs = '?',
                        type = float,
                        default = False,
                        help = "minimum retention time"
                        )
    parser.add_argument('--max-rt',
                        nargs = '?',
                        type = float,
                        default = False,
                        help = "maximum retention time"
                        )
    args = parser.parse_args()
    args.mod = [mod if mod != 'None' else None for mod in args.mod]
    args.nonredend_label = [_ if _ != 'None' else None 
                            for _ in args.nonredend_label]
    args.redend_label = [_ if _ != 'None' else None 
                            for _ in args.redend_label]
    
    default_file_name = 'COS_{m}_DP{dp}_{a}_{n}{s}{q}.csv'.format(
        m = ''.join(args.monomers),
        dp = args.dp[0] if args.dp[0] == args.dp[1] else '{0}-{1}'.format(args.dp[0],
                                                                     args.dp[1]),
        a = '_'.join(sorted(args.adducts)),
#         c = ''.join([str(charge) for charge in args.charges]),
        n = '_{0}-mers'.format('-'.join([str(nm) for nm in args.number_of_molecules])) 
            if args.number_of_molecules != [1] else '',
#         l = '_' + args.redend_label if args.redend_label else '',
#         l = '_' + args.label if args.label else '',
        s = '_{0:.0f}ng'.format(args.standard_mass) if args.standard_mass else '',
        q = '_{0}'.format(args.quant_standard) if args.standard_mass else ''
        )
    
    if os.path.isdir(args.filename):
        filename = os.path.join(args.filename, default_file_name)
    else:
        filename = args.filename
    
    candidates_file = inout.CandidatesFile(filename)
    candidates_file.create_all_possible( 
        monomers = args.monomers,
        min_dp = args.dp[0], 
        max_dp = args.dp[1], 
        adducts = args.adducts, 
#         charges = args.charges, 
        molecules = args.number_of_molecules, 
        modifications = args.mod,
        standard = args.quant_standard if args.standard_mass else None, 
        std_mass = args.standard_mass, 
        nonredend_labels = args.nonredend_label,
        redend_labels = args.redend_label,
#         label = args.label, #TODO: redundant REMOVE
        min_rt = args.min_rt,
        max_rt = args.max_rt, 
        nonredend_types = args.nonredend_ftype,
        redend_types = args.redend_ftype)

    candidates_file.write()
