#!/usr/bin/env python3
# coding: utf-8
"""cos
chitooligosaccharides
"""
import argparse
from cosms.core import chem

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = """
        Chitooligosaccharide elemental composition, mass, etc.
        """,
        formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('oligo',
                        type = str,
                        help = """Short name of oligosaccharide, e.g. 'A2D3'.
                            Available monosaccharides (chemical formula w/o H2O):
                            Hex (C6 H10 O5), 
                            GlcNAc or A (C8 H13 N1 O5), 
                            GlcN or D (C6 H11 N1 O4), 
                            GlcNf (C7 H11 N1 O5), 
                            R (C8 1H10 2H3 N1 O5), 
                            Rstar (12C6 13C2 1H10 2H3 N1 O5), 
                            GalA or G (C6 H8 O6), 
                            GalAm or M (C7 H10 O6), 
                            GalAmd (C7 1H7 2H3 O6), 
                            GalAmdc (12C6 13C1 1H7 2H3 O6), 
                            Rhamnose (C6 H10 O), 
                            RhamnoseS (C6 H10 O7 S1), 
                            Xylose (C5 H8 O4), 
                            XyloseS (C5 H8 16O S1), 
                            GlcA (C6 H8 O6)
                            """
                        )
    parser.add_argument('-a', '--adduct',
                        nargs = '?',
                        type = str,
                        choices = ['4H', '4Na', '4K', '3H', '3Na', '3K', '2H', 
                                   '2Na', '2K', 'H+NH4', 'H+Na', 'H+K', 'H', 
                                   'Na', 'K', 'NH4', '2Na-H', '2K-H', 'ACN+H', 
                                   'ACN+Na', '3H-', '2H-', 'H-', 'Na-2H', 
                                   'FA-H', 'Cl', 'HAc-H', 'TFA-H'],
                        default = 'H',
                        help = """Adduct ion 
                            (overall charge is derived from the ion, 
                            e.g. adduct 2H has charge 2+, 
                            adduct H- has charge 1-)"""
                        )
    parser.add_argument('-f', '--fragment-type',
                        nargs = '?',
                        type = str,
                        choices = ['b', 'c', 'y', 'z'],
                        help = "Fragment ion type"
                        )
    parser.add_argument('-l', '--label',
                        nargs = '?',
                        type = str,
                        choices = ['18O', 'Ox', 'UDP', 'Meth', 'Methd'],
                        help = "Label at reducing end"
                        )
    parser.add_argument('-n', '--number-of-molecules',
                        nargs = '?',
                        type = int,
                        default = 1,
                        choices = [1, 2],
                        help = "e.g. '2' for dimeric ions"
                        )
    parser.add_argument('--mod',
                        nargs = '+',
#                         type = str,
#                         choices = ['None', 'H2Oloss'],
                        default = False,
                        help = "Modification"
                        )
    args = parser.parse_args()
#     args.mod = [mod if mod != 'None' else None for mod in args.mod]
    
    ion = chem.Ion(chem.monomeric_composition(args.oligo),
#               charge = args.charge, 
              adduct = args.adduct, 
              redend_label = args.label,
              modification = args.mod,
              molecules = args.number_of_molecules,
              fragment_type = args.fragment_type)
    
    print("\n{0} {1}\n".format(ion.name, ion.precursor_type))
    print("Isotopic composition\t", ', '.join(['{0}: {1}'.format(k, v) for k, v 
        in sorted(ion.chem_comp.items())]))
    print("Components\t\t", ', '.join(['{0}: {1}'.format(k, v) for k, v 
        in sorted(ion.components.items())]))
    print("Potential sequences\t", ', '.join([''.join(i) for i in sorted(ion.sequences)]))
    print("Mass\t\t\t", round(ion.mass, 3))
    print("m/z\t\t\t", round(ion.mz, 3))
    if not args.fragment_type:
        print("Fragment ions")
        print("Frag. DP  m/z     Oligo  Ion-type")
        for f_type, ions in sorted(ion.fragment_ions.items()):
            if f_type[1] == 1:  
                print('{0}-ions'.format(f_type[0]))
            for f in sorted(ions, key = lambda x: abs(x.mz)):
                print("{frag}{dp} {c: >4.0f} {mz:>8.2f}  {oligo} \t {type}".format(
                    frag = f.fragment_type, 
                    dp = f.dp, 
                    c = f.charge,
                    mz = f.mz, 
                    oligo = f.name,
                    type = f.precursor_type
                    ))
        
