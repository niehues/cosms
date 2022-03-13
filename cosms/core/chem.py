#!/usr/bin/env python3
# coding: utf-8
"""chem
"""
import unittest
import re
import itertools
from collections import defaultdict as ddict
from cosms.core import misc


ATOMIC_MASSES = {
    # atomic masses [u]
    # proton and electron: Mohr et al. 2011, The 2010 CODATA Recommended Values
    # of the Fundamental Physical Constants
    # other: Wang et al. 2012, The Ame2012 atomic mass evaluation
    'p': 1.00727646681,
    'e-': 5.4857990946*10**-4,
    '1n': 1.008664915850,
    '1H': 1.007825032231,
    '2H': 2.014101778120,
    '12C': 12.0,
    '13C': 13.003354835071,
    '14N': 14.003074004426,
    '15N': 15.000108898884,
    '16O': 15.994914619566,
    '17O': 16.999131756500,
    '18O': 17.999159612858,
    '19F': 18.9984031627,
    '23Na': 22.9897692820,
    '28Si': 27.976926534649,
    '31P': 30.973761998417,
    '32S': 31.972071174408,
    '35Cl': 34.96885268,
    '39K': 38.963706486
    }


FORMULAS = {
    # Monosaccharide composition in glycosidic bond:
        'Fru': {'12C': 6, '1H': 10, '16O': 5}, #Fruktose
    'Fructose': {'12C': 6, '1H': 10, '16O': 5}, #Fruktose
        'All': {'12C': 6, '1H': 10, '16O': 5},#Allulose
    'Allulose': {'12C': 6, '1H': 10, '16O': 5},#Allulose
        'Fuc': {'12C':6, '1H': 10, '16O': 4}, #Fucose
    'Fucose': {'12C':6, '1H': 10, '16O': 4}, #Fucose
        'GalN': {'12C': 6, '1H': 11, '14N': 1, '16O': 4}, # Galactosamin
    'Galactosamin': {'12C': 6, '1H': 11, '14N': 1, '16O': 4}, # Galactosamin
        'GlcA': {'12C': 6, '1H': 8, '16O': 6}, # Glucoronic acid
    'Glucuronic acid': {'12C': 6, '1H': 8, '16O': 6}, # Glucoronic acid
        'ManN': {'12C': 6, '1H': 11, '14N': 1, '16O': 4}, # Mannoseamin
    'Mannosamine': {'12C': 6, '1H': 11, '14N': 1, '16O': 4}, # Mannoseamin
        'Rib': {'12C': 5, '1H': 8, '16O': 4}, # Ribose
    'Ribose': {'12C': 5, '1H': 8, '16O': 4}, # Ribose
        'Ara': {'12C': 5, '1H': 8, '16O': 4}, #Arabinose
    'Arabinose': {'12C': 5, '1H': 8, '16O': 4}, #Arabinose
	    'Man': {'12C': 6, '1H': 10, '16O': 5}, #Mannose
    'Mannose': {'12C': 6, '1H': 10, '16O': 5}, #Mannose
	    'Glc': {'12C': 6, '1H': 10, '16O': 5}, # Glucose
    'Glucose': {'12C': 6, '1H': 10, '16O': 5}, # Glucose
	    'Gal': {'12C': 6, '1H': 10, '16O': 5}, # Galcatose
    'Galactose': {'12C': 6, '1H': 10, '16O': 5}, # Galcatose
        'Hex': {'12C': 6, '1H': 10, '16O': 5}, # hexose e.g. glucose
        #Disaccharides
        'Suc': {'12C': 12, '1H': 20, '16O': 10}, # Sucrose
    'Sucrose': {'12C': 12, '1H': 20, '16O': 10}, # Sucrose
        'Koji': {'12C': 12, '1H': 20, '16O': 10}, # Kojibiose
    'Kojibiose': {'12C': 12, '1H': 20, '16O': 10}, # Kojibiose
        # chitosan
        'A': {'12C': 8, '1H': 13, '14N': 1, '16O': 5}, # GlcNAc
    'GlcNAc': {'12C': 8, '1H': 13, '14N': 1, '16O': 5}, # GlcNAc
        'D': {'12C': 6, '1H': 11, '14N': 1, '16O': 4}, # GlcN
    'Dtbdmsone': {'12C': 12, '1H': 25, '14N': 1, '16O': 4, '28Si': 1}, # GlcN
    'Dtbdmstwo': {'12C': 18, '1H': 39, '14N': 1, '16O': 4, '28Si': 2}, # GlcN
#         'Dtbdmsthree': {'12C': 24, '1H': 53, '14N': 1, '16O': 4, '28Si': 3}, # GlcN
        'GlcN': {'12C': 6, '1H': 11, '14N': 1, '16O': 4}, # GlcN
        'F': {'12C': 7, '1H': 11, '14N': 1, '16O': 5}, # GlcN formyl
        'GlcNf': {'12C': 7, '1H': 11, '14N': 1, '16O': 5}, # GlcN formyl
        'DGal':{'12C': 12, '1H': 23, '14N': 1, '16O': 9}, # GlcNGal reduktive Aminierung
        'GlcNGal':{'12C': 12, '1H': 23, '14N': 1, '16O': 9}, # GlcNGal reduktive Aminierung
        'R': {'12C': 8, '1H': 10, '2H': 3, '14N': 1, '16O': 5}, # reac. GlcN
        'Rstar': {'12C': 6, '13C': 2, '1H': 10, '2H': 3, '14N': 1, '16O': 5},
        'GlcNSuc': {'12C': 10, '1H': 15, '14N': 1, '16O': 7}, # succinylated GlcN
        'S': {'12C': 10, '1H': 15, '14N': 1, '16O': 7}, #GlcNSuc
        # pectin
        'GalA': {'12C': 6, '1H': 8, '16O': 6},
    'Galacturonic acid' : {'12C': 6, '1H': 8, '16O': 6},
        'GalAm': {'12C': 7, '1H': 10, '16O': 6},
        'GalAac': {'12C': 8, '1H': 10, '16O': 7}, # GalAm -H +COCH3
        'G': {'12C': 6, '1H': 8, '16O': 6},
        'M': {'12C': 7, '1H': 10, '16O': 6},
        'GalAmd': {'12C': 7, '1H': 7, '2H': 3, '16O': 6},
        'GalAmdc': {'12C': 6, '13C': 1, '1H': 7, '2H': 3, '16O': 6},
        # ulvan
    'Rhamnose': {'12C': 6, '1H': 10, '16O': 4},
        'Rha': {'12C': 6, '1H': 10, '16O': 4},
    'RhamnoseS': {'12C': 6, '1H': 10, '16O': 7, '32S': 1},
    'Xylose': {'12C': 5, '1H': 8, '16O': 4},
        'Xyl': {'12C': 5, '1H': 8, '16O': 4},
    'XyloseS': {'12C': 5, '1H': 8, '16O': 7, '32S': 1},
        'XylS': {'12C': 5, '1H': 8, '16O': 7, '32S': 1},
        # functional groups and other:
        'H2O': {'1H': 2, '16O': 1},
        'SO3': {'32S': 1, '16O': 3},
        'HSO4': {'32S': 1, '16O': 4, '1H': 1},
#         'TBDMS': {'1H': 14, '12C': 6, '28Si': 1},
    # adducts:
        # positive ion mode
        '4H': {'p': 4}, # 4+
        '4Na': {'23Na': 4, 'e-': -4}, # 4+
        '4K': {'39K': 4, 'e-': -4}, # 4+
        '3H': {'p': 3}, # 3+
        '3Na': {'23Na': 3, 'e-': -3}, # 3+
        '3K': {'39K': 3, 'e-': -3}, # 3+
        '2H': {'p': 2}, # 2+
        '2Na': {'23Na': 2, 'e-': -2}, # 2+
        '2K': {'39K': 2, 'e-': -2}, # 2+
        'H+NH4': {'14N': 1, '1H': 3, 'p': 2}, # 2+
        'H+Na': {'23Na': 1, 'e-': -1, 'p': 1}, # 2+
        'H+K': {'39K': 1, 'e-': -1, 'p': 1}, # 2+
        'H': {'p': 1}, # 1+
        'Na': {'23Na': 1, 'e-': -1}, # 1+
        'K': {'39K': 1, 'e-': -1}, # 1+
        'NH4': {'14N': 1, '1H': 3, 'p': 1}, # 1+
        '2Na-H': {'23Na': 2, 'e-': -2, 'p': -1}, # 1+
        '2K-H': {'39K': 2, 'e-': -2, 'p': -1}, # 1+
        'ACN+H': {'12C': 2, '1H': 3, '14N': 1, 'p': 1}, # 1+
        'ACN+Na': {'12C': 2, '1H': 3, '14N': 1, '23Na': 1, 'e-': -1}, # 1+
        # negative ion mode
        '-3H': {'p': -3}, # 3-
        '-2H': {'p': -2}, # 2-
        '-H': {'p': -1}, # 1-
        'Na-2H': {'23Na': 1, 'e-': -1, 'p': -2}, # 1-
        'FA-H': {'12C': 1, '1H': 2, '16O': 2, 'p': -1}, # 1-
#         'ACN+FA-H': {'12C': 3, '1H': 5, '16O': 2, '14N': 1, 'p': -1}, # 1-
        'Cl': {'35Cl': 1, 'e-': -1}, # 1-
        'HAc-H': {'12C': 2, '1H': 4, '16O': 2, 'p': -1}, # 1-
        'TFA-H': {'12C': 2, '1H': 1, '19F': 3, '16O': 2, 'p': -1}, # 1-
    # losses
        'H2Oloss': {'1H': -2, '16O': -1},
#         'CH3loss': {'12C': -1, '1H': -3},
    # precursor_type types (adduct/H2O not included):
        # non reducing end
        'b': {},
        'c': {'1H': 2, '16O': 1},
        # reducing end
        'y': {'1H': 2, '16O': 1},
        'z': {},
        'PMP175': {'12C': 10, '1H': 10, '14N': 2, '16O': 1, 'p': 1},
        'PMPhex': {'12C': 12, '1H': 12, '14N': 2, '16O': 2, 'p': 1},
    # labels at non-reducing end
        'TBDMS': {'1H': 14, '12C': 6, '28Si': 1},
    # labels at reducing end
        'Meth': {'1H': 2, '12C': 1},
        'Methd': {'1H': -1, '12C': 1, '2H': 3},
        '18O': {'16O': -1, '18O': 1},
        'Red': {'1H': 2}, # Reduction of reducing end with NaBH4
        'UDP': {'12C': 9, '1H': 12, '14N': 2, '16O': 11, '31P': 2},
        'PMP': {'12C': 20, '1H': 18, '14N': 4, '16O': 1}, # 2PMP-H2O, derivates ionize as [M+330+H]+ 
        'Ox-2': {'1H': -2},
        'Ox+16': {'16O': +1}
    }


ADDUCT2CHARGE = {
    # positive ion mode
    '4H': 4,
    '4Na': 4,
    '4K': 4,
    '3H': 3,
    '3Na': 3,
    '3K': 3,
    '2H': 2,
    '2Na': 2,
    '2K': 2,
    'H+NH4': 2,
    'H+Na': 2,
    'H+K': 2,
    'H': 1,
    'Na': 1,
    'K': 1,
    'NH4': 1,
    '2Na-H': 1,
    '2K-H': 1,
    'ACN+H': 1,
    'ACN+Na': 1,
    # negative ion mode
    '-3H': -3,
    '-2H': -2,
    '-H': -1,
    'Na-2H': -1,
    'FA-H': -1,
    'Cl': -1,
    'HAc-H': -1,
    'TFA-H': -1
#     'ACN+FA-H': -1
    }


SUBADDUCTS = {
    # Adducts that have to considered during fragmentation
    # positive ion mode
    '4H': ['3H', '2H', 'H'],
    '4Na': ['3Na', '2Na', 'Na'],
    '4K': ['3K', '2K', 'K'],
    '3H': ['2H', 'H'],
    '3Na': ['2Na', 'Na'],
    '3K': ['2K', 'K'],
    '2H': ['H'],
    '2Na': ['Na'],
    '2K': ['K'],
    'H+NH4': ['2H', 'H', 'NH4'],
    'H+Na': ['H', 'Na'],
    'H+K': ['H', 'K'],
    'H': [],
    'Na': [],
    'K': [],
    'NH4': ['H'],
    '2Na-H': [],
    '2K-H': [],
    'ACN+H': ['H'],
    'ACN+Na': ['Na'],
    # negative ion mode
    '-3H': ['-2H', '-H'],
    '-2H': ['-H'],
    '-H': [],
    'Na-2H': [],
    'FA-H': [],
    'Cl': [],
    'HAc-H': [],
    'TFA-H': []
#     'ACN+FA-H': []
    }


STANDARDS = {
    # For quantification
    # standard_name: {standard_monomer: [monomers_to_be_analyzed]}
    'Rstar': {'Rstar': ['A', 'R']},
    'DR': {'D': ['D'], 'R': ['A']},
    'A': {'A': ['A', 'R']},
    'R': {'R': ['A', 'R']},
    'GalAmdc': {'GalAmdc': ['GalAm', 'GalAmd']}
    }


component_masses = {
    # Masses of compounds defined in FORMULAS
    component: sum([
        ATOMIC_MASSES[atom]*count 
        for atom, count in FORMULAS[component].items()
        ])
    for component in FORMULAS
    }
# for k,v in sorted(component_masses.items()):
#     print(k,v)

def monomeric_composition(oligomer_name):
    """Transcribe COS name to monomeric composition
    Write e.g. AD as A1D1
    
    :param str oligomer_name: e.g. 'A2D2'
    
    :return mono_compo: e.g. {'A': 2, 'D': 2}
    :rtype: list
    """
    mono_compo = {}
    for monomer_count in re.findall('\D+\d+', oligomer_name.split('_')[0]):
        split = re.split('(\d+)', monomer_count)
        mono_compo[split[0]] = int(split[1])
    return mono_compo


def oligomer_name(mono_compo):
    """Transcribe monomeric composition to COS name
    
    :param str mono_compo: e.g. {'A': 2, 'D': 2}
    
    :return name: e.g. 'A2D2'
    :rtype: str
    """
    name = ''.join(['{0}{1}'.format(c, count) 
                    for c, count in sorted(mono_compo.items())
                    if count != 0])
    return name


def chemical_composition(mono_compo):
    """Get chemical (isotopic) composition from monomeric composition
    
    :param str mono_compo: e.g. {'A': 2, 'D': 2}; anything from FORMULAS
    
    :return chem_compo: e.g. {'12C': 28, '1H': 50, '14N': 4, '16O': 19}
    :rtype: dict
    """
    chem_compo = ddict(int)
    for c, comp_count in mono_compo.items():
        for atom, count in FORMULAS[c].items():
            chem_compo[atom] += comp_count * count
    return chem_compo


def cos_to_standard(cos, standard):
    try:
        std_composition = {std_monomer: sum([cos[monomer] for monomer in
                                             STANDARDS[standard][std_monomer]])
                           for std_monomer in STANDARDS[standard].keys()}
    except KeyError:
        exit("\n!!! ERROR - Standard '{0}' is not suitable for COS '{1}' !!!\n".format(
            standard, cos.name))
    std = Ion(std_composition,
#               charge = cos.charge, 
              adduct = cos.adduct, 
              redend_label = cos.redend_label,
              molecules = cos.molecules)
    return std

  
def set_mass_lookup(lookup):
    def modify_func(func):
        def molecular_mass(composition):
            mass = sum([lookup[component] * count for component, count in 
                        composition.items()])
            return mass
        return molecular_mass
    return modify_func


@set_mass_lookup(ATOMIC_MASSES)
def molecular_mass_atomic(composition):
    return


@set_mass_lookup(component_masses)
def molecular_mass_components(composition):
    return


def isotopologue_correction(cos, monoiso_intensity):        
    """Correction of peak intensity - deconvolution
    Konstanten fuer Lineargleichungen (y=ax+b, y = DP, x = Intensitaet des 
    monoisotopischen Peaks) zur Berechnung der Anteile des monoisotopischen und 
    +3 peaks an Gesamtintensitaet; Werte von Stefan"""
    monoiso_A_a = -0.0722260488
    monoiso_A_b = 0.9266254621
    monoiso_R_a = -0.0880561333
    monoiso_R_b = 0.9356683
    monoiso_Rstar_a = -0.0801947136
    monoiso_Rstar_b = 0.91495402
    
#     monoiso_A_a = -0.0765106626 # Values determined by script
#     monoiso_A_b = 0.9438030706
#     monoiso_R_a = -0.0553336067
#     monoiso_R_b = 0.8614152754
#     monoiso_Rstar_a = -0.0681137078
#     monoiso_Rstar_b = 0.8709124509
    try:
        a = cos['A']
    except:
        a = 0
    try:
        r = cos['R']
    except:
        r = 0
    try:
        rstar = cos['Rstar']
    except:
        rstar = 0
    corr_quotient = ( a/cos.dp * (monoiso_A_a*cos.dp+monoiso_A_b) +
                      r/cos.dp * (monoiso_R_a*cos.dp+monoiso_R_b) +
                      rstar/cos.dp * (monoiso_Rstar_a*cos.dp+monoiso_Rstar_b) )
    corr_intensity = monoiso_intensity / corr_quotient
    return corr_intensity


def overlap_correction(cos1_a, cos1_r, cos1_intensity, cos2_intensity):
    """Correct for overlapping isotopologue peaks
    cos1 = cos with smaller m/z, its +3 peak intensity will be substracted from cos2
    cos1_intensity = intensity of cos1
    cos2_intensity = intensity of cos with higher m/z, it is to be corrected
    """
    plus3_A_a = 0.010184562
    plus3_A_b = -0.0074847707
    plus3_R_a = 0.0094192546
    plus3_R_b = -0.0130978324
    
#     plus3_A_a = 0.0065444157 # Values determined by script
#     plus3_A_b = -0.0022365069
#     plus3_R_a = 0.0033990976
#     plus3_R_b = 0.0021401249
    cos1_dp = cos1_a + cos1_r
    plus3_proportion = ( cos1_a/cos1_dp * (plus3_A_a*cos1_dp+plus3_A_b) +
                         cos1_r/cos1_dp * (plus3_R_a*cos1_dp+plus3_R_b) )
    corr_intensity = cos2_intensity - (cos1_intensity * plus3_proportion)
    return corr_intensity


class Oligo(dict):
    """Chitooligosaccharide (COS)

    :param dict oligo_dict: Monomeric composition of chitooligosaccharide

    .. attribute:: dp

        Degree of polymerization (DP)
        
    .. attribute:: name
    
        Short name of COS, e.g. 'A2D2'
        
    .. attribute:: sugars
    
        List of monomers 
        
    .. attribute:: components
    
        Monomeric composition and additional components, e.g. 'H2O'
    """
    def __init__(self, oligo_dict):
        super().__init__(oligo_dict)
        self.dp = sum(self.values())
        self.name = oligomer_name(self)
        self.sugars = sorted(self.keys())
        self.components = ddict(int, self)
        self.components['H2O'] += 1
        
    
    @misc.LazyFunction
    def chem_comp(self):
        return chemical_composition(self.components)
    
    
    @misc.LazyFunction
    def mass(self):
        return molecular_mass_components(self.components)
        
 
    @misc.LazyFunction
    def sequences(self):
        #TODO: this should return a list to be usable for pectin 
        """Calculate all possible sequences based on monomeric composition
        
        :return: Possible sequences
        :rtype: set
        """
        _sequences = set()
        num_sugar_types = len(self)
        if num_sugar_types == 1:
#             _sequences.add(self.sugars[0] * self[self.sugars[0]])
            _sequences.add(tuple([self.sugars[0] for _ in range(self[self.sugars[0]])]))
        elif num_sugar_types == 2:
            indices = set(range(self.dp))
            for selected_indices in itertools.combinations(indices, 
                                                           self[self.sugars[0]]):
#                 sequence = ''
                sequence = []
                for index in range(self.dp):
                    if index in selected_indices:
#                         sequence += self.sugars[0]
                        sequence.append(self.sugars[0])
                    else:
#                         sequence += self.sugars[1]
                        sequence.append(self.sugars[1])
                _sequences.add(tuple(sequence))
        else: # NOTE this is expensive - don't use it for high DPs
#             a_sequence = ''.join([sugar * count for sugar, count in self.items()])
            a_sequence = []
            for sugar, count in self.items():
                a_sequence += [sugar for _ in range(count)]
#             _sequences = set([''.join(seq) for seq in itertools.permutations(a_sequence)])
            _sequences = set([tuple(seq) for seq in itertools.permutations(a_sequence)])
        return _sequences
    
    
class Ion(Oligo): 
    """Chitooligosaccharide ion (precursor or fragment ion)
    
    :param dict oligo_dict: Monomeric composition of chitooligosaccharide

    .. attribute:: dp

        Degree of polymerization (DP)
        
    .. attribute:: name
    
        Short name of COS, e.g. 'A2D2'
        
    .. attribute:: sugars
    
        List of monomers 
        
    .. attribute:: components
    
        Monomeric composition and additional components, e.g. 'H2O'
    """
    def __init__(self, 
                 oligo_dict,
#                  charge = None, 
                 adduct = None, 
                 nonredend_label = None,
                 redend_label = None,
                 modification = None,
                 molecules = 1,
                 fragment_type = None,
                 precursor_sequences = None,
                 nonredend_types = ['b', 'c'],
                 redend_types = ['y', 'z'],
                 nonred_maxonly = True):
        Oligo.__init__(self, oligo_dict)
        if adduct.endswith('-'):
            adduct = '-'+adduct[:-1]
        self.charge = ADDUCT2CHARGE[adduct]
        self.adduct = adduct
        self.redend_label = redend_label
        self.nonredend_label = nonredend_label
        if self.nonredend_label:
            self.name += '_{0}'.format(self.nonredend_label)
        self.modification = modification
        self.molecules = molecules
        self.fragment_type = fragment_type
        self.precursor_sequences = precursor_sequences
        self.nonredend_types = nonredend_types
        self.redend_types = redend_types 
        self.nonred_maxonly = nonred_maxonly
        
        for component, count in self.components.items():
            self.components[component] *= self.molecules
        self.components[adduct] += 1 #self.charge*ADDUCT2CHARGE[adduct]
        if self.fragment_type:
            self.components[self.fragment_type] += 1
            self.components['H2O'] -= 1
        if self.redend_label:
            self.components[self.redend_label] += 1
            if self.redend_label == 'UDP':
                self.name = 'UDP-' + self.name
        if self.nonredend_label:
            self.components[self.nonredend_label] += 1 
        if self.modification:
            for mod in self.modification:
                self.components[mod] += 1
        
        self._min_rt = None
        self._max_rt = None
        self._standard = None
        self._ng = None


    def correct(self, value):
        self.mz += value
        self.plus1_mz += value
        self.plus2_mz += value
        self.mass += value
    

    @misc.LazyFunction
    def mz(self):
        return abs(self.mass / self.charge)
    
    
    @misc.LazyFunction
    def plus1_mz(self):
        return self.mz + ((ATOMIC_MASSES['13C']-ATOMIC_MASSES['12C']) / abs(self.charge))
    
    
    @misc.LazyFunction
    def plus2_mz(self):
        return self.mz + (2 * ((ATOMIC_MASSES['13C']-ATOMIC_MASSES['12C']) / abs(self.charge)))


    @misc.LazyFunction
    def precursor_type(self):
        _precursor_type = '[{m}M{mod}{sign1}{adduct}]{charge}{sign2}'.format(
            m = self.molecules if self.molecules > 1 else '',
            mod = self.modification if self.modification else '',
            sign1 = '+' if self.adduct[0] != '-' else '',
            adduct = self.adduct,#.replace('+', ''),
            charge = abs(self.charge) if abs(self.charge) > 1 else '',
            sign2 = '+' if self.charge > 0 else '-')
        return _precursor_type


    @misc.LazyFunction
    def fragment_ions(self, 
#                       redend_types = ['y'],
#                       nonredend_types = ['b'],
                        min_fragment_dp = 1
                      ):
        """All possible fragment ions based on possible sequences
        
#         :param list redend_types: Types of fragment ions, reducing end ['y', 'z']
#         :param list nonredend_types: Types of fragment ions, non-reducing end ['b', 'c']
        :param int min_fragment_dp: Minimum DP of fragment ions
        
        :return: Fragment ions {('y', 2): [<fragment precursor_type #1>, ...]}
        :rtype: dict
        """
        _fragment_ions = ddict(list)
        nonredend_fragments = ddict(set)
        redend_fragments = ddict(set)
        for prec_seq in self.sequences:
            for idx in range(min_fragment_dp, self.dp):
                nonredend_seq = prec_seq[:idx]
                nonredend_oligo = tuple(nonredend_seq.count(sugar) for 
                                        sugar in self.sugars)
                nonredend_fragments[nonredend_oligo].add(prec_seq)
                redend_seq = prec_seq[idx:]
                redend_oligo = tuple(redend_seq.count(sugar) for 
                                     sugar in self.sugars)
                redend_fragments[redend_oligo].add(prec_seq)
        for fragments, types, label in (
            (nonredend_fragments, self.nonredend_types, self.nonredend_label),
            (redend_fragments, self.redend_types, self.redend_label)
            ):
            adducts = [self.adduct] + SUBADDUCTS[self.adduct]
            molecules_s = list(range(1, self.molecules+1)) 
            for frag_oligo, frag_type, adduct, m in itertools.product(
                    fragments, types, adducts, molecules_s):  
                frag_ion = Ion(
                    {sugar: count for sugar, count in 
                        zip(self.sugars, frag_oligo) if count > 0},
                    adduct = adduct,
                    molecules = m,
                    redend_label = label if frag_type in ['x','y','z'] else None,
                    nonredend_label = label if frag_type in ['a','b','c'] else None,
                    fragment_type = frag_type,
                    precursor_sequences = fragments[frag_oligo]
                    )
                _fragment_ions[frag_type, frag_ion.dp].append(frag_ion)
        return _fragment_ions


    @misc.LazyFunction #@property
    def key(self):
        return round(self.mz, 3), self.name 
    
    
    @property
    def min_rt(self):
        return self._min_rt
    @min_rt.setter
    def min_rt(self, value):
        self._min_rt = value
        return 
    
    
    @property
    def max_rt(self):
        return self._max_rt
    @max_rt.setter
    def max_rt(self, value):
        self._max_rt = value
        return 
    
    
    @property
    def standard(self):
        return self._standard
    @standard.setter
    def standard(self, value):
        self._standard = value
        return 
    
    
    @property
    def ng(self):
        return self._ng
    @ng.setter
    def ng(self, value):
        self._ng = value
        return
    
    
class TestChem(unittest.TestCase):
    
    def test_monomeric_composition(self):
        self.assertTrue({'A': 2, 'D': 2}, 
                        monomeric_composition('A2D2'))
        self.assertTrue({'A': 2, 'D': 2, 'Rstar': 1}, 
                        monomeric_composition('A2D2Rstar1'))
        self.assertTrue({'A': 2, 'D': 2, 'Rstar': 1}, 
                        monomeric_composition('A2D2Rstar'))
        
    def test_oligomer_name(self):
        self.assertTrue('A2D2',
                        oligomer_name({'A': 2, 'D': 2}))
        self.assertTrue('A2D2Rstar1',
                        oligomer_name({'A': 2, 'D': 2, 'Rstar': 1}))
        
    def test_chemical_composition(self):
        self.assertTrue({'12C': 28, '1H': 50, '14N': 4, '16O': 19}, 
                        chemical_composition({'A': 2, 'D': 2}))
    
    def test_molecular_mass_atomic(self):
        self.assertTrue(
            746.348919581958, 
            molecular_mass_atomic({'12C': 28, '1H': 50, '14N': 4, '16O': 19})
            )
        
    def test_molecular_mass_components(self):
        self.assertTrue(746.348919581958, 
                        molecular_mass_components({'A': 2, 'D': 2}))


class TestOligo(unittest.TestCase):
    def setUp(self):
        self.oligo = Oligo({'A': 2, 'D': 2})
        
    def test_oligo(self):
        self.assertEqual(2, self.oligo['A'])
        self.assertEqual(2, self.oligo['D'])
        self.assertEqual(4, self.oligo.dp)
        self.assertEqual('A2D2', self.oligo.name)
        self.assertEqual({'12C': 28, '1H': 50, '14N': 4, '16O': 19}, 
                        self.oligo.chem_comp)
        self.assertAlmostEqual(746.3069254, self.oligo.mass)
#         self.assertEqual(set(['AADD', 'ADAD', 'ADDA', 'DAAD', 'DADA', 'DDAA']), 
        self.assertEqual(set([('A', 'D', 'D', 'A'),
                              ('D', 'A', 'A', 'D'),
                              ('A', 'A', 'D', 'D'),
                              ('D', 'A', 'D', 'A'),
                              ('D', 'D', 'A', 'A'),
                              ('A', 'D', 'A', 'D')]),
                        self.oligo.sequences)
        

class TestIon(unittest.TestCase):
    def setUp(self):
        self.precursor = Ion({'A': 2, 'D': 2},
#                        charge = 2, 
                       adduct = '2H', 
                       redend_label = '18O',
                       molecules = 1,
                       fragment_type = None)
 
    def test_oligo(self):
        self.assertEqual(2, self.precursor['A'])
        self.assertEqual(2, self.precursor['D'])
        self.assertEqual(4, self.precursor.dp)
        self.assertEqual('[M+2H]2+', self.precursor.precursor_type)
        self.assertEqual('A2D2', self.precursor.name)
        self.assertAlmostEqual(749.31844686+1.00727646681, 
                               self.precursor.mass)
        self.assertAlmostEqual(375.16286166, self.precursor.mz)
        self.assertEqual(
            {'12C': 28, '1H': 50, '14N': 4, '16O': 18, '18O': 1, 'p': 2}, 
            self.precursor.chem_comp)
        self.assertEqual(set([('A', 'D', 'D', 'A'),
                              ('D', 'A', 'A', 'D'),
                              ('A', 'A', 'D', 'D'),
                              ('D', 'A', 'D', 'A'),
                              ('D', 'D', 'A', 'A'),
                              ('A', 'D', 'A', 'D')]),
                        self.precursor.sequences)
        self.assertEqual(12, len(self.precursor.fragment_ions))
        std = cos_to_standard(self.precursor, 'DR')
        self.assertEqual(std['R'], self.precursor['A'])
        self.assertEqual(std['D'], self.precursor['D'])


class TestIon2(unittest.TestCase):
    def setUp(self):
        self.precursor = Ion({'A': 1, 'D': 1},
#                        charge = 1, 
                       adduct = 'H', 
                       redend_label = None,
                       molecules = 2,
                       fragment_type = None)
 
    def test_oligo(self):
        self.assertEqual(1, self.precursor['A'])
        self.assertEqual(1, self.precursor['D'])
        self.assertEqual(2, self.precursor.dp)
        self.assertEqual('[2M+H]+', self.precursor.precursor_type)
        self.assertEqual('A1D1', self.precursor.name)
        self.assertAlmostEqual(765.32476655, self.precursor.mass)
        self.assertAlmostEqual(765.32476655, self.precursor.mz)
        self.assertEqual({'12C': 28, '1H': 52, '14N': 4, '16O': 20, 'p': 1}, 
                        self.precursor.chem_comp)
        std = cos_to_standard(self.precursor, 'DR')
        self.assertEqual(std['R'], self.precursor['A'])
        self.assertEqual(std['D'], self.precursor['D'])


class TestNames(unittest.TestCase):
    def setUp(self):
        self.oligo1 = Ion({'Hex': 2, 'GalA': 1},
                          adduct = 'FA-H', 
                          molecules = 1)
        self.oligo2 = Ion({'Hex': 1},
                          adduct = 'FA-H', 
                          molecules = 2)
        self.oligo3 = Ion({'Hex': 1},
                          adduct = '3H', 
                          molecules = 1)
        self.oligo4 = Ion({'Hex': 1},
                          adduct = '-2H', 
                          molecules = 1)
        
    def test_oligo(self):
        self.assertEqual('[M+FA-H]-', self.oligo1.precursor_type)
        self.assertEqual('[2M+FA-H]-', self.oligo2.precursor_type)
        self.assertEqual('[M+3H]3+', self.oligo3.precursor_type)
        self.assertEqual('[M-2H]2-', self.oligo4.precursor_type)
 
 
class TestPMPlabeledMonosaccharide(unittest.TestCase):
    def setUp(self):
        self.glc = Ion({'Hex': 1},
                        adduct = 'H', 
                        molecules = 1,
                        redend_label = 'PMP',
                        redend_types = ['PMP175', 'PMPhex'])
        self.glcn = Ion({'GlcN': 1},
                          adduct = 'H', 
                          molecules = 1,
                          redend_label = 'PMP')
        
    def test_oligo(self):
        self.assertAlmostEqual(511.2, self.glc.mz, places = 1)
        self.assertAlmostEqual(510.2, self.glcn.mz, places = 1)
        for k,f in self.glc.fragment_ions.items():
            print(k, f)