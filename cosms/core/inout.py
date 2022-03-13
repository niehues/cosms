#!/usr/bin/env python3
# coding: utf-8
"""inout
"""
import unittest
import os
import csv
import itertools
from collections import defaultdict as ddict
from operator import attrgetter
from cosms.core import chem
from cosms.core import misc
from cosms.core import seq
# from cosms.candidates import candidates_file


def write_csv_file(cls):
    def _write_csv_file(cls):
        try:
            cls.fieldnames += list(sorted(cls.additional_fieldnames))
        except:
            pass
        out_file = csv.DictWriter(open(cls.file_path, 
                                       'w',
                                       newline = ''), 
                                  cls.fieldnames,
                                  dialect = 'excel') #USA/UK format
        out_file.writeheader()
        for rowdict in cls.rows:
            out_file.writerow(rowdict)
        print("Write {0} at \"{1}\"".format(cls.__class__.__name__, 
                                            cls.file_path))
    return _write_csv_file


def check_csv_dialect(cls):
    def _check_csv_dialect(cls):
        csv_file = open(cls.file_path, 'r')
        sample = csv_file.read(1024)
        dialect = csv.Sniffer().sniff(sample, delimiters = [',', ';', '\t'])
        return dialect
    return _check_csv_dialect
    

class CandidatesFile():
    def __init__(self, file_path):
        self.file_path = file_path
        self.fieldnames = ['Oligo', 
                           'Standard', 
                           'Standard mass [ng]', 
                           'M',
                           'Adduct', 
                           'Charge', 
                           'm/z',
#                            'Label',#TODO: redundant REMOVE
                           'Non-red. end label',
                           'Red. end label', 
                           'Modification',
                           'DP',
                           'Min RT', 
                           'Max RT',
                           'Non-red. end fragment ions',
                           'Red. end fragment ions'
                           ]
        self.rows = []
        self.cos = []
        self.standards = {}
        
        
    def create_all_possible(self, 
               monomers = ['A', 'D'],
               min_dp = 1,
               max_dp = 6,
               adducts = ['H', '2H'],
#                charges = [1, 2],
               molecules = [1],
               modifications = None,
               standard = None, #e.g. 'Rstar' or 'DR'
               std_mass = None,
               nonredend_labels = [],
               redend_labels = [],
               label = None, #TODO: redundant REMOVE
               min_rt = None,
               max_rt = None,
               nonredend_types = ['b', 'c'],
               redend_types = ['y', 'z']):
        if not modifications: 
            modifications = [None]
        if label:
            redend_labels = [label]
        nonredend_types_4file = ', '.join(nonredend_types)
        redend_types_4file = ', '.join(redend_types)
        for dp in range(min_dp, max_dp + 1):
            for combo in itertools.combinations_with_replacement(monomers, dp):
                oligo_dict = {m: combo.count(m) for m in monomers}
                for adduct, mol, mod, nonredend_label, redend_label in itertools.product(
                    adducts,
                    molecules,
                    modifications,
                    nonredend_labels,
                    redend_labels):
                    _cos = chem.Ion(
                        oligo_dict,
#                                    charge = charge, 
                        adduct = adduct, 
                        nonredend_label = nonredend_label,
                        redend_label = redend_label,
                        molecules = mol,
                        modification = mod,
                        nonredend_types = nonredend_types,
                        redend_types = redend_types)
                    if standard:
                        std_cos = chem.cos_to_standard(_cos, standard)
                        if _cos.name == std_cos.name:
                            continue  
                        std_cos.ng = float(std_mass)
                        self.standards[std_cos.key] = std_cos
                        _cos.standard = std_cos.key
                    if min_rt:
                        _cos.min_rt = min_rt
                        _cos.max_rt = max_rt
                    self.cos.append(_cos)
                    self.rows.append({
                        'Oligo': _cos.name, 
                        'Standard': std_cos.name 
                        if standard else '', 
                        'Standard mass [ng]': std_mass 
                        if std_mass else '', 
                        'M': _cos.molecules,
                        'Adduct': _cos.adduct 
                            if not _cos.adduct.startswith('-')
                            else _cos.adduct[1:]+'-', 
                        'Charge': _cos.charge, 
#                                       'Label': _cos.redend_label, #TODO: remove, obsolete 
                        'Non-red. end label': _cos.nonredend_label,
                        'Red. end label': _cos.redend_label,    
                        'Modification': mod
                        if mod else '',
                        'DP': _cos.dp,
                        'Min RT': _cos.min_rt, 
                        'Max RT': _cos.max_rt,
                        'Non-red. end fragment ions': nonredend_types_4file,
                        'Red. end fragment ions': redend_types_4file,
                        'm/z': '{0:.3f}'.format(_cos.mz)})
        
        
    def read(self):
        in_file = csv.DictReader(open(self.file_path, 'r'),
                                 dialect = self.check())
        for _row in in_file:
            if _row['Oligo'] == '':
                continue
            oligo_dict = chem.monomeric_composition(_row['Oligo'])
            mod = [_row['Modification']
                   ] if 'Modification' in _row and _row['Modification'
                                                        ] != '' else None
            nonredend_label = _row['Non-red. end label'
                ] if 'Non-red. end label' in _row and _row['Non-red. end label'
                                                           ] != '' else None
            redend_label = _row['Red. end label'
                ] if 'Red. end label' in _row and _row['Red. end label'
                                                       ] != '' else None
            nonred_maxonly = False if 'nonred_maxonly' in _row and _row[
                'nonred_maxonly'].lower() not in ['', 'true', 'yes', 'y'] else True
            nonredend_types = []
            if 'Non-red. end fragment ions' in _row:
                for t in _row['Non-red. end fragment ions'].split(', '):
                    if t:
                        nonredend_types.append(t)
            redend_types = []
            if 'Red. end fragment ions' in _row:
                for t in _row['Red. end fragment ions'].split(', '):
                    if t:
                        redend_types.append(t)
            _cos = chem.Ion(
                oligo_dict,
#                             charge = int(_row['Charge']), 
                adduct = _row['Adduct'], 
                nonredend_label = nonredend_label,
                redend_label = redend_label,
                molecules = int(_row['M']),
                modification = mod,
                nonredend_types = nonredend_types,
                redend_types = redend_types,
                nonred_maxonly = nonred_maxonly)
            if _row['Standard']:
                std_composition = chem.monomeric_composition(_row['Standard'])
                std_cos = chem.Ion(
                                std_composition,
#                                charge = _cos.charge, 
                               adduct = _cos.adduct, 
                               redend_label = _cos.redend_label,
                               molecules = _cos.molecules,
                               modification = mod)
                std_cos.ng = float(_row['Standard mass [ng]'])
                try:
                    std_cos.min_rt = float(_row['Min RT'])
                    std_cos.max_rt = float(_row['Max RT'])
                except ValueError:
                    pass
#                 std_cos.standard = std_cos.key
                self.standards[std_cos.key] = std_cos
                _cos.standard = std_cos.key
            try:
                _cos.min_rt = float(_row['Min RT'])
                _cos.max_rt = float(_row['Max RT'])
            except ValueError:
                pass
            self.rows.append(_row)
            self.cos.append(_cos)
        return self.cos
    
    
    @write_csv_file
    def write(self):
        return
    
    @check_csv_dialect
    def check(self):
        return
    
    
    def clear(self):
        self.rows = []
        self.cos = []
        self.standards = {}
        
    
class MassFile():
    def __init__(self, file_path):
        self.file_path = file_path
        self.fieldnames = ['m/z', 
                           'Oligo', 
                           'Adduct', 
                           'Charge', 
                           'n-mer',
                           'DP',
                           'Non-red. end label',
                           'Red. end label',
                           'Modification']
        self.rows = []
    
    
    def add_candidates_file(self, filepath):
        candidates_file = CandidatesFile(filepath)
        cos_list = candidates_file.read()
#         for cos in sorted(cos_list, key=attrgetter('mz')):
        for cos in sorted(cos_list, key=lambda x: abs(x.mz)):
            self.rows.append({'m/z': '{0:.3f}'.format(cos.mz), 
                           'Oligo': cos.name, 
                           'Adduct': cos.adduct
                                if not cos.adduct.startswith('-')
                                else cos.adduct[1:]+'-', 
                           'Charge': cos.charge, 
                           'n-mer': cos.molecules, 
                           'DP': cos.dp,
                           'Non-red. end label': cos.nonredend_label,
                           'Red. end label': cos.redend_label,
                           'Modification': cos.modification})
    
    
    @write_csv_file
    def write(self):
        return
    
    @check_csv_dialect
    def check(self):
        return
    

class MS1File():
    def __init__(self, file_path):
        self.file_path = file_path
        self.fieldnames = ['File name', 
                           'Oligo', 
                           'Ion type', 
                           'Intensity [arb. unit]', 
                           'Mass [ng]',
                           'Amount of substance [nmol]',
                           'SNR',
                           'Number of spectra',
                           'Peak start time [min]',
                           'Peak end time [min]',
                           'Peak length [min]',
                           'Measured m/z',
                           'Calculated m/z',
                           'Adduct',
                           'Charge', 
                           'n-mer',
                           'DP']
#                            'A', 
#                            'D',
#                            'R',
#                            'Rstar']
        self.rows = []

    
    def add(self, cos, entry):
        entry['Oligo'] = cos.name
        entry['Ion type'] = cos.precursor_type
        entry['Calculated m/z'] = round(cos.mz, 3)
        entry['Adduct'] = cos.adduct if not cos.adduct.startswith('-') else cos.adduct[1:]+'-'
        entry['Charge'] = cos.charge
        entry['n-mer'] = cos.molecules
        entry['DP'] = cos.dp
        if 'Mass [ng]' in entry:
            entry['Amount of substance [nmol]'] = (entry['Mass [ng]'] /
                                                   cos.mass) 
        for monomer, count in cos.items():
            if monomer not in self.fieldnames:
                self.fieldnames.append(monomer)
            entry[monomer] = count
        self.rows.append(entry)


    def overlap_corr(self):
        """Correct for overlapping isotopologue peaks
        """
        lookup = ddict(dict) # DP - Calculated m/z - idx
        for idx, row in enumerate(self.rows):
            lookup[int(row['DP'])][float(row['Calculated m/z'])] = idx
        for dp, mz_dict in lookup.items():
            if len(mz_dict) > 1:
                mz_list = list(sorted(mz_dict))
                for i in range(1, len(mz_list)):
                    if 2.7 < mz_list[i]-mz_list[i-1] < 3.3: 
                        idx2 = lookup[dp][mz_list[i]]
                        cos2_intensity = float(self.rows[idx2]['Intensity [arb. unit]'])
                        idx1 = lookup[dp][mz_list[i-1]]
                        cos1_intensity = float(self.rows[idx1]['Intensity [arb. unit]'])
                        try:
                            cos1_a = float(self.rows[idx1]['A'])
                        except:
                            cos1_a = 0
                        try:
                            cos1_r = float(self.rows[idx1]['R'])
                        except:
                            cos1_r = 0
                        corr_intensity = chem.overlap_correction(cos1_a, 
                                                                 cos1_r, 
                                                                 cos1_intensity, 
                                                                 cos2_intensity)
                        corr_factor = corr_intensity / cos2_intensity
                        self.rows[idx2]['Intensity [arb. unit]'] = '{0:e}'.format(corr_intensity)
                        try:
                            self.rows[idx2]['Mass [ng]'] *= corr_factor
                        except:
                            continue
                        try:
                            self.rows[idx2]['Amount of substance [nmol]'] *= corr_factor 
                        except:
                            continue


    @write_csv_file
    def write(self):
        return
    
    @check_csv_dialect
    def check(self):
        return
    
    
    def ms_precision(self, min_snr = 10):
        in_file = csv.DictReader(open(self.file_path, 'r'),
                                 dialect = 'excel')
        total = ddict(int)
        weighting_sum = ddict(int)
        for row in in_file:
            diff = float(row['Measured m/z']) - float(row['Calculated m/z'])
            total[row['File name']] += diff * float(row['SNR'])
            weighting_sum[row['File name']] += float(row['SNR'])
        rows = []
        for f in sorted(total):
            mean_diff = total[f] / weighting_sum[f]
            rows.append({'File name': f,
                         'meas.-calc. m/z': mean_diff})
        correction_file = CorrectionFile(self.file_path.replace('.csv', 
                                                                '_precision.csv'),
                                         rows = rows)
        correction_file.write()
    
    
    def file_vs_cos(self,
                    min_snr = None,
                    adducts = None,
#                     charges = None,
                    max_peak_length = None,
                    minint = None,
                    minmaxdp = None,
                    candidates_file = None
                    ):
        """Write csv file with files in rows and oligos in columns
        """
        if adducts:
            adducts = [_ if not _.endswith('-') else '-'+_[:-1] 
                       for _ in adducts]
        in_file = csv.DictReader(open(self.file_path, 'r'),
                                 dialect = 'excel')
        dp_cos = set()
        dp_cos_quant = set()
        rel_int_rows = ddict(dict)
        mol_rows = ddict(dict)
        ng_rows = ddict(dict)
        quantitative = False
        if candidates_file:
            # header with all oligos
            cos_candidates = CandidatesFile(candidates_file)
            cos_candidates.read()
            for cos in cos_candidates.cos:
                _a = cos['A'] if 'A' in cos else 0
                dp_cos.add((cos.dp, _a, cos.name))
                dp_cos_quant.add((cos.dp, _a, cos.name))
        for row in in_file:
            if not row['Intensity [arb. unit]']: # for windows
                continue
            if min_snr and float(row['SNR']) < min_snr:
                continue
            if adducts and row['Adduct'] not in adducts:
                continue
#             if charges and int(row['Charge']) not in charges:
#                 continue
            if max_peak_length and float(row['Peak length [min]']) > max_peak_length:
                continue 
            if minmaxdp and not minmaxdp[0] <= int(row['DP']) <= minmaxdp[1]:
                continue
            _a = int(row['A']) if 'A' in row and row['A'] != '' else 0
            dp_cos.add((int(row['DP']), _a, row['Oligo']))
            rel_intensity = float(row['Intensity [arb. unit]'])
            if row['Oligo'] in rel_int_rows[row['File name']].keys():
                rel_int_rows[row['File name']][row['Oligo']].append(rel_intensity)
            else:
                rel_int_rows[row['File name']][row['Oligo']] = [rel_intensity]
            if row['Amount of substance [nmol]'] != '':
                mol = float(row['Amount of substance [nmol]'])
                ng = float(row['Mass [ng]'])
                dp_cos_quant.add((int(row['DP']), _a, row['Oligo']))
                quantitative = True
                if row['Oligo'] in mol_rows[row['File name']].keys():
                    mol_rows[row['File name']][row['Oligo']].append(mol)
                    ng_rows[row['File name']][row['Oligo']].append(ng)
                else:
                    mol_rows[row['File name']][row['Oligo']] = [mol]
                    ng_rows[row['File name']][row['Oligo']] = [ng]
#             except:
#                 pass
        if quantitative:
            file_vs_cos_fieldnames = ['File name'] + [cos for _dp, _a, cos 
                in sorted(dp_cos_quant)]
            out_file_name = self.file_path.replace(
                '.csv', '_file_vs_cos_snr{0}_nmol.csv'.format(min_snr))
            out_file_name_rel = self.file_path.replace(
                '.csv', '_file_vs_cos_snr{0}_moleFraction.csv'.format(min_snr))
            out_file_name_ng = self.file_path.replace(
                '.csv', '_file_vs_cos_snr{0}_ng.csv'.format(min_snr))
            out_file_ng = csv.DictWriter(open(out_file_name_ng, 'w', newline = ''), 
                                  file_vs_cos_fieldnames)
            out_file_ng.writeheader()
        else:
            file_vs_cos_fieldnames = ['File name'] + [cos for _dp, _a, cos 
                in sorted(dp_cos)]
            out_file_name = self.file_path.replace(
                '.csv', '_file_vs_cos_snr{0}_arbInt.csv'.format(min_snr))
            out_file_name_rel = self.file_path.replace(
                '.csv', '_file_vs_cos_snr{0}_relArbInt.csv'.format(min_snr))
        
        out_file = csv.DictWriter(open(out_file_name, 'w', newline = ''), 
                                  file_vs_cos_fieldnames)
        out_file_rel = csv.DictWriter(open(out_file_name_rel, 'w', newline = ''), 
                                  file_vs_cos_fieldnames)
        out_file.writeheader()
        out_file_rel.writeheader()
        if not quantitative:
            for file_name, outdict in sorted(rel_int_rows.items()):
                total = 0
                for cos, intensities in outdict.items():
                    outdict[cos] = sum(intensities)
                    total += outdict[cos] 
                outdict['File name'] = file_name
                out_file.writerow(outdict)
#             for file_name, outdict in rel_int_rows.items():
#                 total = sum([i for k, i in outdict.items() if k != 'File name'])
                for cos, intensity in outdict.items():
                    if cos != 'File name':
                        outdict[cos] = intensity / total
                out_file_rel.writerow(outdict)
        if quantitative:
            for file_name, outdict in sorted(mol_rows.items()):
                total = 0
                for cos, intensities in outdict.items():
                    outdict[cos] = sum(intensities) / len(intensities)
                    total += outdict[cos] 
                    ng_rows[file_name][cos] = (sum(ng_rows[file_name][cos]) /
                                               len(ng_rows[file_name][cos]))
                outdict['File name'] = file_name
                out_file.writerow(outdict)
                ng_rows[file_name]['File name'] = file_name
                out_file_ng.writerow(ng_rows[file_name])
                for cos, intensity in outdict.items():
                    if cos != 'File name':
                        outdict[cos] = intensity / total
                out_file_rel.writerow(outdict)
    
    
    def RT_plot(self,
                min_snr = 3):
        import numpy as np
        in_file = csv.DictReader(open(self.file_path, 'r'),
                                 dialect = 'excel')
        RTstart_dict, RTend_dict = {}, {}
        RTstart_dict_weighted, RTend_dict_weighted = {}, {}
        for row in in_file:
            try:
                RTstart_dict[int(row['DP'])]
            except:
                RTstart_dict[int(row['DP'])] = ddict(list)
                RTend_dict[int(row['DP'])] = ddict(list)
                RTstart_dict_weighted[int(row['DP'])] = ddict(list)
                RTend_dict_weighted[int(row['DP'])] = ddict(list)
            RTstart_dict[int(row['DP'])
                         ][row['Oligo']
                           ].append(float(row['Peak start time [min]']))
            RTend_dict[int(row['DP'])
                       ][row['Oligo']
                         ].append(float(row['Peak end time [min]']))
            snr = float(row['SNR'])
            for _ in range(int(snr)):
                RTstart_dict_weighted[int(row['DP'])
                             ][row['Oligo']
                               ].append(float(row['Peak start time [min]']))
                RTend_dict_weighted[int(row['DP'])
                           ][row['Oligo']
                             ].append(float(row['Peak end time [min]']))
        for dp in sorted(RTstart_dict):
#             print(dp)
            for oligo, RTs in sorted(RTstart_dict[dp].items()):
                start_median = np.median(RTs)
                end_median = np.median(RTend_dict[dp][oligo])
                start = np.average([_ for _ in RTstart_dict_weighted[dp][oligo] if start_median-2 < _ < start_median+2])
                end = np.average([_ for _ in RTend_dict_weighted[dp][oligo] if end_median-2 < _ < end_median+2])
                if np.isnan(start):
                    continue
                print('{1:5.2f}\t{2:5.2f}\t{0}'.format(oligo,
                                                         start,
                                                         end))


class FragmentFile():
    def __init__(self, file_path):
        self.file_path = file_path
        self.fieldnames = ['File name', 
                           'Precursor oligo', 
                           'Precursor ion type', 
                           'Precursor m/z',
                           'Fragment oligo',
                           'Fragment ion type', #e.g. y2, b4++
                           'Fragment intensity [arb. unit]', 
                           'Relative intensity (w.r.t. fragment type)',
                           'Number of spectra',
                           'Peak start time [min]',
                           'Peak end time [min]',
                           'Peak length [min]',
                           'Measured fragment m/z (average)',
                           'Calculated fragment m/z',
                           'Adduct',
                           'Charge',
                           'DP',
#                            'A', 
#                            'D',
#                            'R',
                           'Potential precursor sequences']
        self.additional_fieldnames = set()
        self.rows = []
    
    
    def add(self, cos, entry):
        entry['Fragment oligo'] = cos.name
        entry['Fragment ion type'] = '{0}{1}'.format(cos.fragment_type, cos.dp)
        entry['Fragment ion type'] += cos.charge * '+' if cos.charge > 1 else ''
        entry['Calculated fragment m/z'] = round(cos.mz, 3)
        entry['Adduct'] = cos.adduct
        entry['Charge'] = cos.charge
        entry['DP'] = cos.dp
#         entry['Potential precursor sequences'] = ', '.join(sorted(cos.precursor_sequences))
        entry['Potential precursor sequences'] = ', '.join([_ for _ in 
                                                            [''.join(__) for __ in 
                                                             sorted(cos.precursor_sequences)]])
        for monomer, count in cos.items():
            entry[monomer] = count
        self.rows.append(entry)
        for k in cos.keys():
            self.additional_fieldnames.add(k)


    def ms_precision(self, min_snr = 10):
        in_file = csv.DictReader(open(self.file_path, 'r'),
                                 dialect = 'excel')
        total = ddict(int)
        weighting_sum = ddict(int)
        for row in in_file:
            diff = (float(row['Measured fragment m/z (average)']) - 
                    float(row['Calculated fragment m/z']))
            total[row['File name']] += diff * float(row['Fragment intensity [arb. unit]'])
            weighting_sum[row['File name']] += float(row['Fragment intensity [arb. unit]'])
        rows = []
        for f in sorted(total):
            mean_diff = total[f] / weighting_sum[f]
            rows.append({'File name': f,
                         'meas.-calc. m/z': mean_diff})
        correction_file = CorrectionFile(self.file_path.replace('.csv', 
                                                                '_precision.csv'),
                                         rows = rows)
        correction_file.write()

    
    def sequencing(self):
        in_file = csv.DictReader(open(self.file_path, 'r'),
                                 dialect = 'excel')
        seq_file = SequencesFile(self.file_path.replace('.csv', '_seq.csv'))
        potential_precursor_sequences = {}
        prec_fragments = ddict(dict)
        current_file_cos = (False, False)
        for row in in_file:
#             print(row)
            if current_file_cos != (False, False):
                if current_file_cos != (row['File name'], row['Precursor oligo']): 
#                     prec_compo = chem.monomeric_composition(row['Precursor oligo'])
#                     prec_ion = chem.Ion(prec_compo,
#     #                                     charge = int(row['Charge']), 
#                                         adduct = row['Adduct'])
#                     print(prec_ion)
#                     print(prec_fragments) 
#                     print(potential_precursor_sequences)
                    seq_results = seq.determine_sequences(
                        precursor = prec_ion,
                        fragments = prec_fragments,
                        prec_seq = potential_precursor_sequences)
                    
                    if seq_results:
                        for seqrow in seq_results:
                            seqrow['File name'] = row['File name']
                            seqrow['Number of spectra'] = row['Number of spectra']
                            seqrow['Peak start time [min]'] = row['Peak start time [min]']
                            seqrow['Peak end time [min]'] = row['Peak end time [min]']
                            seqrow['Peak length [min]'] = row['Peak length [min]']
                            seq_file.add(prec_ion, seqrow)
                    
                    potential_precursor_sequences = {}
                    prec_fragments = ddict(dict)
            mono_compo = chem.monomeric_composition(row['Fragment oligo'])
            fragment_ion = chem.Ion(
                mono_compo,
#                                     charge = int(row['Charge']), 
                adduct = 'H',#row['Adduct'], 
#                                     redend_label = row[''],
#                                     modification = None,
#                                     n_mer = 1,
                fragment_type = row['Fragment ion type'][0],
                precursor_sequences = set(row['Potential precursor sequences'].split(', ')))
            prec_fragments[fragment_ion.fragment_type, 
                           fragment_ion.dp][fragment_ion.key] = (float(row['Fragment intensity [arb. unit]']), 
                                                                 float(row['Relative intensity (w.r.t. fragment type)']))
            potential_precursor_sequences[fragment_ion.key] = fragment_ion.precursor_sequences
            current_file_cos = (row['File name'], row['Precursor oligo'])
            prec_compo = chem.monomeric_composition(row['Precursor oligo'])
            prec_ion = chem.Ion(
                prec_compo,
#                                     charge = int(row['Charge']), 
                adduct = row['Adduct'],
                nonredend_types = ['b'],
                 redend_types = ['y'],)
        prec_compo = chem.monomeric_composition(row['Precursor oligo'])
        prec_ion = chem.Ion(prec_compo,
#                             charge = int(row['Charge']), 
                            adduct = row['Adduct'],
                            nonredend_types = ['b'],
                 redend_types = ['y'],)
        seq_results = seq.determine_sequences(precursor = prec_ion,
                                    fragments = prec_fragments,
                                    prec_seq = potential_precursor_sequences)
    
        if seq_results:
            for seqrow in seq_results:
                seqrow['File name'] = row['File name']
                seqrow['Number of spectra'] = row['Number of spectra']
                seqrow['Peak start time [min]'] = row['Peak start time [min]']
                seqrow['Peak end time [min]'] = row['Peak end time [min]']
                seqrow['Peak length [min]'] = row['Peak length [min]']
                seq_file.add(prec_ion, seqrow)
        seq_file.write()


    @write_csv_file
    def write(self):
        return


class SequencesFile():
    def __init__(self, file_path):
        self.file_path = file_path
        self.fieldnames = ['File name', 
                           'Precursor oligo', 
                           'Precursor ion type', 
                           'Precursor m/z',
                           'Sequence',
                           'Relative intensity (w.r.t. oligo)',
                           'Number of spectra',
                           'Peak start time [min]',
                           'Peak end time [min]',
                           'Peak length [min]',
                           'Notes',
                           'Adduct',
                           'Charge',
                           'DP',
                           'All possible sequences'
#                            'A', 
#                            'D',
#                            'R',
                           ]
        self.additional_fieldnames = set()
        self.rows = []
        
        
    def add(self, cos, entry):
        entry['Precursor oligo'] = cos.name
        entry['Precursor ion type'] = cos.precursor_type
        entry['Precursor m/z'] = round(cos.mz, 3)
        entry['Adduct'] = cos.adduct
        entry['Charge'] = cos.charge
        entry['DP'] = cos.dp
        for monomer, count in cos.items():
            entry[monomer] = count
        self.rows.append(entry)
        for k in cos.keys():
            self.additional_fieldnames.add(k)
        

    @write_csv_file
    def write(self):
        return
    

    #todo pa at red nonred ends
class QuantSeqFile():
    def __init__(self, 
                 ms1filevscos, 
                 sequences_file,
                 file_overview):
        self.ms1filevscos = ms1filevscos
        self.sequences_file = sequences_file
        self.file_overview = file_overview
        file_path = sequences_file.replace('.csv', '_incl_' +
                                           os.path.split(ms1filevscos)[1])
        self.file_path = file_path
        self.fieldnames = ['File name MS1', 
                           'File name MS2',
                           'Oligo', 
                           'Sequence']
        if ms1filevscos.endswith('nmol.csv'):
            self.amount_key = 'Amount of substance [mol]'
        else:
            self.amount_key = 'Mole fraction'
        self.fieldnames += [self.amount_key, 'All possible sequences']
        
        self.rows = []
    
        seq_file = csv.DictReader(open(self.sequences_file, 'r'),
                                  dialect = 'excel')
        sequenced = set()
        for row in seq_file:
            if row['File name'] not in self.ms2_to_ms1.keys():
                continue
            outdict = {'File name MS2': row['File name'],
                       'Oligo': row['Precursor oligo'], 
                       'Sequence': row['Sequence'],
                       'All possible sequences': row['All possible sequences']}
            try:
                outdict[self.amount_key] = float(
                    row['Relative intensity (w.r.t. oligo)']
                    ) * self.quant[self.ms2_to_ms1[row['File name']]
                                   ][row['Precursor oligo']]
                outdict['File name MS1'] = self.ms2_to_ms1[row['File name']]
                sequenced.add((self.ms2_to_ms1[row['File name']], 
                               row['Precursor oligo']))
            except: # no quantitative ms1 results
                pass
            self.rows.append(outdict)
        for ms1_file_name, cos_dict in self.quant.items():
            for cos, amount in cos_dict.items():
                if (ms1_file_name, cos) in sequenced:
                    continue
                outdict = {'File name MS1': ms1_file_name,
                           'Oligo': cos,
                           self.amount_key: amount}
                mono_compo = chem.monomeric_composition(cos)
                if len(mono_compo) == 1:
                    sequence = ''
                    for m, count in mono_compo.items():
                        sequence += m * count
                    outdict['Sequence'] = sequence
                    outdict['All possible sequences'] = sequence # todo check if this is required
                else:
                    _cos = chem.Oligo(mono_compo)
                    outdict['Sequence'] = _cos.dp * 'x'
                    outdict['All possible sequences'] = ', '.join(list(sorted(
                        [''.join(_) for _ in _cos.sequences])))
                self.rows.append(outdict)
        self.write()
        
        
    # count PAs todo
    def pa(self, monomers = None, pa_length = None):
        if pa_length == 1:
            patterns = set(monomers)
            idxs = ((-3,), (-2,), (-1,), (0,), (1,), (2,))
        elif pa_length == 2:
            patterns = set()
            for mono1 in monomers:
                for mono2 in monomers:
                    patterns.add(mono1 + mono2)
            idxs = ((-2, -1), (0, 1))
        elif pa_length == 3:
            patterns = set()
            for mono1 in monomers:
                for mono2 in monomers:
                    for mono3 in monomers:
                        patterns.add(mono1 + mono2 + mono3)
            idxs = ((-3, -2, -1), (0, 1, 2))
        
        file_idx_sum = {}
        filecos_idx_xpos_sum = {}
        filecos_idx_xpos_seq = {}
        all_sequences = ddict(set)
        for row in self.rows:
            if 'File name MS1' not in row: # not quantified
                continue
            if row['File name MS1'] not in file_idx_sum: # create entry
                file_idx_sum[row['File name MS1']] = {idx: {p: 0 
                                                            for p in patterns | 
                                                            set(['unknown'])
                                                            }
                                                      for idx in idxs}
            if (row['File name MS1'], row['Oligo']) not in filecos_idx_xpos_sum:
                filecos_idx_xpos_sum[row['File name MS1'], row['Oligo']] = {
                    idx: {p: ddict(int) 
                          for p in patterns | set(['unknown'])
                          } for idx in idxs}
                filecos_idx_xpos_seq[row['File name MS1'], row['Oligo']] = {
                    idx: {p: ddict(set) 
                          for p in patterns | set(['unknown'])
                          } for idx in idxs}
#             if 'Sequence' not in row: # unknown # TODO check if this commenting is problem
#                 for idx in idxs:
#                     file_idx_sum[row['File name MS1']][idx]['unknown'
#                         ] += row[self.amount_key]
#                 continue
            for idx in idxs:
                p_to_check = ''
                for i in idx:
                    try:
                        p_to_check += row['Sequence'][i]
                    except IndexError: 
                        p_to_check = False
                        break
                if not p_to_check:
#                     if idx == (0,):    
#                         print('\t', p_to_check, '\t',idx, '\t', row['Sequence'],'\t',row[self.amount_key])
                    file_idx_sum[row['File name MS1']][idx]['unknown'
                        ] += row[self.amount_key]
                else:
                    if p_to_check not in patterns:
                        p_to_check = 'unknown'
                    xpos = tuple([pos for pos, m in enumerate(row['Sequence']) 
                                  if m == 'x'])
                    seqs = set(row['All possible sequences'].split(', '))
                    if xpos:
#                         if idx == (0,):    
#                             print('\t', '!',p_to_check, '\t',idx, '\t', row['Sequence'],'\t',row[self.amount_key])
                        filecos_idx_xpos_sum[row['File name MS1'], row['Oligo']][
                            idx][p_to_check][xpos] += row[self.amount_key]
                        filecos_idx_xpos_seq[row['File name MS1'], row['Oligo']][
                            idx][p_to_check][xpos] |= seqs
                    else:
                        all_sequences[row['File name MS1'], 
                                      idx] |= seqs
                        file_idx_sum[row['File name MS1']][idx][p_to_check
                            ] += row[self.amount_key]
#                         if idx == (0,):    
#                             print('\t', p_to_check, '\t',idx, '\t', row['Sequence'],'\t',row[self.amount_key])
        for (file_name, cos), idx_dict in filecos_idx_xpos_seq.items():
            for idx, p_dict in idx_dict.items():
                for p, xpos_dict in sorted(p_dict.items()):
                    if not xpos_dict:
                        continue
#                     xpos_order = []
                    for xpos1, seqs1 in xpos_dict.items():
#                         if not seqs1 <= all_sequences[file_name, 
#                                                       idx]:
                        if not len(seqs1 & all_sequences[file_name, idx]) > 0:
#                             if idx == (0,):    
#                                 print(seqs1, all_sequences[file_name, 
#                                                       idx])
                            amount = filecos_idx_xpos_sum[
                                file_name, cos][idx][p][xpos1]
#                             if idx == (0,):    
#                                 print('\t', cos, '\t',idx, '\t',p, '\t',xpos1, '\t',amount)
                            file_idx_sum[file_name][idx][p] += amount
                            all_sequences[file_name, idx] |= seqs1
#                     for xpos1, seqs1 in xpos_dict.items():
#                         if not seqs1 <= all_sequences[file_name, 
#                                                       idx]:
#                         if len(seqs1 - all_sequences[file_name, idx]) > 0:
#                             if idx == (0,):    
#                                 print('''!!!''')
#                                 print(seqs1 - all_sequences[file_name, idx])
#                                 amount = filecos_idx_xpos_sum[
#                                                               file_name, cos][idx][p][xpos1]
#                                 print('\t', cos, '\t',idx, '\t',p, '\t',xpos1, '\t',amount)
                            
        for file_name, idx_dict in file_idx_sum.items():
            print()
            print(file_name)
            for idx, pa_dict in sorted(idx_dict.items()): 
                diff = (sum(self.quant[file_name].values()) -
                        sum(pa_dict.values()))
                if abs(diff) > 0.01:
                    print('WARNING - sth. wrong with PA.', file_name, idx, diff)
                for p, amount in sorted(pa_dict.items()):
                    print(idx, p, amount)
        
        out_file_name = self.file_path.replace('.csv', 
                                               '_PA{0}.csv'.format(pa_length))
        out_file = PaFile(out_file_name,
                          file_idx_sum)
        out_file.write()
        
    
    @misc.LazyFunction
    def quant(self):
        ms1file = csv.DictReader(open(self.ms1filevscos, 'r'),
                                 dialect = 'excel')
        _quant = {}
        for row in ms1file:
            if row['File name'] not in self.ms2_to_ms1.values():
                continue
            _quant[row['File name']] = {}
            for cos, amount in row.items():
                if cos == 'File name' or amount == '':
                    continue
                _quant[row['File name']][cos] = float(amount)
        return _quant
                
    
    @misc.LazyFunction
    def ms2_to_ms1(self):
        fi_ov = csv.DictReader(open(self.file_overview, 'r'),
                                 dialect = 'excel')
        _ms2_to_ms1 = {}
        for row in fi_ov:
            _ms2_to_ms1[row['MS2']] = row['MS1']
        return _ms2_to_ms1    
    
        
    @write_csv_file
    def write(self):
        return
    

class PaFile():
    def __init__(self, file_name, file_idx_sum):
        self.file_path = file_name
        tmp_fieldnames = set()
        tmp2_fieldnames = set()
        out_dicts = ddict(dict)
        for f, idx_dict in file_idx_sum.items():
            for idx, p_dict in idx_dict.items():
                idx_sum = sum([amount for p, amount in p_dict.items()
                               if p != 'unknown'])
                for p, amount in p_dict.items():
                    idxs = ''.join(['{0:+d}'.format(i if i < 0 else i+1)
                                    for i in idx])
                    key = '{0} {1}'.format(idxs, p)
                    tmp_fieldnames.add((idx, key))
                    out_dicts[f][key] = amount
                    if p != 'unknown':
                        key_norm = key + ' norm.'
                        tmp2_fieldnames.add((idx, key_norm))
                        if idx_sum != 0:
                            out_dicts[f][key_norm] = amount / idx_sum
        self.fieldnames = ['File name']
        self.fieldnames += [key for idx, key in sorted(list(tmp_fieldnames))]
        self.fieldnames += [key for idx, key in sorted(list(tmp2_fieldnames))]
        self.rows = []
        for f, out_dict in sorted(out_dicts.items()):
            out_dict['File name'] = f
            self.rows.append(out_dict) 
        
        
    @write_csv_file
    def write(self):
        return


class CorrectionFile():
    def __init__(self, file_name, rows = None):
        self.file_path = file_name
        self.fieldnames = ['File name',
                           'meas.-calc. m/z']
        if rows:
            self.rows = rows
        else:
            self.rows = []
            if os.path.exists(self.file_path):
                in_file = csv.DictReader(open(self.file_path, 'r'),
                                         dialect = self.check())
                for row in in_file:
                    self.rows.append({'File name': row['File name'],
                                      'meas.-calc. m/z': float(row['meas.-calc. m/z'])})
    
    @misc.LazyFunction
    def corr(self):
        f_to_corr = {}
        for row in self.rows:
            f_to_corr[row['File name']] = float(row['meas.-calc. m/z'])
        return f_to_corr
        
    
    @check_csv_dialect
    def check(self):
        return
    
    
    @write_csv_file
    def write(self):
        return


class MS2ScreenFile:
    def __init__(self, file_name):
        self.file_path = file_name
        self.rows = []
        self.fieldnames = ['File name']
        self.additional_fieldnames = set()
        self.data = ddict(dict)
        return
    
    
    def add(self, file_name, key, value):
        self.data[file_name][key] = value
#         entry['Precursor oligo'] = cos.name
#         entry['Precursor ion type'] = cos.precursor_type
#         entry['Precursor m/z'] = round(cos.mz, 3)
#         entry['Adduct'] = cos.adduct
#         entry['Charge'] = cos.charge
#         entry['DP'] = cos.dp
#         for monomer, count in cos.items():
#             entry[monomer] = count
#         self.rows.append(entry)
        self.additional_fieldnames.add(key)
    
    
    def write(self):
        for f, d in sorted(self.data.items()):
            d['File name'] = f
            self.rows.append(d)
        @write_csv_file
        def _write(self):
            return
        _write(self)
        # normalize
        self.file_path = self.file_path.replace('.csv', '_norm.csv') 
        for idx, row in enumerate(self.rows):
            kd = ddict(list)
            kdi = {}
            for k in self.fieldnames[1:]:
                kl = k.split('_')
                try:
                    _ = (kl[-3], kl[-1])
                except:
                    _ = (kl[-1],)
                kd[_].append(k)
                kdi[k] = (_)
            new_row = {'File name': row['File name']}
            for k, v in row.items():
                if k == 'File name':
                    continue
                new_row[k] = v / sum([row[_] for _ in kd[kdi[k]] if _ in row])
            self.rows[idx] = new_row 
        self.additional_fieldnames = set()
        _write(self)
        return


class TestSamplesFile(unittest.TestCase):
    def setUp(self):
        self.f = CandidatesFile('../tests/TestSamplesFile.csv')
        
    def test_write(self):
        self.f.create_all_possible(
            monomers = ['A', 'D', 'Hex'],
            min_dp = 1,
            max_dp = 10,
            adducts = ['H','Na','K','NH4',
                       '2H', '2Na', '2K',
                       '3H', '3Na', '3K'],
            molecules = [1, 2])
        
        self.f.write()
        self.f.clear()
        self.f.read()
        self.f.check()


class TestSamplesFilePectin(unittest.TestCase):
    def setUp(self):
        self.f = CandidatesFile('../tests/TBDMS.csv')
        
    def test_write(self):
        self.f.create_all_possible(monomers = ['D','Dtbdmsone', 'Dtbdmstwo','Dtbdmsthree'],
                                   min_dp = 1,
                                   max_dp = 6,
                                   adducts = ['H', '2H', 'Na', '2Na', 'NH4', 'K', 'H+NH4', 'ACN+H'],
                                   molecules = [1],
                                    nonredend_labels = [None],
                                    redend_labels = [None],
                                    nonredend_types = ['b', 'c'],
                                    redend_types = ['y', 'z']
                        )
        
        self.f.write()
        self.f.clear()
        self.f.read()
        self.f.check()
        
        
class TestMassFile(unittest.TestCase):
    def setUp(self):
        self.cand = '../tests/TBDMS.csv'
        self.m = MassFile('../tests/TBDMS_Masses.csv')
    def test_mass_list(self):
        self.m.add_candidates_file(self.cand)
        self.m.write()
        
        
# class TestRTs(unittest.TestCase):
#     def setUp(self):
#         self.filename = '/mnt/m/Group/fachbereiche/biologie/r0bpmo/user/Anna N/MS/MS Data others/Margareta/17-08-01_results2/COS_GM_DP3-10_H+_charge-1-2_MS1-results.csv'
#     def test_mass_list(self):
#         ms1 = MS1File(self.filename)
#         ms1.RT_plot(20)
