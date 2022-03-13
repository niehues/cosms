#!/usr/bin/env python3
# coding: utf-8
"""mrm
Analysis of LC-MS runs with MRM
"""
import csv
import os
import sys
import pymzml
import math
import bisect
import argparse
from collections import defaultdict as ddict
from core import misc

parser = argparse.ArgumentParser(
description = """mrm_noise.py script for analysis of LC-MS MRM runs""",
formatter_class = argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('mzml',
                type = str,
                help = """mzML file OR folder containing multiple mzML files""")
                
parser.add_argument('candidates',
                type = str,
                help = """CSV file with four columns: 
                            "Precursor m/z", "Fragment m/z", "Min RT", and "Max RT"  
                            that contain m/z values of interest 
                            ("Min RT" and "Max RT" are optional)""")
parser.add_argument('-p', '--precursor',
                nargs = '?',
                type = float,
                default = None,
                help = """precursor ion m/z error""")
        
parser.add_argument('-f', '--fragment',
                nargs = '?',
                type = float,
                default = None,
                help = """fragment ion m/z error""")

parser.add_argument('-hp', '--highest_peaks',
                nargs = '?',
                default = None,
                type = int,
                choices = [i for i in range(0,20)],
                help = """Number of peaks in which the fragment mass need to be found""")

def write_csv_file(cls):
    def _write_csv_file(cls):
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


class MzmlFile():
    def __init__(self, 
                 file_path, candidates_data,
                 precursor_precision = None,
                 fragment_precision = None,
                 highest_peaks = None):
        """Set up pymzml run
        """
        self.file_path = file_path
        self.precursor_precision = (precursor_precision 
                                    if precursor_precision != None else 0.3)
        self.fragment_precision = (fragment_precision
                                   if fragment_precision != None else 0.3)
        self.highest_peaks = (highest_peaks if highest_peaks != None else 15)
        self.high_precision = 5e-6
        self.precursor_noise_level = {}
        self.precursor_temp_dict = {}
        self.file_name = os.path.split(file_path)[1].replace('.mzML', '')
        print('Analyzing "{0}"'.format(self.file_name))
        self.noise_levels(candidates_data)
        return
        
        
    def run(self):
        """Set up pymzml.run.Reader
        """
        _run = pymzml.run.Reader(self.file_path, 
                                 MS1_Precision = self.high_precision,
                                 MSn_Precision = self.high_precision)
        return _run
    
    
    def mrm(self,
            candidates_data):
        p_mzs = {(name, mz): 
                 (mz - self.precursor_precision,
                  mz + self.precursor_precision)
                for name, mz in candidates_data.candidates}
        f_mzs = {(name, mz): 
                 {f: (f - self.fragment_precision,
                      f + self.fragment_precision)
                  for f in fs}
                 for (name, mz), fs in candidates_data.candidates.items()}
        fragment_ions = {(name, mz): 
                         {f: [] for f in fs}
                         for (name, mz), fs 
                         in candidates_data.candidates.items()}
        for spectrum in self.run():
            if spectrum.ms_level == 2:
                precursor_mz = float(spectrum['MS:1000744'])
                tic = float(spectrum['MS:1000285'])
                precursor_intensity = None
                for element in spectrum.element.getiterator():
                    if element.get('accession', default = None
                                   ) == 'MS:1001141':
#                 for element in spectrum.xmlTree:
#                     if element.get('accession') == 'MS:1001141':
                        precursor_intensity = float(element.get('value'))
                        break
#                 if not precursor_intensity: # in case of MRM
#                     continue
                precursors_in_spec = []
                for (name, mz), (min_mz, max_mz) in p_mzs.items():
                    if min_mz < precursor_mz < max_mz:
                            precursors_in_spec.append((name, mz))
                if not precursors_in_spec:
                    continue
                try:
                    highest_peaks = spectrum.highest_peaks(self.highest_peaks)
                except:
                    continue
                for name, mz in precursors_in_spec:
                    for fmz, (min_fmz, max_fmz) in f_mzs[name, mz].items():
                        if (candidates_data.RTs[name, mz, fmz]['min'] and 
                            (candidates_data.RTs[name, mz, fmz]['min'] > 
                            spectrum['scan start time'])):
                            continue
                        if (candidates_data.RTs[name, mz, fmz]['max'] and 
                            (candidates_data.RTs[name, mz, fmz]['max'] < 
                            spectrum['scan start time'])):
                            continue  
                        for _mz, _i in highest_peaks:
                            if min_fmz < _mz < max_fmz:
                                try:
                                    noise_level = self.precursor_noise_level[mz]
                                except:
                                    noise_level = 0
                                fragment_ions[name, mz][fmz].append(
                                    (spectrum['scan start time'],
                                     precursor_intensity,
                                     _i,  
                                     tic,
                                     noise_level))
        return fragment_ions
        
    def noise_levels(self, candidates_data):
        """Determine retention time dependent noise levels in MS2 spectra 
        """
        intensities_for_precursor = {}
        precursor_dict = {}
        
        p_mzs = {(name, mz): 
                 (mz - self.precursor_precision,
                  mz + self.precursor_precision)
                for name, mz in candidates_data.candidates}
        
        for spectrum in self.run():
            if spectrum.ms_level == 2:
                precursor_mz = float(spectrum['MS:1000744'])
                start_scan_time = float(spectrum['MS:1000016'])
        #Creates a dictionary with the different precursor masses as keys and the corresponding retention times and intensities of each spectrum as values
                for (name, mz), (min_mz, max_mz) in p_mzs.items():
                    if min_mz < precursor_mz < max_mz:
                        if mz not in precursor_dict.keys():
                            precursor_dict[mz] = []
                            precursor_dict[mz].append((start_scan_time, list(spectrum.i)))
                        else:
                            precursor_dict[mz].append((start_scan_time, list(spectrum.i)))
                        
        for precursor, RT_intensity in precursor_dict.items():
                ##TODO: faster
            for i, j in RT_intensity:
                try:
                    intensities_for_precursor[precursor] += j
                except:
                    intensities_for_precursor[precursor] = j
            try:
                intensities_for_precursor[precursor][3]
            except IndexError:
                continue
            intensities_for_precursor[precursor].sort()
            #plt.hist(intensities_for_precursor[precursor], bins=2000)
            #plt.show()
            mad1 = misc.median_absolute_deviation(intensities_for_precursor[precursor], 
                                                    sorting = False)
            std_dev_estimate1 = misc.std_dev_from_mad(mad1)
            intensity_threshold = 3*std_dev_estimate1
            idx = bisect.bisect(intensities_for_precursor[precursor], intensity_threshold)
            below_threshold = intensities_for_precursor[precursor][:idx]
            mad2 = misc.median_absolute_deviation(below_threshold, 
                                                    sorting = False)
            std_dev_estimate2 = misc.std_dev_from_mad(mad2)
            self.precursor_noise_level[precursor] = std_dev_estimate2
        return
    
    
class CandidatesFile():
    def __init__(self, file_path):
        """ CSV file that contains the following columns:
        Name    
        Precursor m/z    
        Fragment m/z    
        Min RT    
        Max RT
        """
        self.file_path = file_path
        self.in_file = csv.DictReader(open(file_path, 'r'),
                                      dialect = self.check())
        self.candidates = ddict(list)
        self.RTs = ddict(dict)
        self.precision = {}
        for row in self.in_file:
            pmz = float(row['Precursor m/z'])
            fmz = float(row['Fragment m/z'])
            try:
                name = row['Name']
            except KeyError: # no column "Name" in csv file 
                name = pmz
#             try: 
            self.candidates[name, pmz].append(fmz)
#             except:
#                 self.candidates[name, pmz] = [fmz]
            try:
                min_rt = float(row['Min RT'])
            except:
                min_rt = None
            try:
                max_rt = float(row['Max RT'])
            except:
                max_rt = None
            self.RTs[name, pmz, fmz]['min'] = min_rt
            self.RTs[name, pmz, fmz]['max'] = max_rt
        return
    
    
    @check_csv_dialect
    def check(self):
        return
        

class MRMResultsFile():
    def __init__(self, file_path):
        """ CSV file with results
        """
        self.file_path = file_path
        self.fieldnames = ['File name',
                           'Name',
                           'Precursor m/z', 
                           'Fragment m/z',
                           'Precursor intensity',
                           'Fragment intensity',
                           'TIC', 
                           'Min RT', 
                           'Max RT',
                           'S2N-ratio']
        self.rows = []
        return
    
    
    def add(self, file_name, fragment_ions):
        """ Add row to self.rows
        """
        for (name, pmz), fmz_d in fragment_ions.items():
            for fmz, l in fmz_d.items():
                if not l:
                    continue
                precursor_i, fragment_i, tic, noise = 0, 0, 0, 0
                for _rt, _pi, _fi, _tic, _nl in l:
                    if _pi != None:
                        precursor_i += _pi 
                    fragment_i += _fi
                    tic += _tic
                    noise += _nl
                if noise != 0:
                    s2n = round(fragment_i/noise,2)
                elif noise == 0:
                    s2n = 'could not be calculated'
                entry = {'File name': file_name,
                         'Name': name,
                         'Precursor m/z': pmz, 
                         'Fragment m/z': fmz,
                         'Precursor intensity': '{0:.0f}'.format(precursor_i),
                         'Fragment intensity': '{0:.0f}'.format(fragment_i),
                         'TIC': '{0:.0f}'.format(tic), 
                         'Min RT': l[0][0], 
                         'Max RT': l[-1][0],
                         'S2N-ratio': s2n}
                self.rows.append(entry)
        return
                
                
    def sort_rows(self):
        self.rows.sort(key = lambda x: (x['File name'], x['Name'], 
                                        x['Precursor m/z'], x['Fragment m/z']))
        return
        
        
    @write_csv_file
    def write(self):
        return        


def file_list(path, file_type = 'mzML'):
    """ Returns mzML files found in given directory
    """
    files = []
    if os.path.isdir(path):
        for dirpath, _, filenames in os.walk(path, 
                                                    topdown = True, 
                                                    followlinks = True):
                for f in sorted(filenames):
                    if f.endswith(file_type):
                        file_path = os.path.join(path, dirpath, f)
                        files.append(file_path)
    else:
        files.append(path)
    return files


if __name__ == '__main__':
    args = parser.parse_args()

    mzml_file = args.mzml
    candidates_file = args.candidates
    candidates_data = CandidatesFile(candidates_file)
    precursor_precision = args.precursor
    fragment_precision = args.fragment
    highest_peaks = args.highest_peaks
    
    output_dir, tail = os.path.split(mzml_file)
    if not os.path.exists(output_dir):
        os.mkdir(output_dir) # make directory
    out_file = MRMResultsFile(os.path.join(
        output_dir, '{0}_{1}{2}{3}{4}_results.csv'.format(
            tail.replace('.mzML', ''),
            os.path.split(candidates_file)[1].replace('.csv',''),
            '_{0}'.format(precursor_precision) if precursor_precision else '',
            '_{0}'.format(fragment_precision) if fragment_precision else '',
            '_{0}'.format(highest_peaks) if highest_peaks else '')))
    for f in file_list(mzml_file):     
        infile = MzmlFile(f, candidates_data = candidates_data,
                          precursor_precision = precursor_precision,
                          fragment_precision = fragment_precision,
                          highest_peaks = highest_peaks)
        fragment_ions = infile.mrm(candidates_data)
        out_file.add(infile.file_name, fragment_ions)
    out_file.sort_rows()
    out_file.write() # write csv output file
    
    
