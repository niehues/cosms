#!/usr/bin/env python3
# coding: utf-8
"""mzml
analyze ms1 and ms2 data; requires pymzml [github.com/pymzml/pymzML]
"""
import bisect
import os
import copy
import numpy as np
from collections import defaultdict as ddict
from cosms.core import chem
from cosms.core import misc
from cosms.core import inout
try:
    from cosms.core import plot
except:
    print("Could not import module: plot")
from cosms.core import seq
import pymzml # github.com/pymzml/pymzML


class Chromatogram(): 
    def __init__(self, x_values = None, intensities = None):
        self.x_values = list(x_values) if x_values != None else [] 
        self.intensities = list(intensities) if intensities != None else []

    
    def add(self, x, i):
        self.x_values.append(x)
        self.intensities.append(i)
        
    def sort(self):
        tmp_x = []
        tmp_y = []
        for x,y in sorted(zip(self.x_values, self.intensities)):
            tmp_x.append(x)
            tmp_y.append(y)
        self.x_values = tmp_x
        self.intensities = tmp_y
    
    
    @property
    def xy(self):
        return list(zip(self.x_values, self.intensities))
    
    
    @property
    def length(self):
        return len(self.x_values)


    def smooth(self, x_window = 0.2, grid_size = 0.02):
#     def smooth(self, x_window = 0.1, grid_size = 0.02):
        """Smoothing for peak detection (not peak measuring) by triangular 
        moving average
        """
        self._i_smooth = []
        factor = 1 / grid_size
        self._x_grid = [x / factor 
                        for x in range(int(round(self.x_values[0] * factor)),
                                       int(round(self.x_values[-1] * factor)))]
        for x in self._x_grid:
            min_idx = bisect.bisect_left(self.x_values, x - x_window)
            max_idx = bisect.bisect_right(self.x_values, x + x_window)
            try:
                max_i = max(self.intensities[min_idx:max_idx])
            except:
                self._i_smooth.append(0)
                continue
            window_trifactor_i = [((x_window - abs(x - m)) / x_window, i) 
                for m, i in zip(self.x_values[min_idx:max_idx],
                                self.intensities[min_idx:max_idx]
                                ) if i < max_i]
            self._i_smooth.append(sum([trifactor * i 
                                 for trifactor, i in window_trifactor_i]))
        self.smoothlen = len(self._i_smooth)
        return 
    
    
    @misc.LazyFunction
    def i_smooth(self):
        try:
            self._i_smooth
        except:
            self.smooth()
        return self._i_smooth


    @misc.LazyFunction
    def x_grid(self):
        try:
            self._x_grid
        except:
            self.smooth()
        return self._x_grid


    def elution_peak(self, min_anomer_height = 0.2):
        """Find elution peak in extracted precursor_type chromatogram that was 
        smoothed by triangular moving average
        """
        min_rt, max_rt = None, None
        try:
            max_i = max(self.i_smooth)
        except:
            return min_rt, max_rt
        max_i_idx = self.i_smooth.index(max_i)
        if 0 < max_i_idx < len(self.i_smooth):
            half_max_i = max_i / 2
            for right_idx, i in enumerate(self.i_smooth[max_i_idx:]):
                if i < half_max_i:
                    break
            for left_idx, i in enumerate(reversed(self.i_smooth[:max_i_idx])):
                if i < half_max_i:
                    break
            right_idx += max_i_idx
            left_idx = max_i_idx - left_idx
            min_rt, max_rt = self.x_grid[left_idx], self.x_grid[right_idx]
            rt_length = max_rt - min_rt
            half_rt_length = rt_length / 2
            min_rt -= half_rt_length
            max_rt += half_rt_length
            # test for anomer peak
            min_idx = bisect.bisect_left(self.x_grid, max_rt)
            max_idx = bisect.bisect_right(self.x_grid, max_rt + rt_length) 
            found_anomeric_peak = True
            try:
                max_i_ano = max(self.i_smooth[min_idx:max_idx])
            except:
                found_anomeric_peak = False
            if found_anomeric_peak:
                max_i_idx_ano = min_idx + self.i_smooth[
                    min_idx:max_idx].index(max_i_ano)
                if (max_i_ano > min_anomer_height * max_i and 
                    min_idx < max_i_idx_ano < max_idx):
                    half_max_i_ano = max_i_ano / 2
                    for right_idx_ano, i in enumerate(
                        self.i_smooth[max_i_idx_ano:]):
                        if i < half_max_i_ano:
                            break
                    right_idx = max_i_idx_ano + right_idx_ano
                    max_rt = self.x_grid[right_idx] + half_rt_length
        #TODO: check signal-to-noise 
        return min_rt, max_rt


    def cos_in_smoothed_spec(self, 
                             cos, 
                             precision_da,
                             min_plus1_i = 0.05,
                             min_mono_i = 0.15):
        if cos.dp == 1:
            min_plus1_i = 0.01
#         return_points = []
        mono_mz_range, pl1_mz_range = None, None
        mono_i = []
        pl1_i = []
        pl2_i = []
        min_mz = cos.mz - precision_da
        max_mz = cos.mz + precision_da
        pl1_min_mz = cos.plus1_mz - (precision_da * 1.15)
        pl1_max_mz = cos.plus1_mz + (precision_da * 1.15)
        pl2_min_mz = cos.plus2_mz - (precision_da * 1.3)
        pl2_max_mz = cos.plus2_mz + (precision_da * 1.3)
        for mz, i in zip(self.x_grid, self.i_smooth):
            if min_mz < mz < max_mz:
                mono_i.append(i)
            elif pl1_min_mz < mz < pl1_max_mz:
                pl1_i.append(i)
            elif pl2_min_mz < mz < pl2_max_mz:
                pl2_i.append(i)
        
        if mono_i and pl1_i and pl2_i:
            mono_max = max(mono_i)
            pl1_max = max(pl1_i)
            pl2_max = max(pl2_i)
            mono_idx = self.i_smooth.index(mono_max)
            pl1_idx = self.i_smooth.index(pl1_max)
            pl2_idx = self.i_smooth.index(pl2_max)
            if (pl1_max > min_plus1_i * mono_max and 
                mono_max > min_mono_i * pl1_max):
                if (self.i_smooth[mono_idx-1] < mono_max > 
                    self.i_smooth[mono_idx+1]): 
                    # local maximum for monoisotopic peak
                    if (self.i_smooth[pl1_idx-1] < pl1_max > 
                        self.i_smooth[pl1_idx+1]):
#                         print(cos, self.x_grid)
                        if ((not pl2_idx+2 > self.smoothlen) and 
                            self.i_smooth[pl2_idx-1] < pl2_max > 
                            self.i_smooth[pl2_idx+1]):
                            mono_half_max_i = mono_max / 2
                            for mono_right_idx, i in enumerate(
                                self.i_smooth[mono_idx:]):
                                if i < mono_half_max_i:
                                    break
                            mono_left_idx = None
                            for mono_left_idx, i in enumerate(
                                reversed(self.i_smooth[:mono_idx])):
                                if i < mono_half_max_i:
                                    break
                            mono_right_idx += mono_idx
                            if mono_left_idx:
                                mono_left_idx = mono_idx - mono_left_idx
                            else:
                                mono_left_idx = 0
                            
                            if (self.x_grid[mono_right_idx] - 
                                self.x_grid[mono_left_idx] <
                                3 * precision_da):
                            
                                pl1_half_max_i = pl1_max / 2
                                for pl1_right_idx, i in enumerate(
                                    self.i_smooth[pl1_idx:]):
                                    if i < pl1_half_max_i:
                                        break
                                for pl1_left_idx, i in enumerate(
                                    reversed(self.i_smooth[:pl1_idx])):
                                    if i < pl1_half_max_i:
                                        break
                                pl1_right_idx += pl1_idx
                                pl1_left_idx = pl1_idx - pl1_left_idx
                                
                                if pl1_left_idx > mono_right_idx:
                                    if (self.x_grid[pl1_right_idx] - 
                                        self.x_grid[pl1_left_idx] <
                                        3 * precision_da):
                                        if not (max(self.i_smooth[mono_right_idx:pl1_left_idx])
                                                > (self.i_smooth[mono_right_idx] or 
                                                   self.i_smooth[pl1_left_idx])):
                                            mono_half_fwhm = (self.x_grid[mono_right_idx] - 
                                                              self.x_grid[mono_left_idx]) / 2
                                            pl1_half_fwhm = (self.x_grid[pl1_right_idx] - 
                                                             self.x_grid[pl1_left_idx]) / 2
                                            mono_mz_range = (self.x_grid[mono_left_idx] - mono_half_fwhm,
                                                             self.x_grid[mono_right_idx] + mono_half_fwhm)
                                            pl1_mz_range = (self.x_grid[pl1_left_idx] - pl1_half_fwhm,
                                                            self.x_grid[pl1_right_idx] + pl1_half_fwhm)
    #                                         return_points += [
    #                                             (self.x_grid[mono_idx],
    #                                              self.i_smooth[mono_idx]),
    #                                             (self.x_grid[pl1_idx],
    #                                              self.i_smooth[pl1_idx]),
    #                                             (self.x_grid[pl2_idx],
    #                                              self.i_smooth[pl2_idx])
    #                                                           ]
        
        return mono_mz_range, pl1_mz_range

        
def cos_in_spec(cos,
                peaks,
                precision_da = False,
                mono_only = False,
                check_charge = False,
                return_mz = False,
                min_plus1_i = 0.05,
                min_mono_i = 0.15):
    """Check for oligo precursor_type peak in spectrum
    (monoisotopic, first, and second isotope peaks)
    """
    intensity = None
    measured_mz = None
    min_mz = cos.mz - precision_da
    max_mz = cos.mz + precision_da
    try:
        mono_i = sum([_i for _mz, _i in peaks 
                      if min_mz < _mz < max_mz])
    except IndexError:
        mono_i = False
    if mono_i:
        if mono_only:
            intensity = mono_i
            if return_mz:
                mz_values = [(_mz, _i) for _mz, _i in peaks
                             if min_mz < _mz < max_mz]
                # weighted average
                measured_mz = sum([_mz * _i for _mz, _i in mz_values]
                                  ) / sum([_i for _mz, _i in mz_values]) 
        else:
            pl1_min_mz = cos.plus1_mz - precision_da
            pl1_max_mz = cos.plus1_mz + precision_da
            pl1_i = sum([_i for _mz, _i in peaks
                           if pl1_min_mz < _mz < pl1_max_mz])
            if pl1_i > min_plus1_i * mono_i and mono_i > min_mono_i * pl1_i:
                if return_mz:
                    mz_values = [(_mz, _i) for _mz, _i in peaks
                                 if min_mz < _mz < max_mz]
                    # weighted average
                    measured_mz = sum([_mz * _i for _mz, _i in mz_values]) / sum([_i for _mz, _i in mz_values]) 
                intensity = mono_i + pl1_i
                # check if second peak's intensity is within 5-100% of 
                # first peak's intensity:
                plus2_min_mz = cos.plus2_mz - precision_da
                plus2_max_mz = cos.plus2_mz + precision_da
                pl2_i = sum([_i for _mz, _i in peaks
                               if plus2_min_mz < _mz < plus2_max_mz])
                if pl2_i:
                    intensity += pl2_i
                if check_charge:
                    double_charge = 2 * abs(cos.charge)
                    double_mz = cos.mz + (chem.ATOMIC_MASSES['1n'] / 
                                          double_charge)
                    double_min_mz = double_mz - precision_da
                    double_max_mz = double_mz + precision_da
                    double_i = sum([_i for _mz, _i in peaks 
                                    if double_min_mz < _mz < double_max_mz])
                    if double_i > pl1_i:
                        intensity = None
    return measured_mz, intensity


class MzmlFile():
    def __init__(self, 
                 file_path,
                 precision_da = 0.12):
        """Set up pymzml run and analyze ms1 or ms2 spectra
        """
        self.file_path = file_path
        self.precision_da = precision_da
#         self.precursor_precision_da = precursor_precision_da
        self.ms1_precision = precision_da/1000
        self.msn_precision = precision_da/1000
        self.high_precision = 5e-6
        self.file_name = os.path.split(file_path)[1].replace('.mzML', '')
        self.rough_noise_rt_window = 0.1
        self.smooth_noise_rt_window = 0.5
#         self.noise_mz_window = 50
        print('Analyzing "{0}"'.format(self.file_name))
        
        
    def run(self):
        """Set up pymzml.run.Reader
        """
        _run = pymzml.run.Reader(self.file_path, 
                                MS1_Precision = self.high_precision,
                                MSn_Precision = self.high_precision,
                                )
        return _run
    
    
    @misc.LazyFunction
    def noise_levels(self):
        """Determine retention time dependent noise levels in MS1 spectra 
        """
        ms_level = 1
        rt_to_round_rt = {1: {}, -1: {}}
        round_rt_set = {1: set(), -1: set()}
        rt_intensities = {1: ddict(list), -1: ddict(list)}
        for spectrum in self.run():
            if spectrum.ms_level == ms_level:
                ion_mode = 1
                for element in spectrum.element.getiterator():
                    if element.get('accession', default = None
                                   ) == 'MS:1000129': #i.e. negative scan
                        ion_mode = -1
                round_rt = round(spectrum['scan start time'] / 
                    self.rough_noise_rt_window) * self.rough_noise_rt_window
                rt_to_round_rt[ion_mode][spectrum['scan start time']] = round_rt
                round_rt_set[ion_mode].add(round_rt)
                rt_intensities[ion_mode
                               ][spectrum['scan start time']
                                 ] += list(spectrum.i)
                
        rough_noise_levels = {1: {}, -1: {}}
        for ion_mode in (1, -1):
            for round_rt in round_rt_set[ion_mode]:
                min_rt = round_rt - self.rough_noise_rt_window
                max_rt = round_rt + self.rough_noise_rt_window
                intensities_in_window = []
                #TODO: faster
                for _rt, _i in rt_intensities[ion_mode].items():
                    if min_rt < _rt < max_rt:
                        intensities_in_window += _i
                try:
                    intensities_in_window[3]
                except IndexError:
                    continue
                intensities_in_window.sort()
                mad1 = misc.median_absolute_deviation(intensities_in_window, 
                                                      sorting = False)
                std_dev_estimate1 = misc.std_dev_from_mad(mad1)
                intensity_threshold = 3*std_dev_estimate1
                idx = bisect.bisect(intensities_in_window, intensity_threshold)
                below_threshold = intensities_in_window[:idx]
                mad2 = misc.median_absolute_deviation(below_threshold, 
                                                      sorting = False)
                std_dev_estimate2 = misc.std_dev_from_mad(mad2)
                rough_noise_levels[ion_mode][round_rt] = std_dev_estimate2
            
        noise_levels = {1: {}, -1: {}}
        for ion_mode in (1, -1):
            for round_rt in round_rt_set[ion_mode]:
                min_rt = round_rt - self.smooth_noise_rt_window
                max_rt = round_rt + self.smooth_noise_rt_window
                noise_in_window = [std_dev for _rt, std_dev in 
                                   rough_noise_levels[ion_mode].items() 
                                   if min_rt < _rt < max_rt]
                noise_levels[ion_mode][round_rt] = sum(noise_in_window) / len(noise_in_window)
        return noise_levels
    
    
    @misc.LazyFunction
    def noise_levels_rt_and_mz(self):
        """Determine retention time dependent noise levels in MS1 spectra 
        also m/z dependent
        """
        ms_level = 1
        rt_to_round_rt = {1: {}, -1: {}}
        round_rt_set = {1: set(), -1: set()}
        rt_intensities = {1: ddict(list), -1: ddict(list)}
        round_mz_list = list(range(self.noise_mz_window, 
                                   1800, 
                                   self.noise_mz_window))
        mz_window_size = self.noise_mz_window
        for spectrum in self.run():
            if spectrum.ms_level == ms_level:
                ion_mode = 1
                for element in spectrum.element.getiterator():
                    if element.get('accession', default = None
                                   ) == 'MS:1000129': #i.e. negative scan
                        ion_mode = -1
                round_rt = round(spectrum['scan start time'] / 
                    self.rough_noise_rt_window) * self.rough_noise_rt_window
                rt_to_round_rt[ion_mode][spectrum['scan start time']] = round_rt
                round_rt_set[ion_mode].add(round_rt)
                for mz in round_mz_list:
                    lowidx = bisect.bisect(spectrum.mz, mz-mz_window_size)
                    highidx = bisect.bisect(spectrum.mz, mz+mz_window_size)
                    rt_intensities[ion_mode
                                   ][spectrum['scan start time'], mz
                                     ] += list(spectrum.i)[lowidx:highidx]
                
        rough_noise_levels = {1: {}, -1: {}}
        for ion_mode in (1, -1):
            for round_rt in round_rt_set[ion_mode]:
                min_rt = round_rt - self.rough_noise_rt_window
                max_rt = round_rt + self.rough_noise_rt_window
                intensities_in_window = []
                #TODO: faster
                for mz in round_mz_list:
                    for (_rt, _mz), _i in rt_intensities[ion_mode].items():
                        if (min_rt < _rt < max_rt) and _mz == mz:
                            intensities_in_window += _i
                    try:
                        intensities_in_window[3]
                    except IndexError:
                        continue
                    intensities_in_window.sort()
                    mad1 = misc.median_absolute_deviation(intensities_in_window, 
                                                          sorting = False)
                    std_dev_estimate1 = misc.std_dev_from_mad(mad1)
                    intensity_threshold = 3*std_dev_estimate1
                    idx = bisect.bisect(intensities_in_window, intensity_threshold)
                    below_threshold = intensities_in_window[:idx]
                    mad2 = misc.median_absolute_deviation(below_threshold, 
                                                          sorting = False)
                    std_dev_estimate2 = misc.std_dev_from_mad(mad2)
                    rough_noise_levels[ion_mode][round_rt, mz] = std_dev_estimate2
            
        noise_levels = {1: {}, -1: {}}
        for ion_mode in (1, -1):
            for round_rt in round_rt_set[ion_mode]:
                for mz in round_mz_list:
                    min_rt = round_rt - self.smooth_noise_rt_window
                    max_rt = round_rt + self.smooth_noise_rt_window
                    noise_in_window = [std_dev for (_rt, _mz), std_dev in 
                                       rough_noise_levels[ion_mode].items() 
                                       if (min_rt < _rt < max_rt) and mz == _mz]
                    noise_levels[ion_mode][round_rt, mz] = sum(noise_in_window) / len(noise_in_window)
        return noise_levels
    
             
    def ms1(self, 
            cos_cand,
            min_xic_len = 5, 
            out_file = None, 
            plotting = False,
            min_snr = 3,
            mz_correction = None,
            force_rt = False,
            deconv = False):
        """Create extracted ion chromatograms (XICs) of oligo; find elution peaks;
        determine intensity of signal and amount of oligo if standard is present.
        Write csv output file and plot chromatograms.
        """
        cos_list = [standard for _, standard in 
                    sorted(cos_cand.standards.items())] + cos_cand.cos
        if mz_correction:
            cos_list = copy.deepcopy(cos_list)
            for cos in cos_list:
                cos.correct(mz_correction[self.file_name])
        xi_specs = {cos.key: [] 
                    for cos in cos_list}
        peak_rts = {}
        bpc = Chromatogram() #TODO: usually not required
        for spectrum in self.run():
            if spectrum.ms_level == 1:
                # <cvParam cvRef="MS" accession="MS:1000130" name="positive scan"/>
                # <cvParam cvRef="MS" accession="MS:1000129" name="negative scan"/>
                ion_mode = 1
                for element in spectrum.element.getiterator():
                    if element.get('accession', default = None
                                   ) == 'MS:1000129': #i.e. negative scan
                        ion_mode = -1
                try:
                    bpc.add(spectrum['scan start time'], 
                            spectrum.highest_peaks(1)[0][1])
                except:
                    continue
                for cos in cos_list:
                    if ion_mode * cos.charge < 0:
                        continue
                    if not cos.min_rt or cos.min_rt < spectrum['scan start time'
                                                               ] < cos.max_rt:
                        left = bisect.bisect_left(spectrum.mz, cos.mz - 2) 
                        right = bisect.bisect_right(spectrum.mz, cos.mz + 5) 
                        spec_section = [
                            (mz, i) 
                            for mz, i in
                            zip(spectrum.mz[left:right], spectrum.i[left:right])]
#                             spectrum.peaks('raw')[left:right]]
                        xi_specs[cos.key].append((spectrum['scan start time'],
                                                  spec_section))
        cos_i = {}
        if plotting:
            plot.plot_chromatograms(title = 'Base Peak Chromatogram (BPC)',
                                    mtext = self.file_name,
                                    chromatograms_line = [bpc.xy])
            plot.plot_chromatograms(
                title = 'Noise Levels',
                mtext = self.file_name,
                chromatograms_line = [[(x, y) for x, y in 
                                       sorted(self.noise_levels[ion_mode].items())] 
                                      for ion_mode in sorted(self.noise_levels)]
                )
        for cos in cos_list:
            if not xi_specs[cos.key]:
                continue
            xic = Chromatogram() 
            for rt, spectrum in xi_specs[cos.key]:
                _mz, _intensity = cos_in_spec(cos, 
                                              spectrum,
                                              precision_da = self.precision_da, 
                                              check_charge = True)
                if not _intensity:
                    continue
                xic.add(rt, _intensity)
            if xic.length < min_xic_len:
                continue
            if cos.standard in peak_rts:
                min_peak_rt, max_peak_rt = peak_rts[cos.standard]
            else:
                if not force_rt:
                    min_peak_rt, max_peak_rt = xic.elution_peak()
                else:
                    min_peak_rt, max_peak_rt = cos.min_rt, cos.max_rt
            if (not min_peak_rt):# or 
#                 (cos.min_rt and min_peak_rt < cos.min_rt) or # NOTE this caused problems when elution peaks were very broad (almost not visible)
#                 (cos.max_rt and max_peak_rt > cos.max_rt)):
                continue
            num_spectra = 0
            noise_level = 0
            all_peaks = []
            ion_mode = 1 if cos.charge > 0 else -1
            for rt, spectrum in xi_specs[cos.key]:
                if min_peak_rt < rt < max_peak_rt:
                    num_spectra += 1
                    all_peaks += spectrum
                    noise_level += self.noise_levels[ion_mode][
                        round(rt / self.rough_noise_rt_window) * 
                        self.rough_noise_rt_window]
#                     , 
#                         round(cos.mz / self.noise_mz_window) * self.noise_mz_window]
                        # NOTE TODO summing noise levels - is this valid??
            if not all_peaks:
                continue
            all_peaks.sort()
            
            smoothed_spec = Chromatogram()
            for mz, i in all_peaks:
                smoothed_spec.add(mz, i)
            smoothed_spec.smooth(self.precision_da/2, self.precision_da/10)
            mono_mz_range, pl1_mz_range = smoothed_spec.cos_in_smoothed_spec(cos, 
                                               self.precision_da)
            if not mono_mz_range:
                continue
            intensity = 0
            tmp_mz = 0
            tmp_i = 0
            for mz, i in all_peaks:
                if mono_mz_range[0] < mz < mono_mz_range[1]:
                    intensity += i
                    tmp_i += i
                    tmp_mz += mz*i
                elif pl1_mz_range[0] < mz < pl1_mz_range[1] and not deconv:
                    intensity += i
            if deconv and "D" not in cos:
                intensity = chem.isotopologue_correction(cos, intensity)
            snr = intensity / noise_level
            if snr < min_snr:
                continue
            try:
                measured_mz = tmp_mz / tmp_i
            except: # i.e. plus1 isotpe peak is there but no monoisotopic
                continue
#             if abs(measured_mz - cos.mz) > self.precision_da:
#                 print(cos.name, cos.precursor_type, cos.mz, measured_mz)
#                 continue
            cos.intensity = intensity
            cos_i[cos.key] = intensity
            peak_rts[cos.key] = (min_peak_rt, max_peak_rt)
#             ng = ''
            if out_file:
                out_dict = {'File name': self.file_name, 
                            'Intensity [arb. unit]': '{0:e}'.format(intensity), 
#                             'Mass [ng]': ng,
                            'SNR': round(snr, 1),
                            'Number of spectra': num_spectra,
                            'Peak start time [min]': min_peak_rt,
                            'Peak end time [min]': max_peak_rt,
                            'Peak length [min]': max_peak_rt-min_peak_rt,
                            'Measured m/z': round(measured_mz, 3)}
                if cos.standard and cos.standard in cos_i:
                    ng = intensity / cos_i[cos.standard
                        ] * cos_cand.standards[cos.standard].ng
                    out_dict['Mass [ng]'] = ng
                out_file.add(cos, out_dict)
            
            if plotting and snr > min_snr:
                plot.plot_chromatograms(
                    title = '{0} {1}\nSNR = {2:.1f}'.format(cos.name, 
                                                  cos.precursor_type,
                                                  snr),
                    mtext = self.file_name,
                    chromatograms_area = [[(x, y) for x, y in 
                                           zip(xic.x_values, xic.intensities) 
                                           if min_peak_rt < x < max_peak_rt]],
                    chromatograms_points = [[(x, y) for x, y in 
                                           zip(xic.x_values, xic.intensities)]],
                    labels = ['XIC'],
                    rt_range = [bpc.x_values[0], bpc.x_values[-1]]
                    )
                plot.plot_spectrum([(mz, i) for mz, i in zip(smoothed_spec.x_grid,
                                                             smoothed_spec.i_smooth)], 
                                   cos = cos,
                                   title = '{0} {1}\n{2} spectra'.format(
                                        cos.name,
                                        cos.precursor_type,
                                        num_spectra),
                                   plot_type = 'l',
                                   lwd = 1)
        
        
    def ms2(self, 
            cos_cand,
            out_file = None, 
            seq_file = None,
            plotting = False,
            precursor_precision = 0.3,
            fragments_per_spec_threshold = 1,
            mz_correction = None):
        """Identify fragment ions of precursor oligo in MS2 scans.
        Sequence oligo based on these fragment ions. 
        Write csv output files and plot spectra.
        """
#         precursor_precision = self.precursor_precision_da
        tic_chromatogram = Chromatogram()
        if mz_correction:
            cos_cand = copy.deepcopy(cos_cand)
            for cos in cos_cand.cos:
                cos.correct(mz_correction[self.file_name])
                for (frag_type, dp), fragments in cos.fragment_ions.items():
                    for f in fragments:
                        f.correct(mz_correction[self.file_name])
                    
        precursor_mz_ranges = {p.key: (p.mz - precursor_precision,
                                       p.mz + precursor_precision)
                               for p in cos_cand.cos}
        precursors = {p.key: p for p in cos_cand.cos}
        num_spectra = {p.key: 0 for p in cos_cand.cos}
        retention_times = {p.key: set() for p in cos_cand.cos}
        if plotting:
            merged_spectra = {p.key: pymzml.spec.Spectrum(
                measuredPrecision = self.precision_da / p.mz) 
                          for p in cos_cand.cos}
        merged_fragments = {p.key: {f_type_dp: ddict(int)
                                    for f_type_dp in p.fragment_ions}
                            for p in cos_cand.cos}
        fragments_mz = {p.key: {f_type_dp: ddict(list)
                                for f_type_dp in p.fragment_ions}
                        for p in cos_cand.cos}
        precursor_chromatograms = {p.key: Chromatogram()
                                   for p in cos_cand.cos}
        fragment_ions = {}
        for spectrum in self.run():
            if spectrum.ms_level == 2:
#                 print('\t\t\tMS{0:.0f}\tScan ID {1}     '.format(
#                     spectrum.ms_level, spectrum['id'] ), end = '\r')
                precursor_mz = float(spectrum['MS:1000744'])
                tic = float(spectrum['MS:1000285'])
                tic_chromatogram.add(spectrum['scan start time'], tic)
                precursors_in_spec = []
                for prec_key, (min_mz, max_mz) in precursor_mz_ranges.items():
                    if min_mz < precursor_mz < max_mz:
                        if not precursors[prec_key].min_rt or precursors[
                            prec_key].min_rt < spectrum['scan start time'
                            ] < precursors[prec_key].max_rt:
                            if precursors_in_spec:
                                print("WARNING: more than one precursor mass")
                            precursors_in_spec.append(prec_key)
                if not precursors_in_spec:
                    continue 
                try:
                    highest_peaks = spectrum.highest_peaks(10)
                except AttributeError:
                    if not spectrum.peaks('raw'): # i.e. empty spectrum
                        continue
                    else:
                        print('Error in spectrum #', spectrum['id'])
                        raise
                for prec_key in precursors_in_spec:
                    found_fragments = 0
                    tmp_fragments = {}
                    for (frag_type, dp), fragments in precursors[prec_key
                            ].fragment_ions.items():
                        for f in fragments:
                            _mz, _i = cos_in_spec(f, 
                                                  highest_peaks, 
                                                  self.precision_da, 
                                                  mono_only = True)
                            if _i:
                                found_fragments += 1
                                tmp_fragments[frag_type, dp, f.key] = _i
                    
                    if found_fragments > fragments_per_spec_threshold: 
                        #filtering has impact on sequencing results
                        num_spectra[prec_key] += 1
                        retention_times[prec_key].add(spectrum['scan start time'])
                        precursor_chromatograms[prec_key
                                                ].add(spectrum['scan start time'], 
                                                      tic)
                        for (frag_type, dp), fragments in precursors[prec_key
                            ].fragment_ions.items():
                            for f in fragments:
                                _mz, _i = cos_in_spec(f, 
                                                      spectrum.peaks('raw'), 
                                                      self.precision_da, 
                                                      mono_only = True,
                                                      return_mz = True)
                                if _i:
                                    fragment_ions[prec_key, f.key] = f
                                    merged_fragments[prec_key][frag_type, dp
                                                           ][f.key
                                                             ] += _i
                                    fragments_mz[prec_key][frag_type, dp
                                                           ][f.key
                                                             ].append((_mz, _i))
                        if plotting:
                            merged_spectra[prec_key] += spectrum.deRef()
        
        for prec_key, num in num_spectra.items():
            if num < 1:
                continue
            if plotting and merged_spectra[prec_key].peaks('raw'):
                plot.plot_chromatograms(
                    title = '{0} {1}'.format(precursors[prec_key].name, 
                                             precursors[prec_key].precursor_type),
                    mtext = self.file_name,
                    chromatograms_points = [precursor_chromatograms[prec_key].xy]
#                     labels = ['XIC'],
#                     rt_range = [bpc.x_values[0], bpc.x_values[-1]]
                    )
                plot.plot_spectrum(
                    merged_spectra[prec_key], 
                    cos = precursors[prec_key],
                    title = '{0} {1}\n({2} MS2 scans)'.format(
                        precursors[prec_key].name, 
                        precursors[prec_key].precursor_type,
                        num_spectra[prec_key]))
                    
            min_rt = round(min(retention_times[prec_key]), 2)
            max_rt = round(max(retention_times[prec_key]), 2)
            elution_duration = max_rt - min_rt
            potential_precursor_sequences = {}
            prec_fragments = ddict(dict)
            for (frag_type, dp), frag_intensity in merged_fragments[prec_key].items():
                total_fragtype_intensity = sum(frag_intensity.values())
                for f_key, intensity in frag_intensity.items():
                    relative_intensity = intensity / total_fragtype_intensity
                    mz_i_list = fragments_mz[prec_key][frag_type, dp][f_key]
                    measured_mz = (sum([_mz * _i for _mz, _i in mz_i_list]) /
                                   sum([_i for _mz, _i in mz_i_list]))
                    out_dict = {
                        'File name': self.file_name, 
                        'Precursor oligo': precursors[prec_key].name, 
                        'Precursor ion type': precursors[prec_key].precursor_type, 
                        'Precursor m/z': round(precursors[prec_key].mz, 3),
                        'Fragment intensity [arb. unit]': '{0:e}'.format(intensity), 
                        'Relative intensity (w.r.t. fragment type)': round(relative_intensity, 4),
                        'Number of spectra': num,
                        'Peak start time [min]': min_rt,
                        'Peak end time [min]': max_rt,
                        'Peak length [min]': elution_duration,
                        'Measured fragment m/z (average)': measured_mz
                        }
                    prec_fragments[frag_type, dp][f_key] = (intensity, relative_intensity)
                    out_file.add(fragment_ions[prec_key, f_key], out_dict)
                    potential_precursor_sequences[f_key] = fragment_ions[prec_key, f_key].precursor_sequences 
            #Sequencing:
            seq_results = seq.determine_sequences(precursor = precursors[prec_key],
                                    fragments = prec_fragments,
                                    prec_seq = potential_precursor_sequences)
            if seq_results:
                for row in seq_results:
#                     print(row)
                    row['File name'] = self.file_name
                    row['Number of spectra'] = num
                    row['Peak start time [min]'] = min_rt
                    row['Peak end time [min]'] = max_rt
                    row['Peak length [min]'] = elution_duration
                    seq_file.add(precursors[prec_key], row)
            
        if plotting:
            plot.plot_chromatograms(
                title = 'Total Ion Current (TIC) Chromatogram',
                mtext = self.file_name,
                chromatograms_line = [tic_chromatogram.xy])


    def ms2screen(self, 
            cos_cand,
            out_file = None, 
            screen_file = None,
            precursor_precision = 0.3,
            fragments_per_spec_threshold = 1,
            mz_correction = None):
        """Identify fragment ions of precursor oligo in MS2 scans.
        For non-quantitative data intended for PCA analysis (mutant screening).
        Write csv output files.
        """
        #NOTE: TODO: this is basically a copy of ms2(), clean up duplicate code
#         precursor_precision = self.precursor_precision_da
        tic_chromatogram = Chromatogram()
        if mz_correction:
            cos_cand = copy.deepcopy(cos_cand)
            for cos in cos_cand.cos:
                cos.correct(mz_correction[self.file_name])
                for (frag_type, dp), fragments in cos.fragment_ions.items():
                    for f in fragments:
                        f.correct(mz_correction[self.file_name])
                    
        precursor_mz_ranges = {p.key: (p.mz - precursor_precision,
                                       p.mz + precursor_precision)
                               for p in cos_cand.cos}
        precursors = {p.key: p for p in cos_cand.cos}
        num_spectra = {p.key: 0 for p in cos_cand.cos}
        retention_times = {p.key: set() for p in cos_cand.cos}
        merged_fragments = {p.key: {f_type_dp: ddict(int)
                                    for f_type_dp in p.fragment_ions}
                            for p in cos_cand.cos}
        fragments_mz = {p.key: {f_type_dp: ddict(list)
                                for f_type_dp in p.fragment_ions}
                        for p in cos_cand.cos}
        precursor_chromatograms = {p.key: Chromatogram()
                                   for p in cos_cand.cos}
        fragment_ions = {}
        for spectrum in self.run():
            if spectrum.ms_level == 2:
#                 print('\t\t\tMS{0:.0f}\tScan ID {1}     '.format(
#                     spectrum.ms_level, spectrum['id'] ), end = '\r')
                precursor_mz = float(spectrum['MS:1000744'])
                tic = float(spectrum['MS:1000285'])
                tic_chromatogram.add(spectrum['scan start time'], tic)
                precursors_in_spec = []
                for prec_key, (min_mz, max_mz) in precursor_mz_ranges.items():
                    if min_mz < precursor_mz < max_mz:
                        if not precursors[prec_key].min_rt or precursors[
                            prec_key].min_rt < spectrum['scan start time'
                            ] < precursors[prec_key].max_rt:
                            if precursors_in_spec:
                                print("WARNING: more than one precursor mass")
                            precursors_in_spec.append(prec_key)
                if not precursors_in_spec:
                    continue 
                highest_peaks = spectrum.highest_peaks(10)
                for prec_key in precursors_in_spec:
                    found_fragments = 0
                    tmp_fragments = {}
                    for (frag_type, dp), fragments in precursors[prec_key
                            ].fragment_ions.items():
                        for f in fragments:
                            _mz, _i = cos_in_spec(f, 
                                                  highest_peaks, 
                                                  self.precision_da, 
                                                  mono_only = True)
                            if _i:
                                found_fragments += 1
                                tmp_fragments[frag_type, dp, f.key] = _i
                    
                    if found_fragments > fragments_per_spec_threshold: 
                        #filtering has impact on sequencing results
                        num_spectra[prec_key] += 1
                        retention_times[prec_key].add(spectrum['scan start time'])
                        precursor_chromatograms[prec_key
                                                ].add(spectrum['scan start time'], 
                                                      tic)
                        for (frag_type, dp), fragments in precursors[prec_key
                            ].fragment_ions.items():
                            for f in fragments:
                                _mz, _i = cos_in_spec(f, 
                                                      spectrum.peaks('raw'), 
                                                      self.precision_da, 
                                                      mono_only = True,
                                                      return_mz = True)
                                if _i:
                                    fragment_ions[prec_key, f.key] = f
                                    merged_fragments[prec_key][frag_type, dp
                                                           ][f.key
                                                             ] += _i
                                    fragments_mz[prec_key][frag_type, dp
                                                           ][f.key
                                                             ].append((_mz, _i))
        
        for prec_key, num in num_spectra.items():
            if num < 1:
                continue
            min_rt = min(retention_times[prec_key])
            max_rt = max(retention_times[prec_key])
            elution_duration = max_rt - min_rt
            potential_precursor_sequences = {}
            prec_fragments = ddict(dict)
            for (frag_type, dp), frag_intensity in merged_fragments[prec_key].items():
                total_fragtype_intensity = sum(frag_intensity.values())
                for f_key, intensity in frag_intensity.items():
                    relative_intensity = intensity / total_fragtype_intensity
                    mz_i_list = fragments_mz[prec_key][frag_type, dp][f_key]
                    measured_mz = (sum([_mz * _i for _mz, _i in mz_i_list]) /
                                   sum([_i for _mz, _i in mz_i_list]))
                    out_dict = {
                        'File name': self.file_name, 
                        'Precursor oligo': precursors[prec_key].name, 
                        'Precursor ion type': precursors[prec_key].precursor_type, 
                        'Precursor m/z': round(precursors[prec_key].mz, 3),
                        'Fragment intensity [arb. unit]': '{0:e}'.format(intensity), 
                        'Relative intensity (w.r.t. fragment type)': round(relative_intensity, 4),
                        'Number of spectra': num,
                        'Peak start time [min]': '{0:.2f}'.format(min_rt),
                        'Peak end time [min]': '{0:.2f}'.format(max_rt),
                        'Peak length [min]': elution_duration,
                        'Measured fragment m/z (average)': measured_mz
                        }
                    prec_fragments[frag_type, dp][f_key] = (intensity, relative_intensity)
                    out_file.add(fragment_ions[prec_key, f_key], out_dict)
                    potential_precursor_sequences[f_key] = fragment_ions[prec_key, f_key].precursor_sequences
                    # add data to the screen file
                    cos = fragment_ions[prec_key, f_key]
                    f_ion_type = '{0}{1}'.format(cos.fragment_type, cos.dp)
                    f_ion_type += cos.charge * '+' if cos.charge > 1 else ''
                    entry_key = '{0}_{1}_{2}'.format(
                        precursors[prec_key].name, cos.name, f_ion_type)
                    screen_file.add(self.file_name, entry_key, intensity)
            # add TIC to screen file
            entry_key = '{0}_TIC'.format(
                precursors[prec_key].name)
            intensity = 0
            for x, y in tic_chromatogram.xy:
                if min_rt <= x <= max_rt:
                    intensity += y
            screen_file.add(self.file_name, entry_key, intensity)
            
    
    
def file_list(path, file_type = 'mzML'):
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


def assign_ms_type(out, mzml_analysis):
    def modify_analysis(func):         
        def _analysis(mzml_files = None,
                      cos_candidates_file_name = None,
                      out_file_name = None,
                      precision_da = 0.12,
                      precursor_precision_da = 0.3,
                      plotting = False,
                      mz_correction = False,
                      min_snr = 3,
                      force_rt = False,
                      deconv = False):
            cos_candidates = inout.CandidatesFile(cos_candidates_file_name)
            cos_candidates.read()
            out_file = out(out_file_name)
            kwargs = {'out_file': out_file, 
                     'plotting': plotting}
            if mzml_analysis.__name__ == 'ms2':
                seq_file = inout.SequencesFile(out_file_name.replace('.csv', '_seq.csv'))
                kwargs['seq_file'] = seq_file 
                kwargs['precursor_precision'] = precursor_precision_da
            elif mzml_analysis.__name__ == 'ms2screen':
                screen_file = inout.MS2ScreenFile(out_file_name.replace('.csv', '_screen.csv'))
                kwargs['screen_file'] = screen_file 
                kwargs['precursor_precision'] = precursor_precision_da
                del kwargs['plotting']
            else:
                kwargs['min_snr'] = min_snr
                kwargs['force_rt'] = force_rt
                kwargs['deconv'] = deconv
            if mz_correction:
                correction_file = inout.CorrectionFile(mz_correction)
                kwargs['mz_correction'] = correction_file.corr
            if plotting:
                plot_file = plot.PlotFile(out_file_name.replace('.csv', '.pdf'))
            for num, mzml_file_path in enumerate(file_list(mzml_files, file_type = 'mzML')):
                print('\n{0}\t#{1}'.format(mzml_analysis.__name__, num+1).encode("utf-8").decode("ascii"))
                mzml_file = MzmlFile(mzml_file_path, 
                                     precision_da = precision_da)
                mzml_analysis(mzml_file, 
                              cos_candidates, 
                              **kwargs)
            if mzml_analysis.__name__ == 'ms1' and deconv:
                out_file.overlap_corr()
            out_file.write()
            if mzml_analysis.__name__ == 'ms2':
                seq_file.write()
            if mzml_analysis.__name__ == 'ms2screen':
                screen_file.write()
            if plotting:
                plot_file.close()
        return _analysis
    return modify_analysis


@assign_ms_type(inout.MS1File, MzmlFile.ms1)
def ms1_analysis():
    return


@assign_ms_type(inout.FragmentFile, MzmlFile.ms2)
def ms2_analysis():
    return


@assign_ms_type(inout.FragmentFile, MzmlFile.ms2screen)
def ms2_screen():
    return


# import unittest
# from cosms.tests import testing
#  
# class TestMS1(testing.ProfilingTest):
#     def testMS1(self):
# #     def setUp(self):
#         self.dir = "/home/anna/development/chitosan/cosms/tests/TestMS1/"
#         self.mzml = "DP1-6_AR_Rstern.mzML"
#         self.mzml = "1_pLUX01_RFP_100%Glc_RD1_01_47057.d.mzML"
#         self.cand = "COS_AR_DP1-6_H_.csv"
#         self.cand = "COS_AHex_DP1-10_2H_2Na_H_NH4_Na_.csv"
#         ms1_analysis(
#             mzml_files = os.path.join(self.dir, self.mzml),
#             cos_candidates_file_name = os.path.join(self.dir, self.cand),
#             out_file_name = os.path.join(self.dir, "MS1results.csv"),
#             precision_da = 0.2,
#             plotting = True)
        

# from cosms.tests import testing
# class TestMzmlMS1_AR(testing.ProfilingTest):
#     def test_mzml_ms1_AR(self):
#         self.mzml_file = '../tests/Test Data/DP5_AR.mzML'
#         self.cand = '../tests/Test Data/DP5_AR_candidates.csv'
#         _cand = inout.CandidatesFile(self.cand)
#         _cand.create_all_possible(['A', 'R'], 5, 5, ['H+'], [1], 
#                                   min_rt=1.5, max_rt=3)
#         _cand.write()
#         ms1_analysis(self.mzml_file, 
#                      self.cand, 
#                      '../tests/Test Data/DP5_AR_MS1-results.csv',
#                      plotting = True)
# 
# 
# class TestMzmlMS1_AR_Rstern(cosms.tests.testing.ProfilingTest):
#     def test_mzml_ms1_AR_Rstern(self):
#         self.mzml_file = '../Test Data/DP1-6_AR_Rstern.mzML'
#         self.cand = '../Test Data/DP1-6_AR_Rstern_candidates.csv'
#         _cand = inout.CandidatesFile(self.cand)
#         _cand.create_all_possible(['A', 'R'], 1, 6, ['H+'], [1, 2], [1],
#                                   standard = 'R*', std_mass = 75)
#         _cand.write()
#         ms1_analysis(self.mzml_file, 
#                      self.cand, 
#                      '../Test Data/DP1-6_AR_Rstern_MS1-results.csv',
#                      plotting = True)


# class TestMzmlMS2_AR_18O(cosms.tests.testing.ProfilingTest):
#     def test_mzml_ms2_AR_18O(self):
#         self.mzml_file = '../Test Data/DP1-6_AR_18O_MS2.mzML'
#         self.cand = '../Test Data/DP1-6_AR_18O_candidates.csv'
#         _cand = inout.CandidatesFile(self.cand)
#         _cand.create_all_possible(['A', 'R'], 1, 6, ['H+'], [1], 
#                                   label = '18O')
#         _cand.write()
#         ms2_analysis(self.mzml_file, 
#                      self.cand, 
#                      '../Test Data/DP1-6_AR_18O_MS2-results.csv',
#                      plotting = False)

# class TestMzmlMS2_AR_18O(cosms.tests.testing.ProfilingTest):
#     def test_mzml_ms2_AR_18O(self):
#         self.mzml_file = '/home/anna/Data/Stefan CL/2015-10-26_Proben Pieter MS2/ms2/'
#         self.cand = '/home/anna/Data/Stefan CL/2015-10-26_Proben Pieter MS2/DP4-5_AR_18O_candidates.csv'
# #         _cand = inout.CandidatesFile(self.cand)
# #         _cand.create_all_possible(['A', 'R'], 4, 5, ['H+'], [1], 
# #                                   label = '18O')
# #         _cand.read()
#         ms2_analysis(self.mzml_file, 
#                      self.cand, 
#                      '/home/anna/Data/Stefan CL/2015-10-26_Proben Pieter MS2/DP4-5_AR_18O_MS2-results.csv',
#                      plotting = False)

# class TestMzmlMS1_AR(cosms.tests.testing.ProfilingTest):
#     def test_mzml_ms1_AR(self):
#         self.mzml_file = '/home/anna/Data/Stefan CL/2015-10-26_Proben Pieter MS2/ms1/'
#         self.cand = '/home/anna/Data/Stefan CL/2015-10-26_Proben Pieter MS2/DP4-5_AR_candidates_ms1.csv'
#         _cand = inout.CandidatesFile(self.cand)
#         _cand.create_all_possible(['A', 'R'], 4, 5, ['H+'], [1], label = '18O')
#         _cand.write()
#         ms1_analysis(self.mzml_file, 
#                      self.cand, 
#                      '/home/anna/Data/Stefan CL/2015-10-26_Proben Pieter MS2/DP4-5_AR_MS1-results.csv',
#                      plotting = False)

# class TestMzmlMS2_AR_18O2lea(cosms.tests.testing.ProfilingTest):
#     def test_mzml_ms2_AR_18O(self):
#         self.mzml_file = '/home/anna/Data/Lea H/2015-10-19/MS2_2015-10-26/'#A5 48h A3R2_GC3_01_15499.d.mzML'#A5 144h- A3R2_GC6_01_15510.d.mzML'#A6 144h- A3R3_GD5_01_15542.d.mzML'
#         self.cand = '/home/anna/Data/Lea H/2015-10-19/DP1-6_AR_18O_candidates.csv'
# #         _cand = inout.CandidatesFile(self.cand)
# #         _cand.create_all_possible(['A', 'R'], 1, 6, ['H+'], [1], 
# #                                   label = '18O')
# #         _cand.write()
#         ms2_analysis(self.mzml_file, 
#                      self.cand, 
#                      '/home/anna/Data/Lea H/2015-10-19/DP1-6_AR_18O_MS2-results.csv',
#                      plotting = False)

# class TestMzmlMS1_DP4AR_lea(cosms.testing.ProfilingTest):
#     def test(self):
#         self.mzml_file = '/home/anna/Data/Lea H/151103 Lea/ms1/'
#         self.cand = '/home/anna/Data/Lea H/151103 Lea/DP4_AR_candidates.csv'
#         _cand = inout.CandidatesFile(self.cand)
#         _cand.create_all_possible(['A', 'R'], 4, 4, ['H+'], [1])
#         _cand.write()
#         ms1_analysis(self.mzml_file, 
#                      self.cand, 
#                      '/home/anna/Data/Lea H/151103 Lea/DP4_AR_MS1-results.csv',
#                      plotting = False)

# class TestMzmlMS2_AR_18Opaparao(cosms.tests.testing.ProfilingTest):
#     def test_mzml_ms2_AR_18O(self):
#         self.mzml_file = '/home/anna/Data/Papa Rao/2015-10-29/ms2/'
#         self.cand = '/home/anna/Data/Papa Rao/2015-10-29/DP2-6_AR_18O_candidates.csv'
# #         _cand = inout.CandidatesFile(self.cand)
# #         _cand.create_all_possible(['A', 'R'], 2, 6, ['H+'], [1], 
# #                                   standard = 'R*', std_mass = 112,
# #                                   label = '18O')
# #         _cand.write()
#         ms2_analysis(self.mzml_file, 
#                      self.cand, 
#                      '/home/anna/Data/Papa Rao/2015-10-29/DP2-6_AR_18O_MS2-results.csv',
#                      plotting = False)
