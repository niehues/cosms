#!/usr/bin/env python3
'''
Created on Aug 30, 2017

@author: anna
'''
import csv
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from matplotlib import colors
from scipy.optimize import curve_fit
from scipy.special import lambertw
from scipy.interpolate import interp1d
import seaborn
import sys
from collections import defaultdict as ddict
from matplotlib import rcParams
rcParams['mathtext.default'] = 'regular'


def plot_RTs(infile, all = True):
    in_file = csv.DictReader(open(infile, 'r'))
    plot_dict = {}
    AD = True
    for idx, row in enumerate(in_file):
        try:
            row['D']
            sugar = 'D'
        except:
            try:
                row['R']
            except:
                row['R'] = 0
            sugar = 'R'
            AD = False
        key = (int(row['DP']),
                int(row[sugar]) if row[sugar] != '' else 0,
               row['Oligo'])
        try:
            plot_dict[key]
        except:
            plot_dict[key] = {'Run number': [],
                              'Peak start time [min]': [],
                              'Peak end time [min]': [],
                              'SNR': []}
        if 'Run number' in row.keys():
            plot_dict[key]['Run number'
                           ].append(int(row['Run number']))
        else:
            plot_dict[key]['Run number'
                           ].append(idx)
        plot_dict[key]['Peak start time [min]'
                       ].append(float(row['Peak start time [min]']))
        plot_dict[key]['Peak end time [min]'
                       ].append(float(row['Peak end time [min]']))
        plot_dict[key]['SNR'
                       ].append(float(row['SNR']))               
    figure_count = 1
    subplot_count = 1
    summary_labels = []
    summary_dp = []
    summary_av_rt = []
    all4fit_summary_dp = ddict(list)
    all4fit_summary_av_rt = ddict(list)
    bluemap = plt.get_cmap('Blues')
    redmap = plt.get_cmap('Reds')
    cNorm  = colors.Normalize(vmin = 0, 
                              vmax = 100)
    bluescalarMap = cm.ScalarMappable(norm = cNorm, 
                                   cmap = bluemap)
    redscalarMap = cm.ScalarMappable(norm = cNorm, 
                                   cmap = redmap)
    for idx, ((dp, d, cos), values) in enumerate(sorted(plot_dict.items())):
        if len(values['SNR']) < 4 or dp > 15: #TODO:
            continue
        summary_labels.append((dp, d, cos))
        summary_dp.append(dp)
        meanrt = (np.array(values['Peak start time [min]']) + 
                  np.array(values['Peak end time [min]'])) / 2
        weighted_average = np.average(meanrt, weights = np.log(np.array(values['SNR'])))
        summary_av_rt.append(weighted_average)
        all4fit_summary_dp['A', dp-d] += [dp for _ in range(meanrt.size)]
        all4fit_summary_av_rt['A', dp-d] += list(meanrt)
        all4fit_summary_dp['D', d] += [dp for _ in range(meanrt.size)]
        all4fit_summary_av_rt['D', d] += list(meanrt)
        if all:
            plt.figure(figure_count)
            plt.subplot(int('33{0}'.format(subplot_count)))
            plt.title(cos)
            plt.xlabel('Run number')
            plt.ylabel('Retention time [min]')
            maxsnr = np.log(max(values['SNR']))
            for x,y1,y2,z in zip(values['Run number'],
                                 values['Peak start time [min]'],
                                 values['Peak end time [min]'],
                                 values['SNR']):
                colorvalue = 100*np.log(z)/maxsnr if np.log(z) < maxsnr else 100*maxsnr
                plt.plot(x, y1, '.',
                         color = bluescalarMap.to_rgba(colorvalue))
                plt.plot(x, y2, '.',
                         color = redscalarMap.to_rgba(colorvalue))
            subplot_count += 1
            if subplot_count > 9:
                subplot_count = 1
                figure_count += 1
    
    if AD:
        plt.figure(figure_count+1)
    #     plt.title('Retention times of chito-oligosaccharides')
        plt.plot(summary_av_rt,
                        summary_dp,
                        'w.')
        a_series = ddict(list)
        d_series = ddict(list)
        for idx, (dp, d, cos) in enumerate(summary_labels):
            if d == 0:
                label = r'$A_{{{0}}}$'.format(dp-d)
            elif dp == d:
                label = r'$D_{{{0}}}$'.format(d)
            else:
                label = r'$A_{{{0}}}D_{{{1}}}$'.format(dp-d, d)
            plt.annotate(label,
                        (summary_av_rt[idx],
                         summary_dp[idx]),
                        horizontalalignment = 'center', 
                        verticalalignment = 'center')
            if dp-d > -1:
                a_series[dp-d].append((summary_av_rt[idx], dp))
            if d > -1:
                d_series[d].append((summary_av_rt[idx], dp))
        plt.xlabel('Retention time [min]')
        plt.ylabel('Degree of polymerization')
        if True:
            #### plot fits
            for a, rt_dp in a_series.items():
                x, y = zip(*rt_dp)
                plt.plot(x, y, '-')
                continue
                xnew = np.linspace(x[0], 
                                   x[-1], 
                                   num = int((x[-1]-x[0])*10), 
                                   endpoint = True)
                # cubic splines interpolation
                try:
                    interpolation = interp1d(x, y, kind = 'cubic')
                except:
                    continue
                plt.plot(xnew, interpolation(xnew), '-')
            for d, rt_dp in d_series.items():
                x, y = zip(*rt_dp)
                plt.plot(x, y, '-')
    
    return plt
	

if __name__ == '__main__':
    all = True# if len(sys.argv) < 3 else False
    plot_RTs(sys.argv[1], all)