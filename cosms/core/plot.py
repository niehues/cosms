#!/usr/bin/env python3
# coding: utf-8
"""plotting functions
rpy2 required
"""
import os
import math


#     import rpy2
import rpy2.robjects as robjects
from rpy2.robjects import r
from rpy2.robjects.packages import importr
graphics = importr('graphics')
grdevices = importr('grDevices')


COLORS = {# Masataka Okabe & Kei Ito, 2002
            # http://jfly.iam.u-tokyo.ac.jp/color/#pallet
            'black': r.rgb(0.00, 0.00, 0.00)[0],
            'blue': r.rgb(0.00, 0.45, 0.70)[0], 
            'red': r.rgb(0.80, 0.40, 0.00)[0],
            'green': r.rgb(0.00, 0.60, 0.50)[0],
            'purple': r.rgb(0.80, 0.60, 0.70)[0],
            'orange': r.rgb(0.90, 0.60, 0.00)[0],
            'sky blue': r.rgb(0.35, 0.70, 0.90)[0],
            'yellow': r.rgb(0.95, 0.90, 0.25)[0]
            }

line_colors = [COLORS[col] for col in ['blue', 'red', 'green', 'purple']] * 30
area_colors = [COLORS[col] for col in ['orange', 'sky blue', 'yellow']] * 30


class PlotFile():
def __init__(self, 
                out_file, 
                paper = 'a4', 
                width = 6, 
                height = 8,
                rows = 2,
                columns = 1):
    grdevices.pdf(out_file, paper = paper, width = width, height = height)
    graphics.par(mfrow = r.c(rows, columns))
#         graphics.layout() # TODO

#     def rcode(self, code): # TODO evaluate str or sth similar
#         ri.parse(code) # rinterface requires some update...

def close(self):
    grdevices.dev_off()
    

def empty_plot():
graphics.plot_new()


def plot_chromatograms(title = None,
                    mtext = None,
                    chromatograms_area = None,
                    chromatograms_line = None,
                    chromatograms_points = None,
                    texts = None,
                    labels = None,
                    rt_range = None
#                   spec_mz = None,
#                   spec_i = None
                    ):
"""Plot chromatogram(s) and safe to pdf

:param str title: Title of graphic
:param list chromatograms_area: Lines in background
:param list chromatograms_line: Points in foreground
:param list labels: Names for legend
:param str out_file: Path to output pdf
:param list rt_range: Min. and max. retention time
"""
if not chromatograms_area:
    chromatograms_area = []
if not chromatograms_line:
    chromatograms_line = []
if not chromatograms_points:
    chromatograms_points = []
if not title:
    title = 'Chromatogram'

if rt_range:
    min_rt, max_rt = rt_range[0], rt_range[1]
else:
    all_rt = [rt for chromatogram in chromatograms_area + 
                chromatograms_line + chromatograms_points 
                for rt, _i in chromatogram]
    min_rt = min(all_rt)
    max_rt = max(all_rt)
all_i = [i for chromatogram in chromatograms_area + 
            chromatograms_line + chromatograms_points 
            for rt, i in chromatogram]
min_i = 0
max_i = max(all_i) * 1.05

graphics.plot(0, 0,
                type = 'n',
                main = title,
                xlab = 'Retention Time [min]',
                ylab = 'Intensity [arb. unit]',
                xlim = r.c( min_rt, max_rt ),
                ylim = r.c( min_i, max_i ),
                yaxs = 'i', xaxs = 'i'
                )
graphics.mtext(mtext, adj = 1)

if chromatograms_area:
    for num, chromatogram in enumerate(chromatograms_area):
        x = robjects.FloatVector([chromatogram[0][0]] +
                                    [rt for rt, i in chromatogram ] +
                                    [chromatogram[-1][0]])
        y = robjects.FloatVector([0]+[i for rt, i in chromatogram ]+[0])
        graphics.polygon(x, y, lty = 'blank', col = area_colors[num])
        
if chromatograms_line:
    for num, chromatogram in enumerate(chromatograms_line):
        x = robjects.FloatVector([rt for rt, i in chromatogram])
        y = robjects.FloatVector([i for rt, i in chromatogram])
        graphics.points(x, y, type = 'l', lwd = 1, col = line_colors[num])

if chromatograms_points:
    for num, chromatogram in enumerate(chromatograms_points):
        x = robjects.FloatVector([rt for rt, i in chromatogram])
        y = robjects.FloatVector([i for rt, i in chromatogram])
        graphics.points(x, y, type = 'p')
        
if labels:
    graphics.legend('topright',
                    legend = r.c(labels),
                    lwd = 1,
                    col = robjects.StrVector(line_colors),
                    bty = 'n'
                    )
    
if texts:
    for time, intensity, text in texts:
        graphics.text( time, intensity, text )
        
        
def plot_spectrum(spectrum,
                cos = None,
                title = None,
                plot_type = 'h',
                lwd = 0.5,
                additional = None):
try:
    mz = spectrum.mz
    i = spectrum.i
except:
    mz = [_mz for _mz, _i in spectrum]
    i = [_i for _mz, _i in spectrum]
if title == None:
    title = 'Spectrum'
yax_max = max(i) * 1.05
log_yax = int(math.log(yax_max, 10))
graphics.plot(robjects.FloatVector(mz),
                robjects.FloatVector(i),
                main = title,
                xlab = 'm/z', 
                ylab = 'Intensity [arb.]',
                type = plot_type, yaxs = 'i', xaxs = 'i', lwd = lwd,
                ylim = r.c(0, yax_max),
            xlim = r.c(min(mz), cos.mz+10) if cos != None else r.c(50,1000),
#                   xlim = r.c(cos.mz*0.25, cos.mz*1.05),
                axes=False,
                las = 1
                )
graphics.axis(1, at=r.c(list(range(0,1500,1))), 
                labels=False, tcl=-0.25) #x-axis
graphics.axis(1, at=r.c(list(range(0,1500,1)))) #x-axis

yax_marks_at = list(graphics.axTicks(2))
yax_marks_at_minordiff = int(yax_marks_at[1]/2) 
yax_marks_at_minor = list(range(0, 
                                int(yax_marks_at[-1]+yax_marks_at[1]),
                                yax_marks_at_minordiff
                                ))
graphics.axis(2, las=1, at=r.c(yax_marks_at_minor),
                labels=False, tcl=-0.25)
graphics.axis(2, las=1, at=r.c(yax_marks_at),
                labels=r.c(['{0:.1f}'.format(v/(10**log_yax)) for v in yax_marks_at])
                ) #y-axis
graphics.box()
graphics.mtext(r.parse(text=r.paste("x10", "^{0}".format(log_yax), sep="")),
#                    'Intens.\nx10{0}'.format(log_yax), 
                side = 2, adj = 1, line = 3)#, outer = True)

if cos:
    graphics.points(cos.mz, 0, pch = 1, cex = 1, xpd = True, 
                    col = COLORS['red'])
if additional:
    for x, y in additional:
        graphics.points(x, y, pch = 2, cex = 1, xpd = True, 
                    col = COLORS['blue'])
    
    
