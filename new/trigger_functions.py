import re
import glob

import numpy as np
from __future__ import (division, print_function)

from gwpy.table.lsctables import SnglBurstTable
from gwpy.plotter import HistogramPlot
from glue.ligolw import table
from vco_functions import *


def plot_name(channel, start, end, tag, format='png'):
    """Build T050017-compatible plot name
    """
    c = re.sub('[:_-]', '_', channel).replace('_', '-', 1).upper()
    s = int(start)
    d = int(round(end - start))
    return '%s_%s-%d-%d.%s' % (c, tag, s, d, format)


def get_omicron_files(ifo, segment, trig_dir=None, channel='GDS-CALIB_STRAIN'):
    if trig_dir is not None:
        trig_dir = '/home/detchar/triggers/O1/%s/%s' % (ifo, channel)
    st = segment[0]
    et = segment[1]
    dir1 = str(st)[0:5]
    dirend = str(et)[0:5]
    dirs = np.arange(int(dir1), int(dirend) + 1)
    files = []
    for directory in dirs:
        files_to_add = sorted(glob.glob('%s/%d/*.xml' % (trig_dir, directory)))
        files.append(files_to_add)
    return files


def plot_vco_hist(vco_trigs, segment, channel):
    start = segment[0]
    end = segment[1]
    png = plot_name(channel, start, end, 'VCO-HIST', format='png')
    plot = HistogramPlot(
        amps[amp], weights=1 / abs(segment), bins=100, log=True)
    plot.set_xlabel('VCO Frequency [kHz from 79 MHz]')
    plot.set_ylabel('Rate [Hz]')
    plot.set_title(channel.replace('_', '\_') + ' omicron glitch histogram')
    plot.suptitle(str(spans[span][0]) + '-' + str(spans[span][-1]))
    plot.savefig(png)
    plot.close()
    print('%s written' % png)


def get_vco_trigs(vco, times):
    print('Finding mean values of %s' % channel)
    s = times.size
    amp = numpy.zeros(s)
    for i, t in enumerate(times):
        t = t.seconds + t.nanoseconds / 1e9
        idx1 = int(
            (t - vco.times.value[0] - 0.01) * vco.sample_rate.value)
        idx2 = int(
            (t - vco.times.value[0] + 0.01) * vco.sample_rate.value)
        temp = vco[idx1:idx2]
        amp[i] = (temp.mean().value) - 79e6 / 1.e3
        print('    Processed trigger %d/%d' % (i + 1, s), end='\r')
    print('    Processed trigger %d/%d' % (s, s))
    return amp


def read_omicron_trigs(omicron_files, segment):
    """
    """
    trigs1 = SnglBurstTable.read(omicron_files[0], format='ligolw')
    trigs = table.new_from_template(trigs1)
    for omicron_file in omicron_files:
        trigs2 = table.new_from_template(trigs1)
        trigs2 = SnglBurstTable.read(omicron_file, format='ligolw')
        trigs.extend(filter(lambda row: row.get_peak() in segment, trigs2))
    times = []
    print('%d found' % len(trigs))
    times = trigs.get_peak()
    s = times.size
    return trigs, times, s
