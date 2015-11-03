from __future__ import (division, print_function)
import re
import glob

import numpy as np

from gwpy.table.lsctables import SnglBurstTable
from gwpy.segments import Segment
from gwpy.plotter import HistogramPlot
from glue.ligolw import table
from vco_functions import *
import matplotlib.pyplot as plt
import markup


def plot_name(channel, start, end, tag, format='png'):
    """Build T050017-compatible plot name
    """
    c = re.sub('[:_-]', '_', channel).replace('_', '-', 1).upper()
    s = int(start)
    d = int(round(end - start))
    return '%s_%s-%d-%d.%s' % (c, tag, s, d, format)


def get_all_vco_triggers(ifo, seg, frames=False, fit=True,
                         channel='GDS-CALIB_STRAIN'):
    st = int(seg[0])
    et = int(seg[1])
    st2 = st
    First = 1
    # break into 1000s, run into memory issues otherwise
    amps = []
    while (et > st2):
	dur = min(et - st2, 1000)
        vco = generate_fast_vco(ifo, Segment(st2, st2+dur), frames=frames, fit=fit)
        trigs = SnglBurstTable.fetch("%s:%s" % (ifo, channel), 'omicron', st2, st2+dur,
				 filt=lambda x: st2 <= x.get_peak() < st2+dur)
        if First:
	    trigs2 = table.new_from_template(trigs)
	    First=0
	trigs2.extend(trigs)
	vtrigs = get_vco_trigs(vco, trigs, channel=channel)
	for vtrig in vtrigs:
	    amps.append(vtrig)
	st2 += 1000

    central_freqs = trigs2.get_column('peak_frequency')
    snrs = trigs2.get_column('snr')
    
    amps = np.asarray(amps)
    amps = amps[~(np.isnan(amps))]
    central_freqs = central_freqs[~(np.isnan(amps))]
    snrs = snrs[~(np.isnan(amps))]
    return amps, central_freqs, snrs


def plot_vco_hist(vco_trigs, segment, channel):
    start = segment[0]
    end = segment[1]
    png = plot_name(channel, start, end, 'VCO-HIST', format='png')
    plot = HistogramPlot(
        vco_trigs, weights=1. / abs(int(end) - int(start)), bins=100, log=True)
    plot.set_xlabel('VCO Frequency [kHz from 79 MHz]')
    plot.set_ylabel('Rate [Hz]')
    plot.set_title(channel.replace('_', '\_') + ' omicron glitch histogram', fontsize=14)
    plot.suptitle(str(segment[0]) + '-' + str(segment[-1]))
    plot.savefig(png)
    plot.close()
    print('%s written' % png)

def plot_omicron_scatter(vco_trigs, central_freqs, snrs, segment, channel):
    snrs[np.where(snrs>20)] = 20
    start = segment[0]
    end = segment[1]
    png = plot_name(channel, start, end, 'OMICRON-SCATTER', format='png')
    plt.scatter(vco_trigs, central_freqs, c=snrs)
    plt.xlabel('VCO Frequency [kHz from 79 MHz]', fontsize=12)
    plt.ylabel('Omicron central frequency [Hz]', fontsize=12)
    plt.title('%s Omicron Scatter' % channel.replace('_','\_'), fontsize=14)
    plt.suptitle('%d-%d' % (start, end))
    ax = plt.gca()
    ax.set_ylim(0,8192)
    plt.savefig(png)


def get_vco_trigs(vco, trigs, channel='GDS-CALIB_STRAIN'):
    #print('Finding mean values of %s' % channel)
    times = trigs.get_peak()
    amp = np.zeros(len(trigs))
    for i, t in enumerate(times):
        t = t.seconds + t.nanoseconds / 1e9
        idx1 = int(
            (t - vco.times.value[0] - 0.1) * vco.sample_rate.value)
        idx2 = int(
            (t - vco.times.value[0] + 0.1) * vco.sample_rate.value)
        temp = vco[idx1:idx2]
        amp[i] = ((temp.mean().value) - 79000000) / 1.e3
    return amp
