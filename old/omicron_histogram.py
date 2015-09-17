#!/usr/bin/env python

"""Read Omicron triggers, and plot
originally written by Laura Nuttall 
modified by Pat Meyers 
patrick.meyers@ligo.org

---------
This code takes a channel and a VCO file (produced by get_imc_vco.py)
and finds all of the omicron triggers during that time and histograms them
against the VCO value to see if there are specific VCO values that seem to
correlate to more omicron glitches than others (these are likely whistles)
NOTE: currently you must change the "channel" and "vco_files" variables
in the code itself
---------
"""

from __future__ import (division, print_function)

import re
import glob

import numpy

from gwpy.segments import Segment
import matplotlib.pyplot as plt
from gwpy.table.lsctables import SnglBurstTable
from gwpy.time import to_gps
from gwpy.timeseries import TimeSeries
from gwpy.plotter import HistogramPlot
from gwpy.plotter.tex import label_to_latex
from gwpy.io import nds as ndsio
from glue.ligolw import table

# ----------- INPUTS -----------#
vco_files = [
  '/home/meyers/aDQ/Whistles/get_imc_vco/L1-vcoprediction-1126256591-8100.hdf']

channel = 'L1:GDS-CALIB_STRAIN'
# ------------------------------#


def plot_name(channel, start, end, tag, format='png'):
    """Build T050017-compatible plot name
    """
    c = re.sub('[:_-]', '_', channel).replace('_', '-', 1).upper()
    t = re.sub('[:_\-\s]', '_', tag).upper()
    s = int(start)
    d = int(round(end - start))
    return '%s_%s-%d-%d.%s' % (c, tag, s, d, format)


def _get_vco_data(vco_file, times):
    print('Finding mean values of %s' % channel)
    s = times.size
    amp = numpy.zeros(s)
    vco = TimeSeries.from_hdf5(vco_file)
    for i, t in enumerate(times):
        t = t.seconds + t.nanoseconds / 1e9
        idx1 = int(
            (t - vco.times.value[0] - 0.01) * vco.sample_rate.value)
        idx2 = int(
            (t - vco.times.value[0] + 0.01) * vco.sample_rate.value)
        temp = vco[idx1:idx2]
        amp[i] = (temp.mean().value - 3e6) / 1.e3
        print('    Processed trigger %d/%d' % (i + 1, s), end='\r')
    print('    Processed trigger %d/%d' % (s, s))
    return amp


def _read_omicron_trigs(omicron_files, span):
    trigs1 = SnglBurstTable.read(omicron_files[0], format='ligolw')
    trigs = table.new_from_template(trigs1)
    for omicron_file in omicron_files:
        trigs2 = table.new_from_template(trigs1)
        trigs2 = SnglBurstTable.read(omicron_file, format='ligolw')
        trigs.extend(filter(lambda row: row.get_peak() in span, trigs2))
    times = []

    print('%d found' % len(trigs))
    times = trigs.get_peak()
    s = times.size
    return trigs, times, s


def plot_trig_rates(trigs, spans, channel, files=None, files_to_combine=None):
    if files is not None and files_to_combine is not None:
        final_trigs_to_plot = {}
        durations = {'combined': 0}
        combined_trigs = table.new_from_template(trigs[files[0]])
        for f in files_to_combine:
            combined_trigs.extend(trigs[f])
        for f in trigs.keys():
            if f not in files_to_combine:
                final_trigs_to_plot[f] = trigs[f]
                durations[f] = abs(spans[f])
            else:
                durations['combined'] += abs(spans[f])
        final_trigs_to_plot['combined'] = combined_trigs
        trigs = final_trigs_to_plot
    else:
        durations = {}
        for span in spans.keys():
            durations[span] = abs(spans[span])

    First = 1
    bins = numpy.logspace(0, 2, num=200)
    for key in durations.keys():
        # print(len(trigs[key]))
        try:
            label = str(spans[key][0]) + '-' + str(spans[key][-1])
        except KeyError:
            label = 'combined locks'
        if First:
            plot = trigs[key].hist('snr', logbins=True,
                                   log=True,
                                   histtype='step',
                                   weights=1. / abs(durations[key]),
                                   label=label,
                                   bins=bins,
                                   cumulative=-1
                                   )
            ax = plot.gca()
            First = 0
            continue
        snrs = trigs[key].get_column('snr')
        ax.hist(snrs, logbins=True, log=True, histtype='step',
                weights=1. / abs(durations[key]), label=label,
                bins=bins, cumulative=-1)
    ax.set_yscale('log')
    ax.set_ylabel('Trigger Rate [Hz]')
    ax.set_title(channel.replace('_', '\_'))
    ax.legend()
    ax.set_ylim(1e-4, 1)
    ax.set_xscale('log')
    ax.set_xlim(5, 100)
    plot.save(channel + 'LOCK-TRIG-RATES')


def _hist_triggers(trigs, amp, times, nbins=100, bin_low=1, bin_high=1):
    bin_size = (amp.max() - amp.min()) / float(nbins)
    bins = numpy.arange(amp.min(), amp.max() + bin_size, bin_size)
    hist = numpy.zeros(bins.size)
    hist_indices = {}
    for i in range(bins.size - 1):
        hist_indices[str(i)] = []
        for j in range(times.size):
            if amp[j] > bins[i] and amp[j] < bins[i + 1]:
                hist[i] += 1
                hist_indices[str(i)].append(j)
    idx = numpy.where(hist == hist.max())[0]
    if len(idx) > 1:
        idx = idx[0]
    high = 79 + bins[idx + bin_high] / 1e3
    low = 79 + bins[idx - bin_low] / 1e3
    print('we expect whistles occur between %s and %s MHz'
          % (low, high))

    new_trigs = table.new_from_template(trigs)
    for i in range(len(trigs)):
        if amp[i] > bins[idx - bin_low] and amp[i] < bins[idx + bin_high]:
            new_trigs.extend(trigs[i:i + 1])
    times = []
    snrs = []
    freqs = []
    for trig in new_trigs:
        peak = trig.get_peak()
    snrs.append(trig.snr)
    freqs.append(trig.peak_frequency)
    times.append(peak.seconds + peak.nanoseconds * 1e-9)
    f = open(
        'L1-TRIG-TIMES-' + str(span[0]) + '-' + str(span[-1]) + '.txt', 'w')
    for time, snr, freq in zip(times, snrs, freqs):
        f.write('%f\t%f\t%f\n' % (time, snr, freq))
    f.close()
    return new_trigs, hist


def plot_vco_hist(amps, spans, channel):
    for amp, span in zip(amps.keys(), spans.keys()):
        start = spans[span][0]
        end = spans[span][-1]
        png = plot_name(channel, start, end, 'VCO-HIST', format='png')
        plot = HistogramPlot(amps[amp], weights=1 / abs(spans[span]),
                             bins=100, log=True)
        plot.set_xlabel('Amplitude [kHz from 79 MHz]')
        plot.set_ylabel('Rate [Hz]')
    plot.set_title(channel.replace('_', '\_') +
                   ' omicron glitch histogram')
    plot.suptitle(str(spans[span][0]) + '-' + str(spans[span][-1]))
    ax = plot.gca()
    plot.save(png)
    plot.close()
    print('%s written' % png)
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------

# get Omicron triggers
print('Reading DARM triggers', end=' ')
spans = {}


ifo = channel[0:2]
# omicron channel directory
chan_dir = channel[3:] + '_Omicron'

channel_texname = channel.replace('_', '\_')
xml_dir = '/home/detchar/triggers/O1/%s/%s/*/' % (ifo, chan_dir)
files = sorted(glob.glob('%s/*.xml*' % xml_dir))
for v in vco_files:
    start = int(v.split('-')[-2])
    end = start + int(v.split('-')[-1][:-4])
    spans[v] = Segment(start, end)

trigs = {}
times = {}
s = {}
amp = {}
hist = {}
new_trigs = {}

counter = 0
for v in vco_files:
    start = int(v.split('-')[-2])
    end = start + int(v.split('-')[-1][:-4])
    span = Segment(start, end)
    counter += 1
    trigs[v], times[v], s[v] = _read_omicron_trigs(files, span)
    amp[v] = _get_vco_data(v, times[v])
    new_trigs[v], hist = _hist_triggers(trigs[v], amp[v], times[v], nbins=100,
                                        bin_low=2, bin_high=2)


plot_trig_rates(trigs, spans, channel + '-ALL', files=None)
plot_vco_hist(amp, spans, channel)
plot_trig_rates(new_trigs, spans, channel + '-WHISTLE', files=None)

# -----------------------------------------------------------------------------


# new_trigs, hist = _hist_triggers(trigs, amp, nbins=100)

rate = trigs[vco_files[0]].event_rate(300, 'snr', [5, 6, 8], operator='>=')
plot = rate.plot(label='name', marker='o', linestyle='')
ax = plot.gca()
ax.set_yscale('log')
ax.set_ylabel('Trigger Rate [Hz]')
ax.set_title('%s - %s' % ('5 min Rate Trend', channel_texname))
ax.legend(loc='upper left')
plot.savefig(channel + '-EVENT-RATE-' + str(start) + '-' + str(end - start))
plot.close()
