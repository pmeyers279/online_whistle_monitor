from __future__ import (division, print_function)
import re
import glob

import numpy as np

from gwpy.table.lsctables import SnglBurstTable
from gwpy.segments import Segment
from gwpy.plotter import HistogramPlot
from glue.ligolw import table
from vco_functions import *
import markup


def plot_name(channel, start, end, tag, format='png'):
    """Build T050017-compatible plot name
    """
    c = re.sub('[:_-]', '_', channel).replace('_', '-', 1).upper()
    s = int(start)
    d = int(round(end - start))
    return '%s_%s-%d-%d.%s' % (c, tag, s, d, format)


def get_omicron_files(ifo, segment, trig_dir=None, channel='GDS-CALIB_STRAIN'):
    if trig_dir is not None:
        trig_dir = '/home/detchar/triggers/O1/%s/%s_Omicron' % (ifo, channel)
    st = segment[0]
    et = segment[1]
    dir1 = str(st)[0:5]
    dirend = str(et)[0:5]
    dirs = np.arange(int(dir1), int(dirend) + 1)
    files = []
    for directory in dirs:
        #print('looking for trigs in %s/%d' % (trig_dir, directory))
        files_to_add = sorted(
            glob.glob('%s/%d/*.xml.gz' % (trig_dir, directory)))
	for f in files_to_add:
            files.append(f)
    return files


def get_all_vco_triggers(omicron_files,seg,frames=False, fit=True,
                         channel='GDS-CALIB_STRAIN'):
    amps = []
    for f in omicron_files:
	f2 = f.split('/')[-1]
        # parse file name
        ifo = f2.split('-')[0]
        st = int(f2.split('-')[2])
        dur = int(f2.split('-')[3].split('.')[0])
        et = st + dur
	if st > int(seg[0]) and et < int(seg[1]):
            # generate fast vco
            vco = generate_fast_vco(ifo, Segment(st, et), frames=0, fit=fit)
            trigs, t, s = read_omicron_trigs(f, Segment(st, et))
            vtrigs = get_vco_trigs(vco, t, channel=channel)
            for vtrig in vtrigs:
                amps.append(vtrig)
    return amps


def plot_vco_hist(vco_trigs, segment, channel):
    start = segment[0]
    end = segment[1]
    png = plot_name(channel, start, end, 'VCO-HIST', format='png')
    plot = HistogramPlot(
        vco_trigs, weights=1. / abs(int(end) - int(start)), bins=100, log=True)
    plot.set_xlabel('VCO Frequency [kHz from 79 MHz]')
    plot.set_ylabel('Rate [Hz]')
    plot.set_title(channel.replace('_', '\_') + ' omicron glitch histogram')
    plot.suptitle(str(segment[0]) + '-' + str(segment[-1]))
    plot.savefig(png)
    plot.close()
    #print('%s written' % png)


def get_vco_trigs(vco, times, channel='GDS-CALIB_STRAIN'):
    #print('Finding mean values of %s' % channel)
    s = times.size
    amp = np.zeros(s)
    for i, t in enumerate(times):
        t = t.seconds + t.nanoseconds / 1e9
        idx1 = int(
            (t - vco.times.value[0] - 0.01) * vco.sample_rate.value)
        idx2 = int(
            (t - vco.times.value[0] + 0.01) * vco.sample_rate.value)
        temp = vco[idx1:idx2]
        amp[i] = ((temp.mean().value) - 79000000) / 1.e3
    #    print('    Processed trigger %d/%d' % (i + 1, s), end='\r')
    #print('    Processed trigger %d/%d' % (s, s))
    return amp


def read_omicron_trigs(omicron_files, segment):
    """
    """
    if isinstance(omicron_files, list):
	##print('HI THERE PAT. FUCK YOU!!!')
    	trigs1 = SnglBurstTable.read(omicron_files[0], format='ligolw')
    else:
	# i'm being lazy so i don't have to rewrite code below
	omicron_files = [omicron_files]
    	trigs1 = SnglBurstTable.read(omicron_files, format='ligolw')
    trigs = table.new_from_template(trigs1)
    for omicron_file in omicron_files:
        trigs2 = table.new_from_template(trigs1)
        trigs2 = SnglBurstTable.read(omicron_file, format='ligolw')
        trigs.extend(filter(lambda row: row.get_peak() in segment, trigs2))
    times = []
    #print('%d found' % len(trigs))
    times = trigs.get_peak()
    s = times.size
    return trigs, times, s


def make_html_page(segs, channel='L1:GDS-CALIB_STRAIN'):
    channel = channel.replace(':', '-')
    page = markup.page()
    page.init(title="VCO Histogram to Search for Whistles in %s" % channel)
    page.div()
    for seg in segs:
        page.div()
        page.img(src="%s_VCO_HIST-%d-%d.png" %
                 (channel, int(seg[0]), int(seg[1])))
        page.p("VCO Histogram")
        page.div.close()
    page.div.close()
