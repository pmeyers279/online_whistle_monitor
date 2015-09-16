from gwpy.segments import DataQualityFlag
from gwpy.table.lsctables import SnglBurstTable
from glue.ligolw import table
import numpy as np
import re


def create_name(channel, start, end, tag, format='png'):
    """Build T050017-compatible plot name
    """
    c = re.sub('[:_-]', '_', channel).replace('_', '-', 1).upper()
    t = re.sub('[:_\-\s]', '_', tag).upper()
    s = int(start)
    d = int(round(end - start))
    return '%s_%s-%d-%d.%s' % (c, tag, s, d, format)


def get_segments(ifo, st, et, flag='DMT-ANALYSIS_READY:1'):
    flag = '%s:flag' % ifo
    segments = DataQualityFlag.dqsegdb(
        flag, st, et, url='https://segments.ligo.org')
    return segments


def write_segments(ifo, st, et, flag='DMT-ANALYSIS_READY:1',
                   filename=None, directory=None):
    segments = get_segments(ifo, st, et, flag=flag)
    flag_name = '%s:flag' % ifo
    if filename is None:
        filename = create_name(flag_name, st, et, 'SEGS', format='xml.gz')
    if directory is None:
        segments.write(filename)
    else:
        segments.write('%s/%s' % (directory, filename))


def read_omicron_trigs(omicron_files, span):
    trigs1 = SnglBurstTable.read(omicron_files[0], format='ligolw')
    trigs = table.new_from_template(trigs1)
    for omicron_file in omicron_files:
        trigs2 = table.new_from_template(trigs1)
        trigs2 = SnglBurstTable.read(omicron_file, format='ligolw')
        trigs.extend(filter(lambda row: row.get_peak() in span, trigs2))
    times = []
    print('%d found' % len(trigs))
    times = trigs.get_peak()
    return trigs, times


def get_vco_data(vco_file, times):
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


def hist_triggers(trigs, amp, times, nbins=100, bin_low=1,
                  bin_high=1, lines=None):
    bin_size = (amp.max() - amp.min()) / float(nbins)
    bins = numpy.arange(amp.min(), amp.max() + bin_size, bin_size)
    hist = numpy.zeros(bins.size)

    idxs = []
    if lines is not None:
        for line in lines:
            idxs.append(np.where(line - 15 < bins < line + 15)[0])
        idxs = np.asarray(idxs)
    else:
        idxs = numpy.where(bins < 1e10)[0]

    idx = numpy.where(hist[idxs] == hist[idxs].max())[0]
    if len(idx) > 1:
        idx = idx[0]

    high = 79 + bins[idx + bin_high] / 1e3
    low = 79 + bins[idx - bin_low] / 1e3
    print('we expect whistles occur between %s and %s'
          % (bins[idx - bin_low], bins[idx + bin_high]))

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


def plot_trig_rates(trigs, spans, channel, files=None):
    durations = {}
    for span in spans.keys():
        durations[span] = abs(spans[span])

    First = 1
    bins = np.logspace(0, 2, num=200)
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
    return plot
    # plot.save(channel + 'LOCK-TRIG-RATES')


def plot_vco_hist(amps, spans, channel):
    for amp, span in zip(amps.keys(), spans.keys()):
        start = spans[span][0]
        end = spans[span][-1]
        png = plot_name(channel, start, end, 'VCO-HIST', format='png')
        plot = HistogramPlot(
            amps[amp], weights=1 / abs(spans[span]), bins=100, log=True)
        plot.set_xlabel('Amplitude [kHz from 79 MHz]')
        plot.set_ylabel('Rate [Hz]')
        plot.set_title(channel.replace('_', '\_') + ' omicron glitch histogram')
        plot.suptitle(str(spans[span][0]) + '-' + str(spans[span][-1]))
        ax = plot.gca()
        plot.save(png)
        plot.close()
        print('%s written' % png)
