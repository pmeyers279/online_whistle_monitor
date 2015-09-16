import numpy as np
from gwpy.plotter import HistogramPlot
from glue.ligolw import table


class Triggers(object):

    """docstring for Triggers"""

    def __init__(self, trigs=None, amps=None,
                 files=None, spans=None, channel=None, lines=None):
        super(Triggers, self).__init__()
        if trigs:
            self.trigs = trigs
        else:
            self.trigs = {}
        if amps:
            self.amps = amps
        else:
            self.amps = {}
        if self.lines:
            self.lines = lines
        else:
            self.lines = []
        if self.spans:
            self.spans = spans
        else:
            self.spans = {}
        if channel is not None:
            self.channel = channel
        else:
            self.channel = 'Unknown Channel'

    def extend(self, trigs, amps, span, file):
        self.trigs[file] = trigs
        self.amps[file] = amps
        self.files.append(file)
        self.spans[file] = span

    def get_files(self):
        return self.trigs.keys()

    def plot_vco_hist_all(self):
        # get total duration
        dur = 0
        amps = []
        for file in self.get_files():
            dur += abs(self.spans[file])
            amps.append(self.amps[file])
        plot = HistogramPlot(amps, weights=1 / dur,
                             bins=len(self.get_files()) * 100,
                             log=True)
        plot.set_xlabel('VCO frequency at time of glitch [kHz from 79 MHz]')
        plot.set_ylabel('Rate [Hz]')
        plot.set_title(
            self.channel.replace('_', '\_') + 'omicron glitch histogram')
        return plot

    def plot_vco_hist(self, file):
        plot = HistogramPlot(self.amps[file],
                             weights=1 / abs(self.spans[file]),
                             bins=100, log=True)
        plot.set_xlabel('VCO frequency at time of glitch [kHz from 79 MHz]')
        plot.set_ylabel('Rate [Hz]')
        plot.set_title(
            self.channel.replace('_', '\_') + 'omicron glitch histogram')
        return plot

    def plot_trig_rates(self, bins=None):
        if bins is None:
            bins = np.logspace(0, 2, 200)
        First = True
        for file in self.get_files():
            label = '%d-%d' % (self.spans[file][0], self.spans[file][-1])
            if First:
                First = False
                plot = self.trigs[file].hist('snr', logbins=True,
                                             log=True,
                                             histtype='step',
                                             weights=1. /
                                             abs(self.spans[file]),
                                             bins=bins,
                                             label=label,
                                             cumulative=-1)
                ax = plot.gca()
                continue
            snrs = self.trigs[file].get_column('snr')
            ax.hist(snrs, logbins=True, log=True, histtype='step',
                    weights=1. / abs(self.spans[file]),
                    label=label, bins=bins, cumulative=-1)
        ax.set_yscale('log')
        ax.set_ylabel('Trigger Rate [Hz]')
        ax.set_title(self.channel.replace('_', '\_'))
        return plot

    def hist_triggers(self, file, lines=None, nbins=100,
                      bin_low=1, bin_high=1):
        bin_size = (self.amps[file].max() - self.amps[file].min()
                    / float(nbins))
        bins = np.arange(self.amps[file].max(), self.ams[file].max(),
                         bin_size)
        hist = np.zeros(bins.size)
        idxs = []
        idx = []

        for i, bin in enumerate(bins):
            for amp in self.amps[file]:
                if bin < amp < bin + bin_size:
                    hist[i] += 1

        for line in lines:
            idxs.append(np.where(line - 15 < bins < line + 15)[0])
            idx.append(np.where(hist[idxs] == hist[idxs].max())[0])

        new_trigs = Triggers(channel=self.channel)
        trigs = table.new_from_template(self.trigs[file])
        amps = []

        for i in range(len(self.trigs[file])):
            for ids in idx:
                if self.amps[file][i] > bins[ids - bin_low]\
                        and self.amps[file][i] < bins[ids + bin_high]:
                    amps.append(self.amps[file][i])
                    trigs.extend(self.trigs[file][i:i + 1])
        new_trigs.extend(trigs, amps, self.spans[file], file)
        return new_trigs
