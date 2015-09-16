from gwpy.timeseries import TimeSeries
import os
from scipy.signal import detrend
from glue import datafind
import shutil
import sys

"""
code for creating spectrograms from triggers produced by
omicron_histogram.py

takes 4 arguments

Parameters:
    channel : str
        channel you want to create spectrograms for (probably,
        but not necessarily, channel whose triggers you're plotting)
    file : str
        file whose triggers you want to plot
    dir : str
        base directory inside of which specgrams/%st-%et the
        plots will be saved
    outfile : str
        output zip file that has all of the plots in it.
"""

channel = sys.argv[1]
file = sys.argv[2]
dir = sys.argv[3]
outfile = sys.argv[4]
t = 5
dir2 = dir + '/specgrams'
os.mkdir(dir2)

st = int(file.split('-')[-2])
et = int(file.split('-')[-1][:-4])
extra_dir = '%d-%d' % (st, et)
os.system('mkdir -p %s%s' % (dir2, extra_dir))
f = open(file, 'r')
times = []
snrs = []
freqs = []

for line in f:
    times.append(float(line.split('\t')[0]))
    snrs.append(float(line.split('\t')[1]))
    freqs.append(float(line.split('\t')[2]))

# for now let's just assume you're not looking at a lock
# longer than a day...
for i in range(int(86400 / t)):
    new_times = []
    new_freqs = []
    for time, snr, freq in zip(times, snrs, freqs):
        if time > st and time < st + t:
            new_times.append(time)
            new_freqs.append(freq)
    if len(new_times) == 0:
        st += t
        continue
    connection = datafind.GWDataFindHTTPConnection()
    cache = connection.find_frame_urls(
        'L', 'L1_C', st, st + t, urltype='file')
    data = TimeSeries.read(cache, channel, st, st + t)
    data2 = detrend(data)
    data = TimeSeries(
        data2, dx=data.dx, sample_rate=data.sample_rate, x0=data.x0)
    specgram = data.spectrogram2(fftlength=.1, overlap=0.1 * 0.9)
    specgram = specgram.ratio('median')
    plot = specgram.plot(vmin=1, vmax=10, norm='log')
    plot.add_colorbar(label='amplitude relative to median')
    ax = plot.gca()
    ax.set_ylim(40, 7e4)
    ax.set_yscale('log')
    # ax.scatter(new_times,new_freqs,'x',color='r')
    for time, freq in zip(new_times, new_freqs):
        ax.scatter(time, freq, marker='x', color='r', s=160)
    ax.set_title('expect %d triggers' % (len(new_times)))
    plotfile = '%s/%s/%s/%d/%d' % (dir2, extra_dir,
                                   channel.replace(':', '/'), st, et)
    plot.savefig(plotfile)
    print 'saved specgrams/%s-SPEC-%d-%d' % (channel.replace(':', '-'),
                                             st, et)
    st += t

shutil.make_archive(outfile, 'zip', root_dir=dir, base_dir=dir2)
