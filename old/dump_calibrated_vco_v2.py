import sys
from numpy import *
from gwpy.timeseries import TimeSeries
from glue import datafind

# originally written by Andy Lundgren, modified by Pat Meyers
# patrick.meyers@ligo.org


fname = sys.argv[1]
data = load(fname)

# Figure out the times covered by the file from the filename
# I should start using HDF5 so I can store metadata
temp = fname.split('.')[0]
temp = temp.split('-')
ifo = temp[0]
st, dur = int(temp[-2]), int(temp[-1])
et = st + dur

maxidx = len(data)
width = 45

weights = 1. - ((arange(-width, width) / float(width))**2)

# The VCO frequencies are integers so we could dither them
# to avoid quantization error if we wanted to be fancy
# but it seems to make no differece
if False:
    from numpy.random import triangular
    data[:, 1] += triangular(-1., 0., 1., size=len(data))

# Just fit the whole thing at once, to get a single coefficient
a, b = polyfit(data[:, 0], data[:, 1], 1)
print "%.1f %u" % (a, b)

# Slide through the data fitting PSL to IMC for data around each sample
coeffs = []
for idx in xrange(maxidx):
    idx1 = max(0, idx - width)
    idx2 = min(idx + width, maxidx)
    coeffs.append(polyfit(data[idx1:idx2, 0], data[idx1:idx2, 1], 1,
                          w=weights[idx1 - idx + width:idx2 - idx + width]))
coeffs = array(coeffs)
times = arange(len(coeffs)) + 0.5
connection = datafind.GWDataFindHTTPConnection()
cache = connection.find_frame_urls(
    ifo[0], '%s_R' % ifo, st, et, urltype='file')

imc = TimeSeries.read(cache, "%s:IMC-F_OUT_DQ" % ifo, st, et)
imc = imc[::16384 / 256]
print imc
samp_times = arange(len(imc)) / 256.

coeffs0 = interp(samp_times, times, coeffs[:, 0])
coeffs1 = interp(samp_times, times, coeffs[:, 1]) - 7.6e7

vco_interp = coeffs0 * imc.data + coeffs1

chan = "%s:IMC-VCO_PREDICTION" % (ifo,)
vco_data = TimeSeries(vco_interp, epoch=st, sample_rate=imc.sample_rate.value,
                      name=chan, channel=chan)
vco_data.write("%s-vcoprediction-%u-%u.hdf" % (ifo, st, dur), format='hdf')
