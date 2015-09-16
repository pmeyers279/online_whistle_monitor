import sys
from numpy import *
from gwpy.timeseries import TimeSeries
from scipy.interpolate import interp1d

fname = sys.argv[1]
data = load(fname)

# Figure out the times covered by the file from the filename
# I should start using HDF5 so I can store metadata
temp = fname.split('.')[0]
temp = temp.split('-')
ifo = temp[0]
st, dur = int(temp[-2]), int(temp[-1])
et = st + dur

offset = 7.9e7
fvco = interp1d(arange(len(data)) + 0.5, data[:, 1] - offset, kind='cubic')

chan = "%s:IMC-VCO_INTERPOLATION" % (ifo,)
vco_data = TimeSeries(vco_interp, epoch=st, sample_rate=imc.sample_rate.value,
                      name=chan, channel=chan)
vco_data.write("%s-vcointerpolation-%u-%u.hdf" % (ifo, st, dur), format='hdf')
