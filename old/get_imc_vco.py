from gwpy.timeseries import TimeSeries, TimeSeriesDict
from glue import datafind
import numpy
import scipy.signal
import sys
"""
-originally written by Andy Lundgren
-modified by Pat Meyers
-patrick.meyers@ligo.org

---------------
This code takes a start time and a duration and
calculates a "fast" vco channel by fitting the slow
VCO readback to the faster IMC-F readback on sliding
45 second scales (since the calibration between the channels
is not necessarily consistent)
----------------

takes 3 parameters:

Parameters:
-----------
    ifo : str
        interferometer (e.g. H1)
    st : int
        start time
    dur : int
        duration
"""

chan1_pat = '%s:SYS-TIMING_C_FO_A_PORT_11_SLAVE_CFC_FREQUENCY_5'
chan2_pat = '%s:IMC-F_OUT_DQ'


def calibrate_imc_pslvco(ifo, start_time, dur, cache=None):
    st, et = start_time, start_time + dur
    if cache:
        pslvco = TimeSeries.read(cache, chan1_pat % ifo, start=st, end=et)
        imc = TimeSeries.read(cache, chan2_pat % ifo, start=st, end=et)
    else:
        imc = TimeSeries.fetch(chan2_pat % ifo, st, et)
        pslvco = TimeSeries.fetch(chan1_pat % ifo, st, et)

    arr_psl = pslvco[8::16]
    arr_imc = imc

    tmp1 = (arr_imc[8192::16384])[:-1]
    tmp2 = arr_psl[1:]
    a, b = numpy.polyfit(tmp1, tmp2, 1)

    return a, b


def dump_calibrated_data(fname):
    data = numpy.load(fname)

    # Figure out the times covered by the file from the filename
    # I should start using HDF5 so I can store metadata
    temp = fname.split('.')[0]
    temp = temp.split('-')
    ifo = temp[0]
    st, dur = int(temp[-2]), int(temp[-1])
    et = st + dur

    maxidx = len(data)
    width = 45

    weights = 1. - ((numpy.arange(-width, width) / float(width))**2)

    # The VCO frequencies are integers so we could dither them
    # to avoid quantization error if we wanted to be fancy
    # but it seems to make no differece
    if False:
        from numpy.random import triangular
        data[:, 1] += triangular(-1., 0., 1., size=len(data))

    # Just fit the whole thing at once, to get a single coefficient
    a, b = numpy.polyfit(data[:, 0], data[:, 1], 1)
    print "%.1f %u" % (a, b)

    # Slide through the data fitting PSL to IMC for data around each sample
    coeffs = []
    for idx in xrange(maxidx):
        idx1 = max(0, idx - width)
        idx2 = min(idx + width, maxidx)
        coeffs.append(numpy.polyfit(data[idx1:idx2, 0], data[idx1:idx2, 1], 1,
                                    w=weights[idx1 - idx + width:idx2 - idx + width]))
    coeffs = numpy.array(coeffs)
    times = numpy.arange(len(coeffs)) + 0.5
    connection = datafind.GWDataFindHTTPConnection()
    cache = connection.find_frame_urls(
        ifo[0], '%s_R' % ifo, st, et, urltype='file')

    imc = TimeSeries.read(cache, "%s:IMC-F_OUT_DQ" % ifo, st, et)
    imc = imc[::16384 / 256]
    print imc
    samp_times = numpy.arange(len(imc)) / 256.

    coeffs0 = numpy.interp(samp_times, times, coeffs[:, 0])
    coeffs1 = numpy.interp(samp_times, times, coeffs[:, 1]) - 7.6e7

    vco_interp = coeffs0 * imc.data + coeffs1

    chan = "%s:IMC-VCO_PREDICTION" % (ifo,)
    vco_data = TimeSeries(vco_interp, epoch=st,
                          sample_rate=imc.sample_rate.value,
                          name=chan, channel=chan)
    vco_data.write("%s-vcoprediction-%u-%u.hdf" % (ifo, st, dur), format='hdf')

if __name__ == '__main__':
    ifo = sys.argv[1]
    st = int(sys.argv[2])
    dur = int(sys.argv[3])
    et = st + dur

    connection = datafind.GWDataFindHTTPConnection()
    cache = connection.find_frame_urls(
        ifo[0], '%s_R' % ifo, st, et, urltype='file')

    print "Get IMC"
    imc = TimeSeries.read(cache, chan2_pat % ifo, st, et)
    print "Get psl"
    pslvco = TimeSeries.read(cache, chan1_pat % ifo, st, et + 1)
    print "Downsample psl"
    pslvco = pslvco[16 + 8::16]
    print "Downsample imc"
    imc_srate = int(imc.sample_rate.value)
    imc = imc[imc_srate / 2::imc_srate]

    print "Saving"
    data = numpy.array((imc.data, pslvco.data)).T
    fname = "%s-imc-vco-%u-%u.npy" % (ifo, st, dur)
    numpy.save("%s-imc-vco-%u-%u.npy" % (ifo, st, dur), data)
    dump_calibrated_data(fname)
