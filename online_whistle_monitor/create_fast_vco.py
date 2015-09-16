from gwpy.timeseries import TimeSeries
from glue import datafind
import numpy as np

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
    a, b = np.polyfit(tmp1, tmp2, 1)

    return a, b


def write_fast_vco(ifo, st, dur, directory=None, cache=None):
    """
    Parameters:
    -----------
        ifo : start
            interferometer, e.g. 'L1'
        st : int
            start time
        dur : int
            duration

    Returns:
    --------
        NONE : saves file 'L1:IMC-VCO_PREDICTION-st-dur.hdf'
    """
    et = st + dur

    if cache:
        connection = datafind.GWDataFindHTTPConnection()
        cache = connection.find_frame_urls(
            ifo[0], '%s_R' % ifo, st, et, urltype='file')
        print "Get IMC"
        imc = TimeSeries.read(cache, chan2_pat % ifo, st, et)
        print "Get psl"
        pslvco = TimeSeries.read(cache, chan1_pat % ifo, st, et + 1)
    else:
        print "Get IMC"
        imc = TimeSeries.fetch(chan2_pat % ifo, st, et)
        print "Get psl"
        pslvco = TimeSeries.fetch(chan1_pat % ifo, st, et + 1)

    print "Downsample psl"
    pslvco = pslvco[16 + 8::16]
    print "Downsample imc"
    imc_srate = int(imc.sample_rate.value)
    imc = imc[imc_srate / 2::imc_srate]

    print "Saving"
    data = np.array((imc.data, pslvco.data)).T

    maxidx = len(data)
    width = 45

    weights = 1. - ((np.arange(-width, width) / float(width))**2)

    # The VCO frequencies are integers so we could dither them
    # to avoid quantization error if we wanted to be fancy
    # but it seems to make no differece
    # Just fit the whole thing at once, to get a single coefficient
    a, b = np.polyfit(data[:, 0], data[:, 1], 1)
    print "%.1f %u" % (a, b)

    # Slide through the data fitting PSL to IMC for data around each sample
    coeffs = []
    for idx in xrange(maxidx):
        idx1 = max(0, idx - width)
        idx2 = min(idx + width, maxidx)
        coeffs.append(np.polyfit(data[idx1:idx2, 0], data[idx1:idx2, 1], 1,
                                 w=weights[
                                     idx1 - idx + width:idx2 - idx + width]
                                 ))
    coeffs = np.array(coeffs)
    times = np.arange(len(coeffs)) + 0.5

    imc = imc[::16384 / 256]
    print imc
    samp_times = np.arange(len(imc)) / 256.

    coeffs0 = np.interp(samp_times, times, coeffs[:, 0])
    coeffs1 = np.interp(samp_times, times, coeffs[:, 1]) - 7.6e7

    vco_interp = coeffs0 * imc.data + coeffs1

    chan = "%s:IMC-VCO_PREDICTION" % (ifo,)
    vco_data = TimeSeries(vco_interp, epoch=st,
                          sample_rate=imc.sample_rate.value,
                          name=chan, channel=chan)
    if directory is not None:
        filename = "%s/%s-vcoprediction-%u-%u.hdf" %
                         (ifo, st, dur, directory)
        vco_data.to_hdf5(filename)
                         
    else:
        filename = "%s-vcoprediction-%u-%u.hdf" %
                         (ifo, st, dur)
        vco_data.to_hdf5(filename)
    return filename
