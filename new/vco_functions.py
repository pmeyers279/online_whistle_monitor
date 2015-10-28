from gwpy.timeseries import TimeSeries
from glue import datafind
import numpy as np
from scipy.interpolate import UnivariateSpline


def generate_fast_vco(ifo, segment, frames=False, fit=True):
    """
    Parameters:
    -----------
        ifo : start
            interferometer, e.g. 'L1'
        segment : array like
            time segment. first entry start second entry end
        frames : bool
            read from frames or nds2
        fit : bool
            fit from imc-f (default)
            or spline interpolation

    Returns:
    --------
        vco_data : saves file 'L1:IMC-VCO_PREDICTION-st-dur.hdf'
    """
    st = segment[0]
    et = segment[1]
    chan1_pat = '%s:SYS-TIMING_C_FO_A_PORT_11_SLAVE_CFC_FREQUENCY_5'
    chan2_pat = '%s:IMC-F_OUT_DQ'

    if frames:
        connection = datafind.GWDataFindHTTPConnection()
        cache = connection.find_frame_urls(
            ifo[0], '%s_R' % ifo, st, et, urltype='file')
        if fit:
            print "Get IMC"
            imc = TimeSeries.read(cache, chan2_pat % ifo, st, et)
        print "Get psl"
        pslvco = TimeSeries.read(cache, chan1_pat % ifo, st, et + 1)
    else:
        if fit:
            print "Get IMC"
            imc = TimeSeries.fetch(chan2_pat % ifo, st, et)
        print "Get psl"
        pslvco = TimeSeries.fetch(chan1_pat % ifo, st, et + 1)

    print "Downsample psl"
    pslvco = pslvco[16 + 8::16]

    print "Saving"
    if fit:
        print "Downsample imc"
        imc_srate = int(imc.sample_rate.value)
        imc = imc[imc_srate / 2::imc_srate]
        data = np.array((imc.data, pslvco.data)).T
        vco_interp = fit_with_imc(data, imc)
    else:
        vco_interp = interp_spline(pslvco)

    chan = "%s:IMC-VCO_PREDICTION" % (ifo,)
    vco_data = TimeSeries(vco_interp, epoch=st,
                          sample_rate=imc.sample_rate.value,
                          name=chan, channel=chan)
    return vco_data


def fit_with_imc(data, imc):
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
    samp_times = np.arange(len(imc)) / 256.

    coeffs0 = np.interp(samp_times, times, coeffs[:, 0])
    coeffs1 = np.interp(samp_times, times, coeffs[:, 1])

    vco_interp = coeffs0 * imc.data + coeffs1
    return vco_interp


def interp_spline(pslvco, srate=256):
    fvco = UnivariateSpline(np.arange(len(pslvco))+0.5, pslvco, k=3)
    vco_interp = fvco(np.arange(srate*len(pslvco))/float(srate))
    return vco_interp

