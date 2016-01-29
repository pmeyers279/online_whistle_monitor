import optparse
from gwpy.segments import segments
from gwpy.table.lsctables import SnglBurstTable
import numpy as np
from matplotlib.pyplot import rc
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
rc('text', usetex=True)


# list of expected lines for O1
lines = [(2350, 2450), (2750, 2850), (4200, 4300), (5500, 5700)]
vco_lines = [(13, 21), (54.5, 62.5)]
freq_cutoffs = [8192, 4096, 2048]


def parse_command_line():
    """
    parse command line
    """
    parser = optparse.OptionParser()
    parser.add_option(
        "-f", "--file", help="trigger file", dest="file",
        type="str", default=None)
    parser.add_option(
        "-l", "--flow", help="low omicron central frequency", dest="flow",
        type=float, default=0)
    parser.add_option(
        "-g", "--fhigh", help="high omicron central frequency", dest="fhigh",
        type=float, default=8192)
    parser.add_option(
        "-o", "--snr-cutoff", help="omicron snr cutoff", dest="cutoff",
        type=float, default=5)
    parser.add_option(
        "--cmax", help="max color", dest="cmax", type=float, default=20)
    params, args = parser.parse_args()
    return params


def remove_lines(data, lines):
    for line in lines:
        mask1 = data[:, 1] < line[0]
        mask2 = data[:, 1] > line[1]
        # add masks since they should be
        # mutually exclusive
        data = data[(mask1 + mask2), :]
    return data


def vco_line_mask(VCO_vals, vco_lines, prior_mask=None):
    first = 1
    for line in vco_lines:
        if first:
            mask = (VCO_vals < line[0]) + (VCO_vals > line[1])
            first = 0
        else:
            mask = np.multiply(
                mask, (VCO_vals < line[0]) + (VCO_vals > line[1]))
    if prior_mask is not None:
        mask = np.multiply(mask, prior_mask)
    return VCO_vals[mask]


def gauss(x, *p):
    A, mu, sigma = p
    return A * np.exp(-(x - mu)**2 / (2. * sigma**2))

# p0 is the initial guess for the fitting coefficients (A, mu and sigma above)


params = parse_command_line()

# load data
data = np.load(params.file)
data = remove_lines(data, lines)
VCO_vals = data[:, 0]
central_freqs = data[:, 1]
omicron_snr = data[:, 2]
old = omicron_snr
omicron_snr[np.where(omicron_snr > params.cmax)[0]] = params.cmax

# make cuts based on parameters
mask1 = central_freqs > params.flow
mask2 = central_freqs < params.fhigh
mask3 = omicron_snr > params.cutoff
idxs_final = np.multiply(np.multiply(mask1, mask2), mask3)

# parse channel name
channel = params.file.split('-')[0] + '-' + params.file.split('-')[1]
st = int(params.file.split('-')[-2])
dur = int(params.file.split('-')[-1][:-4])

colors = ['blue', 'yellow', 'red']
expected_rate = {}
fig = plt.figure()
for i in range(len(freq_cutoffs)):
    expected_rate[freq_cutoffs[i]] = {}
    prior_mask = np.multiply(idxs_final, (central_freqs < freq_cutoffs[i]))
    VCO_line_masked = vco_line_mask(VCO_vals, vco_lines, prior_mask=prior_mask)
    n_tot, bins_tot = np.histogram(VCO_line_masked, 75)
    n2, bins_tot2 = np.histogram(VCO_vals[prior_mask], 75)
    bin_centres = (bins_tot[1:] + bins_tot[:-1]) / 2.
    # p0 is the initial guess for the fitting coefficients (A, mu and sigma
    # above)
    p0 = [1., 0., 1.]
    coeff, var_matrix = curve_fit(gauss, bin_centres, n_tot, p0=p0)
    # Get the fitted curve
    hist_fit = gauss(bin_centres, *coeff)
    for vco_line in vco_lines:
        expected_rate[freq_cutoffs[i]][vco_line] = np.interp((vco_line[0] + vco_line[1]) / 2.,
                                                             bin_centres, hist_fit)
        print ((vco_line[0] + vco_line[1]) / 2.)
        print expected_rate[freq_cutoffs[i]][vco_line] / float(dur)
        plt.scatter(
            (vco_line[0] + vco_line[1]) /
            2., expected_rate[freq_cutoffs[i]][vco_line] / float(dur),
            s=64, c=colors[i])
    plt.bar(bins_tot2[:-1], n2 / float(dur),
            width=(bins_tot[1] - bins_tot[0]),
            label='Central Freqs $\leq$ %d' % (freq_cutoffs[i]),
            facecolor=colors[i], edgecolor=colors[i], alpha=0.2)
    plt.plot(bin_centres, hist_fit / float(dur), c=colors[i], linewidth=2)
plt.xlabel('VCO Frequency [kHz from 79 MHz]')
plt.ylabel('Rate [Hz]')
plt.title('Omicron Rate Plot for all triggers (expect gaussian')
ax2 = plt.gca()
handles, labels = ax2.get_legend_handles_labels()
ax2.legend(handles, labels)
ax2.set_yscale('log')
plt.ylim((1e-5, 1))
plt.savefig('GDS-CALIB_STRAIN-VCO-HIST-%d-%d' % (st, dur))
plt.close()


first = 1
# do this for each line
for vco_line in vco_lines:
    whistle_mask = np.multiply(
        np.multiply((VCO_vals > vco_line[0]),
                    VCO_vals < vco_line[1]), idxs_final)
    # get a histogram of all VCO values for glitches...
    n, bins = np.histogram(omicron_snr[whistle_mask], bins=50)
    rate = n / float(dur)
    fig = plt.figure()
    # plt.bar(bins[:-1], rate, width=(bins[1] - bins[0]),
    #        facecolor='blue', alpha=0.3)
    for i in range(len(freq_cutoffs)):
        # histogram ALL triggers in range
        # VCO_line_masked = vco_line_mask(VCO_vals, vco_lines, prior_mask=np.multiply(idxs_final,
        #                                                                             (central_freqs < freq_cutoffs[i])))
        # n_tot, bins_tot = np.histogram(VCO_line_masked, 75)
        # bin_centres = (bins_tot[1:] + bins_tot[:1]) / 2.
        # p0 is the initial guess for the fitting coefficients (A, mu and sigma
        # above)
        # p0 = [1., 0., 1.]
        # coeff, var_matrix = curve_fit(gauss, bin_centres, n_tot, p0=p0)
        # Get the fitted curve
        # hist_fit = gauss(bin_centres, *coeff)
        # expected_rate = np.interp(
        #     (vco_line[0] + vco_line[1]) / 2., bin_centres, hist_fit)
        # print ((vco_line[0] + vco_line[1]) / 2.)
        # print expected_rate / float(dur)
        # Plot expected total rate based on gaussian fit
        print ((vco_line[0] + vco_line[1]) / 2.)
        print expected_rate[freq_cutoffs[i]][vco_line] / float(dur)
        plt.scatter(
            (bins[0] + bins[1]) / 2,
            expected_rate[freq_cutoffs[i]][vco_line] / float(dur),
            s=64, c=colors[i])
        # histogram just whistle triggers
        n, bins = np.histogram(
            omicron_snr[
                np.multiply(whistle_mask, (central_freqs < freq_cutoffs[i]))],
            bins=bins)
        rate = n / float(dur)
        rate = np.cumsum(rate[::-1])[::-1]
        # plot cumulative rate vs. SNR
        plt.bar(bins[:-1], rate, width=(bins[1] - bins[0]),
                edgecolor=colors[i], facecolor=colors[i],
                label='Central Freqs $\leq$ %d' % freq_cutoffs[i], alpha=0.2)
    ax = plt.gca()
    handles, labels = ax.get_legend_handles_labels()
    plt.xlabel('Omicron SNR')
    plt.ylabel('Rate [Hz]')
    plt.title('Omicron glitch rates contributions for VCO Line %4.2f - %4.2f MHz' %
              (79 + vco_line[0] / 100., 79 + vco_line[1] / 100.))
    ax.legend(handles, labels)
    ax.set_yscale('log')
    plt.savefig('GDS-CALIB_STRAIN-RATES-%s-%d-%d' %
                (str(int((vco_line[1] + vco_line[0]) / 2.)), st, dur))
    plt.close()
