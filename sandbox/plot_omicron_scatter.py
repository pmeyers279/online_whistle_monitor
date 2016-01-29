"""
Read Omicron + VCO triggers and plot them
"""
import optparse
from gwpy.segments import segments
from gwpy.table.lsctables import SnglBurstTable
import numpy as np
from matplotlib.pyplot import rc
import matplotlib.pyplot as plt
rc('text', usetex=True)


# list of expected lines for O1
lines = [(2350, 2450), (2750, 2850), (4200, 4300), (5500, 5700)]


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

params = parse_command_line()

# load data
data = np.load(params.file)
data = remove_lines(data, lines)
VCO_vals = data[:, 0]
central_freqs = data[:, 1]
omicron_snr = data[:, 2]
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
print st
print dur

# make scatter plot
fig = plt.figure()
plt.scatter(VCO_vals[idxs_final], central_freqs[idxs_final],
            c=np.log10(omicron_snr[idxs_final]),
            s=16)
plt.xlabel(r'VCO frequency [kHz from 79 MHz]')
plt.ylabel(r'Omicron central frequency')
plt.title(r'%s for %d - %d, color is omicron snr' %
          (channel.replace('_', '\_'), st, int(st + dur)))
#plt.xlim((min(VCO_vals[idxs_final]), max(VCO_vals[idxs_final])))
plt.xlim((50,70))
plt.ylim((params.flow, params.fhigh))
plt.savefig('%s-OMICRON-SCATTER-%d-%d' % (channel, st, dur))
