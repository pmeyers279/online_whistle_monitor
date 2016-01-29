import optparse
from gwpy.segments import segments
from gwpy.table.lsctables import SnglBurstTable
import numpy as np
from matplotlib.pyplot import rc
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import markup
rc('text', usetex=True)


# list of expected lines for O1
lines = [(2350, 2450), (2750, 2850), (4200, 4300), (5500, 5700)]
#vco_lines = [(13, 21), (54.5, 62.5)]
vco_lines = [(-214,-206),(-175,-170), (-142,-136), (-114,-110), (-80,-70), (-50,-30)]
freq_cutoffs = [8192, 4096]#, 2048]


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
    parser.add_option(
	"--line-file", help="file that has expected whistle lines in it", dest="line_file",
	type="str", default="lines.txt")
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

def load_lines(lines_file):
    width = 8.192
    line_data = np.loadtxt(lines_file)
    vco_lines = []
    scalings = []
    for i in range(len(line_data[:,0])):
	# note we do things in kHz
	linewidth=width/float(line_data[i,1])
    	vco_lines.append((int(line_data[i,0])-linewidth, int(line_data[i,0])+linewidth))
	scalings.append(int(line_data[i,1]))
    return vco_lines,scalings

params = parse_command_line()

# load expected vco_lines
vco_lines,scalings = load_lines(params.line_file)
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

# set up colors
colors = ['blue', 'yellow', 'red']
expected_rate = {}
expected_rate_corrected = {}
WHISTLE_SNR = {}
WHISTLE_RATE = {}
# start VCO omicron histogram plot
fig = plt.figure()

# loop over frequency cut offs
for i in range(len(freq_cutoffs)):
    # nested dicts for expected rates and eventually 
    # whistle SNR. whsitle SNR and rate will come later...
    expected_rate[freq_cutoffs[i]] = {}
    expected_rate_corrected[freq_cutoffs[i]] = {}
    WHISTLE_SNR[freq_cutoffs[i]] = {}
    WHISTLE_RATE[freq_cutoffs[i]] = {}
    # generate mask for expected whistle bins
    prior_mask = np.multiply(idxs_final, (central_freqs < freq_cutoffs[i]))
    # create unmasked ans masked masked vco histograms
    VCO_line_masked = vco_line_mask(VCO_vals, vco_lines, prior_mask=prior_mask)
    n_tot, bins_tot = np.histogram(VCO_line_masked, 75)
    n2, bins_tot2 = np.histogram(VCO_vals[prior_mask], 75)
    bin_centres = (bins_tot[1:] + bins_tot[:-1]) / 2.
    # p0 is the initial guess for the fitting coefficients (A, mu and sigma
    # above).
    p0 = [1,np.median(VCO_line_masked), (np.max(VCO_line_masked) - np.min(VCO_line_masked)) / 4.]
    # fit masked vco histogram with gaussian
    coeff, var_matrix = curve_fit(gauss, bin_centres, n_tot, p0=p0)
    hist_fit = gauss(bin_centres, *coeff)
    for vco_line in vco_lines:
	# interpolate for expected rate at line
        expected_rate[freq_cutoffs[i]][vco_line] = np.interp((vco_line[0] + vco_line[1]) / 2.,
                                                             bin_centres, hist_fit) / float(dur)
	# reweight for bin sizes
	expected_rate_corrected[freq_cutoffs[i]][vco_line] = expected_rate[freq_cutoffs[i]][vco_line] * (vco_line[1] - vco_line[0]) / float(bins_tot2[1] - bins_tot2[0]) 
	# plot expected rate points
        plt.scatter(
            (vco_line[0] + vco_line[1]) /
            2., expected_rate[freq_cutoffs[i]][vco_line],
            s=64, c=colors[i])
    # plot unmasked histogram
    plt.bar(bins_tot2[:-1], n2 / float(dur),
            width=(bins_tot[1] - bins_tot[0]),
            label='Central Freqs $\leq$ %d' % (freq_cutoffs[i]),
            facecolor=colors[i], edgecolor=colors[i], alpha=0.2)
    # plot gaussian fit to unmasked histogram (note fit is for MASKED histogram!)
    plt.plot(bin_centres, hist_fit / float(dur), c=colors[i], linewidth=2)
plt.xlabel('VCO Frequency [kHz from 79 MHz]')
plt.ylabel('Rate [Hz]')
plt.title('Omicron Rate Plot for all triggers (expect gaussian')
ax2 = plt.gca()
handles, labels = ax2.get_legend_handles_labels()
ax2.legend(handles, labels)
ax2.set_yscale('log')
plt.ylim((1./(2*dur), 1))
plt.savefig('%s-VCO-HIST-%d-%d' % (channel, st, dur))
plt.close()


# do this for each line
for (vco_line,scaling) in zip(vco_lines,scalings):
    # create a mask to pick out whistles (NOT remove them like above...)
    whistle_mask = np.multiply(
        np.multiply((VCO_vals > vco_line[0]),
                    VCO_vals < vco_line[1]), idxs_final)
    # get a histogram of omicron snr for this whistle frequency
    n, bins = np.histogram(omicron_snr[whistle_mask], bins=50)
    rate = n / float(dur)
    fig = plt.figure()
    fig.subplots_adjust(hspace=0.5)
    # plt.bar(bins[:-1], rate, width=(bins[1] - bins[0]),
    #        facecolor='blue', alpha=0.3)
    plt.subplot(211)
    ax = []
    ax.append(fig.add_subplot(211))
    ax.append(fig.add_subplot(212))
    for i in range(len(freq_cutoffs)):
        ax[0].scatter(
            (bins[0] + bins[1]) / 2,
            expected_rate_corrected[freq_cutoffs[i]][vco_line],
            s=64, c=colors[i])
        # histogram just whistle triggers
        n, bins = np.histogram(
            omicron_snr[
                np.multiply(whistle_mask, (central_freqs < freq_cutoffs[i]))],
            bins=bins)
        rate = n / float(dur)
        rate = np.cumsum(rate[::-1])[::-1]
        WHISTLE_RATE[freq_cutoffs[i]][vco_line] = rate[0] - expected_rate_corrected[freq_cutoffs[i]][vco_line]
        WHISTLE_SNR[freq_cutoffs[i]][vco_line] = rate[0] / float(expected_rate_corrected[freq_cutoffs[i]][vco_line])
        # plot cumulative rate vs. SNR
        ax[0].bar(bins[:-1], rate, width=(bins[1] - bins[0]),
                edgecolor=colors[i], facecolor=colors[i],
                label='Central Freqs $\leq$ %d' % freq_cutoffs[i], alpha=0.2)
    handles, labels = ax[0].get_legend_handles_labels()
    ax[0].set_xlabel('Omicron SNR')
    ax[0].set_ylabel('Rate [Hz]')
    vco_center = (vco_line[0] + vco_line[1])/2.
    ax[0].set_title(r'Omicron glitch rates contributions for VCO Line %4.3f MHz' % (79 + vco_center / 1000.))
    ax[0].legend(handles, labels)
    ax[0].set_yscale('log')
    ax[1].scatter(VCO_vals, central_freqs,c=np.log10(omicron_snr))
    ax[1].plot([vco_center,vco_center + 8.192/scaling],[0, 8192],c='r')
    ax[1].plot([vco_center,vco_center - 8.192/scaling],[0, 8192],c='r',label='Expected V shape')
    ax[1].plot([vco_line[0],vco_line[0]],[0, 8192],c='k')
    ax[1].plot([vco_line[1],vco_line[1]],[0, 8192],c='k',label='Bin Borders for rate')
    ax[1].set_xlim([vco_center - 8.192, vco_center + 8.192])
    ax[1].set_ylim((params.flow, max(central_freqs[np.multiply(VCO_vals >vco_line[0], VCO_vals < vco_line[1])])))
    ax[1].set_title(r'Omicron scatter plot (colors=SNR) around RF Line')
    ax[1].set_xlabel(r'VCO Frequency [kHz from 79 MHz]')
    ax[1].set_ylabel(r'Omicron Central Frequency [Hz]')
    plt.savefig('%s-%s-%d-%d' %
                (channel, str(int((vco_line[1] + vco_line[0]) / 2.)), st, dur))
    plt.close()

############
# GENERATE HTML
############

page = markup.page()
page.init(css=("main.css",
	   "https://maxcdn.bootstrapcdn.com/bootstrap/3.3.2/css/bootstrap.min.css"),
      title="VCO Summary for Whistles for %s %d-%d" % (channel, st, dur))
page.div(class_="container")
page.div(class_="jumbotron")
page.h2("%s"%(channel.split('/')[-1]))
page.div.close()

page.div(class_="panel panel-danger")
page.div(class_="panel-heading")
page.h4("Parameters and Channels")
page.div.close()
page.div(class_="panel-body")
page.table(class_="table")
page.tr()
page.td("Frequency cutoff")
page.td("RF line central frequency [MHz]")
page.td("Assumed scaling factor")
page.td("Whistle rate [Hz]")
page.td("Whistle SNR")
page.tr.close()
for co in freq_cutoffs:
    page.td('%d' % co,rowspan=len(vco_lines) + 1)
    for (vco_line,scale) in zip(vco_lines,scalings):
	page.tr()
	page.td('%4.3f'%(79 + (vco_line[0] + vco_line[1])/2000.))
	page.td('%d'%(scale))
	page.td('%4.3f'%(WHISTLE_RATE[co][vco_line]))
	page.td('%4.3f'%(WHISTLE_SNR[co][vco_line]))
	page.tr.close()
page.table.close()
page.div.close() # panel-body
page.div.close()

page.div(class_="panel panel-primary")
page.div(class_="panel-heading")
page.h3("VCO histogram of all  omicron triggers")
page.div.close()
page.div(class_="panel-body")
page.p("This plot shows a histogram for all triggers, sorted by the frequency of the IMC VCO at the time of the glitch")
page.p("The lines indicate a gaussian fit to the background distribution after the bins assumed to be associated with whistles have been removed")
page.p("The dots indicate the expected background rate associated with that VCO frequency for this lock")
page.div(class_='col-md-6')
page.a(href='%s-%s-%d-%d.png' %
	(channel.split('/')[-1], 'VCO-HIST', st, dur))
page.img(src='%s-%s-%d-%d.png' %
	(channel.split('/')[-1], 'VCO-HIST', st, dur))
page.a.close()
page.div.close() # img
page.div.close()

page.div.close()

for vco_line in vco_lines:
    page.div(class_="panel panel-success")
    page.div(class_="panel-heading")
    # make sure to divide by extra factor of 1000...we're going to MHz now....
    page.h3('Line at %4.3f MHz'%(79 + (vco_line[0] + vco_line[1])/2000.), align='center')
    page.div.close()
    page.div(class_="panel-body")
    page.p("(Top) Cumulative rate plot for the whistle frequency.")
    page.p("(Top) Dots on rate plot indicate the expected background, taken from the top plot on this page (corrected for bin size)")
    page.p("(Bottom) Omicron glitch plot. We expect a 'V' shape for whistles. The red is the assumed shape based on scaling factor. The black indicates the width of the bin used for the rate and background analyses")
    page.div(class_='col-md-6')
    page.a(href='%s-%s-%d-%d.png' %
                (channel.split('/')[-1], str(int((vco_line[1] + vco_line[0]) / 2.)), st, dur))
    page.img(src='%s-%s-%d-%d.png' %
                (channel.split('/')[-1], str(int((vco_line[1] + vco_line[0]) / 2.)), st, dur))
    page.a.close()
    page.div.close() # img
    page.div.close()
    page.div.close()


page.div.close()# container
#outpage = '%d-%d.html' % (params.st, params.et)
f = open('%s-%d-%d.html'%(channel, st, dur),'w')
print >> f, page
f.close()




