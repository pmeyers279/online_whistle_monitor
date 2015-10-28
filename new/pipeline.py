from vco_functions import *
from trigger_functions import *
import optparse
from gwpy.segments import Segment


def parse_command_line():
    """
    parse command line
    """
    parser = optparse.OptionParser()
    parser.add_option(
        "-o", "--ifo", help="interferometer, default is L1", dest="ifo",
        type=str, default="L1")
    parser.add_option(
        "-s", help="start time", dest="st", type=int, default=None)
    parser.add_option(
        "-e", help="end time", dest="et", type=int, default=None)
    parser.add_option(
        "-c", "--channel", help="channel for which user wants histogram.",
        dest="channel", type=str, default="GDS-CALIB_STRAIN")
    parser.add_option(
        "--fit-vco", help="fit vco with imc-f or not (True or False)",
        dest="fit", type=bool, default=False)
    parser.add_option(
        "--frames", help="read from frames", type=bool, default=False)
    parser.add_option(
        "--trig-dir", help="omicron trigger base directory", type=str,
        dest="trig_dir", default='/home/detchar/triggers/O1/')
    params, args = parser.parse_args()
    return params

params = parse_command_line()

if not params.st and params.et:
    raise ValueError('need to specify start and end times')

seg = Segment(params.st, params.et)

vco = generate_fast_vco(params.ifo, seg, frames=params.frames, fit=params.fit)
# get list of files that contain omicron triggers
# for our times
files = get_omicron_files(
    params.ifo, seg, channel=params.channel, trig_dir=params.trig_dir)
# get omicron triggers
# and their times
trigs, times, s = read_omicron_trigs(files, seg)

# get VCO frequencies at times of omicron triggers
vco_trigs = get_vco_trigs(vco, times)

# plot VCO triggers in histogram
plot_vco_hist(vco_trigs, seg, '%s:%s' % (params.ifo, params.channel))
