from utils import *
from create_fast_vco import *
from gwpy.segments import Segment

segments = get_segments(ifo, st, et)
directory = '/home/meyers/online_whistle_monitor/'
write_segments(ifo, st, et, directory=directory)

# get omicron files
channel = 'GDS-CALIB_STRAIN_Omicron'
omicron_directory = '/home/detchar/triggers/ER8/L1/%s' % channel
time_dirs = np.arange(st[0:4], et[0:4] + 1)
omicron_files = []
for od in time_dirs:
    files = sorted(glob.glob('%s/%s/*.xml*' % (omicron_directory, od)))
    omicron_files.append(files)

vco_files = []
spans = []
for segment in segments.active:
    st_seg = segment[0]
    et_seg = segment[-1]
    spans.append(Segment(st_set, et_seg))
    filename = write_fast_vco(
        ifo, st_seg, et_seg - st_seg, directory=directory)
    vco_files.append(filename)


triggers = Triggers(channel=channel)
for v in vco_files:
    start = int(v.split('-')[-2])
    end = start + int(v.split('-')[-1][:-4])
    span = Segment(start, end)
    trigs, times = read_omicron_trigs(omicron_files, span)
    amp = get_vco_data(v, times[v])
    triggers.extend(trigs, amp, span, v)

plot_trig_rates(trigs, spans, channel+'-ALL-', files=None)
plot_vco_hist(amp, spans, channel)
plot_trig_rates(new_trigs, spans, channel+'-WHISTLE-', files=None)
