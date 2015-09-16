import subprocess
import os

# get segments
cmd = "ligolw_segment_query_dqsegdb -t https://dqsegdb5.phy.syr.edu -q -s 1117324816 -e 1118361616 -a 'L1:DMT-ANALYSIS_READY' | ligolw_print -t segment -c start_time -c end_time"
p = subprocess.Popen([cmd], stdout=subprocess.PIPE,shell=True)
out = p.stdout.read()
print out
segs = out.split('\n')

starts = []
ends = []
for seg in segs:
    if seg == '':
        continue
    starts.append(int(seg.split(',')[0]))
    ends.append(int(seg.split(',')[1]))
print starts
for start,end in zip(starts,ends):
#    cmd = "python get_imc_vco.py L1 %d %d" % (start,end-start)
#    print cmd
#    os.system(cmd)
    cmd2 = "python dump_interpolated_vco.py L1-imc-vco-%d-%d.npy"%(start,end-start)
    os.system(cmd2)

