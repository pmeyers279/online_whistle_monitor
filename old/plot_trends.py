from gwpy.timeseries import TimeSeries
from os import listdir
import glob
import matplotlib.pyplot as plt
import numpy as np

files = glob.glob('./L1-vcoprediction*')
idx =1 
means= []
starts = []
for file in files:
    
    data = TimeSeries.from_hdf5(file)
    means.append(np.mean(data.value))
    starts.append(data.times.value[0])

fig = plt.figure()
plt.plot(starts-starts[0],means)
plt.savefig('./mean_vco_value')
means_new = means-means[0]
print means_new








