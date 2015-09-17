Currently only the "old/" directory works. 

Generate "fast" VCO channel with

>>> python get_imc_vco.py %(ifo) %(start time) %(duration)

This produces two files: a numpy array with the VCO and the IMC-F data sets and an hdf file that can be opened with gwpy
TimeSeries.from_hdf5('filename') that is the "fast" VCO channel

---------

Generate Omicron histogram and rate estimate of whistles (for a single line):

1. copy over omciron_histogram.py to a new location
2. open it up, add VCO files you want to plot rates for to vco_files = [] and the channel you want to use for triggers
for the channel variable. 
3. run >>> python omicron_histogram.py


NOTE: currently we have a working example in github and I'd like to keep it that way. Soon I'll add something that doesn't require
you to change the code itself.

