The function you want is `histogram/make_histograms`

This function will query the segment database and generate histograms for each lock in between start and end times. It is configured to search by default for omicron triggers generated for DARM on the detchar account. The simplest use case is on the LLO/LHO cluster

```bash
./make_histogram -o L1 -s %(start) -e %(end)
```
Will produce a file titled `L1-GDS-CALIB_STRAIN_VCO-HIST-%(start)-%(duration).png` in the working directory.

All flags for `make_histograms`:

* `-s`
    * Start time (required)
* `-e`
    * End time (required)
* `-o`, `--ifo`
    * **Default** : L1
    * *recommended that you enter something anyway*
    * Observatory (e.g. H1, L1)
* `-c`, `--channel`
    * **Default** : GDS-CALIB_STRAIN
    * %(ifo):%(channel)_Omicron is the channel 
        used to search for omicron triggers.
* `--frames`
    *  **Default** : 0
    * **1** : read from frames
    * **0** : read from nds2

* `--fit-vco`
    * **Default** : 1
    * **1** : create fast vco by linear fitting with mc-f
    * **0** : create fast vco channel by spline 
            interpolation on slow vco.
* `--flag`
    * **Default** : DMT-ANALYSIS_READY:1
    * DQ flag used to query segment database
* `--trig-dir`
    * **Default** : `/home/detchar/triggers/O1/`
    * Base directory to search for triggers
        * `/%(ifo)/%(ifo):%(channel)/%(str(st)[0:5])` is currently automatically appended
        to query current set-up, so be careful changing this.
