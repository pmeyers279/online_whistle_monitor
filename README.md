The function you want is `histogram/make_histograms`

It is configured to search by default for omicron triggers generated for DARM on the detchar account. The simplest use case is on the LLO/LHO cluster

```bash
./make_histogram -o L1 -s %(start) -e %(end)
```
Will produce a file titled `L1-GDS-CALIB_STRAIN_VCO-HIST-%(start)-%(duration)`. 

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
    *  **Default** : False
    * **True** : read from frames
    * **False** : read from nds2

* `--fit-vco`
    * **Default** : True
    * **True** : create fast vco by linear fitting with mc-f
    * **False** : create fast vco channel by spline 
            interpolation on slow vco.
* `--flag`
    * **Default** : DMT-ANALYSIS_READY:1
    * DQ flag used to query segment database
* `--trig-dir`
    * **Default** : `/home/detchar/triggers/O1/`
    * Base directory to search for triggers
        * `/%(ifo)/%(ifo):%(channel)/%(str(st)[0:5])` is currently automatically appended
        to query current set-up, so be careful changing this.
