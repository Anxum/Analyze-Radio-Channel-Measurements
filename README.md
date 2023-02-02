# Analyze-Radio-Channel-Measurements
A script to quantify radio channel measurements performed with a correlative sounder and a Zadoff-Chu sequence.

The measurement data this script was originally intended for can be found at [insert Link] and is described in [insert paper]. The provided files are intended to perform the following tasks:

## Measurement.py
A collection of methods to provide the fundamental computations to extract coherence parameters RX power levels from the raw data. The names of the files are expected to be in the following format: '{ Name }_{ Date }_{ Time }_{ fc }MHz_{ fs }MSps_{ capture_interval }ms.dat'

## analyze-measurements.py
process the computed coherence parameters into graphs and save them. This script may be executed with the following arguments:

Argument | Description | Default 
--- | --- | --- 
--inpath | Specifies the path to a single measurement or to a folder containing measurements | input 
--zc-path | Specifies the filepath to the zadoff-chu-sequence. | zc_sequence.npy
--outpath | Specifies the path to store the results. | output
--thresh | Specifies the threshold for coherence calculation. Values are accepted from 0 to 1. | 0.9
--corr-method | Specifies the correlation method in use. Options are: 'classic', 'cyclic' or 'both'. | both
--window-size | Specifies the window size in seconds to calculate the coherence time. | 1
--pts-per-window | Specifies the amount of points used per window to calculate the coherence time. | 10
