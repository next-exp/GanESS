import os
import re
import sys
import glob
import linecache

import pandas as pd
import tables as tb

from invisible_cities. io.rwf_io           import rwf_writer 
from invisible_cities. io.run_and_event_io import run_and_event_writer

inpath     = sys.argv[1]
run_number = int(sys.argv[2])
n_pmts     = 7

infiles    = glob.glob(inpath+'/*.txt')[:]
infiles.sort(key=lambda x: os.path.getmtime(x))

for i, infile in enumerate(infiles):
    n_samples_file = int  (re.findall("\d+", linecache.getline(infile, 3))[0])
    tbin           = float(re.findall("\d*\.?\d+", linecache.getline(infile, 4))[1])    

    all_data       = pd.read_csv(infiles[0], delimiter="\t", iterator=True, names=['tbin'] + [str(i) for i in range(8)])
    
    outfile = inpath + f'Run_{run_number:04}_file_{i:04}.h5'
    if os.path.isfile(outfile):
        print(f'File {i:04} ({os.path.basename(infile)}) of run {run_number:04} already exists, skipping.')
        continue

    with tb.open_file(outfile, mode="w", title="Test file") as h5file:
        write_rwf = rwf_writer(h5file, group_name='RD', table_name='pmtrwf', n_sensors=n_pmts, waveform_length=n_samples_file)
        write_evt = run_and_event_writer(h5file)
        
        while True:
            try:
                event_info     = all_data.get_chunk(5)
                event_number   = int(re.findall("\d+", event_info.tbin.values[0])[0])
                timestamp      = int(re.findall("\d+", event_info.tbin.values[1])[0])
                wvfs           = all_data.get_chunk(n_samples_file).values[:, 1:-1].T
        
                write_rwf(wvfs)
                write_evt(run_number, event_number, timestamp)
            except StopIteration:
                break