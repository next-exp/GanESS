import os
import re
import sys
import glob
import linecache
import argparse


import pandas as pd
import tables as tb

from invisible_cities. io.rwf_io           import rwf_writer 
from invisible_cities. io.run_and_event_io import run_and_event_writer


n_pmts     = 7
aux    = '{*}'

parser = argparse.ArgumentParser()
parser.add_argument('--input', '-i',  help="Path to raw file"              , type=str)
parser.add_argument('--run'  , '-r',  help="Run number"                    , type=int)
parser.add_argument('--hdf5' , '-h5', help="Convert to hdf5"               , action='store_false')
parser.add_argument('--split', '-s',  help="Split into individual events"  , action='store_false')
parser.add_argument('--join' , '-j',  help="Join in files with 1000 events", action='store_false')
args   = parser.parse_args()

inpath     = args.input
run_number = args.run
fhdf5      = args.hdf5
fsplit     = args.split
fjoin      = args.join

if fsplit:
    os.system(f"csplit -z {inpath} /Event/ '{aux}' -s -b '%08d.txt' -f '{os.path.dirname(inpath)}/Event_'")

nevts_file = 1000
inpath = os.path.dirname(inpath) + '/'
infiles = glob.glob(inpath+'Event*txt')

if fjoin:
    for i in range(len(infiles)//nevts_file + 1):
        eventfiles = sorted(glob.glob(f'{inpath}Event_{i:05}[0-9][0-9][0-9].txt'))
        for eventfile in eventfiles:
            os.system(f'cat {eventfile} >> {inpath}Run_{run_number}_file_{i}.txt')
            os.system(f'rm {eventfile}')


infiles    = glob.glob(inpath+'/Run*.txt')[:]
infiles.sort(key=lambda x: os.path.getmtime(x))

if fhdf5:
    for i, infile in enumerate(infiles):
        n_samples_file = int  (re.findall("\d+", linecache.getline(infile, 3))[0])
        tbin           = float(re.findall("\d*\.?\d+", linecache.getline(infile, 4))[1])    

        all_data       = pd.read_csv(infile, delimiter="\t", iterator=True, names=['tbin'] + [str(i) for i in range(8)])

        outfile = infile.replace('.txt', '.h5')
        os.system(f'rm {infile}')

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
        linecache.clearcache()
