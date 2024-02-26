import os
import re
import sys
import glob
import linecache
import argparse
import glob
import multiprocessing

import pandas as pd
import tables as tb

from invisible_cities. io.rwf_io           import rwf_writer 
from invisible_cities. io.run_and_event_io import run_and_event_writer

n_pmts     = 7

def run_decode(infile, ifile, run_number):
    n_samples_file = int  (re.findall("\d+", linecache.getline(infile, 3))[0])
    tbin           = float(re.findall("\d*\.?\d+", linecache.getline(infile, 4))[1])    
    all_data       = pd.read_csv(infile, delimiter="\t", iterator=True, names=['tbin'] + [str(i) for i in range(8)])

    outfile = new_dir + f'raw/Run_{run_number}_file_{ifile}_raw.h5'
    os.system(f'rm {infile}')

    if os.path.isfile(outfile):
        print(f'File {outfile} ({os.path.basename(infile)}) of run {run_number:04} already exists, skipping.')
    else:
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

def run_decode_updated_pd(infile, ifile, run_number):
    n_samples_file = int  (re.findall("\d+", linecache.getline(infile, 3))[0])
    tbin           = float(re.findall("\d*\.?\d+", linecache.getline(infile, 4))[1])    
    all_data       = pd.read_csv(infile, delimiter="\t", iterator=True, names=['tbin'] + [str(i) for i in range(9)])

    outfile = new_dir + f'raw/Run_{run_number}_file_{ifile}_raw.h5'
    os.system(f'rm {infile}')

    if os.path.isfile(outfile):
        print(f'File {outfile} ({os.path.basename(infile)}) of run {run_number:04} already exists, skipping.')
    else:
        with tb.open_file(outfile, mode="w", title="Test file") as h5file:
            write_rwf = rwf_writer(h5file, group_name='RD', table_name='pmtrwf', n_sensors=n_pmts, waveform_length=n_samples_file)
            write_evt = run_and_event_writer(h5file)
    
            while True:
                try:
                    event_info     = all_data.get_chunk(5).iloc[: , :-1]
                    event_number   = int(re.findall("\d+", event_info.tbin.values[0])[0])
                    timestamp      = int(re.findall("\d+", event_info.tbin.values[1])[0])
                    wvfs           = all_data.get_chunk(n_samples_file).iloc[: , :-1].values[:, 1:-1].T
    
                    write_rwf(wvfs)
                    write_evt(run_number, event_number, timestamp)
                except StopIteration:
                    break
    linecache.clearcache()

def run_city(origin, city, infile, outfile, config_file):
    os.system(f'city {city} {config_file} -i {infile} -o {outfile}')

parser = argparse.ArgumentParser()
parser.add_argument('--input'    , '-i',  help="Path to raw file", nargs='+')
parser.add_argument('--run'      , '-r',  help="Run number"      , type=int)
parser.add_argument('--cores'    , '-c',  help="# of CPUs"       , type=int, default=4)
parser.add_argument('--decode'   , '-d',  help="Decode data"     , action='store_false')
parser.add_argument('--pmaps'    , '-p',  help="Pmap config file", default="")


args        = parser.parse_args()

infiles     = args.input
infiles.sort(key=os.path.getmtime)

if len(infiles)==0:
    sys.exit(f'No input files.')

inpath      = os.path.dirname(infiles[0])
run_number  = args.run
ncores      = args.cores
fdecode     = args.decode
fpmaps_conf = args.pmaps
fpmaps      = fpmaps_conf != ""

    
new_dir  = inpath + f'/R{run_number}/'
if not os.path.isdir(new_dir):
    os.system(f"mkdir {new_dir}") 

if fdecode:
    decode_dir = new_dir + 'raw/'
    if not os.path.isdir(decode_dir):
        os.system(f"mkdir {decode_dir}") 

    items = [(infile, ifile, run_number) for ifile, infile in enumerate(infiles)]
  
    with multiprocessing.Pool(processes=ncores) as pool:
        pool.starmap(run_decode_updated_pd, items)
    if fpmaps:
        infiles = glob.glob(decode_dir + "*_raw.h5")

if fpmaps:
    origin = "raw"
    city   = "octavia"
    pmaps_dir = new_dir + f'{city}/'
    if not os.path.isdir(pmaps_dir):
        os.system(f"mkdir {pmaps_dir}") 
    os.system(f"cp {fpmaps_conf} {pmaps_dir}octavia.conf") 
    items  = [(origin, city, ifile, ifile.replace(f"_{origin}.h5",f"_{city}.h5").replace(origin, city), fpmaps_conf) 
              for ifile in infiles]
  
    with multiprocessing.Pool(processes=ncores) as pool:
        pool.starmap(run_city, items)