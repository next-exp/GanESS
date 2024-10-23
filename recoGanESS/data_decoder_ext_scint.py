import os
import re
import sys
import glob
import linecache
import argparse
import glob
import multiprocessing
import pywt

import numpy as np
import pandas as pd
import tables as tb

from collections import defaultdict

import gres.database.load_db as db
from   gres.cities.components import calibrate_pmts, calibrate_pmts_wf_data

from invisible_cities. io.rwf_io           import rwf_writer 
from invisible_cities. io.run_and_event_io import run_and_event_writer
import invisible_cities.io.pmaps_io as pmapio
import invisible_cities.io.run_and_event_io as rio
from scipy.signal import butter,filtfilt


parser = argparse.ArgumentParser()
parser.add_argument('--input'    , '-i' , help="Path to raw file", nargs='+')
parser.add_argument('--run'      , '-r' , help="Run number"      , type=int)
parser.add_argument('--cores'    , '-c' , help="# of CPUs"       , type=int, default=4)
parser.add_argument('--decode'   , '-d' , help="Decode data"     , action='store_false')
parser.add_argument('--pmaps'    , '-p' , help="Pmap config file", default="")
parser.add_argument('--pmaps_dst', '-pd', help="Pmap DST flag"   , action='store_true')
parser.add_argument('--detector', '-det', help="Pmap DST flag"   , default="gap")


args        = parser.parse_args()

run_number  = args.run
ncores      = args.cores
fdecode     = args.decode
fpmaps_conf = args.pmaps
fpmaps      = fpmaps_conf != ""
fpmaps_dst  = args.pmaps_dst
fdet_name   = args.detector

data_pmt   = db.DataPMT(fdet_name, run_number)

n_pmts     = len(data_pmt)

n_baseline = 225
timebin    = 8


n_s1_0     = int(9 * 1000 / timebin)
n_s1_1     = int(14 * 1000 / timebin)
n_s1_2     = int(19 * 1000 / timebin)
n_s1_3     = int(24 * 1000 / timebin)
n_s1_4     = int(29 * 1000 / timebin)
n_s1_5     = int(34 * 1000 / timebin)
n_s1       = int(39 * 1000 / timebin)

n_mid      = int(24 * 1000 / timebin)
n_maw      = 250
thr_maw    = 0.2

def denoise_wavelet(signal, wavelet='db1', level=1):
    """
    Denoise a 1D signal using wavelet transform with soft thresholding and universal thresholding.

    Args:
    signal (ndarray): Input signal to be denoised.
    wavelet (str): Name of the wavelet to be used (default is 'db1').
    level (int): Number of wavelet decomposition levels (default is 1).

    Returns:
    ndarray: Denoised signal.
    """
    # Perform wavelet decomposition
    coeffs = pywt.wavedec(signal, wavelet, level=level)

    # Estimate the universal threshold (using median absolute deviation)
    sigma = np.median(np.abs(coeffs[-1])) / 0.6745
    threshold = sigma * np.sqrt(2 * np.log(len(signal)))

    # Thresholding: soft thresholding
    coeffs = [pywt.threshold(c, threshold, mode="soft") for c in coeffs]

    # Reconstruct the signal
    denoised_signal = pywt.waverec(coeffs, wavelet)
    
    return denoised_signal


fs     = 1/8e-9    # sample rate, Hz
cutoff = 5e6       # desired cutoff frequency of the filter, Hz (lower will pass allowed)
nyq    = 0.5 * fs  # Nyquist Frequency
order  = 5      

def butter_lowpass_filter(data, nyq, cutoff, fs, order):
    normal_cutoff = cutoff / nyq
    # Get the filter coefficients 
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    y = filtfilt(b, a, data)
    return y



def run_decode(infile, ifile, run_number):
    n_samples_file = int  (re.findall("\d+", linecache.getline(infile, 3))[0])
    tbin           = float(re.findall("\d*\.?\d+", linecache.getline(infile, 4))[1])    
    all_data       = pd.read_csv(infile, delimiter="\t", iterator=True, names=['tbin'] + [str(i) for i in range(8)])

    outfile = new_dir + f'raw/Run_{run_number}_file_{ifile}_raw.h5'
#    os.system(f'rm {infile}')

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
                    wvfs           = all_data.get_chunk(n_samples_file).values[:, 1:].T
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
#    os.system(f'rm {infile}')

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

def generate_wvf_dataset(infile, outfile):
    rwf_to_cwf      = calibrate_pmts_wf_data(fdet_name, run_number, n_baseline, n_maw, thr_maw)
    
    x_pmt = data_pmt.X.values
    y_pmt = data_pmt.Y.values
    
    pmts_info = {}
    with tb.open_file(infile, 'r') as h5in:
        evts      = np.array([f[0] for f in h5in.root.Run.events[:]])
        timestamp = np.array([f[1] for f in h5in.root.Run.events[:]])
        all_wvfs  = h5in.root.RD.pmtrwf[:]
        all_cwvfs = [rwf_to_cwf(wvfs) for wvfs in all_wvfs]
        cwvfs_pmt = np.array([cwvfs[0] for cwvfs in all_cwvfs])
        cwvfs_sum = np.array([cwvfs[2] for cwvfs in all_cwvfs])

        lim_n_s1   = min(n_s1, cwvfs_sum.shape[1]-1)
        lim_n_s1_0 = min(n_s1_0, cwvfs_sum.shape[1]-1)
        lim_n_s1_1 = min(n_s1_1, cwvfs_sum.shape[1]-1)
        lim_n_s1_2 = min(n_s1_2, cwvfs_sum.shape[1]-1)
        lim_n_s1_3 = min(n_s1_3, cwvfs_sum.shape[1]-1)
        lim_n_s1_4 = min(n_s1_4, cwvfs_sum.shape[1]-1)
        lim_n_s1_5 = min(n_s1_5, cwvfs_sum.shape[1]-1)
        lim_n_mid = min(n_mid, cwvfs_sum.shape[1]-2)
#        d_cwvfs_sum = np.array([butter_lowpass_filter(cwvfs[2], nyq, cutoff, fs, order) for cwvfs in all_cwvfs])
        all_sum    = cwvfs_sum.sum(axis=1)
        all_h      = cwvfs_sum.max(axis=1)
        all_t      = cwvfs_sum.argmax(axis=1)
        all_sum_s2  = cwvfs_sum[:,    lim_n_s1:].sum(axis=1)
        all_h_s2    = cwvfs_sum[:,    lim_n_s1:].max(axis=1)
        all_t_s2   = (cwvfs_sum[:,    lim_n_s1:].argmax(axis=1)+lim_n_s1)*timebin
#        d_all_sum    = d_cwvfs_sum.sum(axis=1)
#        d_all_h      = d_cwvfs_sum.max(axis=1)

        pmts_sum = cwvfs_pmt.sum(axis=2)
        pmts_h   = cwvfs_pmt.max(axis=2)
    
        ped_sum      = cwvfs_sum[:,    :n_baseline].sum(axis=1)
        ped_pmts_sum = cwvfs_pmt[:, :, :n_baseline].sum(axis=2)

        all_sum_s1  = cwvfs_sum[:,    :lim_n_s1].sum(axis=1)
        all_h_s1    = cwvfs_sum[:,    :lim_n_s1].max(axis=1)
        all_t_s1    = cwvfs_sum[:,    :lim_n_s1].argmax(axis=1) * timebin

        all_sum_s1_0  = cwvfs_sum[:, lim_n_s1_0:lim_n_s1_1].sum(axis=1)
        all_sum_s1_1  = cwvfs_sum[:, lim_n_s1_1:lim_n_s1_2].sum(axis=1)
        all_sum_s1_2  = cwvfs_sum[:, lim_n_s1_2:lim_n_s1_3].sum(axis=1)
        all_sum_s1_3  = cwvfs_sum[:, lim_n_s1_3:lim_n_s1_4].sum(axis=1)
        all_sum_s1_4  = cwvfs_sum[:, lim_n_s1_4:lim_n_s1_5].sum(axis=1)
        all_sum_s1_5  = cwvfs_sum[:, lim_n_s1_5:lim_n_s1  ].sum(axis=1)

        #all_t_s2    = cwvfs_sum[:,    lim_n_s1:].argmax(axis=1) * timebin

        all_sum_mid  = cwvfs_sum[:,    lim_n_mid:lim_n_s1].sum(axis=1)
        all_h_mid    = cwvfs_sum[:,    lim_n_mid:lim_n_s1].max(axis=1)
        all_t_mid    = cwvfs_sum[:,    lim_n_mid:lim_n_s1].argmax(axis=1)
        pmts_sum_s1 = cwvfs_pmt[:, :, :lim_n_s1].sum(axis=2)
        pmts_h_s1   = cwvfs_pmt[:, :, :lim_n_s1].max(axis=2)
        pmts_sum_s2 = cwvfs_pmt[:, :, lim_n_s1:].sum(axis=2)
        pmts_h_s2   = cwvfs_pmt[:, :, lim_n_s1:].max(axis=2)
                

        sum_e0  = cwvfs_sum[:, 0]
        pmts_e0 = cwvfs_pmt[:, :, 0]
        
        xb = np.zeros(sum_e0.shape)
        yb = np.zeros(sum_e0.shape)
        
        for i in range(n_pmts):
            xb += x_pmt[i] * pmts_sum[:, i]
            yb += y_pmt[i] * pmts_sum[:, i]
            pmts_info[f'e_pmt{i}']       = pmts_sum    [:, i]
            pmts_info[f'h_pmt{i}']       = pmts_h      [:, i]
            pmts_info[f'e_pmt{i}_s1']    = pmts_sum_s1 [:, i]
            pmts_info[f'h_pmt{i}_s1']    = pmts_h_s1   [:, i]
            pmts_info[f'e_pmt{i}_s2']    = pmts_sum_s2 [:, i]
            pmts_info[f'h_pmt{i}_s2']    = pmts_h_s2   [:, i]
            pmts_info[f'ped_sum_pmt{i}'] = ped_pmts_sum[:, i]
            pmts_info[f'e0_pmt{i}']      = pmts_e0     [:, i]

        xb = xb/all_sum
        yb = yb/all_sum
        rb = np.sqrt(xb**2+yb**2)

        dst = pd.DataFrame({'event':evts, 'timestamp':timestamp, 'x':xb, 'y':yb, 'r':rb,
                            'energy':all_sum, 'height':all_h, #'d_energy':d_all_sum, 'd_height':d_all_h,
                            'ped_sum':ped_sum, 'e0':sum_e0,
                            's1_energy':all_sum_s1, 's1_height':all_h_s1,
                            's1_energy_0':all_sum_s1_0,
                            's1_energy_1':all_sum_s1_1,
                            's1_energy_2':all_sum_s1_2,
                            's1_energy_3':all_sum_s1_3,
                            's1_energy_4':all_sum_s1_4,
                            's1_energy_5':all_sum_s1_5,
                            's2_energy':all_sum_s2, 's2_height':all_h_s2,
                            'mid_energy':all_sum_mid, 'mid_height':all_h_mid,
                            's1_t':all_t_s1, 'mid_t':all_t_mid, 's2_t':all_t_s2,'t_all':all_t,
                            **pmts_info,
                            'file':np.full(sum_e0.shape, os.path.basename(infile))})
        
        dst.to_hdf(outfile, 'WVF')

def generate_pmaps_dst(infile, outfile, fS2=True):
    evt_kr      = []
    t_kr        = []
    peak_kr     = []
    npeak_kr    = []
    ene_kr      = []
    h_kr        = []
    t0_kr       = []
    t1_kr       = []
    filename_kr = []
    
    ene_pmts    = defaultdict(list)
    
    pm         = pmapio.load_pmaps_as_df(infile)[1] if fS2 else pmapio.load_pmaps_as_df(infile)[0]
    timestamp  = rio.read_run_and_event(infile)[1]
    map_table  = timestamp.set_index('evt_number')['timestamp'].to_dict()
    pm['timestamp'] = pm['event'].map(map_table)
    
    g          = pm.groupby(['event', 'peak'])
    evt_kr.extend(g.event.min().to_list())
    t_kr.extend(g.timestamp.min().to_list())
    peak_kr.extend(g.peak.min().to_list())
    ene_kr.extend(g.ene.sum().to_list())
    h_kr  .extend(g.ene.max().to_list())
    t0_kr .extend(g.time.min().to_list())
    t1_kr .extend(g.time.max().to_list())
    
    pm_pmts    = pmapio.load_pmaps_as_df(infile)[4] if fS2 else pmapio.load_pmaps_as_df(infile)[3]
    gpm        = pm_pmts.groupby(['event', 'peak', 'npmt'])

    all_pmts_ene = gpm.ene.sum().reset_index()
    all_pmts_h   = gpm.ene.max().reset_index()
    
    for i in range(n_pmts):
        ene_pmts[f'e_pmt{i}'].extend(all_pmts_ene[all_pmts_ene.npmt==i].ene.to_list())
        ene_pmts[f'h_pmt{i}'].extend(all_pmts_h  [all_pmts_h  .npmt==i].ene.to_list())
        
    filename_kr.extend(np.full(fill_value=os.path.basename(infile),shape=len(g.time.max().to_list())))
    
    kr_data = pd.DataFrame({'event':evt_kr, 'timestamp':t_kr,
                            'peak':peak_kr, 'energy':ene_kr, 'height':h_kr, 
                            't0':t0_kr, 't1':t1_kr, 
                            **ene_pmts,
                            'file':filename_kr})
    kr_data['width'] = kr_data['t1'] - kr_data['t0'] 
    xb = np.zeros(len(kr_data))
    yb = np.zeros(len(kr_data))
    
    for i in range(7):
        xb += x[i] * kr_data[f'e_pmt{i}']
        yb += y[i] * kr_data[f'e_pmt{i}']
    kr_data['x'] = xb/kr_data.energy
    kr_data['y'] = yb/kr_data.energy
    
    kr_data['npeaks'] = kr_data.groupby('event').peak.transform('count')
    kr_data.to_hdf(outfile, 'S2' if fS2 else 'S1')






infiles     = args.input
infiles.sort(key=os.path.getmtime)

if len(infiles)==0:
    sys.exit(f'No input files.')

inpath      = os.path.dirname(infiles[0])





x = data_pmt.X.values
y = data_pmt.Y.values

    
new_dir  = inpath + f'/R{run_number}/'

if not os.path.isdir(new_dir):
    if 'raw' in inpath:
        pass
    elif 'octavia' in inpath:
        pass
    else:
        os.system(f"mkdir {new_dir}") 

if fdecode:
    decode_dir = new_dir + 'raw/'
    if not os.path.isdir(decode_dir):
        os.system(f"mkdir {decode_dir}") 

    items = [(infile, ifile, run_number) for ifile, infile in enumerate(infiles)]
  
    with multiprocessing.Pool(processes=ncores) as pool:
        #pool.starmap(run_decode_updated_pd, items)
        pool.starmap(run_decode, items)
    if fpmaps:
        infiles = glob.glob(decode_dir + "*_raw.h5")

if fpmaps:
    origin = "raw"
    city   = "octavia"
    pmaps_dir = os.path.dirname(infiles[0]).replace(origin, city)#new_dir + f'{city}/'
    if not os.path.isdir(pmaps_dir):
        os.system(f"mkdir {pmaps_dir}") 
    os.system(f"cp {fpmaps_conf} {pmaps_dir}/octavia.conf") 
    items  = [(origin, city, ifile, ifile.replace(f"_{origin}.h5",f"_{city}.h5").replace(origin, city), fpmaps_conf) 
              for ifile in infiles]
    with multiprocessing.Pool(processes=ncores) as pool:
        pool.starmap(run_city, items)

    if fpmaps_dst:
        infiles = glob.glob(pmaps_dir + f"/*_{city}.h5")


if fpmaps_dst:
    origin = "raw"
    city   = "octavia"
    pmaps_dir = os.path.dirname(infiles[0]).replace(origin, city)#new_dir + f'{city}/'    pmaps_dir = os.path.dirname(infiles[0]) if fpmaps else os.path.dirname(infiles[0]).re
    infiles = sorted(glob.glob(pmaps_dir + f"/*_{city}.h5"))
    origin      = 'octavia'
#    infiles     = glob.glob(pmaps_dir + f"/*_{origin}.h5") if not fpmaps else infiles
    items_wvf = [(ifile.replace(f"{origin}", "raw"), ifile.replace(f"_{origin}.h5",f"_{origin}_DST.h5")) for ifile in infiles]
    items_S1  = [(ifile, ifile.replace(f"_{origin}.h5",f"_{origin}_DST.h5"), False) for ifile in infiles]
    items_S2  = [(ifile, ifile.replace(f"_{origin}.h5",f"_{origin}_DST.h5"), True ) for ifile in infiles]

    with multiprocessing.Pool(processes=ncores) as pool:
        pool.starmap(generate_wvf_dataset, items_wvf)
        pool.starmap(generate_pmaps_dst  , items_S1)
        pool.starmap(generate_pmaps_dst  , items_S2)
    
    dsts_wvf = []
    dsts_s1  = []
    dsts_s2  = []
    for ifile in infiles:
        dst_file = ifile.replace(f"_{origin}.h5",f"_{origin}_DST.h5")
        data_wvf, data_s1, data_s2 = pd.read_hdf(dst_file, 'WVF'), pd.read_hdf(dst_file, 'S1'), pd.read_hdf(dst_file, 'S2')
        dsts_wvf.append(data_wvf)
        dsts_s1 .append(data_s1)
        dsts_s2 .append(data_s2)
        os.system(f'rm {dst_file}')
    
    out_dst = pmaps_dir + f'/DST_run_{run_number}_pmaps.h5'
    pd.concat(dsts_wvf).to_hdf(out_dst, 'WVF')
    pd.concat(dsts_s1 ).to_hdf(out_dst, 'S1')
    pd.concat(dsts_s2 ).to_hdf(out_dst, 'S2')


