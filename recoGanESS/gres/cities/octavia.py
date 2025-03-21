"""
-----------------------------------------------------------------------
                                 Octavia
-----------------------------------------------------------------------
This city finds the signal pulses within the waveforms produced by the
detector.
This includes a number of tasks:
    - Removes the PMTs baseline, calibrates them and produced a PMT-summed waveform.
    - Apply a threshold to the PMT-summed waveform.
    - Find pulses in the PMT-summed waveform.
    - Build the PMap object.
"""
import tables as tb

from invisible_cities.core                   import          tbl_functions as tbl
from invisible_cities.core                   import        system_of_units as units
from invisible_cities. io  .run_and_event_io import run_and_event_writer
from invisible_cities. io  .      trigger_io import       trigger_writer

from invisible_cities. dataflow            import dataflow as fl
from invisible_cities. dataflow.dataflow   import push
from invisible_cities. dataflow.dataflow   import pipe
from invisible_cities. dataflow.dataflow   import sink

from invisible_cities.cities.components import city
from invisible_cities.cities.components import print_every
from invisible_cities.cities.components import collect
from invisible_cities.cities.components import copy_mc_info
from invisible_cities.cities.components import zero_suppress_wfs
from invisible_cities.cities.components import WfType

from invisible_cities. core .configure  import             EventRangeType
from invisible_cities. core .configure  import             OneOrManyFiles

from .components import compute_and_write_pmaps
from .components import wf_from_files
from .components import get_number_of_active_pmts
from .components import calibrate_pmts

@city
def octavia( files_in        : OneOrManyFiles
           , file_out        : str
           , compression     : str
           , event_range     : EventRangeType
           , print_mod       : int
           , detector_db     : str
           , run_number      : int
           , n_baseline      : int
           , n_maw           : int
           , thr_maw         : float
           , s1_lmin         : int  , s1_lmax      : int
           , s1_tmin         : float, s1_tmax      : float
           , s1_rebin_stride : int  , s1_stride    : int
           , thr_csum_s1     : float
           , s2_lmin         : int  , s2_lmax      : int
           , s2_tmin         : float, s2_tmax      : float
           , s2_rebin_stride : int  , s2_stride    : int
           , thr_csum_s2     : float 
           , pmt_samp_wid    : float
           ):

    #### Define data transformations
    # WaveForm to Calibrated WaveForm
    rwf_to_cwf      = fl.map(calibrate_pmts(detector_db, run_number, n_baseline, n_maw, thr_maw),
                              args = "pmt",
                              out  = ("ccwfs", "ccwfs_maw", "cwf_sum", "cwf_sum_maw"))

    # Find where waveform is above threshold
    zero_suppress    = fl.map(zero_suppress_wfs(thr_csum_s1, thr_csum_s2),
                              args = ("cwf_sum", "cwf_sum_maw"),
                              out  = ("s1_indices", "s2_indices", "s2_energies"))


    event_count_in  = fl.spy_count()
    event_count_out = fl.spy_count()

    evtnum_collect  = collect()

    with tb.open_file(file_out, "w", filters = tbl.filters(compression)) as h5out:

        # Define writers...
        write_event_info_   = run_and_event_writer(h5out)
        write_trigger_info_ = trigger_writer      (h5out, get_number_of_active_pmts(detector_db, run_number))

        # ... and make them sinks
        write_event_info   = sink(write_event_info_  , args=(   "run_number",     "event_number", "timestamp"   ))
        write_trigger_info = sink(write_trigger_info_, args=( "trigger_type", "trigger_channels"                ))


        compute_pmaps, empty_indices, empty_pmaps = compute_and_write_pmaps(
                                             detector_db, run_number, pmt_samp_wid,
                                             s1_lmax, s1_lmin, s1_rebin_stride, s1_stride, s1_tmax, s1_tmin,
                                             s2_lmax, s2_lmin, s2_rebin_stride, s2_stride, s2_tmax, s2_tmin,
                                             h5out)

        result = push(source = wf_from_files(files_in, WfType.rwf),
                      pipe   = pipe(fl.slice(*event_range, close_all=True),
                                    print_every(print_mod),
                                    event_count_in.spy,
                                    rwf_to_cwf,
                                    zero_suppress,
                                    compute_pmaps,
                                    event_count_out.spy,
                                    fl.branch("event_number", evtnum_collect.sink),
                                    fl.fork(write_event_info,
                                            write_trigger_info)),
                      result = dict(events_in   = event_count_in .future,
                                    events_out  = event_count_out.future,
                                    evtnum_list = evtnum_collect .future,
                                    over_thr    = empty_indices  .future,
                                    full_pmap   = empty_pmaps    .future))

        if run_number <= 0:
            copy_mc_info(files_in, h5out, result.evtnum_list,
                         detector_db, run_number)

        return result