"""
-----------------------------------------------------------------------
                                 Thekla
-----------------------------------------------------------------------

From ancient Greek ‘Υπατια: highest, supreme.

This city reads true waveforms from detsim and compute pmaps from them
without simulating the electronics. This includes:
    - Rebin waveforms
    - Produce a PMT-summed waveform.
    - Apply a threshold to the PMT-summed waveform.
    - Find pulses in the PMT-summed waveform.
    - Build the PMap object.
"""
import numpy  as np
import tables as tb

from functools import partial

from invisible_cities. reco                  import sensor_functions     as sf
from invisible_cities. reco                  import tbl_functions        as tbl
from invisible_cities. reco                  import peak_functions       as pkf
from invisible_cities. core                  import system_of_units      as units
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
from invisible_cities.cities.components import get_number_of_active_pmts

from .components import sensor_data
from .components import compute_and_write_pmaps
from .components import wf_from_files

@city
def thekla(files_in, file_out, compression, event_range, print_mod, detector_db, run_number,
           pmt_wfs_rebin, pmt_pe_rms,
           s1_lmin, s1_lmax, s1_tmin, s1_tmax, s1_rebin_stride, s1_stride, thr_csum_s1,
           s2_lmin, s2_lmax, s2_tmin, s2_tmax, s2_rebin_stride, s2_stride, thr_csum_s2,
           pmt_samp_wid=100*units.ns):
    #### Define data transformations
    sd = sensor_data(files_in[0], WfType.mcrd)

    # Raw WaveForm to Corrected WaveForm
    mcrd_to_rwf      = fl.map(rebin_pmts(pmt_wfs_rebin),
                              args = "pmt",
                              out  = "rwf")

    # Add single pe fluctuation to pmts
    simulate_pmt = fl.map(partial(sf.charge_fluctuation, single_pe_rms=pmt_pe_rms),
                          args = "rwf",
                          out = "ccwfs")

    # Compute pmt sum
    pmt_sum          = fl.map(pmts_sum, args = 'ccwfs',
                              out  = 'pmt')

    # Find where waveform is above threshold
    zero_suppress    = fl.map(zero_suppress_wfs(thr_csum_s1, thr_csum_s2),
                              args = ("pmt", "pmt"),
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

        result = push(source = wf_from_files(files_in, WfType.mcrd),
                      pipe   = pipe(fl.slice(*event_range, close_all=True),
                                    print_every(print_mod),
                                    event_count_in.spy,
                                    mcrd_to_rwf,
                                    simulate_pmt,
                                    pmt_sum,
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


def rebin_pmts(rebin_stride):
    def rebin_pmts(rwf):
        rebinned_wfs = rwf
        if rebin_stride > 1:
            # dummy data for times and widths
            times     = np.zeros(rwf.shape[1])
            widths    = times
            waveforms = rwf
            _, _, rebinned_wfs = pkf.rebin_times_and_waveforms(times, widths, waveforms, rebin_stride=rebin_stride)
        return rebinned_wfs
    return rebin_pmts


def pmts_sum(rwfs):
    return rwfs.sum(axis=0)