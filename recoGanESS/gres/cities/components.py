import math

import numpy  as np
import tables as tb
import pandas as pd

from typing import Callable, List, Generator
from scipy.ndimage import uniform_filter1d

from invisible_cities.dataflow                 import dataflow as  fl

from invisible_cities.  cities.components      import signal_finder
from invisible_cities.  cities.components      import check_nonempty_indices 
from invisible_cities.  cities.components      import check_empty_pmap 
from invisible_cities.  cities.components      import WfType
from invisible_cities.  cities.components      import get_event_info
from invisible_cities.  cities.components      import get_run_number
from invisible_cities.  cities.components      import get_pmt_wfs
from invisible_cities.  cities.components      import get_trigger_info
from invisible_cities.  cities.components      import check_lengths
from invisible_cities.  cities.components      import compute_z_and_dt
from invisible_cities.  cities.components      import create_timestamp

from invisible_cities.      io.event_filter_io import event_filter_writer
from invisible_cities.      io.pmaps_io        import pmap_writer
from invisible_cities.      io                 import mcinfo_io

from invisible_cities.  detsim.sensor_utils    import trigger_times
from invisible_cities.   types.ic_types        import minmax
from invisible_cities.    reco.xy_algorithms   import corona
from invisible_cities.core                     import system_of_units as units
from invisible_cities.reco                     import calib_sensors_functions as csf
from invisible_cities.reco                     import peak_functions       as pkf_ic

from .. detsim              import  buffer_functions as bf
from .. io     .rwf_io      import  buffer_writer
from .. reco                import  peak_functions   as pkf
from .. evm    .containers  import  SensorData
from .. evm    .event_model import  KrEvent
from .. database            import  load_db

def calculate_and_save_buffers(buffer_length    : float        ,
                               max_time         : int          ,
                               pre_trigger      : float        ,
                               pmt_wid          : float        ,
                               trigger_threshold: int          ,
                               h5out            : tb.File      ,
                               run_number       : int          ,
                               npmt             : int          ,
                               nsamp_pmt        : int          ,
                               order_sensors    : Callable=None):

    find_signal       = fl.map(signal_finder(buffer_length, pmt_wid,
                                             trigger_threshold     ),
                               args = "pmt_bin_wfs"                 ,
                               out  = "pulses"                      )

    filter_events_signal = fl.map(lambda x: len(x) > 0,
                                  args= 'pulses',
                                  out = 'passed_signal')
    events_passed_signal = fl.count_filter(bool, args='passed_signal')
    write_signal_filter  = fl.sink(event_filter_writer(h5out, "signal"),
                                   args=('event_number', 'passed_signal'))

    event_times       = fl.map(trigger_times                             ,
                               args = ("pulses", "timestamp", "pmt_bins"),
                               out  = "evt_times"                        )

    calculate_buffers = fl.map(bf.buffer_calculator(buffer_length, pre_trigger, pmt_wid),
                               args = ("pulses",
                                       "pmt_bins" ,  "pmt_bin_wfs")        ,
                               out  = "buffers"                            )

    saved_buffers  = "buffers" if order_sensors is None else "ordered_buffers"
    max_subevt     =  math.ceil(max_time / buffer_length)
    buffer_writer_ = fl.sink(buffer_writer( h5out
                                          , run_number = run_number
                                          , n_sens_eng = npmt
                                          , length_eng = nsamp_pmt
                                          , max_subevt = max_subevt),
                             args = ("event_number", "evt_times"  ,
                                     saved_buffers                ))

    find_signal_and_write_buffers = ( find_signal
                                    , filter_events_signal
                                    , fl.branch(write_signal_filter)
                                    , events_passed_signal.filter
                                    , event_times
                                    , calculate_buffers
                                    , order_sensors
                                    , fl.branch(buffer_writer_))

    # Filter out order_sensors if it is not set
    buffer_definition = fl.pipe(*filter(None, find_signal_and_write_buffers))
    return buffer_definition

def rebin_pmts(rebin_stride):
    def rebin_pmts(rwf):
        rebinned_wfs = rwf
        if rebin_stride > 1:
            # dummy data for times and widths
            times     = np.zeros(rwf.shape[1])
            widths    = times
            waveforms = rwf
            _, _, rebinned_wfs = pkf_ic.rebin_times_and_waveforms(times, widths, waveforms, rebin_stride=rebin_stride)
        return rebinned_wfs
    return rebin_pmts


def pmts_sum(rwfs):
    return rwfs.sum(axis=0)

def calibrate_pmts(dbfile, run_number, n_baseline, n_MAU, thr_MAU, pedestal_function=np.mean):
    DataPMT    = load_db.DataPMT(dbfile, run_number = run_number)
    adc_to_pes = np.abs(DataPMT.adc_to_pes.values)
    adc_to_pes = adc_to_pes[adc_to_pes > 0]
    def calibrate_pmts(wf):# -> CCwfs:
        cwf = pedestal_function(wf[:, :n_baseline]) - wf ### Change pulse polarity
        return calibrate_pmts_maw(cwf,
                                  adc_to_pes = adc_to_pes,
                                  n_maw      = n_MAU,
                                  thr_maw    = thr_MAU)
    return calibrate_pmts

#def calibrate_pmts(dbfile, run_number, n_baseline, n_MAU, thr_MAU, pedestal_function=np.mean):
#    DataPMT    = load_db.DataPMT(dbfile, run_number = run_number)
#    adc_to_pes = np.abs(DataPMT.adc_to_pes.values)
#    adc_to_pes = adc_to_pes[adc_to_pes > 0]
#    def calibrate_pmts(wf):# -> CCwfs:
#        cwf = pedestal_function(wf[:, :n_baseline]) - wf ### Change pulse polarity
#        return csf.calibrate_pmts(cwf,
#                                  adc_to_pes = adc_to_pes,
#                                  n_MAU      = n_MAU,
#                                  thr_MAU    = thr_MAU)
#    return calibrate_pmts

def calibrate_pmts_maw(cwfs, adc_to_pes, n_maw=100, thr_maw=3):
    """
    This function is called for PMT waveforms that have
    already been baseline restored and pedestal subtracted.
    It computes the calibrated waveforms and its sensor sum.
    It also computes the calibrated waveforms and sensor
    sum for elements of the waveforms above some value
    (thr_maw) over a MAW that follows the waveform. These
    are useful to suppress oscillatory noise and thus can
    be applied for S1 searches (the calibrated version
    without the MAW should be applied for S2 searches).
    """

    # ccwfs stands for calibrated corrected waveforms
    maw         = uniform_filter1d(cwfs, n_maw, axis=1, mode='constant')
    ccwfs       = csf.calibrate_wfs(cwfs, adc_to_pes)
    ccwfs_maw   = np.where(maw >= thr_maw, ccwfs, 0)

    cwf_sum     = np.sum(ccwfs    , axis=0)
    cwf_sum_maw = np.sum(ccwfs_maw, axis=0)
    return ccwfs, ccwfs_maw, cwf_sum, cwf_sum_maw, maw

def compute_and_write_pmaps(detector_db, run_number, pmt_samp_wid,
                            s1_lmax, s1_lmin, s1_rebin_stride, s1_stride, s1_tmax, s1_tmin,
                            s2_lmax, s2_lmin, s2_rebin_stride, s2_stride, s2_tmax, s2_tmin,
                            h5out):

    # Filter events without signal over threshold
    indices_pass    = fl.map(check_nonempty_indices,
                             args = ("s1_indices", "s2_indices"),
                             out = "indices_pass")
    empty_indices   = fl.count_filter(bool, args = "indices_pass")

    # Build the PMap
    compute_pmap     = fl.map(build_pmap(detector_db, run_number, pmt_samp_wid,
                                         s1_lmax, s1_lmin, s1_rebin_stride, s1_stride, s1_tmax, s1_tmin,
                                         s2_lmax, s2_lmin, s2_rebin_stride, s2_stride, s2_tmax, s2_tmin),
                              args = ("cwfs", "s1_indices", "s2_indices"),
                              out  = "pmap")

    # Filter events with zero peaks
    pmaps_pass      = fl.map(check_empty_pmap, args = "pmap", out = "pmaps_pass")
    empty_pmaps     = fl.count_filter(bool, args = "pmaps_pass")

    # Define writers...
    write_pmap_         = pmap_writer        (h5out,              )
    write_indx_filter_  = event_filter_writer(h5out, "s12_indices")
    write_pmap_filter_  = event_filter_writer(h5out, "empty_pmap" )

    # ... and make them sinks
    write_pmap         = fl.sink(write_pmap_        , args=(        "pmap", "event_number"))
    write_indx_filter  = fl.sink(write_indx_filter_ , args=("event_number", "indices_pass"))
    write_pmap_filter  = fl.sink(write_pmap_filter_ , args=("event_number",   "pmaps_pass"))

    fn_list = (indices_pass,
               fl.branch(write_indx_filter),
               empty_indices.filter,
               compute_pmap,
               pmaps_pass,
               fl.branch(write_pmap_filter),
               empty_pmaps.filter,
               fl.branch(write_pmap))

    # Filter out simp_rwf_to_cal if it is not set
    compute_pmaps = fl.pipe(*filter(None, fn_list))

    return compute_pmaps, empty_indices, empty_pmaps

def build_pmap(detector_db, run_number, pmt_samp_wid,
               s1_lmax, s1_lmin, s1_rebin_stride, s1_stride, s1_tmax, s1_tmin,
               s2_lmax, s2_lmin, s2_rebin_stride, s2_stride, s2_tmax, s2_tmin):
    s1_params = dict(time        = minmax(min = s1_tmin,
                                          max = s1_tmax),
                    length       = minmax(min = s1_lmin,
                                          max = s1_lmax),
                    stride       = s1_stride,
                    rebin_stride = s1_rebin_stride)

    s2_params = dict(time        = minmax(min = s2_tmin,
                                          max = s2_tmax),
                    length       = minmax(min = s2_lmin,
                                          max = s2_lmax),
                    stride       = s2_stride,
                    rebin_stride = s2_rebin_stride)

    datapmt = load_db.DataPMT(detector_db, run_number)
    pmt_ids = datapmt.SensorID[datapmt.Active.astype(bool)].values

    def build_pmap(ccwf, s1_indx, s2_indx): # -> PMap
        return pkf.get_pmap(ccwf, s1_indx, s2_indx,
                            s1_params, s2_params, pmt_ids,
                            pmt_samp_wid)

    return build_pmap

def wf_from_files(paths, wf_type):
    for path in paths:
        with tb.open_file(path, "r") as h5in:
            try:
                event_info  = get_event_info  (h5in)
                run_number  = get_run_number  (h5in)
                pmt_wfs     = get_pmt_wfs     (h5in, wf_type)
                (trg_type ,
                 trg_chann) = get_trigger_info(h5in)
            except tb.exceptions.NoSuchNodeError:
                continue

            check_lengths(pmt_wfs, event_info, trg_type, trg_chann)

            for pmt, evtinfo, trtype, trchann in zip(pmt_wfs, event_info, trg_type, trg_chann):
                event_number, timestamp         = evtinfo.fetch_all_fields()
                if trtype  is not None: trtype  = trtype .fetch_all_fields()[0]

                yield dict(pmt=pmt, run_number=run_number,
                           event_number=event_number, timestamp=timestamp,
                           trigger_type=trtype, trigger_channels=trchann)


def build_pointlike_event(dbfile, run_number, drift_v, reco):
    datapmt   = load_db.DataPMT(dbfile, run_number)
    pmt_xs    = datapmt.X.values
    pmt_ys    = datapmt.Y.values
    pmt_xys   = np.stack((pmt_xs, pmt_ys), axis=1)

    def build_pointlike_event(pmap, selector_output, event_number, timestamp):
        evt = KrEvent(event_number, timestamp * 1e-3)

        evt.nS1 = 0
        for passed, peak in zip(selector_output.s1_peaks, pmap.s1s):
            if not passed: continue

            evt.nS1 += 1
            evt.S1w.append(peak.width)
            evt.S1h.append(peak.height)
            evt.S1e.append(peak.total_energy)
            evt.S1t.append(peak.time_at_max_energy)

        evt.nS2 = 0

        for passed, peak in zip(selector_output.s2_peaks, pmap.s2s):
            if not passed: continue

            evt.nS2 += 1
            evt.S2w.append(peak.width / units.mus)
            evt.S2h.append(peak.height)
            evt.S2e.append(peak.total_energy)
            evt.S2t.append(peak.time_at_max_energy)

            xys = pmt_xys[peak.pmts.ids]
            qs  = peak.pmts.sum_over_times

            try:
                clusters = reco(xys, qs)
            except XYRecoFail:
                c    = NNN()
                Z    = tuple(NN for _ in range(0, evt.nS1))
                DT   = tuple(NN for _ in range(0, evt.nS1))
                Zrms = NN
            else:
                c = clusters[0]
                Z, DT = compute_z_and_dt(evt.S2t[-1], evt.S1t, drift_v)
                Zrms  = peak.rms / units.mus

            evt.X    .append(c.X)  
            evt.Y    .append(c.Y)  
            evt.Xrms .append(c.Xrms)
            evt.Yrms .append(c.Yrms)
            evt.R    .append(c.R)
            evt.Phi  .append(c.Phi)
            evt.DT   .append(DT)
            evt.Z    .append(Z)
            evt.Zrms .append(Zrms)
            evt.qmax .append(max(qs))

        return evt

    return build_pointlike_event

def compute_xy_position(dbfile, run_number, **reco_params):
    # `reco_params` is the set of parameters for the corona
    # algorithm either for the full corona or for barycenter
    datapmt = load_db.DataPMT(dbfile, run_number)

    def compute_xy_position(xys, qs):
        return corona(xys, qs, datapmt, **reco_params)
    return compute_xy_position

def get_number_of_active_pmts(detector_db, run_number):
    datapmt = load_db.DataPMT(detector_db, run_number)
    return np.count_nonzero(datapmt.Active.values.astype(bool))


def mcsensors_from_file(paths     : List[str],
                        db_file   :      str ,
                        run_number:      int ,
                        rate      :    float ) -> Generator:
    """
    Loads the nexus MC sensor information into
    a pandas DataFrame using the IC function
    load_mcsensor_response_df.
    Returns info event by event as a
    generator in the structure expected by
    the dataflow.
    paths      : List of strings
                 List of input file names to be read
    db_file    : string
                 Name of detector database to be used
    run_number : int
                 Run number for database
    rate       : float
                 Rate value in base unit (ns^-1) to generate timestamps
    """

    timestamp = create_timestamp(rate)

    pmt_ids  = load_db.DataPMT(db_file, run_number).SensorID

    for file_name in paths:
        sns_resp = mcinfo_io.load_mcsensor_response_df(file_name              ,
                                                       return_raw = False     ,
                                                       db_file    = db_file   ,
                                                       run_no     = run_number)

        for evt in mcinfo_io.get_event_numbers_in_file(file_name):

            try:
                pmt_indx  = sns_resp.loc[evt].index.isin(pmt_ids)
                pmt_resp  = sns_resp.loc[evt][ pmt_indx]
            except KeyError:
                pmt_resp  = pd.DataFrame(columns=sns_resp.columns)

            yield dict(event_number = evt      ,
                       timestamp    = timestamp(evt),
                       pmt_resp     = pmt_resp)