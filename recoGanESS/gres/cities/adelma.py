"""
-----------------------------------------------------------------------
                                Adelma
-----------------------------------------------------------------------

 -- Sorts MC sensor info into buffers (only PMTs)

This city reads nexus Monte Carlo Sensor information and sorts the
information into data-like buffers without the addition of
electronics noise.

Uses a configured buffer length, pretrigger and threshold for the
positioning of trigger like signal. If more than one trigger is found
separated from each other by more than a buffer width the nexus event
can be split into multiple data-like triggers.
"""

import pandas as pd
import tables as tb

import warnings

from typing import Callable
from typing import     List

from invisible_cities.core                   import        system_of_units as units

from invisible_cities.io    .event_filter_io import    event_filter_writer
from invisible_cities.reco                   import          tbl_functions as tbl

from invisible_cities. dataflow   import                   dataflow as fl

from invisible_cities.cities.components import                       city
from invisible_cities.cities.components import                    collect
from invisible_cities.cities.components import               copy_mc_info
from invisible_cities.cities.components import        mcsensors_from_file
from invisible_cities.cities.components import                print_every
from invisible_cities.cities.components import                  wf_binner
from invisible_cities.cities.components import             check_max_time

from .. detsim.sensor_utils    import        pmt_bin_width
from .. detsim.sensor_utils    import first_and_last_times
from .. detsim.sensor_utils    import        get_n_sensors
from .. detsim.sensor_utils    import         sensor_order

from .  components             import calculate_and_save_buffers

@city
def adelma(files_in     , file_out   , compression      , event_range,
          print_mod    , detector_db, run_number       , max_time   ,
          buffer_length, pre_trigger, trigger_threshold, rate       ):

    max_time    = check_max_time(max_time, buffer_length)
    npmt        = get_n_sensors(detector_db, run_number)
    pmt_wid     = pmt_bin_width_safe_(files_in)
    nsamp_pmt   = int(buffer_length /  pmt_wid)

    extract_tminmax   = fl.map(first_and_last_times_(pmt_wid) ,
                               args = ("pmt_resp")            ,
                               out  = ("min_time", "max_time"))

    bin_calculation   = wf_binner(max_time)
    bin_pmt_wf        = fl.map(binning_set_width(bin_calculation, pmt_wid),
                               args = ("pmt_resp", "min_time", "max_time"),
                               out  = ("pmt_bins", "pmt_bin_wfs")         )

    order_sensors     = fl.map(sensor_order(detector_db, run_number, nsamp_pmt) ,
                               args = ("pmt_bin_wfs", "buffers") ,
                               out  = "ordered_buffers"              )

    filter_events     = fl.map(lambda x : not x.empty,
                               args = ('pmt_resp')   ,
                               out  = 'event_passed' )
    events_with_resp  = fl.count_filter(bool, args="event_passed")

    with tb.open_file(file_out, "w", filters=tbl.filters(compression)) as h5out:
        buffer_calculation = calculate_and_save_buffers( buffer_length
                                                       , max_time
                                                       , pre_trigger
                                                       , pmt_wid
                                                       , trigger_threshold
                                                       , h5out
                                                       , run_number
                                                       , npmt
                                                       , nsamp_pmt
                                                       , order_sensors)

        write_filter   = fl.sink(event_filter_writer(h5out, "detected_events"),
                                 args=("event_number", "event_passed")       )

        event_count_in = fl.spy_count()

        evtnum_collect = collect()

        result = fl.push(source = mcsensors_from_file(files_in   ,
                                                      detector_db,
                                                      run_number ,
                                                      rate       ),
                         pipe   = fl.pipe( fl.slice(*event_range, close_all=True)
                                         , event_count_in.spy
                                         , print_every(print_mod)
                                         , filter_events
                                         , fl.branch(write_filter)
                                         , events_with_resp.filter
                                         , extract_tminmax
                                         , bin_pmt_wf
                                         , buffer_calculation
                                         , "event_number"
                                         , evtnum_collect.sink),
                         result = dict(events_in   = event_count_in.future  ,
                                       events_resp = events_with_resp.future,
                                       evtnum_list = evtnum_collect.future  ))

        copy_mc_info(files_in, h5out, result.evtnum_list,
                     detector_db, run_number)

        return result


def first_and_last_times_(pmt_binwid: float):
    def get_event_time_extremes(pmt_resp : pd.DataFrame):
        return first_and_last_times(pmt_resp, pmt_binwid)
    return get_event_time_extremes


def binning_set_width(binning_fnc: Callable, bin_width: float):
    def bin_calculation_(sensors: pd.DataFrame,
                         t_min  : float       ,
                         t_max  : float       ):
        return binning_fnc(sensors, bin_width, t_min, t_max)
    return bin_calculation_


def pmt_bin_width_safe_(files_in: List[str]):
    for fn in files_in:
        try:
            pmt_wid = pmt_bin_width(fn)
            return pmt_wid
        except (tb.HDF5ExtError, tb.NoSuchNodeError) as e:
            warnings.warn(f' no useful bin widths: {0}'.format(e), UserWarning)
            continue
    return 25 * units.ns