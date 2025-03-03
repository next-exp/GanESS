"""
-----------------------------------------------------------------------
                                Theodora
-----------------------------------------------------------------------
This city processes each combination of S1 and S2 signals previously
selected as pmaps by irene within an event to produce reconstructed
pointlike energy depositions. They consist of three dimensional
coordinates with some associated global S1 and S2 properties. The city
contains a peak/event filter, which can be configured to find events
with a certain number of S1/S2 signals that satisfy certain properties.
The tasks performed are:
    - Classify peaks according to the filter.
    - Filter out events that do not satisfy the selector conditions.
    - Compute S1 properties
    - Compute S2 properties (for each possible S1)
The properties of the S1 signals are:
    - Width  (time over threshold)
    - Energy (PMT- and time-summed amplitude over threshold)
    - Height (maximum amplitude)
    - Time   (waveform time at maximum amplitude)
The properties of the S2 signals are those of the S1 signal plus:
    - X      (reconstructed x position using the barycenter algorithm)
    - Y      (reconstructed y position using the barycenter algorithm)
    - DT     (drift time, computed as S2 time - S1 time)
    - Z      (reconstructed z position using DT / drift_v)
    - Xrms   (standard dev. of the PMT signal in the x coordinate)
    - Yrms   (standard dev. of the PMT signal in the y coordinate)
    - Zrms   (standard dev. of the PMT signal in the z coordinate)
    - R      (radial coordinate from X and Y)
    - Phi    (azimuthal coordinate from X and Y)
"""
from operator import attrgetter

import tables as tb

from invisible_cities.core                   import          tbl_functions as tbl
from invisible_cities.io.         kdst_io import            kr_writer
from invisible_cities.io.run_and_event_io import run_and_event_writer
from invisible_cities.io. event_filter_io import  event_filter_writer

from invisible_cities.dataflow          import dataflow      as fl

from invisible_cities.cities.components import city
from invisible_cities.cities.components import print_every
from invisible_cities.cities.components import pmap_from_files
from invisible_cities.cities.components import peak_classifier
from invisible_cities. core .configure  import             EventRangeType
from invisible_cities. core .configure  import             OneOrManyFiles

from invisible_cities. types.symbols       import               XYReco


from .components import compute_xy_position
from .components import build_pointlike_event  as build_pointlike_event_

from typing import Optional

@city
def theodora( files_in         :  OneOrManyFiles
            , file_out         :  str
            , compression      :  str
            , event_range      :  EventRangeType
            , print_mod        :  int
            , detector_db      :  str
            , run_number       :  int
            , drift_v          :  float
            , s1_nmin          :    int, s1_nmax           :   int
            , s1_emin          :  float, s1_emax           : float
            , s1_wmin          :  float, s1_wmax           : float
            , s1_hmin          :  float, s1_hmax           : float
            , s1_ethr          :  float
            , s2_nmin          :    int, s2_nmax           :   int
            , s2_emin          :  float, s2_emax           : float
            , s2_wmin          :  float, s2_wmax           : float
            , s2_hmin          :  float, s2_hmax           : float
            , s2_ethr          :  float
            , s2_nsipmmin      :    int, s2_nsipmmax       :   int ### SiPM patch to re-use IC code
            , global_reco_algo : XYReco, global_reco_params:  dict
            , include_mc       : Optional[bool] = False
):



#(files_in, file_out, compression, event_range, print_mod, detector_db, run_number,
#             drift_v,
#             s1_nmin, s1_nmax, s1_emin, s1_emax, s1_wmin, s1_wmax, s1_hmin, s1_hmax, s1_ethr,
#             s2_nmin, s2_nmax, s2_emin, s2_emax, s2_wmin, s2_wmax, s2_hmin, s2_hmax, s2_ethr, s2_nsipmmin, s2_nsipmmax, 
#             global_reco_params=dict()):
    # global_reco_params are qth, qlm, lm_radius, new_lm_radius, msipm
    # qlm           =  0 * pes every Cluster must contain at least one PMT with charge >= qlm
    # lm_radius     = -1 * mm  by default, use overall barycenter for KrCity
    # new_lm_radius = -1 * mm  find a new cluster by calling barycenter() on pos/qs of PMTs within
    #                          new_lm_radius of new_local_maximum
    # msipm         =  1       minimum number of PMTs in a Cluster

    classify_peaks        = fl.map(peak_classifier(**locals()),
                                   args = "pmap",
                                   out  = "selector_output")

    pmap_passed           = fl.map(attrgetter("passed"), args="selector_output", out="pmap_passed")
    pmap_select           = fl.count_filter(bool, args="pmap_passed")

    reco_algo             = compute_xy_position(detector_db, run_number, global_reco_algo, **global_reco_params)
    build_pointlike_event = fl.map(build_pointlike_event_(detector_db, run_number, drift_v, reco_algo),
                                   args = ("pmap", "selector_output", "event_number", "timestamp"),
                                   out  = "pointlike_event"                                       )

    event_count_in        = fl.spy_count()
    event_count_out       = fl.spy_count()

    with tb.open_file(file_out, "w", filters = tbl.filters(compression)) as h5out:

        # Define writers...
        write_event_info      = fl.sink(run_and_event_writer(h5out                ), args=("run_number", "event_number", "timestamp"))
        write_pointlike_event = fl.sink(           kr_writer(h5out                ), args="pointlike_event")
        write_pmap_filter     = fl.sink( event_filter_writer(h5out, "s12_selector"), args=("event_number", "pmap_passed"))

        res = fl.push(source = pmap_from_files(files_in),
                    pipe   = fl.pipe(
                        fl.slice(*event_range, close_all=True),
                        print_every(print_mod)                ,
                        event_count_in       .spy             ,
                        classify_peaks                        ,
                        pmap_passed                           ,
                        fl.branch(write_pmap_filter)          ,
                        pmap_select          .filter          ,
                        event_count_out      .spy             ,
                        build_pointlike_event                 ,
                        fl.fork(write_pointlike_event         ,
                                write_event_info              )),
                    result = dict(events_in  = event_count_in .future,
                                  events_out = event_count_out.future,
                                  selection  = pmap_select    .future))

        if run_number <= 0 and include_mc:
            copy_mc_info(files_in, h5out, result.evtnum_list,
                         detector_db, run_number)

        return res
