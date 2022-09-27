import numpy  as np
import tables as tb
import pandas as pd
import warnings

from typing import Callable
from typing import     List
from typing import    Tuple

from               .. database.load_db   import DataPMT
from invisible_cities.io      .mcinfo_io import get_sensor_binning

def pmt_bin_width(file_name: str) -> float:
    """
    Returns PMT bin widths as set in simGanESS.
    """
    sns_bins = get_sensor_binning(file_name)
    if sns_bins.empty or np.any(sns_bins.bin_width <= 0):
        raise tb.NoSuchNodeError('No useful binning info found')
    pmt_wid  = sns_bins.bin_width[sns_bins.index.str.contains( 'Pmt')].iloc[0]
    return pmt_wid

def first_and_last_times(pmt_wfs    : pd.DataFrame,
                         pmt_binwid : float       ) -> Tuple[float, float]:
    """
    Returns the maximum and minimum time of an
    event given the type of detector.
    """
    min_time  = pmt_wfs.time.min()
    max_time  = pmt_wfs.time.max()
    max_time += pmt_binwid
    return min_time, max_time

def sensor_order(detector_db: str, run_number : int, length_pmt : int) -> Callable:
    """
    Casts the event sensor info into the correct order
    adding zeros for sensors which didn't see any signal.
    """
    pmt_ids    = DataPMT (detector_db, run_number).SensorID##[:1] ### Patch until a proper database is done
    n_pmt      = get_n_sensors(detector_db, run_number)
    pmt_shape  = (n_pmt , length_pmt )

    def ordering(sensor_order : pd.Int64Index  ,
                 sensor_resp  : np.ndarray     ,
                 sensor_shape : Tuple[int, int]) -> np.ndarray:
        sensors = np.zeros(sensor_shape, int)
        sensors[sensor_order] = sensor_resp
        return sensors

    def order_and_pad(pmt_resp    : pd.Series                          ,
                      evt_buffers : List[np.ndarray]
                      ) -> List[np.ndarray]:
        pmt_ord  = pmt_ids [ pmt_ids.isin( pmt_resp.index)].index
        return [ordering(pmt_ord , pmts , pmt_shape ) for pmts in evt_buffers]

    return order_and_pad

def get_n_sensors(detector_db: str, run_number: int) -> int:
    """Get the number of sensors for this run"""
    npmt  = DataPMT (detector_db, run_number).shape[0]  ### Patch until a proper database is done
    return npmt