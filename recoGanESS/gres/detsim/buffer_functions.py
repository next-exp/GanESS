import numpy  as np
import pandas as pd

from typing import Callable
from typing import     List
from typing import    Union
from typing import    Tuple

from invisible_cities. detsim .buffer_functions import pad_safe


def buffer_calculator(buffer_len: float, pre_trigger: float, pmt_binwid: float) -> Callable:
    """
    Calculates the output buffers for all sensors
    based on a configured buffer length and pretrigger.

    Parameters
    ----------
    buffer_len  : float
                  Length of buffer expected in mus
    pre_trigger : float
                  Time in buffer before identified signal in mus
    pmt_binwid  : float
                  Width in mus of PMT sample integration
    """
    pmt_buffer_samples  = int(buffer_len  //  pmt_binwid)
    pmt_pretrg_base     = int(pre_trigger //  pmt_binwid)
    pmt_postrg_base     = pmt_buffer_samples - pmt_pretrg_base

    def generate_slice(trigger    : int       ,
                       pmt_bins   : np.ndarray,
                       pmt_charge : np.ndarray) -> np.ndarray:

        npmt_bin  = len(pmt_bins)

        bin_corr   = pmt_bins[trigger] // pmt_binwid
        pmt_pretrg = pmt_pretrg_base + int(bin_corr)
        pmt_postrg = pmt_postrg_base - int(bin_corr)

        pmt_pre    = 0       , trigger - pmt_pretrg
        pmt_pos    = npmt_bin, trigger + pmt_postrg
        pmt_slice  = slice(max(pmt_pre), min(pmt_pos))
        pmt_pad    = -min(pmt_pre), max(0, pmt_pos[1] - npmt_bin)

        return pad_safe( pmt_charge[:,  pmt_slice],  pmt_pad)

    def position_signal(triggers   :       List                  ,
                        pmt_bins   : np.ndarray                  ,
                        pmt_charge : Union[pd.Series, np.ndarray]
                        ) -> List[np.ndarray]:
        """
        Synchronises the SiPMs and PMTs for each identified
        trigger and calls the padding function to fill with
        zeros where necessary.
        """
        if isinstance(pmt_charge, pd.Series):
            pmt_q  = np.asarray(pmt_charge.tolist())
        else:
            pmt_q  =  pmt_charge
        return [generate_slice(trigger, pmt_bins, pmt_q) for trigger in triggers]
    return position_signal