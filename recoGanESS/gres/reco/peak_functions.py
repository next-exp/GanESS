
from invisible_cities.reco.peak_functions import find_peaks
from invisible_cities. evm.pmaps          import S1
from invisible_cities. evm.pmaps          import S2
from invisible_cities. evm.pmaps          import PMap


def get_pmap(ccwf, s1_indx, s2_indx,
             s1_params, s2_params, pmt_ids,
             pmt_samp_wid):
    return PMap(find_peaks(ccwf, s1_indx, Pk=S1, pmt_ids=pmt_ids,
                           pmt_samp_wid=pmt_samp_wid,
                           **s1_params),
                find_peaks(ccwf, s2_indx, Pk=S2, pmt_ids=pmt_ids,
                           pmt_samp_wid  = pmt_samp_wid,
                           **s2_params))