import tables as tb
import numpy  as np
from functools import  partial
from typing    import Callable
from typing    import     List
from typing    import Optional
from typing    import    Tuple

from invisible_cities.evm.nh5         import           MCEventMap

from invisible_cities.io.table_io         import make_table
from invisible_cities.io.run_and_event_io import run_and_event_writer
from invisible_cities.io.rwf_io           import rwf_writer
from invisible_cities.io.rwf_io           import ic_event_number_base

def buffer_writer(h5out, *,
                  run_number :          int           ,
                  n_sens_eng :          int           ,
                  length_eng :          int           ,
                  group_name : Optional[str] =    None,
                  compression: Optional[str] = 'ZLIB4',
                  max_subevt : Optional[int] =      10
                  ) -> Callable[[int, List, List], None]:
    """
    Generalised buffer writer which defines a raw waveform writer
    for each type of sensor as well as an event info writer.
    Each call gives a list of 'triggers' to be written as
    separate events in the output.
    parameters
    ----------
    run_number  : int
                  Run number to be saved in runInfo.
    n_sens_eng  : int
                  Number of sensors in the energy plane.
    length_eng  : int
                  Number of samples per waveform for energy plane.
    group_name  : Optional[str] default None
                  Group name within root where waveforms to be saved.
                  Default directly in root
    compression : Optional[str] default 'ZLIB4'
                  Compression level for output file.
    returns
    -------
    write_buffers : Callable
                    A function which takes event information
                    for the tracking and energy planes and
                    the event timestamps and saves to file.
    """

    eng_writer = rwf_writer(h5out,
                            group_name      =  group_name,
                            compression     = compression,
                            table_name      =     'pmtrd',
                            n_sensors       =  n_sens_eng,
                            waveform_length =  length_eng)

    run_and_event = partial(run_and_event_writer(h5out                  ,
                                                 compression=compression),
                            run_number = run_number                      )

    nexus_map = make_table(h5out, 'Run', 'eventMap', MCEventMap,
                           "event & nexus evt for each index", compression)

    evt_number_generator = ic_event_number_base(max_subevt)
    def write_buffers(nexus_evt :        int ,
                      timestamps: List[  int],
                      events    : List[np.ndarray]) -> None:
        """
        Write run info and event waveforms to file.
        parameters
        ----------
        nexus_evt  :  int
                     Event number from MC output file.
        timestamps : List[int]
                     List of event times
        events     : List[np.ndarray]
                     List of np.ndarray containing the energy info for each identified 'trigger'.
        """
        event_number_base = evt_number_generator(nexus_evt)
        for i, (t_stamp, eng) in enumerate(zip(timestamps, events)):
            ## Save event number and log nexus event number.
            event_number = event_number_base + i
            run_and_event(event_number=event_number, timestamp=t_stamp)
            mrow = nexus_map.row
            mrow["evt_number"] = event_number
            mrow[ "nexus_evt"] = nexus_evt
            mrow.append()
            ##
            eng_writer(eng)

    return write_buffers