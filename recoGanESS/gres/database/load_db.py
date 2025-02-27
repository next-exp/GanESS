import sqlite3
import numpy as np
import pandas as pd
import os
from operator  import itemgetter
from functools import lru_cache


class DetDB:
    gap     = os.environ['GRESDIR'] + '/database/localdb.GaP.sqlite3'
    gap_ext = os.environ['GRESDIR'] + '/database/localdb.GaP_ext.sqlite3'

def tmap(*args):
    return tuple(map(*args))

def get_db(db):
    return getattr(DetDB, db, db)

# Run to take always the same calibration constant, etc for MC files
runNumberForMC = 0

@lru_cache(maxsize=10)
def DataPMT(db_file, run_number=1e5):
    if run_number == 0:
        run_number = runNumberForMC
    conn = sqlite3.connect(get_db(db_file))
    sql = '''select pos.SensorID, map.ElecID "ChannelID", Label "PmtID",
case when msk.SensorID is NULL then 1 else 0 end "Active",
X, Y, coeff_blr, coeff_c, abs(Centroid) "adc_to_pes", noise_rms, Sigma
from ChannelPosition as pos INNER JOIN ChannelMapping
as map ON pos.SensorID = map.SensorID LEFT JOIN
(select * from PmtNoiseRms where MinRun <= {0} and (MaxRun >= {0} or MaxRun is NULL))
as noise on map.ElecID = noise.ElecID LEFT JOIN
(select * from ChannelMask where MinRun <= {0} and {0} <= MaxRun)
as msk ON pos.SensorID = msk.SensorID LEFT JOIN
(select * from ChannelGain where  MinRun <= {0} and {0} <= MaxRun)
as gain ON pos.SensorID = gain.SensorID LEFT JOIN
(select * from PmtBlr where MinRun <= {0} and (MaxRun >= {0} or MaxRun is NULL))
as blr ON map.ElecID = blr.ElecID
where pos.SensorID < 100
and pos.MinRun <= {0} and {0} <= pos.MaxRun
and map.MinRun <= {0} and {0} <= map.MaxRun
and pos.Label LIKE 'PMT%'
order by Active desc, pos.SensorID
'''.format(abs(run_number))
    data = pd.read_sql_query(sql, conn)
    data.fillna(0, inplace=True)
    conn.close()
    return data


@lru_cache(maxsize=10)
def DetectorGeo(db_file):
    conn = sqlite3.connect(get_db(db_file))
    sql = 'select * from DetectorGeo'
    data = pd.read_sql_query(sql, conn)
    conn.close()
    return data


@lru_cache(maxsize=10)
def PMTLowFrequencyNoise(db_file, run_number=1e5):
    conn = sqlite3.connect(get_db(db_file))
    cursor = conn.cursor()

    sqlmapping = '''select SensorID, FEBox from PMTFEMapping
    where MinRun <= {0} and (MaxRun >= {0} or MaxRun is NULL)
    order by SensorID;'''.format(abs(run_number))
    mapping = pd.read_sql_query(sqlmapping, conn)

    ## Now get the frequencies and magnitudes (and ?) for each box
    ## Number of boxes can be different for different detectors, so we need to
    ## find out how many columns are there in the table
    sql = '''PRAGMA table_info('PMTFELowFrequencyNoise');'''
    schema = pd.read_sql_query(sql, conn)
    colnames = schema.name[schema.name.str.contains("FE")].values
    colnames = ', '.join(colnames)

    sqlmagnitudes = '''select Frequency, {0}
    from PMTFELowFrequencyNoise where MinRun <= {1}
    and (MaxRun >= {1} or MaxRun is NULL)'''.format(colnames, abs(run_number))
    frequencies = pd.read_sql_query(sqlmagnitudes, conn)

    return mapping, frequencies.values
