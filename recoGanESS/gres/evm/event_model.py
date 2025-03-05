from invisible_cities.evm.event_model import Event
import numpy as np
class KrEvent(Event):
    """Represents a point-like (Krypton) event."""
    def __init__(self, event_number, event_time):
        Event.__init__(self, event_number, event_time)
        self.nS1   = -1 # number of S1 in the event
        self.S1w   = [] # widht
        self.S1h   = [] # heigth
        self.S1e   = [] # energy
        self.S1t   = [] # time

        self.nS2   = -1 # number of S2s in the event
        self.S2w   = []
        self.S2h   = []
        self.S2e   = []
        self.S2t   = [] # time
        self.qmax  = []

        self.DT    = [] # drift time
        self.Z     = [] # Position (x,y,z,R,phi)
        self.X     = []
        self.Y     = []
        self.R     = []
        self.Phi   = []
        self.Xrms  = [] # error in position
        self.Yrms  = []
        self.Zrms  = []

    def fill_defaults(self):
        if self.nS1 == 0:
            for attribute in ["w", "h", "e", "t"]:
                setattr(self, "S1" + attribute, [np.nan])

        if self.nS2 == 0:
            for attribute in ["w", "h", "e", "t"]:
                setattr(self, "S2" + attribute, [np.nan])

            self.qmax  = 0
            for attribute in ["X", "Y", "R", "Phi", "Xrms", "Yrms", "Zrms"]:
                setattr(self, attribute, [np.nan])

        if not self.nS1 or not self.nS2:
            for attribute in ["Z", "DT"]:
                setattr(self, attribute, [[np.nan] * max(self.nS1, 1)] * max(self.nS2, 1))

    def store(self, table):
        row = table.row

        s1_peaks = range(int(self.nS1)) if self.nS1 else [-1]
        s2_peaks = range(int(self.nS2)) if self.nS2 else [-1]
        self.fill_defaults()

        for i in s1_peaks:
            for j in s2_peaks:
                row["event"  ] = self.event
                row["time"   ] = self.time
                row["s1_peak"] = i
                row["s2_peak"] = j
                row["nS1"    ] = self.nS1
                row["nS2"    ] = self.nS2

                row["S1w"    ] = self.S1w  [i]
                row["S1h"    ] = self.S1h  [i]
                row["S1e"    ] = self.S1e  [i]
                row["S1t"    ] = self.S1t  [i]

                row["S2w"    ] = self.S2w  [j]
                row["S2h"    ] = self.S2h  [j]
                row["S2e"    ] = self.S2e  [j]
                row["S2t"    ] = self.S2t  [j]
                row["qmax"   ] = self.qmax [j]

                row["DT"     ] = self.DT   [j][i]
                row["Z"      ] = self.Z    [j][i]
                row["Zrms"   ] = self.Zrms [j]
                row["X"      ] = self.X    [j]
                row["Y"      ] = self.Y    [j]
                row["R"      ] = self.R    [j]
                row["Phi"    ] = self.Phi  [j]
                row["Xrms"   ] = self.Xrms [j]
                row["Yrms"   ] = self.Yrms [j]
                row.append()

    def __str__(self):
        s = "{0}Event\n{0}".format("#"*20 + "\n")
        for attr in self.__dict__:
            s += "{}: {}\n".format(attr, getattr(self, attr))
        return s

    __repr__ =     __str__