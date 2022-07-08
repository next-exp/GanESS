class S12SelectorOutput:
    """
    Class devoted to hold the output of the S12Selector.
    It contains:
        - passed  : a boolean flag indicating whether the event as
                    a whole has passed the filter.
        - s1_peaks: a sequence with a boolean flag for each peak
                    indicating whether the peak has been selected
                    or not.
        - s2_peaks: same as s1_peaks.
    The __and__ (&) and __or__ (|) methods, allow the user to combine
    two outputs.
    """
    def __init__(self,
                 passed   : bool,
                 s1_peaks : Sequence[bool],
                 s2_peaks : Sequence[bool]):
        self.passed   = passed
        self.s1_peaks = s1_peaks
        self.s2_peaks = s2_peaks

    def __and__(self, other : "S12SelectorOutput") -> "S12SelectorOutput":
        if (len(self.s1_peaks) != len(other.s1_peaks) or
            len(self.s2_peaks) != len(other.s2_peaks)):
            raise ValueError("Cannot and lists of different length")
        s1_peaks = tuple(map(np.logical_and, self.s1_peaks, other.s1_peaks))
        s2_peaks = tuple(map(np.logical_and, self.s2_peaks, other.s2_peaks))
        passed   = self.passed and other.passed
        return S12SelectorOutput(passed, s1_peaks, s2_peaks)


    def __or__(self, other : "S12SelectorOutput") -> "S12SelectorOutput":
        if (len(self.s1_peaks) != len(other.s1_peaks) or
            len(self.s2_peaks) != len(other.s2_peaks)):
            raise ValueError("Cannot or lists of different length")
        s1_peaks = tuple(map(np.logical_or, self.s1_peaks, other.s1_peaks))
        s2_peaks = tuple(map(np.logical_or, self.s2_peaks, other.s2_peaks))
        passed   = self.passed or other.passed
        return S12SelectorOutput(passed, s1_peaks, s2_peaks)

    def __str__(self):
        return dedent("""
                  S12SelectorOutput:
                      Passed  : {self.passed}
                      s1_peaks: {self.s1_peaks}
                      s2_peaks: {self.s2_peaks}
                      """.format(self = self))

    __repr__ = __str__


class S12Selector:
    def __init__(self, **kwds):
        conf = Namespace(**kwds)
        self.s1n = minmax(conf.s1_nmin, conf.s1_nmax)
        self.s1e = minmax(conf.s1_emin, conf.s1_emax)
        self.s1w = minmax(conf.s1_wmin, conf.s1_wmax)
        self.s1h = minmax(conf.s1_hmin, conf.s1_hmax)
        self.s1_ethr = conf.s1_ethr

        self.s2n = minmax(conf.s2_nmin, conf.s2_nmax)
        self.s2e = minmax(conf.s2_emin, conf.s2_emax)
        self.s2w = minmax(conf.s2_wmin, conf.s2_wmax)
        self.s2h = minmax(conf.s2_hmin, conf.s2_hmax)
        self.nsi = minmax(conf.s2_nsipmmin, conf.s2_nsipmmax)
        self.s2_ethr = conf.s2_ethr

    @staticmethod
    def valid_peak(peak   : _Peak,
                   thr    : float,
                   energy : minmax,
                   width  : minmax,
                   height : minmax,
                   nsipm  : minmax = None) -> bool:
        """Returns True if the peak energy, width and height
        is contained in the minmax defined by energy, width and
        height."""
        f1 = energy.contains(peak.energy_above_threshold(thr))
        f2 = width .contains(peak. width_above_threshold(thr))
        f3 = height.contains(peak.height)
        f4 = True
        if nsipm:
            f4 = nsipm.contains(peak.sipms.ids.size)
        return f1 and f2 and f3 and f4

    @staticmethod
    def select_valid_peaks(peaks  : Sequence[_Peak],
                           thr    : float,
                           energy : minmax,
                           width  : minmax,
                           height : minmax,
                           nsipm  : minmax = None) ->Sequence[bool]:
        """
        Takes a sequence of peaks and returns a sequence
        with the outcome of the filter for each peak
        """
        peak_is_valid = partial(S12Selector.valid_peak,
                                thr    = thr,
                                energy = energy,
                                width  = width,
                                height = height,
                                nsipm  = nsipm)
        valid_peaks   = tuple(map(peak_is_valid, peaks))
        return valid_peaks

    def select_s1(self, s1s : Sequence[_Peak]) -> Sequence[bool]:
        """
        Takes a sequence of S1s and returns a sequence with the
        outcome of the filter for each peak.
        """
        passed = self.select_valid_peaks(s1s,
                                         self.s1_ethr,
                                         self.s1e,
                                         self.s1w,
                                         self.s1h)
        return passed

    def select_s2(self, s2s : Sequence[_Peak]) -> Sequence[bool]:
        """
        Takes a sequence of S2s and returns a sequence with the
        outcome of the filter for each peak
        """
        passed = self.select_valid_peaks(s2s,
                                         self.s2_ethr,
                                         self.s2e,
                                         self.s2w,
                                         self.s2h)
        return passed


    def __str__(self):
        return dedent("""
                   S12_selector:
                       s1n     = {self.s1n}
                       s1e     = {self.s1e} pes
                       s1w     = {self.s1w} ns
                       s1h     = {self.s1h} pes
                       s1_ethr = {self.s1_ethr} pes
                       s2n     = {self.s2n}
                       s2e     = {self.s2e} pes
                       s2w     = {self.s2w} ns
                       s2h     = {self.s2h} pes
                       s2_ethr = {self.s2_ethr} pes""".format(self = self))

    __repr__ = __str__