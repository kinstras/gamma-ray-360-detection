# gamma-ray-360-detection
Several modern sensitive particle detectors utilize the technique proposed by H.O.
Anger described in detail in [1], in which the scintillations resulting from a gamma-
ray photon are detected by means of a two-dimensional silicon photomultiplier array
(2d Silicon PhotoMultiplier array), placed at the top or bottom of a [1].
In this thesis, the study of spatial directional emission of radioactive radiation
is analyzed, using a NaI(Ti) detector for isotopes of different energy strikes. In
particular, 8x8 and 12x12 SiPM arrays (Silicon PhotoMultiplier) are placed at the
bottom of the cylindrical NaI crystal, on which gamma-ray photons were incident via
ANTS2 software [2]. Knowledge of quantum and classical physics was required for
a deeper understanding of gamma-radiation and matter interaction and radioactive
decays.

Moving forward, by applying the analytical method of weighted vector to silicon
photomultipliers (SiPM), the direction of the source - isotope that radiates isotrop-
ically in space, can be determined easily, reliably and with great accuracy. Based
on the aforementioned analytical method, a python program was implemented to
calculate the maximum source directivity for each isotope, using functions that vi-
sualize the distributions of visible photons in the NaI crystal and provide details on
each step of implementation.

The isotopes studied and finally simulated are: 57Co with characteristic gamma
photopeaks at 122 keV and 136 keV , 137Cs with characteristic gamma photopeak
of 662 keV and 60Co with characteristic gamma photocurves of 1173 keV and 1332
keV . The isotopes were each placed separately, at a distance of one meter (1 m)
from the detector array (the SiP M array and N aI crystal arrays) and thirty-six (36)
measurements were simulated for each isotope at 10 degrees in azimuth, spatially
covering an entire circle.
