# Radioactive source emission direction estimation using a cylindrical NaI detector and MVA techniques

## Abstract
In this thesis, the study of spatial directional emission of radioactive radiation
is analyzed, using a NaI(Tl) detector for isotopes of different energy strikes. In
particular, 8x8 and 12x12 SiPM arrays (Silicon PhotoMultiplier) are placed at the
bottom of the cylindrical NaI crystal, on which gamma-ray photons were incident via
ANTS2 software.

Moving forward, by applying the analytical method of weighted vector to silicon
photomultipliers (SiPM), the direction of the source - isotope that radiates isotropically
in space, can be determined easily, reliably and with great accuracy. Based
on the aforementioned analytical method, a python program was implemented to
calculate the maximum source directivity for each isotope, using functions that visualize
the distributions of visible photons in the NaI crystal and provide details on
each step of implementation.

The isotopes studied and finally simulated are: <sup>57</sup>Co
 with characteristic gamma photopeaks at 122 keV and 136 keV , <sup>137</sup>Cs with characteristic gamma photopeak
of 662 keV and <sup>60</sup>Co with characteristic gamma photocurves of 1173 keV and 1332
keV. The isotopes were each placed separately, at a distance of one meter (1 m)
from the detector array (the SiPM array and NaI crystal arrays) and thirty-six (36)
measurements were simulated for each isotope at 10 degrees in azimuth, spatially
covering an entire circle.

Radioactive isotopes with characteristic low and high photopeaks energy were
deliberately used to observe it:

1. how the accuracy of determining the azimuthal emission angle is affected as
the energy of the emitted g-photons increases,
2. the optimal accuracy of determining the azimuthal emission angle that can be
achieved with these geometric configurations for all isotopes, especially those
with the highest energy such as 60Co,
3. which layout of the two (8x8 SiPM layout or 12x12 SiPM layout) achieves the
best accuracy of determining the azimuthal emission angle for a given number
of events, equivalent to a given exposure time

## ANTS2 Simulation


ANTS2 is a software for simulation and data processing for Anger Camera-type detection sensors. 
It is based on the ROOT data processing software, which was developed at the European Organization for Nuclear Research (CERN). 
The software can simulate radioactive sources and their interactions with scintillation detectors. Additionally, it has the capability 
to simulate signals from photosensitive sensors, such as photomultiplier tubes (PMTs) and silicon photomultipliers (SiPMs).

### Data input into ANTS2

The configuration of the detector's characteristics can be saved, loaded, or transmitted to the primary 
structure using JSON (JavaScript Object Notation) files.The simulation module performs Monte Carlo simulations
for each gamma photon in the detector, based on the provided settings. By selecting the advanced options, we 
choose the type of radioactive isotope we want to study, assuming it emits radiation isotropically.

### Characteristics of the Detector

Initially, we define the scintillation detector, its geometry, and its type. We have selected a NaI detector with the following geometric characteristics:

- Outer diameter: 51.5 mm (2”)
- Cylinder height: 76.5 mm (3”)

This specific type and geometry is a very common type of scintillation detector that can be found easily and affordably in the market. 
Next, we place the SiPM detector at the bottom of the cylinder, forming a 12x12 two-dimensional array with the following geometric characteristics:

- Square dimensions: 3 mm x 3 mm
- Thickness: 0.1 mm
- Center-to-center distance: 4.2 mm
- PDE (Photon Detection Efficiency): 0.41

Another type of SiPM is the 8x8 two-dimensional array with the following geometric characteristics:

- Square dimensions: 6.07 mm x 6.07 mm
- Thickness: 0.1 mm
- Center-to-center distance: 6.13 mm
- PDE: 0.5

<p align="center">
  <img src="images/12x12NaI.png" alt="Diagram of the project" width="500" />
  <br>
  <i>Visualization of a cylindrical NaI detector in the ANTS2 program. The silicon photomultipliers are positioned at the bottom of the cylinder in a 12x12 array, covering the diameter of the detector</i>
</p>

### Simulation Results of the Isotope <sup>57</sup>Co

The distribution of scintillation photons resulting from the interaction of a single gamma photon of the 57Co isotope with the crystal, with phi = 0°.

The following observations are made:

- Most photons are detected on the right side of the detector. However, a certain number of scintillation photons are distributed across all SiPMs in the detector,
  making it challenging to precisely determine the direction of the point source in these configurations. This highlights the difficulties and the algorithms that need to be developed to accurately determine the direction of the point source in these arrangements.
- The SiPM photomultipliers located at the corners do not detect photons due to the cylindrical geometry of the scintillation detector and the square geometry of the SiPM array
  
<p align="center">
  <img src="images/12x12NaI.png" alt="Diagram of the project" width="500" />
  <br>
  <i>Visualization of a cylindrical NaI detector in the ANTS2 program. The silicon photomultipliers are positioned at the bottom of the cylinder in a 12x12 array, covering the diameter of the detector</i>
</p>







