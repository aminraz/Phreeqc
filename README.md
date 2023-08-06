
The Python scripts presented in this repository can simulate Kinetics
of microbially induced calcium carbonate precipitation (MICP) in a
batch system. This script was used in my article [1]. A geochemical
solver called [Phreeqc](https://www.usgs.gov/software/phreeqc-version-3) is used for geochemical calculations. [Phreeqpy](https://www.phreeqpy.com/)
is a python interface for Phreeqc that is used here. Phreeqc handles
various geochemical calculations using different keyword blocks. A
keyword block is a code segment in the input file processed by
Phreeqc. The full description of keyword blocks is provided in the
[Phreeqc manual](https://pubs.usgs.gov/tm/06/a43/).

(1) Razbani, M. A.; Jettestuen, E.; Røyne, A. Direct Pore-Scale
Numerical Simulation of Microbially Induced Calcium Carbonate
Precipitation. Water Resources Research 2023, 59
(1). <https://doi.org/10.1029/2022WR032988>.
