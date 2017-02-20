from __future__ import print_function
import fatslimlib
from fatslimlib.core_base import Trajectory
from fatslimlib.core_datareading import GroReaderTopol, GroReaderTopolMDAnalysis

def load_topology_fatslim():
    print("Loading topology using old FATSLiM routines")

print("Disclaimer: The goal of this script is to compared the performances of the FATSLiM I/O with the MDAnalysis-based I/O")
