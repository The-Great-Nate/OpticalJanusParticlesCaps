import sys
import os
import numpy as np
import time
import h5py
import pyoptics.forces as pyf
import matplotlib.pyplot as plt
from pyoptics import dipoles
from pyoptics import readyaml
from pyoptics import hydro
from pyoptics import constants as ct
from scipy.spatial.transform import Rotation
import numpy as np


filestem = sys.argv[1]
filename_vtf = filestem+".vtf"
filename_xl = filestem+".xlsx"
filename_hdf = filestem+".h5"
filename_yaml = filestem+".yml"

options = readyaml.Options(filestem)

display = options.display

beam_collection = options.beam_collection

fig,ax = display.plot_intensity(beam_collection)

fig.tight_layout()
plt.show()