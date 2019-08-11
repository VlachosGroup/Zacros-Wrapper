import matplotlib as mat
import platform
if platform.system() == 'Linux':
	mat.use('Agg')

import matplotlib.pyplot as plt
from zacros_wrapper.utils import *
from zacros_wrapper.IO_data import *
from zacros_wrapper.Lattice import Lattice
from zacros_wrapper.KMC_Run import *
from zacros_wrapper.Replicates import Replicates
from zacros_wrapper.RateRescaling import *
