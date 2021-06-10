import os
import sys

# Insert relevant path
initPath = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, initPath)

# Import project libraries
from mmodpy.beams.Beams import MAT_HYSTERETIC_BEAM, MAT_ELASTIC
