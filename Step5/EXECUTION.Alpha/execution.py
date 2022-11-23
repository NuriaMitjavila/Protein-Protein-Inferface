import os
import sys
from pathlib import Path

home = str(Path.home())
myPath = os.path.join(home, "Downloads", "EXECUTION.Alpha", "Scripts")
sys.path.insert(1, myPath)

import structure_checking
import interaction_energy
import ala_scanning
