from systems.three_dim_non_dimensionalised import NonDimensionalThreeDim
from systems.two_dim_non_dimensionalised import NonDimensionalTwoDim
import sys
from systems.system import System
from typing import Dict

systems: Dict[str, System] = {
    "2d": NonDimensionalTwoDim,
    "3d": NonDimensionalThreeDim,
}

if __name__ == "__main__":
    system_name = sys.argv[1]
    if len(sys.argv) >  2 and (sys.argv[2] == 'move'):
        systems[system_name].move()
    else:
        systems[system_name].run()
