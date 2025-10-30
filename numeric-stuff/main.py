from systems.two_dimensional_system import TwoDimensionalSystem
from systems.three_dimensional_system import ThreeDimensionalSystem
from systems.three_dim_non_dimensionalised import NonDimensionalThreeDim
from systems.two_dim_non_dimensionalised import NonDimensionalTwoDim
import sys
from systems.system import System
from typing import Dict

systems: Dict[str, System] = {
    "2d": TwoDimensionalSystem,
    "3d": ThreeDimensionalSystem,
    "3dnon": NonDimensionalThreeDim,
    "2dnon": NonDimensionalTwoDim
}

if __name__ == "__main__":
    system_name = sys.argv[1]
    choice = sys.argv[2]
    if choice == "normal":
        systems[system_name].run()
    if choice == "move":
        systems[system_name].move()
