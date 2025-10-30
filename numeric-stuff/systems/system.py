from abc import ABC, abstractmethod

class System(ABC):
    @abstractmethod
    def run() -> None:
        pass