import torch
from dataclasses import dataclass, field
from typing import List, Dict
import yaml


@dataclass
class TBF:
    """A Trajectory Basis Function (TBF)"""
    positions: torch.Tensor
    momenta: torch.Tensor
    forces: torch.Tensor
    phase: complex = 1.0 + 0j
    amplitude: complex = 1.0 + 0j

    def __post_init__(self):
        """Place Frozen Gaussian Wavepacket in the center of positions and momenta in phase space"""
        self.positions_com = torch.zeros_like(self.positions)
        self.momenta_com = torch.zeros_like(self.momenta)


@dataclass
class Bundle:
    config: Dict
    tbfs: List[TBF] = field(default_factory=list)
    device: torch.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    def __post_init__(self):
        self.parse_config()
        self.initialize_bundle()

    def parse_config(self):
        with open(self.config, "r") as f:
            self.config = yaml.safe_load(f)

    def initialize_bundle(self):
        initial_positions = torch.tensor(
            self.config["initial_positions"], device=self.device
        )
        initial_momenta = torch.tensor(
            self.config["initial_momenta"], device=self.device
        )
        initial_forces = torch.zeros_like(
            initial_positions
        )  # Initialize forces to zero

        initial_tbf = TBF(initial_positions, initial_momenta, initial_forces)
        self.tbfs.append(initial_tbf)

    def add_tbf(
        self,
        positions: torch.Tensor,
        momenta: torch.Tensor,
        forces: torch.Tensor,
        amplitude: complex = 1.0 + 0j,
    ):
        new_tbf = TBF(positions, momenta, forces, amplitude)
        self.tbfs.append(new_tbf)

    def remove_tbf(self, index: int):
        if 0 <= index < len(self.tbfs):
            del self.tbfs[index]
        else:
            raise IndexError("TBF index out of range")

    def get_num_tbfs(self):
        return len(self.tbfs)

    def update_tbf(
        self,
        index: int,
        positions: torch.Tensor,
        momenta: torch.Tensor,
        forces: torch.Tensor,
    ):
        if 0 <= index < len(self.tbfs):
            self.tbfs[index].positions = positions
            self.tbfs[index].momenta = momenta
            self.tbfs[index].forces = forces
        else:
            raise IndexError("TBF index out of range")

    def get_tbf(self, index: int) -> TBF:
        if 0 <= index < len(self.tbfs):
            return self.tbfs[index]
        else:
            raise IndexError("TBF index out of range")

    def propagate_tbfs(self, dt: float):
        for tbf in self.tbfs:
            # Implement your propagation method here
            # This is a placeholder for velocity Verlet or any other propagation method
            tbf.positions += tbf.momenta * dt + 0.5 * tbf.forces * dt**2
            tbf.momenta += 0.5 * tbf.forces * dt
            # You would typically update forces here based on the new positions
            # tbf.forces = self.calculate_forces(tbf.positions)
            tbf.momenta += 0.5 * tbf.forces * dt

    def calculate_forces(self, positions: torch.Tensor) -> torch.Tensor:
        # Implement your force calculation method here
        # This is a placeholder
        return torch.zeros_like(positions)
