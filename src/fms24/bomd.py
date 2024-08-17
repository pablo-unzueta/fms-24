from pathlib import Path
from typing import List

import torch
from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io import read, write


class BOMD(torch.nn.Module):
    def __init__(self, xyz_path: Path, dt: float, timesteps: int):
        super().__init__()
        self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        self.xyz_path = xyz_path
        self.dt = dt
        self.timesteps = timesteps
        self.atoms = None
        self.positions: List[torch.Tensor] = []
        self.momenta: List[torch.Tensor] = []
        self.forces: List[torch.Tensor] = []
        self.dts: List[float] = []
        self.atoms_list = None
        self.initialize()

    def initialize(self):
        if not self.xyz_path.exists():
            raise ValueError(f"XYZ file not found: {self.xyz_path}")

        atoms = read(self.xyz_path)
        self.positions.append(torch.tensor(atoms.get_positions(), device=self.device))
        self.momenta.append(torch.zeros_like(self.positions[-1]))
        self.forces.append(torch.zeros_like(self.positions[-1]))
        self.atoms_list = atoms.get_chemical_symbols()
        self.dts.append(self.dt)

    def update_timestep_info(
        self,
        positions: torch.Tensor,
        momenta: torch.Tensor,
        forces: torch.Tensor,
        dt: float,
    ):
        self.positions.append(positions.clone().detach())
        self.momenta.append(momenta.clone().detach())
        self.forces.append(forces.clone().detach())
        self.dts.append(dt)

    def get_trajectory_tensor(self):
        return (
            torch.stack(self.positions),
            torch.stack(self.momenta),
            torch.stack(self.forces),
            torch.tensor(self.dts),
        )

    @property
    def trajectory(self):
        return self.get_trajectory_tensor()

    def lj_potential(
        self,
        positions: torch.Tensor,
        atoms_list: List[str],
        sigma: float = 1.0,
        epsilon: float = 1.0,
    ):
        r = positions[:, None, :] - positions[None, :, :]
        r_mag = torch.norm(r, dim=-1)
        energy = 4 * epsilon * ((sigma / r_mag) ** 12 - (sigma / r_mag) ** 6).sum()
        forces = (
            -24
            * epsilon
            * (
                (2 * (sigma / r_mag) ** 12 - (sigma / r_mag) ** 6)[:, :, None]
                * r
                / r_mag[:, :, None] ** 2
            ).sum(dim=1)
        )
        return energy, forces

    def call_potential(self):
        """
        Call the potential energy calculator

        !!! warning
            Using dummy Lennard-Jones potential for now
        """
        energy, forces = self.lj_potential(self.positions[-1])
        return energy, forces

    def prop_velocity_verlet(self, dt: float):
        current_v = self.velocities[-1] + 0.5 * dt * self.forces[-1]
        new_pos = self.positions[-1] + dt * current_v

        energy, new_forces = self.call_potential()

        new_v = current_v + 0.5 * dt * new_forces

        self.update_timestep_info(new_pos, new_v, new_forces, dt)

    def write_traj_dump(self, output_path: Path):
        # Create an ASE Atoms object directly with PyTorch tensors
        curr_atoms = Atoms(
            positions=self.positions[-1].cpu().numpy(),
            symbols=self.atoms_list,
            pbc=False,
        )

        # Set additional properties
        calculator = SinglePointCalculator(
            curr_atoms,
            energy=self.energies[-1].item(),
            forces=self.forces[-1].cpu().numpy(),
        )
        curr_atoms.calc = calculator

        # Write the extxyz file
        write(
            output_path / "trajectory.extxyz",
            curr_atoms,
            format="extxyz",
        )
