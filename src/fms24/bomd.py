import torch
from ase.io import write
from ase import Atoms
from ase.calculators.calculator import SinglePointCalculator
import yaml


class BOMD:
    def __init__(self, config="config.yaml"):
        self.config = config
        self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        self.atoms = None
        self.positions = []
        self.velocities = []
        self.forces = []
        self.timesteps = []
        self.atoms_list = None

    def initialize(self):
        with open(self.config, "r") as f:
            config_data = yaml.safe_load(f)

        atoms = Atoms(config_data["xyz"])
        self.positions = torch.tensor(atoms.get_positions(), device=self.device)
        self.velocities = torch.tensor(atoms.get_velocities(), device=self.device)
        self.atoms_list = atoms.get_chemical_symbols()
        self.dt = config_data["dt"]
        self.timesteps = [0.0]

    def update_timestep_info(
        self,
        positions: torch.Tensor,
        velocities: torch.Tensor,
        forces: torch.Tensor,
        dt: float,
    ):
        self.positions.append(positions)
        self.velocities.append(velocities)
        self.forces.append(forces)
        self.timesteps.append(self.timesteps[-1] + dt if self.timesteps else dt)

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

    def write_traj_dump(self):
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
            self.config["output"]["path"] / f"trajectory.extxyz",
            curr_atoms,
            format="extxyz",
        )
