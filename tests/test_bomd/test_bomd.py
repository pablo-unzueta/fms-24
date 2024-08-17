from pathlib import Path

import torch

from fms24 import BOMD
from fms24.utils import parse_config


def sample_timestep_info(n_atoms: int):
    return {
        "positions": torch.zeros((n_atoms, 3)),
        "momenta": torch.zeros((n_atoms, 3)),
        "forces": torch.zeros((n_atoms, 3)),
        "dt": 20,
    }


class TestBOMD:
    def test_initialization(self, sample_config="config.yaml"):
        current_dir = Path(__file__).parent
        config = parse_config(current_dir / sample_config)
        bomd = BOMD(
            xyz_path=current_dir / config["xyz"],
            dt=config["dt"],
            timesteps=config["timesteps"],
        )
        assert bomd.positions[-1].shape == (2, 3)
        assert bomd.atoms_list == ["H", "H"]
        assert bomd.timesteps == 100

    # def test_update_timestep_info(self, sample_config="config.yaml"):
    #     current_dir = Path(__file__).parent
    #     bomd = BOMD(current_dir / sample_config)

    #     timestep_info = sample_timestep_info(2)
    #     bomd.update_timestep_info(
    #         positions=timestep_info["positions"],
    #         momenta=timestep_info["momenta"],
    #         forces=timestep_info["forces"],
    #         dt=timestep_info["dt"],
    #     )
    #     assert torch.all(bomd.positions[-1] == torch.zeros((2, 3)))
    # assert torch.all(bomd.momenta == torch.zeros((2, 3)))
    # assert torch.all(bomd.forces == torch.zeros((2, 3)))
    # assert bomd.dt == 20
