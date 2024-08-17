from pathlib import Path

from fms24 import BOMD


class TestBOMD:
    def test_initialization(self, sample_config="config.yaml"):
        current_dir = Path(__file__).parent
        bomd = BOMD(current_dir / sample_config)
        assert bomd.positions.shape == (2, 3)
        assert bomd.atoms_list == ["H", "H"]
        assert bomd.timesteps == 100
