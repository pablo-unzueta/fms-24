from pathlib import Path

import pytest
import yaml

from fms24.utils import parse_config


@pytest.fixture
def temp_config_file(tmp_path):
    def _create_config(content):
        config_path = tmp_path / "test_config.yaml"
        with open(config_path, "w") as f:
            yaml.dump(content, f)
        return config_path

    return _create_config


def test_parse_config_success(temp_config_file):
    config_content = {"xyz": "path/to/xyz/file.xyz", "dt": 0.1, "timesteps": 100}
    config_path = temp_config_file(config_content)

    result = parse_config(config_path)

    assert result == config_content
    assert isinstance(result, dict)
    assert "xyz" in result
    assert "dt" in result
    assert "timesteps" in result


def test_parse_config_missing_xyz(temp_config_file):
    config_content = {"dt": 0.1, "timesteps": 100}
    config_path = temp_config_file(config_content)

    with pytest.raises(ValueError, match="No xyz file provided"):
        parse_config(config_path)


def test_parse_config_missing_dt(temp_config_file):
    config_content = {"xyz": "path/to/xyz/file.xyz", "timesteps": 100}
    config_path = temp_config_file(config_content)

    with pytest.raises(ValueError, match="No dt provided"):
        parse_config(config_path)


def test_parse_config_missing_timesteps(temp_config_file):
    config_content = {"xyz": "path/to/xyz/file.xyz", "dt": 0.1}
    config_path = temp_config_file(config_content)

    with pytest.raises(ValueError, match="No timesteps provided"):
        parse_config(config_path)


def test_parse_config_file_not_found():
    non_existent_path = Path("non_existent_config.yaml")

    with pytest.raises(FileNotFoundError):
        parse_config(non_existent_path)
