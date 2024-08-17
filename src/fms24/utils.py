from pathlib import Path

import yaml


def parse_config(config_path: Path) -> dict:
    with open(config_path, "r") as f:
        config = yaml.safe_load(f)

    # check if config has xyz key
    if "xyz" not in config:
        raise ValueError("No xyz file provided")
    if "dt" not in config:
        raise ValueError("No dt provided")
    if "timesteps" not in config:
        raise ValueError("No timesteps provided")

    return config
