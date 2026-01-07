"""
Configuration loader utility.

Provides unified config loading with local override support.
"""

from pathlib import Path
from typing import Any

import yaml


def deep_merge(base: dict, override: dict) -> dict:
    """Recursively merge override dict into base dict."""
    for key, value in override.items():
        if key in base and isinstance(base[key], dict) and isinstance(value, dict):
            deep_merge(base[key], value)
        else:
            base[key] = value
    return base


def load_config(config_path: Path | str = "config.yaml") -> dict[str, Any]:
    """
    Load configuration with local overrides.

    Args:
        config_path: Path to main config file (default: config.yaml)

    Returns:
        Merged configuration dictionary

    The function looks for a config.local.yaml in the same directory
    and merges it over the base config if present.
    """
    config_path = Path(config_path)

    with open(config_path) as f:
        config = yaml.safe_load(f)

    # Merge local overrides if exists
    local_path = config_path.parent / "config.local.yaml"
    if local_path.exists():
        with open(local_path) as f:
            local_config = yaml.safe_load(f) or {}
        deep_merge(config, local_config)

    return config


# Make functions available for direct import from file path
__all__ = ["load_config", "deep_merge"]
