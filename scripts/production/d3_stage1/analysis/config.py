"""Configuration objects for DES Y6 chain ingestion."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path


@dataclass(frozen=True)
class ChainConfig:
    path: Path
    parameter_map: dict[str, str] = field(default_factory=lambda: {"w0": "w0", "wa": "wa"})
    weight_column: str | None = "weight"
    log_weight_column: str | None = None
    log_posterior_column: str | None = None
    delimiter: str | None = None
    comment: str = "#"
    skip_rows: int = 0
    burn_fraction: float = 0.0
    thin: int = 1
    filters: dict[str, tuple[float | None, float | None]] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not 0 <= self.burn_fraction < 1:
            raise ValueError("burn_fraction must be in [0, 1)")
        if self.thin < 1:
            raise ValueError("thin must be >= 1")
        if self.weight_column and self.log_weight_column:
            raise ValueError("specify either weight_column or log_weight_column, not both")


DATASET_ROLES = {
    "des_sn_bao": "Geometry-dominated DES baseline",
    "des_3x2pt": "Growth-sensitive mixed geometry+growth LSS likelihood",
    "des_sn_bao_3x2pt": "Primary DES geometry-plus-LSS comparison",
}

PUBLISHED_FOM = {
    "DES SN+DES BAO": 11.0,
    "DES SN+DES BAO+DES Y6 3x2pt": 48.0,
    "DES SN+DESI BAO": 61.0,
    "DES SN+DESI BAO+DES Y6 3x2pt": 110.0,
    "DES SN+DESI BAO+CMB": 202.0,
    "DES SN+DESI BAO+CMB+DES Y6 3x2pt": 222.0,
}
