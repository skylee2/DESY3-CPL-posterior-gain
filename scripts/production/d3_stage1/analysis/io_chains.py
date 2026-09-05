"""Chain readers. Third-party posterior data are intentionally not bundled."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re

import numpy as np
import pandas as pd
from numpy.typing import NDArray

from .config import ChainConfig
from .weighted_stats import effective_sample_size


@dataclass
class ChainData:
    samples: NDArray
    weights: NDArray
    names: tuple[str, ...]
    log_posterior: NDArray | None
    metadata: dict[str, object]


def _cosmosis_header(path: Path, config: ChainConfig) -> list[str] | None:
    """Return names from a CosmoSIS commented first-line header, if present."""
    if not config.comment:
        return None
    with path.open("r", encoding="utf-8") as stream:
        first = stream.readline().strip()
    if not first.startswith(config.comment):
        return None
    names = first[len(config.comment):].strip().split()
    required = set(config.parameter_map.values())
    if config.weight_column:
        required.add(config.weight_column)
    if config.log_weight_column:
        required.add(config.log_weight_column)
    return names if names and required.issubset(names) else None


def _read_text_metadata(path: Path, comment: str) -> dict[str, object]:
    """Read compact scalar metadata from commented CosmoSIS output lines."""
    metadata: dict[str, object] = {}
    if not comment:
        return metadata
    pattern = re.compile(rf"^{re.escape(comment)}\s*([A-Za-z0-9_]+)\s*=\s*(.*?)\s*$")
    wanted = {"sampler", "nsample", "log_z", "log_z_error", "n_varied"}
    with path.open("r", encoding="utf-8") as stream:
        for line in stream:
            match = pattern.match(line)
            if not match or match.group(1) not in wanted:
                continue
            key, raw = match.groups()
            try:
                value: object = float(raw)
                if key in {"nsample", "n_varied"} and float(value).is_integer():
                    value = int(value)
            except ValueError:
                value = raw
            metadata[key] = value
    return metadata


def _read_frame(config: ChainConfig) -> pd.DataFrame:
    path = Path(config.path)
    if not path.exists():
        raise FileNotFoundError(path)
    suffix = path.suffix.lower()
    if suffix == ".csv":
        return pd.read_csv(path, comment=config.comment, skiprows=config.skip_rows)
    if suffix in {".txt", ".dat", ".tsv"}:
        sep = config.delimiter if config.delimiter is not None else r"\s+"
        names = _cosmosis_header(path, config)
        if names is not None:
            return pd.read_csv(
                path,
                sep=sep,
                comment=config.comment,
                skiprows=config.skip_rows,
                header=None,
                names=names,
            )
        return pd.read_csv(path, sep=sep, comment=config.comment, skiprows=config.skip_rows)
    if suffix == ".npy":
        array = np.load(path, allow_pickle=False)
        if array.dtype.names:
            return pd.DataFrame.from_records(array)
        raise ValueError("plain .npy arrays lack parameter names; use .npz or a labeled text format")
    if suffix == ".npz":
        archive = np.load(path, allow_pickle=False)
        return pd.DataFrame({key: archive[key] for key in archive.files})
    raise ValueError(f"unsupported chain format: {suffix}")


def _weights_from_frame(frame: pd.DataFrame, config: ChainConfig) -> NDArray:
    if config.log_weight_column:
        logw = frame[config.log_weight_column].to_numpy(float)
        # Nautilus and other nested samplers commonly store log posterior weights.
        return np.exp(logw - np.nanmax(logw))
    if config.weight_column:
        return frame[config.weight_column].to_numpy(float)
    return np.ones(len(frame), dtype=float)


def read_chain(config: ChainConfig) -> ChainData:
    frame = _read_frame(config)
    original_size = len(frame)
    start = int(np.floor(config.burn_fraction * original_size))
    frame = frame.iloc[start:: config.thin].copy()
    for canonical, (lower, upper) in config.filters.items():
        source = config.parameter_map.get(canonical, canonical)
        if lower is not None:
            frame = frame.loc[frame[source] >= lower]
        if upper is not None:
            frame = frame.loc[frame[source] <= upper]
    source_names = list(config.parameter_map.values())
    missing = [name for name in source_names if name not in frame]
    if missing:
        raise KeyError(f"missing parameter columns: {missing}; available: {list(frame.columns)}")
    samples = frame[source_names].to_numpy(float)
    weights = _weights_from_frame(frame, config)
    logp = None
    if config.log_posterior_column:
        logp = frame[config.log_posterior_column].to_numpy(float)
    finite = np.all(np.isfinite(samples), axis=1) & np.isfinite(weights) & (weights >= 0)
    if logp is not None:
        finite &= np.isfinite(logp)
    samples, weights = samples[finite], weights[finite]
    logp = None if logp is None else logp[finite]
    if len(samples) < 2 or np.sum(weights) <= 0:
        raise ValueError("chain has insufficient valid weighted samples")
    weights = weights / np.sum(weights)
    metadata: dict[str, object] = {
        "path": str(config.path),
        "original_samples": original_size,
        "retained_samples": len(samples),
        "effective_sample_size": effective_sample_size(weights),
        "burn_fraction": config.burn_fraction,
        "thin": config.thin,
    }
    if Path(config.path).suffix.lower() in {".txt", ".dat", ".tsv"}:
        metadata.update(_read_text_metadata(Path(config.path), config.comment))
    return ChainData(
        samples=samples,
        weights=weights,
        names=tuple(config.parameter_map.keys()),
        log_posterior=logp,
        metadata=metadata,
    )
