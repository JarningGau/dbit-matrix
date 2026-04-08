"""Shared path existence checks for workflow drivers (make_cmd entrypoints)."""

from __future__ import annotations

from pathlib import Path


def resolve_config_path(s: str) -> Path:
    p = Path(s).expanduser()
    if p.is_absolute():
        return p
    return (Path.cwd() / p).resolve()


def require_file(label: str, path: Path) -> None:
    if not path.is_file():
        raise ValueError(f"missing required file ({label}): {path}")


def require_dir(label: str, path: Path) -> None:
    if not path.is_dir():
        raise ValueError(f"missing required directory ({label}): {path}")


def looks_like_filesystem_path(s: str) -> bool:
    return "/" in s or s.startswith(".") or s.startswith("~")


def require_script_path(label: str, raw: str) -> None:
    require_file(label, resolve_config_path(raw))


def require_optional_executable_path(label: str, raw: str) -> None:
    if looks_like_filesystem_path(raw):
        require_file(label, resolve_config_path(raw))


def bwa_bwt_marker_path(prefix: str) -> Path:
    return resolve_config_path(str(prefix) + ".bwt")


def require_bwa_index_prefix(prefix: str) -> None:
    path = bwa_bwt_marker_path(prefix)
    require_file("bwa index (.bwt marker next to prefix)", path)


def chunk_names(number_of_split_parts: int) -> list[str]:
    if number_of_split_parts <= 0:
        raise ValueError("number_of_split_parts must be > 0")
    width = max(4, len(str(number_of_split_parts)))
    return [f"{index:0{width}d}" for index in range(1, number_of_split_parts + 1)]


def spike_rhs_paths(spike_in_index: list[str]) -> list[tuple[str, Path]]:
    out: list[tuple[str, Path]] = []
    for item in spike_in_index:
        name, sep, rhs = item.partition("=")
        if not sep or not rhs.strip():
            continue
        out.append((name.strip(), resolve_config_path(rhs.strip())))
    return out


def require_spike_index_paths(spike_in_index: list[str]) -> None:
    for name, path in spike_rhs_paths(spike_in_index):
        require_file(f"spike_in_index[{name}]", path)


def validate_pool_align_shards(sample_work: Path, spike_in_index: list[str]) -> None:
    align_dir = sample_work / "align_shards"
    if not align_dir.is_dir():
        raise ValueError(f"missing align_shards directory: {align_dir}")
    host = sorted(align_dir.glob("*.cb.bam"))
    if not host:
        raise ValueError(f"no host shard BAMs (*.cb.bam) under {align_dir}")
    spike_names: list[str] = []
    for item in spike_in_index:
        name, sep, _ = item.partition("=")
        if sep and name.strip() and name.strip() not in spike_names:
            spike_names.append(name.strip())
    for spike_name in spike_names:
        spike_shards = sorted(align_dir.glob(f"*.{spike_name}.bam"))
        if not spike_shards:
            raise ValueError(
                f"no spike shard BAMs (*.{spike_name}.bam) under {align_dir}"
            )


def validate_optional_path_file(label: str, raw: str | None) -> None:
    if raw in (None, ""):
        return
    require_file(label, resolve_config_path(raw))


def validate_optional_path_dir(label: str, raw: str | None) -> None:
    if raw in (None, ""):
        return
    require_dir(label, resolve_config_path(raw))
