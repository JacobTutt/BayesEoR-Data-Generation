#!/usr/bin/env python3
"""Turn the generated component sky maps into pyuvsim visibilities."""

from __future__ import annotations

import argparse
from pathlib import Path

import yaml


REPOSITORY = Path(__file__).resolve().parent


def simulate_visibilities(
    input_file: Path,
    overwrite: bool = False,
) -> list[Path]:
    """Simulate one UVH5 file for every generated component/FoV sky map.

    Parameters
    ----------
    input_file
        The same observation YAML passed to ``generate_skies.py``. Telescope,
        time and frequency values are written directly into pyuvsim inputs, so
        the visibility calculation cannot silently use a different observation
        from the sky generation.
    overwrite
        Replace existing observation YAML and UVH5 files when true.

    Returns
    -------
    list[pathlib.Path]
        Expected UVH5 output paths, one per component and FoV.

    Raises
    ------
    FileNotFoundError
        If an antenna layout or required skyh5 file is missing.

    Notes
    -----
    Under MPI, rank zero writes the shared telescope and observation files.
    After one barrier, all ranks enter pyuvsim together.
    """
    from mpi4py import MPI
    from pyuvsim.uvsim import run_uvsim

    communicator = MPI.COMM_WORLD
    rank = communicator.Get_rank()
    with input_file.open() as stream:
        inputs = yaml.safe_load(stream)

    output_root = Path(inputs["output_directory"])
    if not output_root.is_absolute():
        output_root = REPOSITORY / output_root
    telescope = inputs["telescope"]
    time_input = inputs["time"]
    frequency_input = inputs["frequency"]
    sky_input = inputs["sky"]

    antenna_layout = Path(telescope["antenna_layout"])
    if not antenna_layout.is_absolute():
        antenna_layout = REPOSITORY / antenna_layout
    if not antenna_layout.exists():
        raise FileNotFoundError(f"Missing antenna layout: {antenna_layout}")

    observation_root = output_root / "observation_files"
    visibility_root = output_root / "visibilities"
    telescope_file = observation_root / "telescope.yaml"
    if rank == 0:
        observation_root.mkdir(parents=True, exist_ok=True)
        visibility_root.mkdir(parents=True, exist_ok=True)
        if overwrite or not telescope_file.exists():
            telescope_parameters = {
                "beam_paths": {
                    0: {
                        "type": telescope["beam"]["type"],
                        "diameter": telescope["beam"]["diameter_m"],
                    }
                },
                "telescope_location": [
                    telescope["latitude_deg"],
                    telescope["longitude_deg"],
                    telescope["altitude_m"],
                ],
                "telescope_name": telescope["name"],
            }
            with telescope_file.open("w") as stream:
                yaml.safe_dump(telescope_parameters, stream, sort_keys=False)

    products: list[tuple[Path, Path, Path]] = []
    for component in ("eor", "gsm", "gleam_ateam"):
        for fov in sky_input[component]["fovs_deg"]:
            label = f"{fov:.1f}".rstrip("0").rstrip(".")
            sky_file = output_root / "sky_models" / component / f"fov-{label}.skyh5"
            observation_file = observation_root / component / f"fov-{label}.yaml"
            visibility_file = visibility_root / component / f"fov-{label}.uvh5"
            if not sky_file.exists():
                raise FileNotFoundError(
                    f"Missing {sky_file}; run generate_skies.py first."
                )
            products.append((observation_file, visibility_file, sky_file))
            if rank != 0:
                continue
            observation_file.parent.mkdir(parents=True, exist_ok=True)
            visibility_file.parent.mkdir(parents=True, exist_ok=True)
            if observation_file.exists() and not overwrite:
                continue
            parameters = {
                "filing": {
                    "outdir": str(visibility_file.parent.resolve()),
                    "outfile_name": visibility_file.stem,
                    "output_format": "uvh5",
                    "clobber": bool(overwrite),
                },
                "freq": {
                    "Nfreqs": frequency_input["n_channels"],
                    "channel_width": frequency_input["channel_width_hz"],
                    "start_freq": frequency_input["start_hz"],
                },
                "sources": {"catalog": str(sky_file.resolve())},
                "telescope": {
                    "array_layout": str(antenna_layout.resolve()),
                    "telescope_config_name": str(telescope_file.resolve()),
                },
                "time": {
                    "integration_time": time_input["integration_time_s"],
                    "start_time": time_input["start_jd"],
                    "end_time": time_input["end_jd"],
                },
                "select": {
                    "redundant_threshold": telescope["redundant_threshold_m"]
                },
            }
            print(f"Writing pyuvsim input: {observation_file}")
            with observation_file.open("w") as stream:
                yaml.safe_dump(parameters, stream, sort_keys=False)

    communicator.Barrier()
    outputs = []
    for observation_file, visibility_file, _ in products:
        outputs.append(visibility_file)
        if visibility_file.exists() and not overwrite:
            if rank == 0:
                print(f"Exists, keeping: {visibility_file}")
            continue
        if rank == 0:
            print(f"Running pyuvsim: {observation_file}", flush=True)
        run_uvsim(str(observation_file))
    return outputs


def main() -> None:
    """Run visibility simulation for all configured component sky files."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "input_file",
        type=Path,
        nargs="?",
        default=REPOSITORY / "inputs/h1c_band2.yaml",
    )
    parser.add_argument("--overwrite", action="store_true")
    arguments = parser.parse_args()
    simulate_visibilities(arguments.input_file.resolve(), arguments.overwrite)


if __name__ == "__main__":
    main()
