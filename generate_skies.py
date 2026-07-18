#!/usr/bin/env python3
"""Read one observation YAML and generate its three sky components."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import yaml

from eor import generate_eor_maps
from gleam_ateam import generate_gleam_ateam_maps
from gsm import generate_gsm_maps


REPOSITORY = Path(__file__).resolve().parent


def generate_skies(input_file: Path, overwrite: bool = False) -> dict[str, list[Path]]:
    """Generate EoR, GSM and GLEAM-like+A-team maps for one observation.

    The function reads the telescope location, time axis, frequency axis and
    component-specific FoVs once. It computes the central pointing and shared
    HEALPix selections, then calls the three independent physical sky
    generators.

    Parameters
    ----------
    input_file
        Observation YAML; see ``inputs/h1c_band2.yaml``.
    overwrite
        Replace existing skyh5 products when true.

    Returns
    -------
    dict
        Map from component name to its generated skyh5 paths.
    """
    from astropy import units as u
    from astropy.coordinates import AltAz, EarthLocation, SkyCoord
    from astropy.time import Time
    from astropy_healpix import healpy as hp

    with input_file.open() as stream:
        inputs = yaml.safe_load(stream)
    output_root = Path(inputs["output_directory"])
    if not output_root.is_absolute():
        output_root = REPOSITORY / output_root

    telescope = inputs["telescope"]
    time_input = inputs["time"]
    frequency_input = inputs["frequency"]
    sky_input = inputs["sky"]
    central_jd = (time_input["start_jd"] + time_input["end_jd"]) / 2
    frequencies_hz = frequency_input["start_hz"] + np.arange(
        frequency_input["n_channels"]
    ) * frequency_input["channel_width_hz"]
    location = EarthLocation.from_geodetic(
        telescope["longitude_deg"] * u.deg,
        telescope["latitude_deg"] * u.deg,
        telescope["altitude_m"] * u.m,
    )

    nside = sky_input["nside"]
    pixels = np.arange(12 * nside**2)
    ra, dec = hp.pix2ang(nside, pixels, lonlat=True)
    horizontal = SkyCoord(ra * u.deg, dec * u.deg, frame="icrs").transform_to(
        AltAz(obstime=Time(central_jd, format="jd"), location=location)
    )
    zenith_angle = np.pi / 2 - horizontal.alt.rad
    healpix_fovs = set(sky_input["eor"]["fovs_deg"] + sky_input["gsm"]["fovs_deg"])
    pixels_by_fov = {
        fov: pixels[zenith_angle <= np.deg2rad(fov / 2)]
        for fov in healpix_fovs
    }

    print(f"Central sky time: JD {central_jd:.12f}")
    return {
        "eor": generate_eor_maps(
            sky_input,
            frequencies_hz,
            pixels_by_fov,
            central_jd,
            output_root / "sky_models/eor",
            overwrite,
        ),
        "gsm": generate_gsm_maps(
            sky_input,
            frequencies_hz,
            pixels_by_fov,
            central_jd,
            output_root / "sky_models/gsm",
            overwrite,
        ),
        "gleam_ateam": generate_gleam_ateam_maps(
            sky_input,
            central_jd,
            location,
            output_root / "sky_models/gleam_ateam",
            overwrite,
        ),
    }


def main() -> None:
    """Run the three-component sky-generation workflow."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "input_file",
        type=Path,
        nargs="?",
        default=REPOSITORY / "inputs/h1c_band2.yaml",
    )
    parser.add_argument("--overwrite", action="store_true")
    arguments = parser.parse_args()
    generate_skies(arguments.input_file.resolve(), arguments.overwrite)


if __name__ == "__main__":
    main()
