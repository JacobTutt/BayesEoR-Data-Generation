"""Generate the EoR, GSM and GLEAM-like+A-team sky components."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping

import numpy as np

from .eor import generate_eor_maps
from .gleam_ateam import generate_gleam_ateam_maps
from .gsm import generate_gsm_maps


def generate_skies(
    config: Mapping[str, Any], overwrite: bool = False
) -> dict[str, list[Path]]:
    """Generate EoR, GSM and GLEAM-like+A-team maps for one observation.

    The function reads the telescope location, time axis, frequency axis and
    component-specific FoVs once. It computes the central pointing and shared
    HEALPix selections, then calls the three independent physical sky
    generators.

    Parameters
    ----------
    config
        Observation configuration dictionary. See the repository-level
        ``example_config.py`` for every required field.
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

    output_root = Path(config["output_directory"]).expanduser().resolve()

    telescope = config["telescope"]
    time_input = config["time"]
    frequency_input = config["frequency"]
    sky_input = config["sky"]
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
            Path(sky_input["gleam_ateam"]["a_team_catalogue"])
            .expanduser()
            .resolve(),
            output_root / "sky_models/gleam_ateam",
            overwrite,
        ),
    }
