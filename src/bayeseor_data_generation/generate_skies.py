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
    component-specific FoVs once. It computes the requested central-time or
    drift-scan-union HEALPix selections, then calls the three independent
    physical sky generators.

    Parameters
    ----------
    config
        Observation configuration dictionary. See the repository-level
        ``example_config.yaml`` for every required field.
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
    import healpy as hp

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

    sky_mask_time_mode = str(
        sky_input.get("sky_mask_time_mode", "center")
    ).lower()
    if sky_mask_time_mode not in {"center", "union"}:
        raise ValueError(
            "sky.sky_mask_time_mode must be 'center' or 'union', got "
            f"{sky_mask_time_mode!r}."
        )

    number_of_times = int(
        round(
            (time_input["end_jd"] - time_input["start_jd"])
            * 86400
            / time_input["integration_time_s"]
        )
    ) + 1
    observation_times_jd = np.linspace(
        time_input["start_jd"], time_input["end_jd"], number_of_times
    )
    mask_times_jd = (
        np.array([central_jd])
        if sky_mask_time_mode == "center"
        else observation_times_jd
    )

    nside = sky_input["nside"]
    healpix_fovs = set(
        sky_input["eor"]["fovs_deg"] + sky_input["gsm"]["fovs_deg"]
    )
    # Use the same ICRS zenith centres and RING-ordered, pixel-centre spherical
    # query discs as the BayesEoR-v2 inference mask. ``mask_times_jd`` contains
    # either the central time alone or the full drift-scan time axis.
    pixel_chunks_by_fov = {fov: [] for fov in healpix_fovs}
    for jd in mask_times_jd:
        zenith = SkyCoord(
            az=0 * u.deg,
            alt=90 * u.deg,
            frame=AltAz(obstime=Time(jd, format="jd"), location=location),
        ).transform_to("icrs")
        zenith_vector = hp.ang2vec(
            zenith.ra.deg, zenith.dec.deg, lonlat=True
        )
        for fov in healpix_fovs:
            pixel_chunks_by_fov[fov].append(
                hp.query_disc(
                    nside,
                    zenith_vector,
                    np.deg2rad(fov / 2),
                    inclusive=False,
                    nest=False,
                )
            )
    pixels_by_fov = {
        fov: np.unique(np.concatenate(pixel_chunks))
        for fov, pixel_chunks in pixel_chunks_by_fov.items()
    }

    print(
        f"Sky mask: {sky_mask_time_mode} mode using {mask_times_jd.size} "
        f"time(s); central JD {central_jd:.12f}"
    )
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
            mask_times_jd,
            location,
            Path(sky_input["gleam_ateam"]["a_team_catalogue"])
            .expanduser()
            .resolve(),
            output_root / "sky_models/gleam_ateam",
            overwrite,
        ),
    }
