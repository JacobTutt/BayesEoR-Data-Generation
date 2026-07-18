"""Generate EoR-like HEALPix sky maps for requested fields of view."""

from pathlib import Path

import numpy as np


def generate_eor_maps(
    sky_inputs: dict,
    frequencies_hz: np.ndarray,
    pixels_by_fov: dict[float, np.ndarray],
    central_jd: float,
    output_directory: Path,
    overwrite: bool = False,
) -> list[Path]:
    """Generate and save the configured EoR-like sky maps.

    Parameters
    ----------
    sky_inputs
        The ``sky`` section of the observation configuration. The EoR subsection
        supplies NSIDE, RMS temperature, random seed and FoV diameters.
    frequencies_hz
        Complete observing frequency axis in Hz.
    pixels_by_fov
        HEALPix indices selected using the configured time-mask mode for every
        requested EoR FoV.
    central_jd
        Midpoint of the observation in Julian days; recorded in file history.
    output_directory
        Directory that will contain ``fov-X.skyh5`` files.
    overwrite
        Replace existing maps when true.

    Returns
    -------
    list[pathlib.Path]
        One skyh5 path for every configured EoR FoV.

    Notes
    -----
    The full-sky random cube is generated before selecting pixels. This
    preserves the historical random-number sequence. The cube is not saved.
    """
    from astropy import units as u
    from pyradiosky import SkyModel

    nside = sky_inputs["nside"]
    eor_inputs = sky_inputs["eor"]
    number_of_pixels = 12 * nside**2
    output_directory.mkdir(parents=True, exist_ok=True)
    outputs = [
        output_directory
        / (f"fov-{fov:.1f}".rstrip("0").rstrip(".") + ".skyh5")
        for fov in eor_inputs["fovs_deg"]
    ]
    if not overwrite and all(path.exists() for path in outputs):
        return outputs

    np.random.seed(eor_inputs["random_seed"])
    intensity = np.random.normal(
        0.0,
        eor_inputs["rms_k"],
        (frequencies_hz.size, number_of_pixels),
    )
    intensity -= intensity.mean(axis=1)[:, None]
    stokes = np.zeros((4, frequencies_hz.size, number_of_pixels)) * u.K
    stokes[0] = intensity * u.K
    full_sky = SkyModel(
        nside=nside,
        hpx_inds=np.arange(number_of_pixels),
        stokes=stokes,
        spectral_type="full",
        freq_array=frequencies_hz * u.Hz,
        frame="icrs",
        hpx_order="ring",
        history=(
            "Gaussian EoR-like sky generated from scratch; "
            f"RMS={eor_inputs['rms_k']} K; seed={eor_inputs['random_seed']}; "
            "each frequency channel shifted to zero spatial mean."
        ),
    )

    mask_mode = str(sky_inputs.get("sky_mask_time_mode", "center")).lower()
    for fov, output in zip(eor_inputs["fovs_deg"], outputs):
        if output.exists() and not overwrite:
            continue
        selected = full_sky.select(
            component_inds=pixels_by_fov[fov], inplace=False
        )
        if mask_mode == "union":
            selected.history += (
                f"\nFoV diameter {fov} deg selected using the drift-scan "
                f"union mask around central JD {central_jd}."
            )
        else:
            selected.history += (
                f"\nFoV diameter {fov} deg selected at central JD {central_jd}."
            )
        print(f"Writing EoR map: {output}")
        selected.write_skyh5(output, clobber=overwrite)
    return outputs
