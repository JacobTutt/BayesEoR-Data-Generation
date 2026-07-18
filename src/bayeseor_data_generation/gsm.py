"""Generate pygdsm diffuse-foreground maps for requested fields of view."""

from pathlib import Path

import numpy as np


def generate_gsm_maps(
    sky_inputs: dict,
    frequencies_hz: np.ndarray,
    pixels_by_fov: dict[float, np.ndarray],
    central_jd: float,
    output_directory: Path,
    overwrite: bool = False,
) -> list[Path]:
    """Generate and save GSM2008 maps at every observing frequency.

    Parameters
    ----------
    sky_inputs
        The ``sky`` section of the observation configuration. The GSM subsection
        supplies FoVs, CMB inclusion and the memory-saving chunk size.
    frequencies_hz
        Complete observing frequency axis in Hz.
    pixels_by_fov
        HEALPix indices selected at the central time for every requested FoV.
    central_jd
        Midpoint of the observation in Julian days; recorded in file history.
    output_directory
        Directory that will contain ``fov-X.skyh5`` files.
    overwrite
        Replace existing maps when true.

    Returns
    -------
    list[pathlib.Path]
        One skyh5 path for every configured GSM FoV.

    Notes
    -----
    pygdsm supplies the spatially varying spectral structure. Evaluating the
    frequencies in chunks changes memory use but not the physical model.
    """
    import healpy
    from astropy import units as u
    from pygdsm import GlobalSkyModel
    from pyradiosky import SkyModel

    nside = sky_inputs["nside"]
    gsm_inputs = sky_inputs["gsm"]
    number_of_pixels = 12 * nside**2
    output_directory.mkdir(parents=True, exist_ok=True)
    outputs = [
        output_directory
        / (f"fov-{fov:.1f}".rstrip("0").rstrip(".") + ".skyh5")
        for fov in gsm_inputs["fovs_deg"]
    ]
    if not overwrite and all(path.exists() for path in outputs):
        return outputs

    intensity = np.empty((frequencies_hz.size, number_of_pixels))
    model = GlobalSkyModel(include_cmb=gsm_inputs["include_cmb"])
    chunk_size = gsm_inputs["frequency_chunk_size"]
    for start in range(0, frequencies_hz.size, chunk_size):
        stop = min(start + chunk_size, frequencies_hz.size)
        native_maps = np.asarray(model.generate(frequencies_hz[start:stop] / 1e6))
        if native_maps.ndim == 1:
            native_maps = native_maps[None, :]
        for local_index, native_map in enumerate(native_maps):
            intensity[start + local_index] = healpy.ud_grade(native_map, nside)

    stokes = np.zeros((4, frequencies_hz.size, number_of_pixels)) * u.K
    stokes[0] = intensity * u.K
    full_sky = SkyModel(
        nside=nside,
        hpx_inds=np.arange(number_of_pixels),
        stokes=stokes,
        spectral_type="full",
        freq_array=frequencies_hz * u.Hz,
        frame="galactic",
        hpx_order="ring",
        history=(
            "GSM2008 generated from scratch with pygdsm, downgraded to "
            f"NSIDE={nside}, and transformed from Galactic to ICRS."
        ),
    )
    full_sky.healpix_interp_transform(frame="icrs", full_sky=True)

    for fov, output in zip(gsm_inputs["fovs_deg"], outputs):
        if output.exists() and not overwrite:
            continue
        selected = full_sky.select(
            component_inds=pixels_by_fov[fov], inplace=False
        )
        selected.history += (
            f"\nFoV diameter {fov} deg selected at central JD {central_jd}."
        )
        print(f"Writing GSM map: {output}")
        selected.write_skyh5(output, clobber=overwrite)
    return outputs
