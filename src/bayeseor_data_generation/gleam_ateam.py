"""Generate the statistical GLEAM-like catalogue plus ten A-team sources."""

import csv
from functools import cached_property
from pathlib import Path

import numpy as np


def _franzen_differential_source_count(flux_jy: float | np.ndarray):
    """Evaluate the Franzen et al. (2019) differential source-count model."""
    coefficients = [3.52, 0.307, -0.388, -0.0404, 0.0351, 0.00600]
    log_flux = np.log10(flux_jy)
    euclidean_normalised = 10 ** sum(
        coefficient * log_flux**power
        for power, coefficient in enumerate(coefficients)
    )
    return flux_jy**-2.5 * euclidean_normalised


class FranzenSourceCounts:
    """Numerically sample the historical GLEAM-like source distribution.

    The expensive cumulative distribution is constructed once and reused after
    imposing the bright-source and confusion limits.
    """

    def __init__(
        self,
        minimum_flux_jy: float,
        maximum_flux_jy: float,
        base_fluxes: np.ndarray | None = None,
        base_cdf: np.ndarray | None = None,
    ) -> None:
        """Construct the cumulative source-count distribution."""
        from scipy.integrate import quad

        self.minimum_flux_jy = minimum_flux_jy
        self.maximum_flux_jy = maximum_flux_jy
        self.base_fluxes = (
            np.logspace(4, -10, 10_000) if base_fluxes is None else base_fluxes
        )
        if base_cdf is None:
            self.base_cdf = np.zeros(self.base_fluxes.size)
            for index, lower_flux in enumerate(self.base_fluxes[1:], start=1):
                self.base_cdf[index] = self.base_cdf[index - 1] + quad(
                    _franzen_differential_source_count,
                    lower_flux,
                    self.base_fluxes[index - 1],
                )[0]
        else:
            self.base_cdf = base_cdf

        mask = (self.base_fluxes > minimum_flux_jy) & (
            self.base_fluxes < maximum_flux_jy
        )
        self.flux_grid = self.base_fluxes[mask]
        self.cdf = self.base_cdf[mask].copy()
        self.cdf -= self.cdf.min()

    @cached_property
    def source_density_per_steradian(self) -> float:
        """Return the integrated source density over the selected flux range."""
        return float(self.cdf.max())

    @cached_property
    def inverse_cdf(self):
        """Return a spline converting uniform probabilities into fluxes."""
        from scipy.interpolate import InterpolatedUnivariateSpline

        return InterpolatedUnivariateSpline(
            self.cdf / self.source_density_per_steradian,
            self.flux_grid,
        )

    def with_maximum_flux(self, maximum_flux_jy: float) -> "FranzenSourceCounts":
        """Return the distribution below a new bright-source limit."""
        return FranzenSourceCounts(
            self.minimum_flux_jy,
            maximum_flux_jy,
            self.base_fluxes,
            self.base_cdf,
        )

    def above_confusion_limit(self, nside: int) -> "FranzenSourceCounts":
        """Retain fluxes giving no more than one source per HEALPix pixel."""
        import healpy

        pixel_area = healpy.nside2pixarea(nside)
        lower_index = np.argwhere(self.cdf * pixel_area > 1)[0, 0]
        return FranzenSourceCounts(
            self.flux_grid[lower_index],
            self.maximum_flux_jy,
            self.base_fluxes,
            self.base_cdf,
        )

    def draw_full_sky_fluxes(self) -> np.ndarray:
        """Draw a Poisson full-sky source count and its flux densities in Jy."""
        number_of_sources = np.random.poisson(
            4 * np.pi * self.source_density_per_steradian
        )
        return self.inverse_cdf(np.random.uniform(size=number_of_sources))


def generate_gleam_ateam_maps(
    sky_inputs: dict,
    central_jd: float,
    mask_times_jd: np.ndarray,
    telescope_location,
    a_team_catalogue: Path,
    output_directory: Path,
    overwrite: bool = False,
) -> list[Path]:
    """Generate and save GLEAM-like+A-team maps for all requested FoVs.

    Parameters
    ----------
    sky_inputs
        The ``sky`` section of the observation configuration. The
        ``gleam_ateam``
        subsection supplies FoVs, random seed, confusion NSIDE and reference
        frequency.
    central_jd
        Observation midpoint stored with the selected sky model.
    mask_times_jd
        Julian dates used for source selection. This contains only the
        observation midpoint for a central mask, or every integration for a
        drift-scan union mask.
    telescope_location
        Astropy EarthLocation constructed from the observation configuration.
    a_team_catalogue
        CSV containing source name, RA, Dec, Stokes-I flux, reference frequency
        and spectral index for the fixed bright-source catalogue.
    output_directory
        Directory that will contain ``fov-X.skyh5`` files.
    overwrite
        Replace existing maps when true.

    Returns
    -------
    list[pathlib.Path]
        One skyh5 path for every configured GLEAM-like+A-team FoV.

    Notes
    -----
    The catalogue is generated over the full sky before the FoV selection.
    This reproduces the HERA validation-sim statistical recipe rather than
    reading an existing GLEAM catalogue.
    """
    from astropy import units as u
    from astropy.coordinates import AltAz, Latitude, Longitude, SkyCoord
    from astropy.time import Time
    from pyradiosky import SkyModel

    point_inputs = sky_inputs["gleam_ateam"]
    if not a_team_catalogue.exists():
        raise FileNotFoundError(f"Missing A-team catalogue: {a_team_catalogue}")
    output_directory.mkdir(parents=True, exist_ok=True)
    outputs = [
        output_directory
        / (f"fov-{fov:.1f}".rstrip("0").rstrip(".") + ".skyh5")
        for fov in point_inputs["fovs_deg"]
    ]
    if not overwrite and all(path.exists() for path in outputs):
        return outputs

    with a_team_catalogue.open(newline="") as stream:
        ateam_rows = list(csv.DictReader(stream))
    ateam_names = np.array([row["name"] for row in ateam_rows])
    ateam_stokes = np.zeros((4, 1, ateam_names.size)) * u.Jy
    ateam_stokes[0, 0] = np.array(
        [float(row["stokes_i_jy"]) for row in ateam_rows]
    ) * u.Jy
    ateam = SkyModel(
        name=ateam_names,
        ra=Longitude([row["ra"] for row in ateam_rows]),
        dec=Latitude([row["dec"] for row in ateam_rows]),
        stokes=ateam_stokes,
        spectral_type="spectral_index",
        spectral_index=np.array(
            [float(row["spectral_index"]) for row in ateam_rows]
        ),
        reference_frequency=np.array(
            [float(row["reference_frequency_hz"]) for row in ateam_rows]
        ) * u.Hz,
        frame="icrs",
        history="Ten-source A-team model used by HERA validation-sim.",
    )

    np.random.seed(point_inputs["random_seed"])
    source_counts = (
        FranzenSourceCounts(1e-10, 1e4)
        .with_maximum_flux(ateam.stokes[0].min().to_value("Jy"))
        .above_confusion_limit(point_inputs["confusion_nside"])
    )
    fluxes = source_counts.draw_full_sky_fluxes()
    right_ascension = np.random.random(fluxes.size) * 2 * np.pi
    colatitude = np.arccos(np.random.random(fluxes.size) * -2 + 1)
    spectral_index = -np.random.normal(0.8, 0.05, size=fluxes.size)
    gleam_stokes = np.zeros((4, 1, fluxes.size)) * u.Jy
    gleam_stokes[0, 0] = fluxes * u.Jy
    gleam_like = SkyModel(
        name=np.array([f"gl{index:06d}" for index in range(fluxes.size)]),
        ra=Longitude(right_ascension * u.rad),
        dec=Latitude((np.pi / 2 - colatitude) * u.rad),
        stokes=gleam_stokes,
        spectral_type="spectral_index",
        spectral_index=spectral_index,
        reference_frequency=np.full(
            fluxes.size, point_inputs["reference_frequency_hz"]
        ) * u.Hz,
        frame="icrs",
        history=(
            "GLEAM-like statistical catalogue generated from Franzen et al. "
            f"source counts with seed {point_inputs['random_seed']}."
        ),
    )

    full_sky = gleam_like.concat(ateam, inplace=False)
    selected_by_fov = {
        fov: np.zeros(full_sky.name.size, dtype=bool)
        for fov in point_inputs["fovs_deg"]
    }
    for jd in mask_times_jd:
        zenith = SkyCoord(
            az=0 * u.deg,
            alt=90 * u.deg,
            frame=AltAz(
                obstime=Time(jd, format="jd"), location=telescope_location
            ),
        ).transform_to("icrs")
        zenith_angle_deg = full_sky.skycoord.separation(zenith).deg
        for fov in point_inputs["fovs_deg"]:
            selected_by_fov[fov] |= zenith_angle_deg <= fov / 2

    # Store apparent coordinates at the central time for both mask modes.
    full_sky.update_positions(Time(central_jd, format="jd"), telescope_location)
    for fov, output in zip(point_inputs["fovs_deg"], outputs):
        if output.exists() and not overwrite:
            continue
        selected = full_sky.select(
            component_inds=np.where(selected_by_fov[fov])[0],
            inplace=False,
        )
        mask_description = (
            f"central JD {central_jd}"
            if mask_times_jd.size == 1
            else (
                f"the union of {mask_times_jd.size} integrations from JD "
                f"{mask_times_jd[0]} to {mask_times_jd[-1]}"
            )
        )
        selected.history += (
            f"\nFoV diameter {fov} deg selected using {mask_description}."
        )
        print(f"Writing GLEAM-like+A-team map: {output}")
        selected.write_skyh5(output, clobber=overwrite)
    return outputs
