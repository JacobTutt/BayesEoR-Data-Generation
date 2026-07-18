#!/usr/bin/env python3
"""Reproduce Jacob Burba's BayesEoR large-FoV simulation products from scratch."""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from functools import cached_property
from pathlib import Path

import numpy as np


FREQ_MIN_HZ = 150_292_968.75
CHANNEL_WIDTH_HZ = 97_656.25
N_FREQS = 180
JD_START = 2_459_999.067040411
JD_END = 2_459_999.127292119
N_TIMES = 243
INTEGRATION_TIME_S = 21.511353650861533
JD_CENTRE = (JD_START + JD_END) / 2
NSIDE = 128
WHITE_NOISE_SIGMA_K = 6.4e-3
WHITE_NOISE_SEED = 92_381_923
GLEAM_SEED = 42
HERA_LAT_DEG = -30.72152777777791
HERA_LON_DEG = 21.428305555555557
HERA_HEIGHT_M = 1073.0000000093132
HISTORICAL_EOR_FOV_DEG = 12.9080728652


def repo_root() -> Path:
    return Path(__file__).resolve().parent


def default_output_root() -> Path:
    return repo_root() / "data" / "large-fov"


def fov_map_label(fov: float) -> str:
    return f"{fov:.1f}"


def fov_short_label(fov: float) -> str:
    return fov_map_label(fov).rstrip("0").rstrip(".")


def frequency_axis_hz() -> np.ndarray:
    return FREQ_MIN_HZ + np.arange(N_FREQS) * CHANNEL_WIDTH_HZ


def hera_location():
    from astropy import units as u
    from astropy.coordinates import EarthLocation

    return EarthLocation.from_geodetic(
        HERA_LON_DEG * u.deg,
        HERA_LAT_DEG * u.deg,
        HERA_HEIGHT_M * u.m,
    )


def _sky_model_from_file(path: Path):
    from pyradiosky import SkyModel

    if hasattr(SkyModel, "from_file"):
        return SkyModel.from_file(path)
    sky = SkyModel()
    sky.read_skyh5(path)
    return sky


def _write_sky_model(sky, path: Path, overwrite: bool) -> None:
    if path.exists() and not overwrite:
        print(f"Exists, keeping: {path}")
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    print(f"Writing: {path}")
    sky.write_skyh5(path, clobber=overwrite)


def white_noise_parent_path(root: Path) -> Path:
    return (
        root
        / "sims/sky-models/white-noise/nside-128/h1c-idr2/band2"
        / "sigma-6.40mK-seed-92381923-all-sky.skyh5"
    )


def gsm_parent_path(root: Path) -> Path:
    return root / "sims/sky-models/gsm/nside-128/h1c-idr2/band2/all-sky.skyh5"


def ptsrc_parent_dir(root: Path) -> Path:
    return root / "sims/sky-models/ptsrc"


def make_white_noise_parent(root: Path, overwrite: bool) -> Path:
    """Generate the 180-channel, NSIDE-128 white-noise EoR-like parent."""
    from astropy import units as u
    from pyradiosky import SkyModel

    path = white_noise_parent_path(root)
    if path.exists() and not overwrite:
        print(f"Exists, keeping: {path}")
        return path

    npix = 12 * NSIDE**2
    # Use NumPy's legacy RNG deliberately: this is Burba's original recipe.
    np.random.seed(WHITE_NOISE_SEED)
    stokes_i = np.random.normal(
        0.0, WHITE_NOISE_SIGMA_K, (N_FREQS, npix)
    )
    stokes_i -= stokes_i.mean(axis=1)[:, None]

    stokes = np.zeros((4, N_FREQS, npix)) * u.K
    stokes[0] = stokes_i * u.K
    history = (
        f"NSIDE={NSIDE} white-noise sky drawn with numpy.random.normal; "
        f"sigma={WHITE_NOISE_SIGMA_K} K; seed={WHITE_NOISE_SEED}; each "
        "frequency map was shifted to zero spatial mean."
    )
    sky = SkyModel(
        nside=NSIDE,
        hpx_inds=np.arange(npix),
        stokes=stokes,
        spectral_type="full",
        freq_array=frequency_axis_hz() * u.Hz,
        history=history,
        frame="icrs",
        hpx_order="ring",
    )
    _write_sky_model(sky, path, overwrite=True)
    return path


def make_gsm_parent(root: Path, overwrite: bool) -> Path:
    """Generate GSM2008 with pygdsm and reproduce Burba's coordinate transform."""
    import healpy
    from astropy import units as u
    from pygdsm import GlobalSkyModel
    from pyradiosky import SkyModel

    path = gsm_parent_path(root)
    if path.exists() and not overwrite:
        print(f"Exists, keeping: {path}")
        return path

    freqs_hz = frequency_axis_hz()
    gsm_native = np.asarray(
        GlobalSkyModel(include_cmb=True).generate(freqs_hz / 1e6)
    )
    if gsm_native.ndim == 1:
        gsm_native = gsm_native[None, :]

    npix = 12 * NSIDE**2
    gsm_nside = healpy.npix2nside(gsm_native.shape[1])
    if gsm_nside != NSIDE:
        gsm = np.empty((N_FREQS, npix), dtype=gsm_native.dtype)
        for index in range(N_FREQS):
            gsm[index] = healpy.ud_grade(gsm_native[index], NSIDE)
    else:
        gsm = gsm_native

    stokes = np.zeros((4, N_FREQS, npix)) * u.K
    stokes[0] = gsm * u.K
    sky = SkyModel(
        nside=NSIDE,
        hpx_inds=np.arange(npix),
        stokes=stokes,
        spectral_type="full",
        freq_array=freqs_hz * u.Hz,
        history=(
            "GSM2008 generated with pygdsm.GlobalSkyModel(include_cmb=True), "
            f"downgraded to NSIDE={NSIDE}, then interpolated from Galactic to ICRS."
        ),
        frame="galactic",
        hpx_order="ring",
    )
    sky.healpix_interp_transform(frame="icrs", full_sky=True)
    _write_sky_model(sky, path, overwrite=True)
    return path


def _random_sphere(count: int) -> tuple[np.ndarray, np.ndarray]:
    longitude = np.random.random(count) * (2 * np.pi)
    uniform_cos_colatitude = np.random.random(count) * -2.0 + 1.0
    colatitude = np.arccos(uniform_cos_colatitude)
    return colatitude, longitude


def _franzen_dnds(flux_jy: float | np.ndarray) -> float | np.ndarray:
    coefficients = [3.52, 0.307, -0.388, -0.0404, 0.0351, 0.00600]
    log_flux = np.log10(flux_jy)
    euclidean_normalised = 10 ** sum(
        coefficient * log_flux**power
        for power, coefficient in enumerate(coefficients)
    )
    return flux_jy**-2.5 * euclidean_normalised


class FranzenSourceCounts:
    """Numerical Franzen et al. (2019) source-count sampler."""

    def __init__(
        self,
        smin: float,
        smax: float,
        base_fluxes: np.ndarray | None = None,
        base_cdf: np.ndarray | None = None,
    ) -> None:
        from scipy.integrate import quad

        self.smin = smin
        self.smax = smax
        self._base_fluxes = (
            np.logspace(4, -10, 10_000) if base_fluxes is None else base_fluxes
        )
        if base_cdf is None:
            self._base_cdf = np.zeros(self._base_fluxes.size)
            for index, lower_flux in enumerate(self._base_fluxes[1:], start=1):
                integral = quad(
                    _franzen_dnds,
                    lower_flux,
                    self._base_fluxes[index - 1],
                )[0]
                self._base_cdf[index] = integral + self._base_cdf[index - 1]
        else:
            self._base_cdf = base_cdf

        self._mask = (self._base_fluxes > smin) & (self._base_fluxes < smax)
        self.flux_grid = self._base_fluxes[self._mask]
        self.cdf = self._base_cdf[self._mask].copy()
        self.cdf -= self.cdf.min()

    @cached_property
    def source_density(self) -> float:
        return float(self.cdf.max())

    @cached_property
    def inverse_cdf(self):
        from scipy.interpolate import InterpolatedUnivariateSpline

        return InterpolatedUnivariateSpline(self.cdf / self.source_density, self.flux_grid)

    def with_bounds(self, smax: float) -> "FranzenSourceCounts":
        return FranzenSourceCounts(
            self.smin, smax, self._base_fluxes, self._base_cdf
        )

    def with_one_source_per_pixel(self, nside: int) -> "FranzenSourceCounts":
        import healpy

        pixel_area = healpy.nside2pixarea(nside)
        index = np.argwhere(self.cdf * pixel_area > 1)[0, 0]
        return FranzenSourceCounts(
            self.flux_grid[index],
            self.smax,
            self._base_fluxes,
            self._base_cdf,
        )

    def sample_full_sky_fluxes(self) -> np.ndarray:
        count = np.random.poisson(4 * np.pi * self.source_density)
        return self.inverse_cdf(np.random.uniform(size=count))


def make_ateam_model():
    from astropy import units as u
    from astropy.coordinates import Latitude, Longitude
    from pyradiosky import SkyModel

    names = np.array(
        [
            "3C444", "CentaurusA", "HydraA", "PictorA", "HerculesA",
            "VirgoA", "Crab", "CygnusA", "CassiopeiaA", "FornaxA",
        ]
    )
    ra = Longitude(
        [
            "22h14m16s", "13h25m28s", "09h18m06s", "05h19m50s",
            "16h51m08s", "12h30m49s", "05h34m32s", "19h59m28s",
            "23h23m28s", "3h22m42s",
        ]
    )
    dec = Latitude(
        [
            "-17d01m36s", "-43d01m09s", "-12d05m44s", "-45d46m44s",
            "04d59m33s", "12d23m28s", "22d00m52s", "40d44m02s",
            "58d48m42s", "-37d12m2s",
        ]
    )
    stokes = np.zeros((4, 1, names.size)) * u.Jy
    stokes[0, 0] = np.array(
        [60, 1370, 280, 390, 377, 861, 1340, 7920, 11900, 750]
    ) * u.Jy
    reference_frequency = np.array([200] * 9 + [154]) * 1e6 * u.Hz
    spectral_index = np.array(
        [-0.96, -0.50, -0.96, -0.99, -1.07, -0.86, -0.22, -0.78, -0.41, -0.825]
    )
    return SkyModel(
        name=names,
        ra=ra,
        dec=dec,
        stokes=stokes,
        spectral_type="spectral_index",
        spectral_index=spectral_index,
        reference_frequency=reference_frequency,
        history="Ten-source A-team model used by HERA validation-sim.",
        frame="icrs",
    )


def make_gleam_like_model(max_flux_density_jy: float):
    from astropy import units as u
    from astropy.coordinates import Latitude, Longitude
    from pyradiosky import SkyModel

    np.random.seed(GLEAM_SEED)
    source_counts = (
        FranzenSourceCounts(1e-10, 1e4)
        .with_bounds(max_flux_density_jy)
        .with_one_source_per_pixel(nside=256)
    )
    flux = source_counts.sample_full_sky_fluxes()
    colatitude, longitude = _random_sphere(flux.size)
    # The historical helper drew +0.8 and Burba's notebook then negated it.
    spectral_index = -np.random.normal(0.8, 0.05, size=flux.size)
    names = np.array([f"gl{index:06d}" for index in range(flux.size)])
    stokes = np.zeros((4, 1, flux.size)) * u.Jy
    stokes[0, 0] = flux * u.Jy
    return SkyModel(
        name=names,
        ra=Longitude(longitude * u.rad),
        dec=Latitude((np.pi / 2 - colatitude) * u.rad),
        stokes=stokes,
        spectral_type="spectral_index",
        reference_frequency=np.full(flux.size, 154e6) * u.Hz,
        spectral_index=spectral_index,
        history=(
            "Synthetic GLEAM-like catalogue from Franzen et al. (2019) source "
            f"counts; seed={GLEAM_SEED}; NSIDE=256 confusion threshold."
        ),
        frame="icrs",
    )


def make_point_source_parents(root: Path, overwrite: bool) -> tuple[Path, Path, Path]:
    directory = ptsrc_parent_dir(root)
    ateam_path = directory / "a-team.skyh5"
    gleam_path = directory / "gleam-like.skyh5"
    combined_path = directory / "gleam-like-with-ateam.skyh5"
    if all(path.exists() for path in (ateam_path, gleam_path, combined_path)) and not overwrite:
        print(f"Point-source parents exist under: {directory}")
        return ateam_path, gleam_path, combined_path

    ateam = make_ateam_model()
    gleam = make_gleam_like_model(ateam.stokes[0].min().to_value("Jy"))
    combined = gleam.concat(ateam, inplace=False)
    _write_sky_model(ateam, ateam_path, overwrite=True)
    _write_sky_model(gleam, gleam_path, overwrite=True)
    _write_sky_model(combined, combined_path, overwrite=True)
    return ateam_path, gleam_path, combined_path


def make_parents(root: Path, overwrite: bool) -> None:
    root.mkdir(parents=True, exist_ok=True)
    make_white_noise_parent(root, overwrite)
    make_gsm_parent(root, overwrite)
    make_point_source_parents(root, overwrite)


def healpix_fov_indices(nside: int, fov_deg: float) -> np.ndarray:
    from astropy import units as u
    from astropy.coordinates import AltAz, SkyCoord
    from astropy.time import Time
    from astropy_healpix import healpy as hp

    all_pixels = np.arange(12 * nside**2)
    ra, dec = hp.pix2ang(nside, all_pixels, lonlat=True)
    horizontal = SkyCoord(ra * u.deg, dec * u.deg, frame="icrs").transform_to(
        AltAz(obstime=Time(JD_CENTRE, format="jd"), location=hera_location())
    )
    zenith_angle = np.pi / 2 - horizontal.alt.rad
    return all_pixels[zenith_angle <= np.deg2rad(fov_deg / 2)]


def _cut_healpix(parent_path: Path, output_path: Path, fov: float, overwrite: bool) -> None:
    parent = _sky_model_from_file(parent_path)
    indices = healpix_fov_indices(int(parent.nside), fov)
    subset = parent.select(component_inds=indices, inplace=False)
    subset.history += (
        f"\nSelected pixel centres with zenith angle <= {fov / 2:.12g} deg "
        f"at JD {JD_CENTRE}."
    )
    _write_sky_model(subset, output_path, overwrite)


def _cut_point_sources(parent_path: Path, output_path: Path, fov: float, overwrite: bool) -> None:
    from astropy.time import Time

    parent = _sky_model_from_file(parent_path)
    parent.update_positions(Time(JD_CENTRE, format="jd"), hera_location())
    zenith_angle_deg = 90 - np.rad2deg(parent.alt_az[0])
    indices = np.where(zenith_angle_deg <= fov / 2)[0]
    subset = parent.select(component_inds=indices, inplace=False)
    subset.history += (
        f"\nSelected sources with zenith angle <= {fov / 2:.12g} deg "
        f"at JD {JD_CENTRE}."
    )
    _write_sky_model(subset, output_path, overwrite)


def make_cutouts(root: Path, eor_fov: float, foreground_fovs: list[float], overwrite: bool) -> None:
    eor_output = (
        root
        / "sims/sky-models/white-noise/nside-128/h1c-idr2/band2"
        / f"sigma-6.40mK-seed-92381923-fov-{fov_map_label(eor_fov)}deg.skyh5"
    )
    _cut_healpix(white_noise_parent_path(root), eor_output, eor_fov, overwrite)

    for fov in foreground_fovs:
        gsm_output = (
            root / "sims/sky-models/gsm/nside-128/h1c-idr2/band2/field1"
            / f"fov-{fov_map_label(fov)}deg.skyh5"
        )
        ptsrc_output = (
            root / "sims/sky-models/ptsrc/h1c-idr2/field1"
            / f"fov-{fov_map_label(fov)}deg.skyh5"
        )
        _cut_healpix(gsm_parent_path(root), gsm_output, fov, overwrite)
        _cut_point_sources(
            ptsrc_parent_dir(root) / "gleam-like-with-ateam.skyh5",
            ptsrc_output,
            fov,
            overwrite,
        )


def _write_obsparam(path: Path, catalog: Path, output_dir: Path, output_name: str, overwrite: bool) -> None:
    import yaml

    if path.exists() and not overwrite:
        print(f"Exists, keeping: {path}")
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)
    instrument = repo_root() / "instrument"
    params = {
        "filing": {
            "outdir": str(output_dir.resolve()),
            "outfile_name": output_name,
            "output_format": "uvh5",
            "clobber": bool(overwrite),
        },
        "freq": {
            "Nfreqs": N_FREQS,
            "channel_width": CHANNEL_WIDTH_HZ,
            "start_freq": FREQ_MIN_HZ,
        },
        "sources": {"catalog": str(catalog.resolve())},
        "telescope": {
            "array_layout": str((instrument / "hex-37-14.6m.csv").resolve()),
            "telescope_config_name": str((instrument / "hex-37-14.6m-airy-diam14.2m.yml").resolve()),
        },
        "time": {
            "integration_time": INTEGRATION_TIME_S,
            "start_time": JD_START,
            "end_time": JD_END,
        },
        "select": {"redundant_threshold": 1.0},
    }
    print(f"Writing: {path}")
    with path.open("w") as stream:
        yaml.safe_dump(params, stream, sort_keys=False)


def make_obsparams(root: Path, eor_fov: float, foreground_fovs: list[float], overwrite: bool) -> None:
    eor_label = fov_short_label(eor_fov)
    eor_catalog = (
        root / "sims/sky-models/white-noise/nside-128/h1c-idr2/band2"
        / f"sigma-6.40mK-seed-92381923-fov-{fov_map_label(eor_fov)}deg.skyh5"
    )
    _write_obsparam(
        root / "sims/airy/white-noise/h1c-idr2/band2/obsparams" / f"fov-{eor_label}deg.yaml",
        eor_catalog,
        root / "sims/airy/white-noise/h1c-idr2/band2/vis",
        f"fov-{eor_label}",
        overwrite,
    )
    for component in ("gsm", "ptsrc"):
        for fov in foreground_fovs:
            label = fov_short_label(fov)
            if component == "gsm":
                catalog = root / "sims/sky-models/gsm/nside-128/h1c-idr2/band2/field1" / f"fov-{fov_map_label(fov)}deg.skyh5"
            else:
                catalog = root / "sims/sky-models/ptsrc/h1c-idr2/field1" / f"fov-{fov_map_label(fov)}deg.skyh5"
            base = root / f"sims/airy/{component}/h1c-idr2/band2/field1"
            _write_obsparam(
                base / "obsparams" / f"fov-{label}deg.yaml",
                catalog,
                base / "vis",
                f"fov-{label}",
                overwrite,
            )


def obsparam_paths(root: Path, eor_fov: float, foreground_fovs: list[float], components: list[str]) -> list[Path]:
    paths: list[Path] = []
    if "white-noise" in components:
        paths.append(root / "sims/airy/white-noise/h1c-idr2/band2/obsparams" / f"fov-{fov_short_label(eor_fov)}deg.yaml")
    for component in ("gsm", "ptsrc"):
        if component in components:
            paths.extend(
                root / f"sims/airy/{component}/h1c-idr2/band2/field1/obsparams" / f"fov-{fov_short_label(fov)}deg.yaml"
                for fov in foreground_fovs
            )
    return paths


def simulate(root: Path, eor_fov: float, foreground_fovs: list[float], components: list[str]) -> None:
    from pyuvsim.uvsim import run_uvsim

    for obsparam in obsparam_paths(root, eor_fov, foreground_fovs, components):
        if not obsparam.exists():
            raise FileNotFoundError(f"Run the obsparams stage first: {obsparam}")
        print(f"Running pyuvsim: {obsparam}", flush=True)
        run_uvsim(str(obsparam))


def _visibility_path(root: Path, component: str, fov: float) -> Path:
    label = fov_short_label(fov)
    if component == "white-noise":
        return root / "sims/airy/white-noise/h1c-idr2/band2/vis" / f"fov-{label}.uvh5"
    return root / f"sims/airy/{component}/h1c-idr2/band2/field1/vis" / f"fov-{label}.uvh5"


def sum_visibilities(root: Path, eor_fov: float, foreground_fovs: list[float], overwrite: bool) -> None:
    from pyuvdata import UVData

    eor_path = _visibility_path(root, "white-noise", eor_fov)
    for fov in foreground_fovs:
        inputs = [eor_path, _visibility_path(root, "gsm", fov), _visibility_path(root, "ptsrc", fov)]
        for path in inputs:
            if not path.exists():
                raise FileNotFoundError(path)
        skies = [UVData.from_file(path) for path in inputs]
        reference = skies[0]
        for other in skies[1:]:
            for field in ("time_array", "freq_array", "ant_1_array", "ant_2_array", "polarization_array"):
                if not np.array_equal(getattr(reference, field), getattr(other, field)):
                    raise ValueError(f"Visibility axes differ for {field}: {inputs}")
        reference.data_array = sum(sky.data_array for sky in skies)
        reference.flag_array = np.logical_or.reduce([sky.flag_array for sky in skies])
        reference.nsample_array = np.minimum.reduce([sky.nsample_array for sky in skies])
        reference.history += "\nSum of white-noise EoR, GSM, and GLEAM-like+A-team simulations."
        eor_label = fov_short_label(eor_fov)
        fg_label = fov_short_label(fov)
        if fg_label == eor_label:
            name = f"wn-gsm-ptsrc-fov-{fg_label}.uvh5"
        else:
            name = f"wn-fov-{eor_label}-gsm-ptsrc-fov-{fg_label}.uvh5"
        output = root / "sims/airy/sum/h1c-idr2/band2/field1" / name
        if output.exists() and not overwrite:
            print(f"Exists, keeping: {output}")
            continue
        output.parent.mkdir(parents=True, exist_ok=True)
        print(f"Writing: {output}")
        reference.write_uvh5(output, clobber=overwrite)


def preprocess(root: Path, eor_fov: float, foreground_fovs: list[float], overwrite: bool) -> None:
    """Run the historical BayesEoR preprocessor from its pinned external checkout."""
    historical_repo = repo_root() / "external" / "BayesEoR"
    script = historical_repo / "scripts" / "data_preprocessing.py"
    if not script.exists():
        raise FileNotFoundError(
            f"Missing {script}. Run ./bootstrap-historical-dependencies.sh first."
        )
    data_dir = root / "sims/airy/sum/h1c-idr2/band2/field1"
    proc_dir = data_dir / "proc"
    proc_dir.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    shim = repo_root() / "historical_shim"
    env["PYTHONPATH"] = str(shim) + os.pathsep + env.get("PYTHONPATH", "")
    instrument_root = root / "inst-models"
    for index, fov in enumerate(foreground_fovs):
        eor_label = fov_short_label(eor_fov)
        fg_label = fov_short_label(fov)
        filename = (
            f"wn-gsm-ptsrc-fov-{fg_label}.uvh5"
            if fg_label == eor_label
            else f"wn-fov-{eor_label}-gsm-ptsrc-fov-{fg_label}.uvh5"
        )
        command = [
            sys.executable, str(script),
            "--data_path", str(data_dir),
            "--filename", filename,
            "--save_dir", str(proc_dir),
            "--inst_model_dir", str(instrument_root),
            "--telescope_name", "hex-37-14.6m",
            "--bl_cutoff_m", "40",
            "--start_freq_MHz", "157.2265625",
            "--nf", "39",
            "--nt", "47",
            "--start_int", "98",
            "--form_pI",
        ]
        if index == 0:
            # The first case creates the common UVW/redundancy model from scratch.
            command.append("--save_model")
        if overwrite:
            command.append("--clobber")
        print("Running:", " ".join(command), flush=True)
        subprocess.run(command, check=True, env=env)


def verify(root: Path, reference_root: Path, eor_fov: float, foreground_fovs: list[float]) -> None:
    """Compare regenerated products with references; references are never inputs."""
    pairs: list[tuple[Path, Path]] = [
        (white_noise_parent_path(root), white_noise_parent_path(reference_root)),
        (gsm_parent_path(root), gsm_parent_path(reference_root)),
        (ptsrc_parent_dir(root) / "gleam-like.skyh5", ptsrc_parent_dir(reference_root) / "gleam-like.skyh5"),
        (ptsrc_parent_dir(root) / "a-team.skyh5", ptsrc_parent_dir(reference_root) / "a-team.skyh5"),
    ]
    for fov in [eor_fov, *foreground_fovs]:
        if fov == eor_fov:
            rel = Path("sims/sky-models/white-noise/nside-128/h1c-idr2/band2") / f"sigma-6.40mK-seed-92381923-fov-{fov_map_label(fov)}deg.skyh5"
            pairs.append((root / rel, reference_root / rel))
        for relbase in (
            Path("sims/sky-models/gsm/nside-128/h1c-idr2/band2/field1"),
            Path("sims/sky-models/ptsrc/h1c-idr2/field1"),
        ):
            rel = relbase / f"fov-{fov_map_label(fov)}deg.skyh5"
            if (reference_root / rel).exists():
                pairs.append((root / rel, reference_root / rel))

    failures = 0
    for generated_path, reference_path in pairs:
        if not generated_path.exists() or not reference_path.exists():
            print(f"MISSING {generated_path} or {reference_path}")
            failures += 1
            continue
        generated = _sky_model_from_file(generated_path)
        reference = _sky_model_from_file(reference_path)
        same_components = generated.Ncomponents == reference.Ncomponents
        same_freqs = np.allclose(generated.freq_array.value, reference.freq_array.value, rtol=0, atol=1e-8)
        same_stokes = np.allclose(generated.stokes.value, reference.stokes.value, rtol=1e-12, atol=1e-10)
        result = same_components and same_freqs and same_stokes
        failures += int(not result)
        print(f"{'PASS' if result else 'FAIL'} {generated_path.relative_to(root)}")
    if failures:
        raise SystemExit(f"{failures} regression comparison(s) failed")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", type=Path, default=default_output_root())
    parser.add_argument("--eor-fov", type=float, default=HISTORICAL_EOR_FOV_DEG)
    parser.add_argument(
        "--foreground-fovs", type=float, nargs="+",
        default=[HISTORICAL_EOR_FOV_DEG, 30.0, 60.0, 90.0, 120.0],
    )
    parser.add_argument(
        "--components", nargs="+", choices=["white-noise", "gsm", "ptsrc"],
        default=["white-noise", "gsm", "ptsrc"],
        help="Components to run during the simulate stage.",
    )
    parser.add_argument("--overwrite", action="store_true")
    subparsers = parser.add_subparsers(dest="stage", required=True)
    subparsers.add_parser("parents", help="Generate all three full-sky parents from scratch.")
    subparsers.add_parser("cutouts", help="Select requested FoVs from newly generated parents.")
    subparsers.add_parser("obsparams", help="Write pyuvsim observation YAML files.")
    subparsers.add_parser("simulate", help="Run pyuvsim for selected components.")
    subparsers.add_parser("sum", help="Sum EoR, GSM, and point-source visibilities.")
    subparsers.add_parser("preprocess", help="Create the historical 39x47 inference vectors.")
    subparsers.add_parser("prepare", help="Run parents, cutouts, and obsparams from scratch.")
    check = subparsers.add_parser("verify", help="Regression-check against reference files.")
    check.add_argument("--reference-root", type=Path, required=True)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    root = args.output_root.resolve()
    if args.stage in ("parents", "prepare"):
        make_parents(root, args.overwrite)
    if args.stage in ("cutouts", "prepare"):
        make_cutouts(root, args.eor_fov, args.foreground_fovs, args.overwrite)
    if args.stage in ("obsparams", "prepare"):
        make_obsparams(root, args.eor_fov, args.foreground_fovs, args.overwrite)
    if args.stage == "simulate":
        simulate(root, args.eor_fov, args.foreground_fovs, args.components)
    elif args.stage == "sum":
        sum_visibilities(root, args.eor_fov, args.foreground_fovs, args.overwrite)
    elif args.stage == "preprocess":
        preprocess(root, args.eor_fov, args.foreground_fovs, args.overwrite)
    elif args.stage == "verify":
        verify(root, args.reference_root.resolve(), args.eor_fov, args.foreground_fovs)


if __name__ == "__main__":
    main()
