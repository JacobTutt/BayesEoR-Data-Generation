# BayesEoR Data Generation

Simulated sky models and interferometric visibilities for
[BayesEoR](https://github.com/JacobTutt/BayesEoR-JAX), generated from a single
YAML configuration.

The package builds three sky components, each for one or more fields of view,
and simulates their visibilities with
[pyuvsim](https://github.com/RadioAstronomySoftwareGroup/pyuvsim) using a
shared telescope, time and frequency setup:

- **21-cm signal:** a Gaussian white-noise sky, with each channel shifted to
  zero mean.
- **Diffuse foregrounds:** the Global Sky Model (GSM2008), through pygdsm.
- **Point-source foregrounds:** a GLEAM-like population drawn from the
  Franzen et al. (2019) source counts, plus fixed bright A-team sources.

Each component is written to its own visibility file, so any combination can be
summed for an analysis.

## Installation

Python 3.10 or newer and an MPI implementation are required. From the
repository root:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e .
```

## Getting started

Set `output_directory` in `example_config.yaml`, then run:

```bash
python example_run.py
```

This generates every sky model in the configuration and, when
`run.simulate_visibilities` is `true`, their visibilities. The example is large,
with 180 channels and five fields of view, so reduce `frequency.n_channels`, the
time range or the field-of-view lists for a quick test.

Existing outputs are kept by default. Set `run.overwrite_existing_files: true`
to replace them.

## Configuration

`example_config.yaml` lists every input. Relative paths are resolved from the
directory containing it.

### Run and output

- **`run.simulate_visibilities`:** also runs pyuvsim after generating the sky
  models.
- **`run.overwrite_existing_files`:** replaces existing sky, observation and
  visibility files. Filenames depend only on the component and field of view,
  so when changing the observation, use a new `output_directory` or enable
  this.
- **`output_directory`:** the root for all outputs.

### Telescope

- **`telescope.name`:** written into the telescope metadata.
- **`telescope.antenna_layout`:** a pyuvsim antenna-layout file, described
  below.
- **`latitude_deg`, `longitude_deg`, `altitude_m`:** the array location. This
  sets both the observation and the zenith-centred sky mask, so changing it
  requires regenerating the sky models.
- **`beam.type`, `beam.diameter_m`:** the pyuvsim analytic beam, shared by every
  antenna, and its aperture diameter in metres. The example uses `airy`.
- **`redundant_threshold_m`:** the position tolerance pyuvsim uses to identify
  redundant baselines.

The antenna layout is a whitespace-separated table with one row per antenna:

```text
Name   Number   BeamID   E      N      U
ANT0   0        0        0.0    0.0    0.0
ANT1   1        0        14.6   0.0    0.0
ANT2   2        0        0.0    14.6   0.0
```

`Name` and `Number` must be unique, and `E`, `N` and `U` are positions in metres
relative to the array centre. `BeamID` is 0 for every antenna, since one beam is
shared. `default_configs/hex-37-14.6m.csv` is a complete 37-antenna example.

### Time and frequency

- **`time.start_jd`, `time.end_jd`:** the first and last integration times, in
  Julian days.
- **`time.integration_time_s`:** the integration length in seconds. This also
  sets the time sampling of the `union` sky mask.
- **`frequency.start_hz`, `channel_width_hz`, `n_channels`:** the first channel,
  the channel spacing and the number of channels, used for both the sky spectra
  and the visibilities.

### Sky models

Each component has its own `fovs_deg` list, so several fields of view can be
generated in one run. A field of view is a full angular diameter around zenith:
30° includes everything within 15° of zenith.

- **`sky.nside`:** the HEALPix NSIDE of the 21-cm and GSM maps. Higher values
  increase the resolution, memory use and runtime.
- **`sky.sky_mask_time_mode`:** the sky kept by the field-of-view mask.
  - `center` (default): the sky visible at the midpoint of the observation.
  - `union`: every pixel or source that enters the field of view at any
    integration.

  Both modes write the same filenames, so use a new output directory or enable
  overwriting when switching.
- **`sky.eor.rms_k`, `random_seed`:** the RMS in kelvin, before each channel is
  shifted to zero mean, and the seed of the realisation.
- **`sky.gsm.include_cmb`:** passed to GSM2008. `frequency_chunk_size` sets how
  many channels are evaluated together, which affects memory use but not the
  sky.
- **`sky.gleam_ateam.a_team_catalogue`:** the bright-source catalogue, with
  columns `name`, `ra`, `dec`, `stokes_i_jy`, `reference_frequency_hz` and
  `spectral_index`. `default_configs/a_team.csv` is the supplied example.
- **`sky.gleam_ateam.random_seed`, `confusion_nside`, `reference_frequency_hz`:**
  the seed of the faint-source population, the HEALPix scale that sets its
  confusion flux limit, and the reference frequency of its sources.

## Outputs

Products are grouped by component under `output_directory`:

```text
sky_models/{eor,gsm,gleam_ateam}/fov-*.skyh5
observation_files/{eor,gsm,gleam_ateam}/fov-*.yaml
visibilities/{eor,gsm,gleam_ateam}/fov-*.uvh5
```

The observation YAML files are the generated pyuvsim inputs. Field-of-view
labels are rounded to one decimal place, for example `fov-12.9`.

## Python interface

The workflow can also be run from Python, with a configuration dictionary in the
same form as `example_config.yaml`. Relative paths are then resolved from the
current directory.

```python
import yaml
from bayeseor_data_generation import generate_skies, simulate_visibilities

with open("example_config.yaml") as stream:
    config = yaml.safe_load(stream)
config.pop("run")

sky_paths = generate_skies(config, overwrite=False)
visibility_paths = simulate_visibilities(config, overwrite=False)
```

### Summing components

Sum the components needed for an analysis with `sum_visibility_files`:

```python
from bayeseor_data_generation import sum_visibility_files

sum_visibility_files(
    input_files=[
        "data/h1c_band2/visibilities/eor/fov-12.9.uvh5",
        "data/h1c_band2/visibilities/gsm/fov-30.uvh5",
        "data/h1c_band2/visibilities/gleam_ateam/fov-30.uvh5",
    ],
    output_file="data/h1c_band2/visibilities/total/eor-12.9-fgs-30.uvh5",
)
```

The inputs must share the same time, frequency, baseline and polarisation
sampling. Pass `overwrite=True` to replace an existing output.

## Using with BayesEoR

The `.uvh5` files can be passed directly to BayesEoR as `data_path`, with its
frequency, time and telescope settings matching the configuration used here.
pyuvsim stores half of Stokes I in each of XX and YY, which matches BayesEoR's
default `pi_norm: 1`.
