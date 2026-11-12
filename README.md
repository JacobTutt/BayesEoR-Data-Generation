# BayesEoR Data Generation

Generate sky models and interferometric visibilities for BayesEoR-style
simulations from one YAML configuration.

The package creates three sky components:

- an EoR-like Gaussian white-noise sky;
- diffuse foreground emission from GSM2008;
- a statistical GLEAM-like point-source population plus fixed bright A-team
  sources.

Each component can be generated for one or more fields of view. The resulting
sky models can then be passed through `pyuvsim` using the same telescope, time
and frequency settings.

## Installation

Python 3.10 or later and an MPI implementation are required. Create an
environment and install the package with:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e .
```

## Quick start

Edit `example_config.yaml`, particularly `output_directory`, then run:

```bash
python example_run.py
```

The example generates all requested sky models and, when
`run.simulate_visibilities` is `true`, their visibility simulations. It is a
large example, so reduce the number of frequency channels, time range or fields
of view for a quicker test.

Existing files are kept by default. Set
`run.overwrite_existing_files: true` to replace them.

## Configuration

`example_config.yaml` contains every supported input. Relative paths are
resolved from the directory containing that file.

### Run and output

- `run.simulate_visibilities` controls whether the workflow stops after making
  the sky models or also runs `pyuvsim`.
- `run.overwrite_existing_files` replaces existing sky, observation and
  visibility files when `true`. When changing an observation, use a new
  `output_directory` or enable this option; otherwise matching filenames are
  deliberately kept.
- `output_directory` is the root for all generated products.

### Telescope and antenna layout

- `telescope.name` is written into the generated telescope metadata.
- `telescope.antenna_layout` points to a `pyuvsim` antenna-layout file. Any
  number or arrangement of antennas can be used.
- `latitude_deg`, `longitude_deg` and `altitude_m` set the array location. The
  location affects both the `pyuvsim` observation and the zenith-centred sky
  mask, so changing it requires regenerating the sky models.
- `beam.type` selects the `pyuvsim` analytic beam. The supplied configuration
  uses `airy`; `beam.diameter_m` is its aperture diameter in metres. The current
  workflow assigns this one beam to every antenna.
- `redundant_threshold_m` is the positional tolerance passed to `pyuvsim` when
  identifying redundant baselines.

The antenna layout is a whitespace-separated table with one row per antenna:

```text
Name   Number   BeamID   E      N      U
ANT0   0        0        0.0    0.0    0.0
ANT1   1        0        14.6   0.0    0.0
ANT2   2        0        0.0    14.6   0.0
```

`Name` and `Number` must identify each antenna uniquely. `E`, `N` and `U` are
positions in metres relative to the array centre. Use `BeamID` 0 for every
antenna because this workflow currently creates one shared beam definition.
`default_configs/hex-37-14.6m.csv` is a complete example.

### Time and frequency

- `time.start_jd` and `time.end_jd` are the first and last integration times in
  Julian days.
- `time.integration_time_s` is the integration length in seconds. It also sets
  the time sampling used by the `union` sky mask.
- `frequency.start_hz` is the first channel frequency, while
  `channel_width_hz` is the channel spacing and `n_channels` is the total
  number of channels. These values define both the generated sky spectra and
  the visibility simulation.

### Sky models

- `sky.nside` is the HEALPix NSIDE used for the EoR and GSM maps. Increasing it
  increases angular resolution, memory use and runtime.
- `sky.sky_mask_time_mode` is either `center` or `union`, as described below.
- Each component's `fovs_deg` is a list, so multiple fields of view can be
  generated in one run and the three components do not need identical lists.
- `sky.eor.rms_k` is the Gaussian RMS in kelvin before each full-sky frequency
  slice is shifted to zero mean. `random_seed` makes the realization
  repeatable.
- `sky.gsm.include_cmb` is passed to GSM2008. `frequency_chunk_size` controls
  how many channels are evaluated together; it changes memory use, not the sky
  model.
- `sky.gleam_ateam.a_team_catalogue` is the fixed bright-source catalogue. Its
  columns are `name`, `ra`, `dec`, `stokes_i_jy`, `reference_frequency_hz` and
  `spectral_index`; `default_configs/a_team.csv` is the supplied example.
- `sky.gleam_ateam.random_seed` makes the statistical faint-source catalogue
  repeatable. `confusion_nside` sets the HEALPix scale used to determine its
  confusion flux limit, and `reference_frequency_hz` is the reference
  frequency assigned to those generated sources.

Fields of view are full angular diameters around zenith. For example, a
30-degree field of view includes sky positions within 15 degrees of zenith.

The `sky.sky_mask_time_mode` option controls how the field-of-view mask is
applied:

- `center` selects the sky visible at the midpoint of the observation. This is
  the default when the option is omitted.
- `union` selects every pixel or point source that enters the field of view at
  any integration during the observation.

Use a new output directory or enable overwriting when changing the mask mode,
because both modes use the same filenames.

## What is generated

### EoR-like sky

This component creates a Gaussian white-noise sky at the configured RMS and
random seed. Each frequency channel is shifted to zero spatial mean before the
field-of-view mask is applied.

### GSM diffuse foreground

This component evaluates `pygdsm.GlobalSkyModel` (GSM2008) at each configured
frequency, converts the maps to the requested HEALPix resolution and ICRS
coordinates, and applies the field-of-view mask.

### GLEAM-like and A-team foreground

This component draws a statistical faint-source population from its
source-count model and adds the bright sources listed in
`default_configs/a_team.csv`.

### Visibilities

For every generated component and field of view, `pyuvsim` uses the configured
array, beam, time and frequency sampling to create a UVH5 visibility file.

## Outputs

Products are grouped by component under the configured output directory:

```text
sky_models/{eor,gsm,gleam_ateam}/fov-*.skyh5
observation_files/{eor,gsm,gleam_ateam}/fov-*.yaml
visibilities/{eor,gsm,gleam_ateam}/fov-*.uvh5
```

The observation YAML files are generated inputs for `pyuvsim`.

## Using the library

The main functions can also be called directly with a configuration dictionary:

```python
from bayeseor_data_generation import generate_skies, simulate_visibilities

sky_paths = generate_skies(config)
visibility_paths = simulate_visibilities(config)
```

## Summing visibility components

Visibility files are kept separate during generation. Combine the exact
components needed for an analysis explicitly with `sum_visibility_files`:

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

The inputs must have compatible time, frequency, baseline and polarization
sampling. Pass `overwrite=True` to replace an existing summed output.
