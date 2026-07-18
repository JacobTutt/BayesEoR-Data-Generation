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

`example_config.yaml` contains the complete input configuration:

- `run` controls whether visibilities are simulated and whether existing files
  are replaced.
- `output_directory` sets where generated products are written.
- `telescope` defines the antenna layout, location, beam and redundancy
  threshold.
- `time` defines the observation start, end and integration time.
- `frequency` defines the first channel, channel width and number of channels.
- `sky` defines the HEALPix resolution, sky-mask mode and settings for each sky
  component.

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

The three component sections provide their own useful controls:

- `sky.eor` sets the fields of view, RMS and random seed for the EoR-like sky.
- `sky.gsm` sets the fields of view, whether to include the CMB and the
  frequency chunk size used while evaluating GSM2008.
- `sky.gleam_ateam` sets the fields of view, A-team catalogue, random seed,
  confusion resolution and reference frequency for the point-source sky.

Relative paths in the example are resolved from the directory containing
`example_config.yaml`.

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
