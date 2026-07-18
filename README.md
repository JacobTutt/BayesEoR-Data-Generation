# BayesEoR Data Generation

A standalone, reusable pipeline for generating the three sky components used in
BayesEoR-style simulations and propagating them through a pyuvsim visibility
calculation. It does not import or modify the BayesEoR inference code and does
not require pre-generated sky or visibility products.

## Repository layout

```text
src/bayeseor_data_generation/
  eor.py                 EoR-like white-noise sky generation
  gsm.py                 diffuse GSM2008 sky generation
  gleam_ateam.py         statistical GLEAM-like plus A-team generation
  generate_skies.py      shared pointing and sky-generation workflow
  simulate.py            pyuvsim observation and visibility workflow
default_configs/
  hex-37-14.6m.csv       example antenna layout
  a_team.csv             fixed bright-source catalogue
example_config.yaml      complete editable example configuration
example_run.py           complete example workflow
```

Generated data, observation files, visibilities, machine environments, batch
scripts and scheduler logs are deliberately not part of the repository.

## Installation

Clone the repository, create an isolated Python environment and install the
package in editable mode:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e .
```

The dependency constraints in `pyproject.toml` retain the historical pyuvsim,
pyuvdata and pyradiosky APIs used while developing this pipeline. Visibility
simulation uses MPI through `mpi4py` and pyuvsim; an MPI implementation must
therefore be available on the machine.

## Running the example

Review `example_config.yaml`, especially the output directory, time and frequency
axes, component FoVs and `run.simulate_visibilities` switch. Then run:

```bash
python example_run.py
```

The Python runner loads the YAML and first calls `generate_skies(config)`. If
`run.simulate_visibilities` is true, it then passes the same dictionary to
`simulate_visibilities(config)`. This guarantees that the sky selection and
visibility calculation use one observation definition.

The supplied configuration is the recovered HERA Band-2 example. It is a real,
large simulation: 180 frequency channels, 243 time integrations, a 37-element
layout, and foreground FoVs out to 120 degrees. For a quick trial, reduce the
channel count, time range and FoV lists before running it.

Existing outputs are preserved. Set `run.overwrite_existing_files: true` in
`example_config.yaml` only when existing products should be replaced.

## Configuration

The YAML has a small `run` section followed by the physical configuration:

- `run`: whether to simulate visibilities and overwrite existing products.
- `output_directory`: destination for all generated products.
- `telescope`: antenna layout, geodetic location, beam and redundancy tolerance.
- `time`: start/end Julian dates and integration time in seconds.
- `frequency`: first channel, channel width and number of channels.
- `sky`: common HEALPix NSIDE plus settings for each physical sky component.

FoVs are full angular diameters centred on zenith at the midpoint of the
observation. For example, a 30-degree FoV retains pixels or point sources within
15 degrees of zenith.

The runner interprets relative paths from the directory containing
`example_config.yaml`, so the example works from any current working directory.
A downstream project can load its own YAML or construct a dictionary and call
the two library functions directly:

```python
from bayeseor_data_generation import generate_skies, simulate_visibilities

sky_paths = generate_skies(my_config)
visibility_paths = simulate_visibilities(my_config)
```

## Required physical inputs

### Telescope and observation

The antenna-layout CSV supplies antenna numbers, beam IDs and East/North/Up
positions in metres. The telescope configuration also defines latitude,
longitude, altitude, pyuvsim beam type and aperture diameter. The time and
frequency sections define the complete visibility sampling axes.

The included layout is a 37-element HERA-like hexagon with 14.6 m antenna
separation. This separation is distinct from the example's 14.2 m Airy-aperture
diameter.

### EoR-like white noise

The EoR component is a Gaussian statistical sky with configurable NSIDE, RMS,
random seed and FoV. The full-sky cube is drawn before the FoV is selected, and
each frequency channel is shifted to zero spatial mean. The supplied values are
6.4 mK RMS and random seed 92381923.

### GSM diffuse foreground

`pygdsm.GlobalSkyModel` evaluates GSM2008 at every requested frequency. It
provides the diffuse foreground's spatially varying spectral structure; the
frequency maps are not generated from one global spectral index. The maps are
converted to the configured NSIDE, transformed from Galactic coordinates to
ICRS and selected around the central-time zenith. Frequency chunking reduces
memory use without changing the physical model.

### GLEAM-like plus A-team foreground

The faint-source population is generated statistically from the Franzen et al.
source-count model using configurable random seed, confusion NSIDE and reference
frequency. It is not read from a pre-generated GLEAM catalogue. Ten fixed bright
A-team sources are read from `default_configs/a_team.csv`, where positions,
fluxes, reference frequencies and spectral indices remain visible and editable.

## Outputs

The example creates:

```text
data/h1c_band2/
  sky_models/{eor,gsm,gleam_ateam}/fov-*.skyh5
  observation_files/{eor,gsm,gleam_ateam}/fov-*.yaml
  visibilities/{eor,gsm,gleam_ateam}/fov-*.uvh5
```

The observation YAML files are generated runtime inputs for pyuvsim. They are
derived from the loaded configuration; users do not maintain a second
observation configuration.

This repository stops at component visibility generation. Combining components
into a mock observation and selecting a BayesEoR inference data vector belong in
the downstream analysis.
