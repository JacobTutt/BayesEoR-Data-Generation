# BayesEoR Data Generation

A small, standalone guide to generating BayesEoR simulation data from scratch.
It starts from a physical description of an observation and produces skyh5 sky
models plus pyuvsim UVH5 visibilities. It does not import or modify the BayesEoR
inference code, and copied historical simulations are not inputs.

## The workflow

One YAML file describes the observation. The code is split only at genuine
physical boundaries:

- `eor.py` generates EoR maps.
- `gsm.py` generates diffuse GSM maps.
- `gleam_ateam.py` generates the statistical point-source+A-team maps.
- `generate_skies.py` computes the central pointing once and calls those three
  component generators.
- `simulate.py` takes the generated maps and calculates their visibilities.

The workflow is:

```text
observation YAML
    |
    +-- telescope location, layout and beam
    +-- start/end time and integration length
    +-- frequency start, spacing and channel count
    +-- HEALPix NSIDE
    +-- EoR FoV(s)
    +-- GSM FoV(s)
    +-- GLEAM-like+A-team FoV(s)
    |
    v
central observation time and zenith direction
    |
    v
EoR, GSM and GLEAM-like+A-team skyh5 files
    |
    v
simulate.py
    |
    v
one pyuvsim observation YAML per sky file
    |
    v
one UVH5 visibility file per component and FoV
```

FoV values are full angular diameters. For example, a 30-degree FoV contains
HEALPix pixel centres or point sources within 15 degrees of zenith. Zenith is
calculated at the midpoint between `start_jd` and `end_jd`.

## Required inputs

All required inputs are visible in `inputs/h1c_band2.yaml`.

### Telescope

| Input | Meaning |
| --- | --- |
| `antenna_layout` | CSV containing antenna numbers, beam IDs and ENU positions in metres |
| `latitude_deg`, `longitude_deg`, `altitude_m` | Telescope location used both for the sky cut and pyuvsim |
| `beam.type` | pyuvsim beam model; the historical example uses `airy` |
| `beam.diameter_m` | Airy aperture diameter; the recovered historical value is 14.2 m |
| `redundant_threshold_m` | pyuvsim tolerance for grouping redundant baselines |

The included `instrument/hex-37-14.6m.csv` is a 37-element hexagonal HERA-like
layout with 14.6 m antenna separation. The separation and Airy aperture diameter
are different quantities.

### Time axis

| Input | Historical example |
| --- | ---: |
| Start time | JD 2459999.067040411 |
| End time | JD 2459999.127292119 |
| Integration time | 21.511353650861533 s |
| Result | 243 integrations |

The central JD is used to decide which sky pixels and point sources fall inside
each FoV. The full start/end range is passed to pyuvsim for the visibility
simulation.

### Frequency axis

| Input | Historical example |
| --- | ---: |
| First channel | 150.29296875 MHz |
| Channel spacing | 97.65625 kHz |
| Number of channels | 180 |
| Final channel | 167.7734375 MHz |

The same frequency array is used by the EoR and GSM sky cubes and by pyuvsim.
Point-source fluxes are represented by reference fluxes and spectral indices;
pyuvsim evaluates them at these frequencies.

### Sky fields

Each sky component has its own `fovs_deg` list. The supplied example requests:

- EoR: 12.9080728652 degrees.
- GSM: 12.9080728652, 30, 60, 90 and 120 degrees.
- GLEAM-like+A-team: 12.9080728652, 30, 60, 90 and 120 degrees.

Changing or adding an FoV requires editing only this list. There are no FoVs
hard-coded into the workflow.

## How each sky is generated

### EoR-like white noise

The historical example is a Gaussian, frequency-independent statistical sky:

- NSIDE 128.
- RMS 6.4 mK.
- NumPy random seed 92381923.
- Every frequency channel is shifted to zero spatial mean.

The random cube is drawn over the full sky in memory before the requested FoV
is selected. This preserves the historical random-number ordering.

### GSM

`pygdsm.GlobalSkyModel` generates GSM2008 at every requested frequency. This is
where the diffuse foreground's spatially varying spectrum comes from: the
frequency maps are not produced using one global spectral index. The maps are
downgraded to NSIDE 128, transformed from Galactic coordinates to ICRS, and
then cut around zenith at the central observation time.

The GSM is evaluated in small frequency chunks to reduce laptop memory use.
This changes memory consumption, not the physical model.

### GLEAM-like plus A-team

The point-source model is generated rather than read from a catalogue file:

- Franzen et al. statistical source counts.
- Random seed 42.
- NSIDE-256 confusion threshold.
- Spectral indices drawn around -0.8.
- Ten explicit A-team sources with their positions, fluxes and spectral indices.

The full statistical population is generated in memory. Sources are then
selected by their zenith angle at the same central time used for the HEALPix
maps.

## Running it

Generate the three sets of skyh5 files:

```bash
python generate_skies.py inputs/h1c_band2.yaml
```

Then take those sky files and generate their pyuvsim YAML and UVH5 files:

```bash
python simulate.py inputs/h1c_band2.yaml
```

Existing files are retained. To regenerate them:

```bash
python generate_skies.py inputs/h1c_band2.yaml --overwrite
python simulate.py inputs/h1c_band2.yaml --overwrite
```

## Output layout

The example writes:

```text
data/h1c_band2/
  sky_models/
    eor/
    gsm/
    gleam_ateam/
  observation_files/
    telescope.yaml
    eor/
    gsm/
    gleam_ateam/
  visibilities/
    eor/
    gsm/
    gleam_ateam/
```

For every `sky_models/<component>/fov-X.skyh5`, the workflow creates a matching
`observation_files/<component>/fov-X.yaml` and, with `--simulate`, a matching
`visibilities/<component>/fov-X.uvh5`.

## Isambard

Install the historical environment once, then run the same command through
Slurm:

```bash
source activate-isambard.sh
sbatch jobs/generate_skies.sbatch inputs/h1c_band2.yaml
```

After inspecting the generated sky and observation files:

```bash
sbatch jobs/simulate.sbatch inputs/h1c_band2.yaml
```

Both jobs request one GH200 superchip division: 72 CPU cores, 120 GB CPU memory
and one GPU. Sky generation uses one process; visibility simulation uses four
MPI tasks with 18 cores each. The GPU is part of the allocation unit; this
historical pyuvsim workflow is primarily CPU/MPI work.

## Scope

This repository stops at component visibility generation. Combining components
into a particular mock dataset and selecting the final BayesEoR inference vector
are analysis choices and should live in the downstream analysis repository.
