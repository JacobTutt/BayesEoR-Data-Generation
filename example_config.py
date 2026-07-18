"""Example observation and sky configuration.

Edit the values in ``CONFIG`` to define a new simulation. Paths are constructed
relative to this file, so the example works regardless of the current working
directory.
"""

from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parent

# The runner imports these two controls. Keeping them outside CONFIG separates
# execution choices from the physical definition of the observation.
OVERWRITE_EXISTING_FILES = False
SIMULATE_VISIBILITIES = True

CONFIG = {
    "output_directory": REPOSITORY_ROOT / "data" / "h1c_band2",
    "telescope": {
        "name": "Hex37-14.6m",
        "antenna_layout": (
            REPOSITORY_ROOT / "default_configs" / "hex-37-14.6m.csv"
        ),
        "latitude_deg": -30.72152777777791,
        "longitude_deg": 21.428305555555557,
        "altitude_m": 1073.0000000093132,
        "beam": {"type": "airy", "diameter_m": 14.2},
        "redundant_threshold_m": 1.0,
    },
    "time": {
        "start_jd": 2459999.067040411,
        "end_jd": 2459999.127292119,
        "integration_time_s": 21.511353650861533,
    },
    "frequency": {
        "start_hz": 150292968.75,
        "channel_width_hz": 97656.25,
        "n_channels": 180,
    },
    "sky": {
        "nside": 128,
        "eor": {
            "fovs_deg": [12.9080728652],
            "rms_k": 0.0064,
            "random_seed": 92381923,
        },
        "gsm": {
            "fovs_deg": [12.9080728652, 30.0, 60.0, 90.0, 120.0],
            "include_cmb": True,
            "frequency_chunk_size": 10,
        },
        "gleam_ateam": {
            "fovs_deg": [12.9080728652, 30.0, 60.0, 90.0, 120.0],
            "a_team_catalogue": (
                REPOSITORY_ROOT / "default_configs" / "a_team.csv"
            ),
            "random_seed": 42,
            "confusion_nside": 256,
            "reference_frequency_hz": 154000000.0,
        },
    },
}
