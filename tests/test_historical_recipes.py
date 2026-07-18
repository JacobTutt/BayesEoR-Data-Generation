from pathlib import Path

import numpy as np

import generate


def test_historical_fov_pixel_counts():
    expected = {
        generate.HISTORICAL_EOR_FOV_DEG: 623,
        30.0: 3352,
        60.0: 13179,
        90.0: 28791,
        120.0: 49153,
    }
    for fov, count in expected.items():
        assert generate.healpix_fov_indices(generate.NSIDE, fov).size == count


def test_frequency_axis():
    frequencies = generate.frequency_axis_hz()
    assert frequencies.size == 180
    assert frequencies[0] == 150_292_968.75
    assert frequencies[-1] == 167_773_437.5


def test_no_reference_path_in_generation_defaults():
    assert "BayesEoR/large-fov" not in str(generate.default_output_root())
