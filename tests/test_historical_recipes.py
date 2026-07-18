"""Fast checks for the documented example observation."""

import csv
from pathlib import Path

import numpy as np
import yaml

import eor
import generate_skies
import gleam_ateam
import gsm
import simulate


INPUT_FILE = Path(__file__).parents[1] / "inputs/h1c_band2.yaml"


def test_example_frequency_axis():
    """The example must describe the historical 180-channel band-2 axis."""
    inputs = yaml.safe_load(INPUT_FILE.read_text())
    frequency = inputs["frequency"]
    values = frequency["start_hz"] + np.arange(
        frequency["n_channels"]
    ) * frequency["channel_width_hz"]
    assert values.size == 180
    assert values[0] == 150_292_968.75
    assert values[-1] == 167_773_437.5


def test_example_observation_has_all_physical_inputs():
    """The user-facing YAML must contain every required observation section."""
    inputs = yaml.safe_load(INPUT_FILE.read_text())
    assert set(inputs) == {
        "output_directory", "telescope", "time", "frequency", "sky"
    }
    assert set(inputs["sky"]) == {"nside", "eor", "gsm", "gleam_ateam"}
    assert inputs["sky"]["eor"]["fovs_deg"] == [12.9080728652]
    assert inputs["sky"]["gsm"]["fovs_deg"][-1] == 120.0
    assert inputs["sky"]["gleam_ateam"]["fovs_deg"][-1] == 120.0


def test_generation_has_no_reference_data_input():
    """The normal workflow must not know where copied Burba files live."""
    input_text = INPUT_FILE.read_text()
    assert "BayesEoR/large-fov" not in input_text
    assert "/projects/u6" not in input_text


def test_source_count_model_is_finite():
    """The embedded GLEAM-like source-count polynomial must remain numerical."""
    values = gleam_ateam._franzen_differential_source_count(
        np.array([1e-3, 1.0, 10.0])
    )
    assert np.all(np.isfinite(values))
    assert np.all(values > 0)


def test_a_team_is_a_documented_input_catalogue():
    """The ten fixed bright sources must live in CSV rather than Python code."""
    inputs = yaml.safe_load(INPUT_FILE.read_text())
    catalogue = INPUT_FILE.parents[1] / inputs["sky"]["gleam_ateam"][
        "a_team_catalogue"
    ]
    with catalogue.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == 10
    assert rows[0]["name"] == "3C444"
    assert rows[-1]["name"] == "FornaxA"
    assert float(rows[-1]["reference_frequency_hz"]) == 154_000_000


def test_public_workflow_functions_are_documented():
    """Every independently runnable physical stage must explain its inputs."""
    functions = [
        eor.generate_eor_maps,
        gsm.generate_gsm_maps,
        gleam_ateam.generate_gleam_ateam_maps,
        generate_skies.generate_skies,
        simulate.simulate_visibilities,
    ]
    assert all(function.__doc__ and "Parameters" in function.__doc__ for function in functions)
