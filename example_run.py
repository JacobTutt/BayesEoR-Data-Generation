#!/usr/bin/env python3
"""Run the complete example sky and visibility-generation workflow."""

from pathlib import Path

import yaml

from bayeseor_data_generation import generate_skies, simulate_visibilities


def main() -> None:
    """Load the YAML and generate its skies and optional visibilities."""
    config_file = Path(__file__).resolve().with_name("example_config.yaml")
    with config_file.open() as stream:
        config = yaml.safe_load(stream)

    run_options = config.pop("run")
    config_directory = config_file.parent

    output_directory = Path(config["output_directory"]).expanduser()
    if not output_directory.is_absolute():
        output_directory = config_directory / output_directory
    config["output_directory"] = output_directory.resolve()

    antenna_layout = Path(config["telescope"]["antenna_layout"]).expanduser()
    if not antenna_layout.is_absolute():
        antenna_layout = config_directory / antenna_layout
    config["telescope"]["antenna_layout"] = antenna_layout.resolve()

    a_team_catalogue = Path(
        config["sky"]["gleam_ateam"]["a_team_catalogue"]
    ).expanduser()
    if not a_team_catalogue.is_absolute():
        a_team_catalogue = config_directory / a_team_catalogue
    config["sky"]["gleam_ateam"]["a_team_catalogue"] = (
        a_team_catalogue.resolve()
    )

    overwrite = run_options["overwrite_existing_files"]
    products = generate_skies(config, overwrite=overwrite)
    for component, paths in products.items():
        print(f"Generated {len(paths)} {component} sky model(s).")

    if run_options["simulate_visibilities"]:
        visibility_paths = simulate_visibilities(config, overwrite=overwrite)
        print(f"Generated {len(visibility_paths)} visibility file(s).")


if __name__ == "__main__":
    main()
