#!/usr/bin/env python3
"""Run the complete example sky and visibility-generation workflow."""

from bayeseor_data_generation import generate_skies, simulate_visibilities

from example_config import CONFIG, OVERWRITE_EXISTING_FILES, SIMULATE_VISIBILITIES


def main() -> None:
    """Generate the configured skies and, when enabled, their visibilities."""
    products = generate_skies(CONFIG, overwrite=OVERWRITE_EXISTING_FILES)
    for component, paths in products.items():
        print(f"Generated {len(paths)} {component} sky model(s).")

    if SIMULATE_VISIBILITIES:
        visibility_paths = simulate_visibilities(
            CONFIG, overwrite=OVERWRITE_EXISTING_FILES
        )
        print(f"Generated {len(visibility_paths)} visibility file(s).")


if __name__ == "__main__":
    main()
