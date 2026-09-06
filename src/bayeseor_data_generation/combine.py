"""Combine explicitly selected visibility simulations."""

from __future__ import annotations

from collections.abc import Sequence
from copy import deepcopy
from pathlib import Path


def _phase_centres_match_except_name(
    reference_catalog: dict, component_catalog: dict
) -> bool:
    """Return whether two phase-centre catalogues differ only in ``cat_name``."""
    import numpy as np

    if reference_catalog.keys() != component_catalog.keys():
        return False
    for identifier in reference_catalog:
        reference = reference_catalog[identifier]
        component = component_catalog[identifier]
        keys = set(reference) | set(component)
        for key in keys - {"cat_name"}:
            reference_value = reference.get(key)
            component_value = component.get(key)
            if reference_value is None or component_value is None:
                if reference_value is not component_value:
                    return False
                continue
            reference_array = np.asarray(reference_value)
            component_array = np.asarray(component_value)
            numeric = np.issubdtype(
                reference_array.dtype, np.number
            ) and np.issubdtype(component_array.dtype, np.number)
            if not np.array_equal(
                reference_array,
                component_array,
                equal_nan=numeric,
            ):
                return False
    return True


def sum_visibility_files(
    input_files: Sequence[str | Path],
    output_file: str | Path,
    overwrite: bool = False,
) -> Path:
    """Sum compatible UVH5 visibility products into one UVH5 file.

    Parameters
    ----------
    input_files
        Two or more UVH5 files to sum. Their visibility sampling and metadata
        must be compatible according to :meth:`pyuvdata.UVData.sum_vis`.
    output_file
        Destination for the summed UVH5 product. It must not be one of the
        input files.
    overwrite
        Replace an existing output file when true. The default is false.

    Returns
    -------
    pathlib.Path
        Absolute path of the summed UVH5 file.

    Raises
    ------
    TypeError
        If ``input_files`` is passed as one path rather than a sequence.
    ValueError
        If fewer than two input files are supplied or the output is also an
        input. PyUVData also raises this when visibility metadata are
        incompatible.
    FileNotFoundError
        If an input file does not exist.
    FileExistsError
        If the output already exists and ``overwrite`` is false.

    Notes
    -----
    Pyuvsim writes visibility products as :class:`pyuvdata.UVData` objects.
    PyUVData's ``sum_vis`` operation performs the compatibility checks and
    elementwise complex-visibility addition used here.
    """
    from pyuvdata import UVData

    if isinstance(input_files, (str, Path)):
        raise TypeError("input_files must be a sequence of at least two paths")

    inputs = [Path(path).expanduser().resolve() for path in input_files]
    if len(inputs) < 2:
        raise ValueError("At least two visibility files are required")
    for path in inputs:
        if not path.is_file():
            raise FileNotFoundError(f"Missing visibility file: {path}")

    output = Path(output_file).expanduser().resolve()
    if output in inputs:
        raise ValueError("output_file must not overwrite an input visibility file")
    if output.exists() and not overwrite:
        raise FileExistsError(f"Output already exists: {output}")

    total = UVData()
    total.read(str(inputs[0]))
    for path in inputs[1:]:
        component = UVData()
        component.read(str(path))

        # Pyuvsim derives these bookkeeping labels from each component's input
        # filename. They therefore differ even when the time/frequency/baseline
        # sampling and physical phase centre are identical. Confirm that the
        # phase centres differ only by name before using the reference labels.
        if not _phase_centres_match_except_name(
            total.phase_center_catalog, component.phase_center_catalog
        ):
            raise ValueError(
                f"Physical phase centres do not match for {inputs[0]} and {path}"
            )
        component.phase_center_catalog = deepcopy(total.phase_center_catalog)
        component.filename = deepcopy(total.filename)
        total = total.sum_vis(component, inplace=False)

    total.history += " Input visibility files: " + ", ".join(
        str(path) for path in inputs
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    total.write_uvh5(str(output), clobber=overwrite)
    return output
