"""Combine explicitly selected visibility simulations."""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path


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
        total = total.sum_vis(component, inplace=False)

    total.history += " Input visibility files: " + ", ".join(
        str(path) for path in inputs
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    total.write_uvh5(str(output), clobber=overwrite)
    return output
