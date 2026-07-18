"""Minimal version shim for the pinned historical preprocessing script.

The official script only imports ``bayeseor`` to record this version string in
its output dictionary. Keeping the shim local avoids importing the inference
package and its unrelated MultiNest runtime during data preprocessing.
"""

__version__ = "1.0.1.dev75+g86d5f92"
