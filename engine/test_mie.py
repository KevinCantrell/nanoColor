"""
The tests in this file all directly compare the results of
the FORTRAN-implemented functions to the Python implemented
functions and ensure they are identical.
"""

# The hypothesis library generates "random" input for our functions
# to ensure that the tests are well covered.
from hypothesis import given
from hypothesis.strategies import floats

import mie  # python
import nmie  # fortran


@given(floats(allow_nan=False, allow_infinity=False, min_value=0, max_value=1E6))
def test_nm(x):
    assert mie.nm(x) == nmie.nm(x)
