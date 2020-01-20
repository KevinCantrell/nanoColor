"""
The tests in this file all directly compare the results of
the FORTRAN-implemented functions to the Python implemented
functions and ensure they are identical.
"""
import math

# The hypothesis library generates "random" input for our functions
# to ensure that the tests are well covered.
from hypothesis import assume, given, settings
from hypothesis.strategies import complex_numbers, floats, integers, lists
import numpy

import mie  # python
import nmie  # fortran


# Make hypothesis not worry about the first iteration of functions
# running long due to the jit compiler.
settings.register_profile("jit", deadline=1E5)
settings.load_profile("jit")

# We will be repeating optins for many of the hypothesis strategies,
# so they are defined once here.
rational = dict(allow_nan=False, allow_infinity=False)
int_opts = dict(min_value=1, max_value=1000)
list_opts = dict(min_size=1, max_size=1000)
float_opts = dict(min_value=-1e20, max_value=1e20, **rational,)
complex_opts = dict(min_magnitude=0, max_magnitude=1e20, **rational,)


@given(floats(min_value=0, max_value=1e6, **rational))
def test_nm(x):
    assert mie.nm(x) == nmie.nm(x)


@given(
    complex_numbers(**complex_opts).filter(lambda x: x != 0), integers(**int_opts),
)
def test_aa1(rx, num):
    a = mie.aa1(rx, num)
    b = nmie.aa1(rx, num)
    assert numpy.all(numpy.equal(a.real, b.real))
    assert numpy.all(numpy.equal(a.imag, b.imag))


@given(
    floats(**float_opts).filter(lambda x: x != 0), integers(**int_opts),
)
def test_aax(a, num):
    assert numpy.all(numpy.equal(mie.aax(a, num), nmie.aax(a, num)))


@given(
    floats(**float_opts).filter(lambda x: x != 0),
    lists(floats(**float_opts).filter(lambda x: x != 0), **list_opts,),
)
def test_cd3x(x, d1x):
    d1x = numpy.asarray(d1x, dtype=numpy.float64)
    a, b = mie.cd3x(x, d1x)
    c, d = nmie.cd3x(x, d1x)
    assert numpy.all(numpy.equal(a.real, c.real))
    assert numpy.all(numpy.equal(a.imag, c.imag))
    assert numpy.all(numpy.equal(b.real, d.real))
    assert numpy.all(numpy.equal(b.imag, d.imag))


@given(
    complex_numbers(min_magnitude=0, max_magnitude=500, **rational)
    .filter(lambda x: x != 0)
    .filter(lambda x: abs(x.imag) < 300),
    integers(**int_opts),
)
def test_bcd(rx, num):
    a, b, c = mie.bcd(rx, num)
    d, e, f = nmie.bcd(rx, num)
    assert numpy.all(numpy.equal(a.real, d.real))
    assert numpy.all(numpy.equal(a.imag, d.imag))
    assert numpy.all(numpy.equal(b.real, e.real))
    assert numpy.all(numpy.equal(b.imag, e.imag))
    assert numpy.all(numpy.equal(c.real, f.real))
    assert numpy.all(numpy.equal(c.imag, f.imag))


@given(
    floats(**float_opts).filter(lambda x: x != 0),
    lists(complex_numbers(**complex_opts).filter(lambda x: x != 0), **list_opts,),
    lists(complex_numbers(**complex_opts).filter(lambda x: x != 0), **list_opts,),
)
def test_qq1(a, ra, rb):
    # Ensure complex arrays are the same length
    n = min(len(ra), len(rb))
    ra = numpy.asarray(ra[:n], dtype=numpy.complex128)
    rb = numpy.asarray(rb[:n], dtype=numpy.complex128)
    assert all(
        math.isclose(x, y, rel_tol=1e-11)
        for x, y in zip(mie.qq1(a, ra, rb), nmie.qq1(a, ra, rb))
    )


def create_rel_rad(n):
    """
    The relative radii (aa) must sum to 1 and each be > 0.
    """
    rel_rad = numpy.zeros(n, dtype=numpy.float64)
    remainder = 1
    for i in range(len(rel_rad) - 1):
        value = numpy.random.random_sample()
        while value == 0.0:  # must be non-zero
            value = numpy.random.random_sample()
        value *= remainder
        rel_rad[i] = value
        remainder -= value
    rel_rad[-1] = remainder
    return rel_rad


# Note that the function abn1 is not explictly tested because it is dependent on
# values generated in the execution of shexqn1. The following tests implictly
# validate abn1.


@given(
    lists(
        complex_numbers(
            allow_nan=False, allow_infinity=False, min_magnitude=0, max_magnitude=10
        ).filter(lambda x: x != 0),
        **list_opts,
    ),
    floats(allow_nan=False, allow_infinity=False, min_value=1e-6, max_value=1),
)
def test_shexqn1(ri_n, x):
    ri_n = numpy.asarray(ri_n)
    assume(max(numpy.abs(ri_n)) * x < 1000)
    # The relative radii (aa) must sum to 1 and be > 0. Easier
    # to construct this explicitly than use hypothesis.
    aa = create_rel_rad(len(ri_n))
    assert all(
        math.isclose(a, b, rel_tol=1e-11)
        for a, b in zip(mie.shexqn1(ri_n, aa, x), nmie.shexqn1(ri_n, aa, x))
    )
