"""
All function and variable names have been taken directly
from nmie-reference/nmie3a.for. I am just as unhappy as
you are about how unmeaningful the names are, and am
equally unhappy about the lack of comments (any clarifying
comments were added as part of porting).
"""
import math

import numba
import numpy

# This module uses the library "numba" to jit-compile
# (jit stands for "just in time") the python code into
# machine code using the LLVM compiler. This makes the
# code ultra-fast. The only price to pay is a small complier
# time on the first use of the function.
#
# The nice thing is that numba figures out the types for us
# so there is no need to specify a function signature (though
# we could if we wanted to).


@numba.jit(nogil=True, cache=True)
def shexqn1(ri_n, aa, x):
    """
    *********   shexqn1 - Spheres: n-layers
                          Theory: exact
                          Results: efficiency factors
    March 1999, AI SPbU
    """
    n_l = len(ri_n)
    xx = numpy.empty(n_l, dtype=numpy.float64)

    ax = 1.0 / x
    xx[0] = x * aa[0]
    xx[-1] = x

    for i in range(1, n_l - 1):
        xx[i] = 0.0
        for j in range(i + 1):
            xx[i] += aa[j]
        xx[i] *= x

    # d1(x), rd3(x), rc(x)
    num = nm(x)
    d1x = aax(ax, num)
    rd3x, rcx = cd3x(x, d1x)

    ari = abs(ri_n[0])
    for i in range(1, n_l):
        ari1 = abs(ri_n[i])
        if ari1 > ari:
            ari = ari1

    num2 = nm(ari * x)

    # rd11(m_1*x_1)
    if (ri_n[0] * xx[0]).imag > 20.0:
        print("k*x > 20", (ri_n[0] * xx[0]).imag)
    rd11 = aa1(ri_n[0] * xx[0], num2)

    # In the original FORTRAN code, the following arrays were
    # over-allocated and then initiallized with zeros, so if num
    # were greater than num2 there was no problem. In the Python
    # implementation, we only allocate what we need, so we need to
    # make sure that the arrays can handle being indexed up to num
    # even though we may only calculate vales up to num2.
    num3 = max(num, num2)
    if num2 < num3:
        rd11 = numpy.append(rd11, numpy.zeros(num3 - num2, dtype=numpy.complex128))
    rrbb = numpy.zeros((n_l, num3), dtype=numpy.complex128)
    rrd1 = numpy.zeros((n_l, num3), dtype=numpy.complex128)
    rrd2 = numpy.zeros((n_l, num3), dtype=numpy.complex128)
    srbb = numpy.zeros((n_l, num3), dtype=numpy.complex128)
    srd1 = numpy.zeros((n_l, num3), dtype=numpy.complex128)
    srd2 = numpy.zeros((n_l, num3), dtype=numpy.complex128)
    for i in range(1, n_l):

        # rd1(m_i*x_i-1), rd2(m_i*x_i-1), rbb(m_i*x_i-1), rcc(m_i*x_i-1),
        if (ri_n[i] * xx[i - 1]).imag > 20.0:
            print("k*x > 20", (ri_n[i] * xx[i - 1]).imag)
        rd1, rd2, rbb = bcd(ri_n[i] * xx[i - 1], num2)
        for j in range(num2):
            rrbb[i, j] = rbb[j]
            rrd1[i, j] = rd1[j]
            rrd2[i, j] = rd2[j]

        # rd1(m_i*x_i), rd2(m_i*x_i), rbb(m_i*x_i), rcc(m_i*x_i),
        if (ri_n[i] * xx[i]).imag > 20.0:
            print("k*x > 20", (ri_n[i] * xx[i]).imag)

        rd1, rd2, rbb = bcd(ri_n[i] * xx[i], num2)
        for j in range(num2):
            srbb[i, j] = rbb[j]
            srd1[i, j] = rd1[j]
            srd2[i, j] = rd2[j]

    num1, ra, rb = abn1(
        ri_n, num, rrbb, rrd1, rrd2, srbb, srd1, srd2, rd11, rd3x, rcx, d1x
    )
    qext, qsca, qbk, qpr = qq1(ax, ra[:num1], rb[:num1])

    qabs = qext - qsca
    alb = qsca / qext
    g = (qext - qpr) / qsca

    return qext, qsca, qabs, qbk, qpr, alb, g


@numba.jit(nogil=True, cache=True)
def nm(x):
    """
    NM-auxiliary function for AA1 & BESSEL
    (number NM is calculated using X)
    see: Trudy Astronom. Observ. LGU V.28,P.14,1971
    for X>1 value of NM was raised
    August 1989, AO LGU
    """
    if x < 1:
        return int(7.5 * x + 9)
    elif x > 100:
        return int(1.0625 * x + 28.5)
    else:
        return int(1.25 * x + 15.5)


@numba.jit(nogil=True, cache=True)
def aa1(rx, num):
    """
    AA1-subroutine for calculations of the ratio of the derivative
       to the function for Bessel functions of half order with
       the complex argument: J'(N)/J(N).
       The calculations are given by the recursive expression
       ``from top to bottom'' beginning from N=NUM.
       RU-array of results.
       A=1/X (X=2*PI*A(particle radius)/LAMBDA - size parameter).
       RI - complex refractive index.
    August 1989, AO LGU
    """
    s = 1.0 / rx
    ru = numpy.empty(num, dtype=numpy.complex128)
    ru[-1] = (num + 1.0) * s
    for j in range(num - 1):
        i = (num - 1) - (j + 1)
        i1 = i + 1
        s1 = (i1 + 1) * s
        ru[i] = s1 - 1.0 / (ru[i1] + s1)
    return ru


@numba.jit(nogil=True, cache=True)
def aax(a, num):
    """
    AAx-subroutine for calculations of the ratio of the derivative
       to the function for Bessel functions of half order with
       the real argument: J'(N)/J(N).
       The calculations are given by the recursive expression
       ``from top to bottom'' beginning from N=NUM.
       RU-array of results.
       A=1/X (X=2*PI*A(particle radius)/LAMBDA - size parameter).
    March 1999, AI SPbU
    """
    ru = numpy.empty(num, dtype=numpy.float64)
    ru[-1] = (num + 1.0) * a
    for j in range(num - 1):
        i = (num - 1) - (j + 1)
        i1 = i + 1
        s1 = (i1 + 1) * a
        ru[i] = s1 - 1.0 / (ru[i1] + s1)
    return ru


@numba.jit(nogil=True, cache=True)
def cd3x(x, d1x):
    """
    CD3X-subroutine for calculations of the ratio of the derivative
       to the function for Riccati-Bessel functions of half order with
       the real argument: zeta'(N)/zeta(N)
       and the ratio of functions: psi(N)/zeta(N).
       The calculations are given by the recursive expression
       ``from bottom to top'' beginning from N=0.
       rd3x, rcx-arrays of results.
       X - size parameter
    March 1999, AI SPbU
    """
    num = len(d1x)
    s1 = 1j
    ax = 1.0 / x

    rd3x = numpy.empty(num, dtype=numpy.complex128)
    rcx = numpy.empty(num, dtype=numpy.complex128)

    rd30 = s1
    rxy = math.cos(2.0 * x) + s1 * math.sin(2.0 * x)
    rc0 = -(1.0 - rxy) / (2.0 * rxy)
    rd3x[0] = -ax + 1.0 / (ax - rd30)
    rcx[0] = rc0 * (ax + rd3x[0]) / (ax + d1x[0])

    for i in range(1, num):
        a1 = (i + 1) * ax
        rd3x[i] = -a1 + 1.0 / (a1 - rd3x[i - 1])
        rcx[i] = rcx[i - 1] * (a1 + rd3x[i]) / (a1 + d1x[i])

    return rd3x, rcx


@numba.jit(nogil=True, cache=True)
def bcd(rx, num):
    """
     BCD-subroutine for calculations of the ratios of the derivative
        to the function for Riccati-Bessel functions of half order with
        the complex argument: psi'(N)/psi(N) and khi'(N)/khi(N)
        and the ratios of functions: psi(N)/khi(N).
        The calculations are given by the recursive expression
        ``from bottom to top'' beginning from N=0.
        rd1, rd2, rbb, rcc-arrays of results.
        rx - (refr. index) * (size parameter)
     March 1999, AI SPbU
    """
    s1 = 1j
    x = rx.real
    y = rx.imag
    rx1 = 1.0 / rx

    rd1 = aa1(rx, num)
    rd2 = numpy.empty(num, dtype=numpy.complex128)
    rbb = numpy.empty(num, dtype=numpy.complex128)
    rd3 = numpy.empty(num, dtype=numpy.complex128)
    rcc = numpy.empty(num, dtype=numpy.complex128)

    # n = 0
    rd30 = s1
    rxy = (math.cos(2.0 * x) + s1 * math.sin(2.0 * x)) * math.exp(-2.0 * y)
    rc0 = -(1.0 - rxy) / (2.0 * rxy)
    rb0 = s1 * (1.0 - rxy) / (1.0 + rxy)

    # n = 1
    rd3[0] = -rx1 + 1.0 / (rx1 - rd30)
    rcc[0] = rc0 * (rx1 + rd3[0]) / (rx1 + rd1[0])
    rd2[0] = (rcc[0] * rd1[0] - rd3[0]) / (rcc[0] - 1.0)
    rbb[0] = rb0 * (rx1 + rd2[0]) / (rx1 + rd1[0])

    for i in range(1, num):
        r1 = (i + 1) * rx1
        rd3[i] = -r1 + 1.0 / (r1 - rd3[i - 1])
        rcc[i] = rcc[i - 1] * (r1 + rd3[i]) / (r1 + rd1[i])
        rd2[i] = (rcc[i] * rd1[i] - rd3[i]) / (rcc[i] - 1.0)
        rbb[i] = rbb[i - 1] * (r1 + rd2[i]) / (r1 + rd1[i])

    return rd1, rd2, rbb


@numba.jit(nogil=True, cache=True)
def qq1(a, ra, rb):
    """
    QQ1-subroutine for calculations of the efficiency factors for
        extinction (QEXT), scattering (QSCA), backscattering (QBK)
        and radiation pressure (QPR) for spherical particles.
    August 1989, AO LGU
    """
    num = len(ra)
    b = 2.0 * a ** 2
    c = 0.0
    d = 0.0
    s = 0j
    r = 0j
    n = 1
    for i in range(num - 1):
        j = i + 1
        n += 2
        r += (j + 0.5) * ((-1) ** j) * (ra[i] - rb[i])
        s += j * (j + 2.0) / (j + 1.0) * (
            ra[i] * ra[i + 1].conjugate() + rb[i] * rb[i + 1].conjugate()
        ) + n // j / (j + 1.0) * (ra[i] * rb[i].conjugate())
        c += (n * (ra[i] + rb[i])).real
        d += (n * (ra[i] * ra[i].conjugate() + rb[i] * rb[i].conjugate())).real

    qext = b * c
    qsca = b * d
    qbk = 2.0 * b * r * r.conjugate()
    qpr = qext - 2.0 * b * s
    return qext, qsca, qbk.real, qpr.real


@numba.jit(nogil=True, cache=True)
def abn1(ri_n, num, rrbb, rrd1, rrd2, srbb, srd1, srd2, rd11, rd3x, rcx, d1x):
    """
    ABn1-subroutine for calculations of the complex coefficients
       A(N), B(N) for n-layered spheres.
       n_l - number of layers
       RI_n(i) - complex refractive indices for innermost layer (1),
       layer2, ... (i = 1, n_l)
       The coefficients are calculated up to the number NUM1.LE.NUM,
       for which |A(N)**2+B(N)**2|.LE.10**(-40)
       RA-array of coefficients A(N), RB-array of coefficients B(N)
    March 1999, AI SPbU
    """

    n_l = len(ri_n)
    sa = numpy.empty(n_l, dtype=numpy.complex128)
    sha = numpy.empty(n_l, dtype=numpy.complex128)
    sb = numpy.empty(n_l, dtype=numpy.complex128)
    shb = numpy.empty(n_l, dtype=numpy.complex128)
    ra = numpy.empty(num, dtype=numpy.complex128)
    rb = numpy.empty(num, dtype=numpy.complex128)

    for i in range(num):

        sa[0] = 0j
        sha[0] = rd11[i]
        sb[0] = 0j
        shb[0] = rd11[i]

        for j in range(1, n_l):
            if abs(ri_n[j] * sha[j - 1] - ri_n[j - 1] * rrd2[j, i]) == 0.0:
                sa[j] = (
                    rrbb[j, i]
                    * (ri_n[j] * sha[j - 1] - ri_n[j - 1] * rrd1[j, i])
                    / (ri_n[j] * sha[j - 1] - ri_n[j - 1] * rrd2[j, i] + 1e-30)
                )
            else:
                sa[j] = (
                    rrbb[j, i]
                    * (ri_n[j] * sha[j - 1] - ri_n[j - 1] * rrd1[j, i])
                    / (ri_n[j] * sha[j - 1] - ri_n[j - 1] * rrd2[j, i])
                )

            if abs(ri_n[j] * shb[j - 1] - ri_n[j - 1] * rrd2[j, i]) == 0.0:
                sb[j] = (
                    rrbb[j, i]
                    * (ri_n[j - 1] * shb[j - 1] - ri_n[j] * rrd1[j, i])
                    / (ri_n[j - 1] * shb[j - 1] - ri_n[j] * rrd2[j, i] + 1e-30)
                )
            else:
                sb[j] = (
                    rrbb[j, i]
                    * (ri_n[j - 1] * shb[j - 1] - ri_n[j] * rrd1[j, i])
                    / (ri_n[j - 1] * shb[j - 1] - ri_n[j] * rrd2[j, i])
                )

            if abs(srbb[j, i] - sa[j]) == 0.0:
                sha[j] = srbb[j, i] * srd1[j, i] / (srbb[j, i] - sa[j] + 1e-30) - sa[
                    j
                ] * srd2[j, i] / (srbb[j, i] - sa[j] + 1e-30)
            else:
                sha[j] = srbb[j, i] * srd1[j, i] / (srbb[j, i] - sa[j]) - sa[j] * srd2[
                    j, i
                ] / (srbb[j, i] - sa[j])

            if abs(srbb[j, i] - sb[j]) == 0.0:
                shb[j] = srbb[j, i] * srd1[j, i] / (srbb[j, i] - sb[j] + 1e-30) - sb[
                    j
                ] * srd2[j, i] / (srbb[j, i] - sb[j] + 1e-30)
            else:
                shb[j] = srbb[j, i] * srd1[j, i] / (srbb[j, i] - sb[j]) - sb[j] * srd2[
                    j, i
                ] / (srbb[j, i] - sb[j])

        # calculations of a(n), b(n)
        ra[i] = rcx[i] * (sha[-1] - ri_n[-1] * d1x[i]) / (sha[-1] - ri_n[-1] * rd3x[i])
        rb[i] = rcx[i] * (ri_n[-1] * shb[-1] - d1x[i]) / (ri_n[-1] * shb[-1] - rd3x[i])

        if abs(ra[i]) + abs(rb[i]) <= 1e-40:
            break

    return i + 1, ra, rb
