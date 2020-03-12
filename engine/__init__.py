"""
Nanoparticle efficiency engine front-end functions

These functions are intended to be user-friendly and
documented.
"""
import math

import numba
import numpy as np

from .mie import shexqn1


def mie_efficiency(diameter, relative_radii, wavelength, refractive_index, medium=1.0):
    """
    Calculates the Mie efficiencies for the given nanosphere at the
    given wavelength in the given surrounding medium.

    Parameters
    ----------
    diameter : float
        The diameter in nm of the nanoparticle to simulate.
    relative_radii : array of float
        The relative radii of the nanoparticle layers, where the first
        element is the innermost layer. All relative radii must sum to
        unity.
    wavelength : float
        The wavelength in nm of light at which to simulate the efficiency.
    refractive_index : array of complex
        The refractive indexes of the nanoparticle layers for the given
        wavelength. The legnth of the array must equal that of
        the relative_radii argument.
    medium : float, optional
        The refracive index of the medium in which to simulate the
        nanoparticle, by default 1.0

    Returns
    -------
    extinction, float
    scattering, float
    absorption, float

    """
    ri = np.asarray(refractive_index)
    rel_rad = np.asarray(relative_radii)
    if len(ri) != len(rel_rad):
        raise ValueError(
            "len(refractive_index) != len(relative_radii) "
            f"({len(ri)} != {len(rel_rad)})"
        )
    _validate_rel_rad(rel_rad)
    return _mie_efficiency(diameter, rel_rad, wavelength, ri, medium)


def mie_efficiencies_by_wavelength(
    diameter, relative_radii, wavelengths, refractive_index, medium=1.0
):
    """
    Calculates the Mie efficiencies for the given nanosphere at the
    given wavelength range in the given surrounding medium.

    Parameters
    ----------
    diameter : float
        The diameter in nm of the nanoparticle to simulate.
    relative_radii : array of float
        The relative radii of the nanoparticle layers, where the first
        element is the innermost layer. All relative radii must sum to
        unity.
    wavelength : array of float
        The wavelengths in nm of light at which to simulate the efficiency.
    refractive_index : 2D array of complex
        The refractive indexes of the nanoparticle layers for each given
        wavelengths. The shape of the array must the length of the
        relative_radii array by the length of the wavelength array.
    medium : float, optional
        The refracive index of the medium in which to simulate the
        nanoparticle, by default 1.0

    Returns
    -------
    extinction, array of float
    scattering, array of float
    absorption, array of float

    """
    ri = np.asarray(refractive_index)
    rel_rad = np.asarray(relative_radii)
    waves = np.asarray(wavelengths)
    if ri.shape != (len(rel_rad), len(waves)):
        raise ValueError(
            "refractive_index.shape != (len(relative_radii), len(wavelengths)) "
            f"({ri.shape} != ({len(rel_rad)}, {len(waves)}))"
        )
    _validate_rel_rad(rel_rad)
    return _mie_efficiencies_by_wavelength(diameter, rel_rad, waves, ri, medium)


def _validate_rel_rad(relative_radii):
    """Ensure the sum of the relative radii is unity"""
    rel_rad_sum = relative_radii.sum()
    close_to_unity = abs(rel_rad_sum - 1.0) <= 1e-9 * max(abs(rel_rad_sum), 1.0)
    if not close_to_unity:
        raise ValueError(
            f"The sum of relative_radii must sum to unity (got {rel_rad_sum})."
        )


@numba.jit(nogil=True, cache=True)
def _mie_efficiency(diameter, rel_rad, wavelength, ri, medium):
    """
    Work that actually calculates the mie efficiencies.

    Input is assumed to be correct.

    Parameters
    ----------
    diameter : float
    relative_radii : array of float
    wavelength : float
    refractive_index : array of complex
    medium : float, optional

    Returns
    -------
    extinction, float
    scattering, float
    absorption, float

    """
    # The size parameter encodes the nanoparticle size as well as the
    # wavelength of light
    size_parameter = math.pi * diameter * medium / wavelength
    return shexqn1(ri, rel_rad, size_parameter)[:3]


@numba.jit(nogil=True, parallel=True, cache=True)
def _mie_efficiencies_by_wavelength(
    diameter, relative_radii, wavelengths, refractive_index, medium
):
    """
    Work that actually calculates the mie efficiencies across wavelengths.

    Input is assumed to be correct.

    Parameters
    ----------
    diameter : float
    relative_radii : array of float
    wavelength : float
    refractive_index : 2D array of complex
    medium : float, optional

    Returns
    -------
    extinction, array of float
    scattering, array of float
    absorption, array of float

    """
    qext = np.empty(len(wavelengths), dtype=np.float64)
    qabs = np.empty(len(wavelengths), dtype=np.float64)
    qsca = np.empty(len(wavelengths), dtype=np.float64)
    for i in numba.prange(len(wavelengths)):
        (
            qext[i],
            qabs[i],
            qsca[i],
        ) = _mie_efficiency(
             diameter,
             relative_radii,
             wavelengths[i],
             refractive_index[:, i],
             medium,
        )
    return qext, qabs, qsca