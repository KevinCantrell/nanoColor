"""
All function and variable names have been taken directly
from nmie-reference/nmie3a.for. I am just as unhappy as
you are about how unmeaningful the names are.
"""


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
