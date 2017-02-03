"""
Based on Blanton et al. 2007: 2007AJ....133..734B
The UBRVI filters are those of Bessell (1990). The ugriz filters are those
determined by M. Doi, D. Eisenstein, and J. Gunn and are available on the
SDSS DR4 Web site ( http://www.sdss.org/dr4/ ). The JHKs filters are those from
Cohen et al. (2003).
"""


## offset is calculated as m_AB - m_vega
offset = {
    "U" : 0.79,
    "B" : 0.09,
    "V" : 0.02,
    "R" : 0.21,
    "I" : 0.45,
    "u" : 0.91,
    "g" : 0.08,
    "r" : 0.16,
    "i" : 0.37,
    "z" : 0.54,
    "J" : 0.91,
    "H" : 1.39,
    "Ks" : 1.85 ,
    "0:1u" : 1.25 ,
    "0:1g" : 0.01 ,
    "0:1r" : 0.04 ,
    "0:1i" : 0.27 ,
    "0:1z" : 0.46,
}

def convert_Vega_to_AB(phot_table):
    """
    Parameters
    ----------
    Returns
    -------
    """

    return phot_table


def convert_AB_to_Vega():
    """
    Parameters
    ----------
    Returns
    -------
    """

    return phot_table
