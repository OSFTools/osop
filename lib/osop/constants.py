"""
(C) Crown Copyright, Met Office. All rights reserved.

This file is part of osop and is released under the BSD 3-Clause license.
See LICENSE in the root of the repository for full licensing details.
"""

"""
Collection of constants used in the OSOP library
"""

# mapping of centre to system based on latest
# systems March 2023. More flexibility may be
# needed on this longer term
SYSTEMS = {
    "ecmwf": "51",
    "meteo_france": "8",
    "dwd": "21",
    "cmcc": "35",
    "ncep": "2",
    "jma": "3",
    "eccc_can": "2",
    "eccc_gem5": "3",
    "ukmo": "602",
}

SYSTEMSFC = {
    "ecmwf": "51",
    "meteo_france": "9",
    "dwd": "22",
    "cmcc": "35",
    "ncep": "2",
    "jma": "3",
    "eccc_can": "4",
    "eccc_gem5": "5",
    "ukmo": "604",
}

# named domains corners - x = longitude, y = latitude
DOMAINS = {
    "ArabCOF": {"x0": -30.0, "x1": 60.0, "y0": -2.5, "y1": 45.0},
    "MedCOF": {"x0": -18.5, "x1": 52.0, "y0": 15.0, "y1": 25.0},
}
