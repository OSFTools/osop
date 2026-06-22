# Introduction

```{eval-rst}
.. toctree::
   :maxdepth: 2
   :caption: Contents:
```

This document provides guidance on using the Objective Seasonal
Outlook Package (OSOP) toolkit for producing seasonal based
forecasting and verification.

This user guide includes three main sections:

* [Initial set_up](set_up)
* [Running the toolkit](run)
* [Developing the toolkit](devel)

This user guide contains instructions on how to get started with
the package using the currently available functions. The code is
able to generate:

* verification plots of GCM skill against observations
* direct tercile forecasts
* CCA based calibrated forecasts

The code can currently download ERa5 data and seasonal forecast models
available from [Copernicus Climate Data Store](https://cds.climate.copernicus.eu/).

Beyond this, the user can modify or change the OSOP toolkit
as needed to suit their requirements.

The package is designed so that there are a small number of top level
bash scripts which the user modifies to set options and the locations
to store data and images. These scripts call more moduler python scripts
which do different stages of the workflow - for example
download of data or creation of skill scores.

The goal is that the user does not need to be experienced in coding to
use the OSOP toolkit at a basic level. Users more advanced in coding
may find that they can edit further down in the software to get more
highly specific use cases. As outlined in the section
[development](devel), we encourage you to share your results with the
community by submitting a pull request on the main repository so that
others can benefit from the same changes. Similarly, if a user feels
they have identified a need that is not available within the product
but does not feel they have the skillset to handle them, then these
issues should be posted on GitHub to seek advice or for developers to
consider.

This documentation is based on the latest version of the OSOP toolkit. Any
future changes to the toolkit by the producers may affect this
guidance document. It is recommended to regularly update the toolkit
locally to ensure that it is in sync with the main remote branch and
includes all recent changes. Guidance on how to do this can be found
in the section [Development](devel).
