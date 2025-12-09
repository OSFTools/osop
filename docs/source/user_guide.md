# OSOP User Guide

```{eval-rst}
.. toctree::
   :maxdepth: 2
   :caption: Contents:
```

## Introduction

This document provides guidance on using the Objective Seasonal
Outlook Package (OSOP) toolkit for producing seasonal based
forecasting and verification. 

This user guide includes three main sections:

* [Initial set_up](set_up)
* [Running the toolkit](run)
* [Developing the toolkit](devel)

The text, as laid out here, is a
baseline for how to get started with the package. As such, the user
can modify or change the OSOP toolkit away from what is suggested here
as needed to suit their requirements. This document is based on the
latest version of the OSOP toolkit. Any
future changes to the toolkit by the producers may affect this
guidance document. It is recommended to regularly update the toolkit
locally to ensure that it is in sync with the main remote branch and
includes all recent changes. Guidance on how to do this can be found
in Appendix 1.

The package consists of a "driver script" structure. That is to say a
few top-level entry point scripts that orchestrate smaller pipeline
scripts. These smaller scripts have modular functions that can be used
in turn to produce hindcast verification, single service forecasts and
functions are written in python, with this holding up the base of the
toolkit. The top-level entry is written in a combination of bash
scripting and yaml to allow for simple execution and variable
management.

With all of this in consideration the user does not need to be
"fluent" in coding or these specific languages to get coherent use of
the OSOP toolkit. It is designed that after set-up, its base usage is
easy and comprehensible as well as versatile enough to cover the main
directives of the product. Users more advanced in coding may find that
they can edit further down in the software to get more highly specific
use cases. As outlined in the appendix, it is encouraged that if this
is done the ideas or solved selection of code should be pushed back up
to the main repository so that others can benefit from the same
changes. Similarly, if a user feels they have identified a need that
is not quite encapsulated within the product but does not feel they
have the skillset to handle them, then these issues should be posted
on the GitHub board to seek advice or for developers to take on.
