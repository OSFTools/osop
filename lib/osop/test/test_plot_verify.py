# (C) Crown Copyright, Met Office. All rights reserved.

# This file is part of osop and is released under the BSD 3-Clause license.
# See LICENSE in the root of the repository for full licensing details.

"""Tests for plot_verify module."""

import cartopy.feature
import pytest

import osop.plot_verify as plot_verify


@pytest.fixture
def config():
    """Fixture that provides an empty config dictionary."""
    return {}


def test_location_country(config):
    """Test that the correct object is returned when a country name is provided in the config."""
    config["border"] = "Morocco"
    local = plot_verify.location(config)

    assert local.name == "admin_0_countries_mar"
    assert isinstance(local, cartopy.feature.NaturalEarthFeature)


def test_location_none(config):
    """Test that False is returned when no location is provided in the config."""
    config["border"] = "None"
    local = plot_verify.location(config)

    assert local == "False"


def test_location_invalid(config):
    """Test that a KeyError is raised when an invalid location is provided in the config."""
    config["border"] = "Spam"

    with pytest.raises(KeyError, match="Location Name does not exist in dictionary"):
        plot_verify.location(config)
