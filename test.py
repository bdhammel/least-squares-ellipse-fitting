from unittest import mock

import numpy as np
import numpy.testing as nptest
import pytest

from ellipse import LsqEllipse


def make_dataset(center, width, height, phi, n_points):
    """Generate Elliptical data with noise"""

    t = np.linspace(0, 1.8 * np.pi, n_points)

    x = (center[0] + width * np.cos(t) * np.cos(phi) - height * np.sin(t) * np.sin(phi))
    y = (center[1] + width * np.cos(t) * np.sin(phi) + height * np.sin(t) * np.cos(phi))

    return np.c_[x, y]


def _normalize_result(major, minor, phi):
    """Map parameters back to original orientation

    If phi > np.pi/2 then we are actually measuring a tall ellipse
    chaneg this to measure angle of the waist
    """

    if phi >= np.pi/2:
        width, height = minor, major
        phi -= np.pi/2
    else:
        width, height = major, minor

    return width, height, phi


# phi needs to be < (1/4 * pi) and width != height or test is degenerate
@pytest.mark.parametrize('center', [[1, 1], [0, 1]])
@pytest.mark.parametrize('width', [.4, 10])
@pytest.mark.parametrize('height', [.2, 3])
@pytest.mark.parametrize('phi', [np.pi/5, np.pi/13])
def test_ellipse_fit(center, width, height, phi):
    X = make_dataset(
        center=center,
        width=width,
        height=height,
        phi=phi,
        n_points=10
    )
    elp = LsqEllipse()
    elp.fit(X)
    _center, _major, _minor, _phi = elp.as_parameters()

    _width, _height, _phi = _normalize_result(_major, _minor, _phi)

    nptest.assert_array_almost_equal(_center, center)
    nptest.assert_almost_equal(_width, width)
    nptest.assert_almost_equal(_height, height)
    nptest.assert_almost_equal(_phi, phi)


def test_minimum_data_points():
    X = make_dataset(
        center=[0, 0],
        width=1,
        height=.5,
        phi=np.pi/20,
        n_points=5
    )
    elp = LsqEllipse()
    elp.fit(X)
    _center, _major, _minor, _phi = elp.as_parameters()

    _width, _height, _phi = _normalize_result(_major, _minor, _phi)

    nptest.assert_array_almost_equal(_center, [0, 0])
    nptest.assert_almost_equal(_width, 1)
    nptest.assert_almost_equal(_height, .5)
    nptest.assert_almost_equal(_phi, np.pi/20)


def test_less_than_minimum_data_points_raises_err():
    X = make_dataset(
        center=[0, 0],
        width=1,
        height=.5,
        phi=0,
        n_points=4
    )
    elp = LsqEllipse()
    with pytest.raises(ValueError):
        elp.fit(X)


def test_cannot_get_coef_without_fitting():
    elp = LsqEllipse()
    with pytest.raises(ValueError):
        elp.coefficients


@pytest.mark.parametrize('n_points', [6, 100])
def test_return_fit_returns_correct_ellipse(n_points):
    X = make_dataset(
        center=[0, 0],
        width=1,
        height=.5,
        phi=0,
        n_points=n_points
    )

    elp = LsqEllipse().fit(X)
    t = np.linspace(0, 1.8 * np.pi, n_points)
    center, major, minor, phi = elp.as_parameters()
    with mock.patch.object(LsqEllipse, 'as_parameters') as as_parameters:
        as_parameters.return_value = (center, *_normalize_result(major, minor, phi))
        x = elp.return_fit(n_points, t=t)

    nptest.assert_array_almost_equal(x, X)


def test_if_perfect_circle():
    X = make_dataset(
        center=[0, 0],
        width=1,
        height=1,
        phi=0,
        n_points=50
    )

    elp = LsqEllipse().fit(X)
    _center, _width, _height, _phi = elp.as_parameters()

    nptest.assert_array_almost_equal(_center, [0, 0])
    nptest.assert_almost_equal(_width, 1)
    nptest.assert_almost_equal(_height, 1)
    # nptest.assert_almost_equal(_phi, 0)


@pytest.mark.xfail()
def test_if_no_ellipse_found():
    """
    This data causes a divide by zero error
    TODO: add in check of this
    """
    X = make_dataset(
        center=[0, 0],
        width=1,
        height=1,
        phi=0,
        n_points=5
    )

    elp = LsqEllipse().fit(X)
    _center, _width, _height, _phi = elp.as_parameters()

    nptest.assert_array_almost_equal(_center, [0, 0])
    nptest.assert_almost_equal(_width, 1)
    nptest.assert_almost_equal(_height, 1)
    # nptest.assert_almost_equal(_phi, 0)
