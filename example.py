import argparse

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Ellipse

from ellipse import LsqEllipse


def get_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('--center', nargs=2, type=float, default=(1, 1), metavar=('X', 'Y'))
    parser.add_argument('--width', default=1, type=float)
    parser.add_argument('--height', default=.6, type=float)
    parser.add_argument('--phi', default=np.pi/5, type=float)
    parser.add_argument('--noise', default=.1, type=float)

    return parser.parse_args()


def make_test_ellipse(center, width, height, phi, noise):
    """Generate Elliptical data with noise

    Parameters
    ----------
    center: list:float
        (<x_location>, <y_location>)
    width: float
        semimajor axis. Horizontal dimension of the ellipse (**)
    height: float
        semiminor axis. Vertical dimension of the ellipse (**)
    phi: float:radians
        tilt of the ellipse, the angle the semimajor axis
        makes with the x-axis

    Returns
    -------
    data:  list:list:float
        list of two lists containing the x and y data of the ellipse.
        of the form [[x1, x2, ..., xi],[y1, y2, ..., yi]]
    """
    t = np.linspace(0, 2*np.pi, 1000)
    x_noise, y_noise = np.random.randn(2, len(t))

    ellipse_x = center[0] + width*np.cos(t)*np.cos(phi)-height*np.sin(t)*np.sin(phi) + noise*x_noise  # noqa: E501
    ellipse_y = center[1] + width*np.cos(t)*np.sin(phi)+height*np.sin(t)*np.cos(phi) + noise*y_noise  # noqa: E501

    return [ellipse_x, ellipse_y]


def report_results(args, center, width, height, phi):

    given_phi = args.phi/np.pi
    found_phi = phi/np.pi

    given_center_msg = f"({args.center[0]:.3f}, {args.center[1]:.3f})"
    found_center_msg = f"({center[0]:.3f}, {center[1]:.3f})"
    given_width_msg = f"{args.width:.3f}"
    found_width_msg = f"{width:.3f}"
    given_height_msg = f"{args.height:.3f}"
    found_height_msg = f"{height:.3f}"
    given_phi_msg = f"{given_phi:.3f}\N{GREEK SMALL LETTER PI}"
    found_phi_msg = f"{found_phi:.3f}\N{GREEK SMALL LETTER PI}"

    print(f"{'Parameter':_<20}{'Given':_<20}{'Found':_<20}")
    print(f"{'center':_<20}{given_center_msg:_<20}{found_center_msg:_<20}")
    print(f"{'width':_<20}{given_width_msg:_<20}{found_width_msg:_<20}")
    print(f"{'height':_<20}{given_height_msg:_<20}{found_height_msg:_<20}")
    print(f"{'phi':_<20}{given_phi_msg:_<20}{found_phi_msg:_<20}")


if __name__ == '__main__':
    args = get_args()

    X1, X2 = make_test_ellipse(
        center=args.center, width=args.width, height=args.height, phi=args.phi, noise=args.noise
    )

    X = np.array(list(zip(X1, X2)))
    reg = LsqEllipse().fit(X)
    center, width, height, phi = reg.as_parameters()

    report_results(args, center, width, height, phi)

    fig = plt.figure(figsize=(6, 6))
    ax = plt.subplot()
    ax.axis('equal')
    ax.plot(X1, X2, 'ro', zorder=1)
    ellipse = Ellipse(
        xy=center, width=2*width, height=2*height, angle=np.rad2deg(phi),
        edgecolor='b', fc='None', lw=2, label='Fit', zorder=2
    )
    ax.add_patch(ellipse)

    plt.xlabel('$X_1$')
    plt.ylabel('$X_2$')

    plt.legend()
    plt.show()
