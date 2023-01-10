#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Utility for plotting verticle profiles, to evaluate decay of boundary layer profile over the distance downwind"""

from __future__ import unicode_literals

import os
from matplotlib import pyplot as plt
import numpy as np



def read_U(filename):
    """ read z and U from vertical profile."""
    U = np.genfromtxt(filename, dtype=np.float)
    return (U[:, 0], U[:, 1:])


def read_p_k_epsilon(filename):
    """read p, k, epsilon from vertical profile."""

    p_k_epsilon = np.genfromtxt(filename, dtype=np.float)
    z = p_k_epsilon[:, 0]
    p = p_k_epsilon[:, 1]
    k = p_k_epsilon[:, 2]
    epsilon = p_k_epsilon[:, 3]
    return z, p, k, epsilon


def main():

    times = map(int, os.listdir('../postProcessing/singleGraph'))
    latestTime = max(times)
    plt.rc('figure', figsize=(18, 10))  # set the figure size (in inches)
    fig1, (ax1, ax2, ax3) = plt.subplots(1, 3)

    filename = '../postProcessing/singleGraph/{time}/{profile}_U.xy'.format(
        time=latestTime,
        profile='line0x504y'
    )
    z, U_inlet = read_U(filename)

    filename = '../postProcessing/singleGraph/{time}/{profile}_U.xy'.format(
        time=latestTime,
        profile='line1710x504y'
    )
    z, U_outlet = read_U(filename)

    filename = '../postProcessing/singleGraph/{time}/{profile}_p_k_epsilon.xy'.format(
        time=latestTime,
        profile='line0x504y'
    )
    z, p_inlet, k_inlet, eps_inlet = read_p_k_epsilon(filename)

    filename = '../postProcessing/singleGraph/{time}/{profile}_p_k_epsilon.xy'.format(
        time=latestTime,
        profile='line1710x504y'
    )
    z, p_outlet, k_outlet, eps_outlet = read_p_k_epsilon(filename)

    ax1.plot(U_inlet[:, 0], z, label='inlet_profile')
    ax1.plot(U_outlet[:, 0], z, label='outlet_profile')

    ax2.plot(k_inlet, z, label='inlet_profile')
    ax2.plot(k_outlet, z, label='outlet_profile')

    ax3.plot(eps_inlet, z, label='inlet_profile')
    ax3.plot(eps_outlet, z, label='outlet_profile')

    fig1.tight_layout(rect=[0.02,0.02,0.98, 0.98])
    ax1.minorticks_on()
    ax2.minorticks_on()
    ax3.minorticks_on()
    
    ax1.legend(loc='best', shadow=True)
    ax2.legend(loc='best', shadow=True)
    ax3.legend(loc='best', shadow=True)

    ax1.set_ylabel('z [m]')
    ax2.set_ylabel('z [m]')
    ax3.set_ylabel('z [m]')

    ax1.set_title('Wind speed')
    ax2.set_title('Turbulent kinetic energy')
    ax3.set_title('Dissipation of turbulent kinetic energy')

    ax1.grid(True)
    ax1.set_xlabel('u [m/s]')

    ax2.grid(b=True)
    ax2.set_xlabel('k [m²/s²]')

    ax3.grid(True)
    ax3.set_xlabel('$\epsilon$ [m²/s³]')
    
    plt.savefig('../zz_pictures/a_az_inletOutletProfiles.png')

    plt.show()


if __name__ == '__main__':
    main()
