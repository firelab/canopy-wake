#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Utility for evaluation of horizontal homogenity."""

from __future__ import division
from __future__ import unicode_literals

import os
from matplotlib import pyplot as plt
import numpy as np
from math import sqrt, log

PROFILE_NAMES = ('mast1', 'mast2', 'mast3', 'mast4')
TREE_HEIGHT = 7.5
REF_HEIGHT = 2 * TREE_HEIGHT
Z0 = 0.06
KARMAN = 0.40

UREF = 6.17
ZREF = 15
Z0 = 0.06

Cmu = 0.09


def ustar(uref, zref, z0):
    """Calculate friction velocity for a neutral ABL logartithmic profile."""
    return KARMAN * uref / log((zref + z0) / z0)


def k_inlet_profile(z, z0, zref, uref):
    """Calculate k inlet profile."""
    k = np.power(ustar(uref, zref, z0), 2) / sqrt(Cmu)
    return k * np.ones(z.shape)


def read_U(filename):
    """ read z and U from vertical profile."""
    U = np.genfromtxt(filename, dtype=np.float)
    return (U[:, 0], U[:, 1:])


def read_epsilon_k_p(filename):
    """read epsilon, k p from vertical profile."""

    epsilon_k_p = np.genfromtxt(filename, dtype=np.float)
    z = epsilon_k_p[:, 0]
    epsilon = epsilon_k_p[:, 1]
    k = epsilon_k_p[:, 2]
    p = epsilon_k_p[:, 3]
    return z, epsilon, k, p


def main():
    
    times = map(int, os.listdir('../postProcessing/verticalProfiles'))
    latestTime = max(times)
    plt.rc('figure', figsize=(18, 10))  # set the figure size (in inches)
    fig1, (ax1, ax2, ax3) = plt.subplots(1, 3)
    
    filename = '../postProcessing/verticalProfiles/{time}/{profile}_U.xy'.format(
        time=latestTime,
        profile='inlet_profile'
    )
    z, U_inlet = read_U(filename)
    
    filename = '../postProcessing/verticalProfiles/{time}/{profile}_U.xy'.format(
        time=latestTime,
        profile='outlet_profile'
    )
    z, U_outlet = read_U(filename)
    
    filename = '../postProcessing/verticalProfiles/{time}/{profile}_p_k_epsilon.xy'.format(
        time=latestTime,
        profile='inlet_profile'
    )
    z, eps_inlet, k_inlet, p_inlet = read_epsilon_k_p(filename)
    
    filename = '../postProcessing/verticalProfiles/{time}/{profile}_p_k_epsilon.xy'.format(
        time=latestTime,
        profile='outlet_profile'
    )
    z, eps_outlet, k_outlet, p_outlet = read_epsilon_k_p(filename)
    
    ax1.plot(U_inlet[:, 0], z, label='inlet_profile', linewidth=1.75)
    ax1.plot(U_outlet[:, 0], z, label='outlet_profile', linewidth=1.75)
    
    ax2.plot(k_inlet, z, label='inlet_profile', linewidth=1.75)
    ax2.plot(k_outlet, z, label='outlet_profile', linewidth=1.75)
    ax2.plot(k_inlet_profile(z, Z0, ZREF, UREF), z, label='Richards & Hoxey', linewidth=1.75)
    
    ax3.plot(eps_inlet, z, label='inlet_profile', linewidth=1.75)
    ax3.plot(eps_outlet, z, label='outlet_profile', linewidth=1.75)
    
    # U_axis[i].minorticks_on()
    ax1.set_ylim((0, 80))
    ax2.set_ylim((0, 80))
    ax3.set_ylim((0, 80))
    
    ax1.tick_params(axis='x', labelsize=18)
    ax1.tick_params(axis='y', labelsize=18)
    ax2.tick_params(axis='x', labelsize=18)
    ax2.tick_params(axis='y', labelsize=18)
    ax3.tick_params(axis='x', labelsize=18)
    ax3.tick_params(axis='y', labelsize=18)
    
    ax1.tick_params(axis="x", length=4.25, width=1.3)
    ax1.tick_params(axis="y", length=4.25, width=1.3)
    ax2.tick_params(axis="x", length=4.25, width=1.3)
    ax2.tick_params(axis="y", length=4.25, width=1.3)
    ax3.tick_params(axis="x", length=4.25, width=1.3)
    ax3.tick_params(axis="y", length=4.25, width=1.3)
    
    ax1.legend(loc='best', shadow=True, fontsize=18)
    ax2.legend(loc='best', shadow=True, fontsize=18)
    ax3.legend(loc='best', shadow=True, fontsize=18)
    
    ax1.set_ylabel('z (m)', fontsize=20)
    ax2.set_ylabel('z (m)', fontsize=20)
    ax3.set_ylabel('z (m)', fontsize=20)
    
    ax1.set_title('Wind speed', fontsize=18)
    ax2.set_title('Turbulent kinetic energy', fontsize=18)
    ax3.set_title('Dissipation of turbulent kinetic energy', fontsize=18)
    
    ax1.grid(True)
    ax1.set_xlabel('u (m/s)', fontsize=20)
    
    ax2.grid(b=True)
    ax2.set_xlabel('k (m²/s²)', fontsize=20)
    
    ax3.grid(True)
    ax3.set_xlabel('$\epsilon$ (m²/s³)', fontsize=20)
    
    fig1.tight_layout(rect=[0.01,0.01,0.99, 0.99])
    
    
    plt.savefig('../zz_pictures/a_inletOutletProfiles.png')
    plt.show()
    

if __name__ == '__main__':
    main()
