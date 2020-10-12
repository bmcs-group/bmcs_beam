import matplotlib.pyplot as plt
from moment_curvature import MomentCurvature

import numpy as np
import sympy as sp


def run_example_with_default_params():
    mc = MomentCurvature(idx=25, n_m=100)
    # mc.h = 600
    # mc.b = 200
    mc.kappa_range = (-0.00002, 0.00002, 100)
    fig, ((ax1, ax2)) = plt.subplots(1, 2, figsize=(10, 5))
    mc.plot(ax1, ax2)
    plt.show()


def run_example_with_t_section_and_custom_params():
    mc = MomentCurvature(idx=25, n_m=100)

    # Material parameters [mm], [N/mm2]
    mc.h = 666
    mc.E_ct = 24000
    mc.E_cc = 24000
    mc.eps_cr = 0.000125
    mc.eps_cy = 0.0010625  # 8.5 * eps_cr_
    mc.eps_cu = 0.0035
    mc.eps_tu = 0.02
    mc.mu = 0.33

    # 2 layers reinforcement details
    mc.A_j = np.array([250, 0])  # A_j[0] for tension steel / A_j[1] for compression steel
    mc.z_j = np.array([0.1 * mc.h, 0.9 * mc.h])
    mc.E_j = np.array([210000, 210000])
    mc.eps_sy_j = np.array([0.002, 0.002])

    # Defining a variable width (T-section as an example)
    z = sp.Symbol('z')
    b_w = 50
    b_f = 500
    h_w = 0.85 * mc.h
    # Beam width b as a function of the height z (the sympy z symbol in MomentCurvatureSymbolic is used)
    mc.b = sp.Piecewise((b_w, z < h_w), (b_f, z >= h_w))

    # If plot_norm is used, use the following:
    # mc.kappa_range = (0, mc.kappa_cr * 100, 100
    mc.kappa_range = (-0.00002, 0.00002, 100)

    # Plotting
    fig, ((ax1, ax2)) = plt.subplots(1, 2, figsize=(10, 5))
    mc.plot(ax1, ax2)
    plt.show()


if __name__ == '__main__':
    run_example_with_default_params()
    # run_example_with_t_section_and_custom_params()