import matplotlib.patches as mpatches
import numpy as np
import sympy as sp
import traits.api as tr
from bmcs_beam.moment_curvature.moment_curvature import MomentCurvature
from bmcs_utils.api import InteractiveModel, Item, View
from matplotlib.patches import PathPatch
from matplotlib.path import Path
from scipy.integrate import cumtrapz
from sympy.physics.continuum_mechanics.beam import Beam

class MomentCurvatureProfile(InteractiveModel):
    name = 'M_K Profile'
    mc = MomentCurvature(kappa_range=(-0.0002, 0.0002, 100),
        idx=25, n_m=100)

    ipw_view = View(

    )

    def subplots(self, fig):
        return fig.subplots(1, 1)

    def update_plot(self, axes):
        ax1 = axes

        ax1.plot(self.mc.kappa_t, self.mc.M_t)
        ax1.set_xlabel('Kappa')
        ax1.set_ylabel('Moment')
        ax1.grid(True)