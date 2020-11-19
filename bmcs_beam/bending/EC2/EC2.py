import matplotlib.patches as mpatches
import numpy as np
import sympy as sp
import traits.api as tr
from bmcs_cross_section.mkappa import MKappa
from bmcs_beam.beam_config.boundary_conditions import BoundaryConditions
from bmcs_utils.api import InteractiveModel, \
    Item, View, Float, Int, FloatEditor, FloatRangeEditor, mpl_align_yaxis
from matplotlib.patches import PathPatch
from matplotlib.path import Path
from sympy.physics.continuum_mechanics.beam import Beam

class UltimateLimitState(InteractiveModel):

    F_max = Float(1, desc='maximum load value', BC=True)
    theta_F = Float(1.0, desc='load factor', BC=True)

    F = tr.Property(depends_on='+BC')

    ipw_view = View(
        Item('theta_F', latex=r'\theta [-]',
             editor=FloatRangeEditor(low=0, high=1)),
        Item('F_max', latex='F_\mathrm{max} [\mathrm{N}]',
             readonly=True),
        Item('G_adj', param=True, latex='G_{adj} [\mathrm{-}]', minmax=(1e-3, 1e-1)),
        Item('n_x', param=True, latex='n_x [\mathrm{-}]'),
    )

    x = tr.Property(depends_on='_GEO')

    @tr.cached_property
    def _get_x(self):
        pass

    def get_M_x(self):
        pass

    def plot_geo(self, ax1):
        pass

    def subplots(self, fig):
        pass

    def update_plot(self, axes):
        ax1, ax2, ax3 = axes
        self.plot_geo(ax1)
        self.plot_MQ(ax2, ax3)
        ax1.axis('equal')
        ax1.autoscale(tight=True)

class ServiceabilityLimitState(InteractiveModel):

    F_max = Float(1, desc='maximum load value', BC=True)
    theta_F = Float(1.0, desc='load factor', BC=True)

    F = tr.Property(depends_on='+BC')

    ipw_view = View(
        Item('theta_F', latex=r'\theta [-]',
             editor=FloatRangeEditor(low=0, high=1)),
        Item('F_max', latex='F_\mathrm{max} [\mathrm{N}]',
             readonly=True),
        Item('G_adj', param=True, latex='G_{adj} [\mathrm{-}]', minmax=(1e-3, 1e-1)),
        Item('n_x', param=True, latex='n_x [\mathrm{-}]'),
    )

    x = tr.Property(depends_on='_GEO')

    @tr.cached_property
    def _get_x(self):
        pass

    def get_M_x(self):
        pass

    def plot_geo(self, ax1):
        pass

    def subplots(self, fig):
        ax1, ax2 = fig.subplots(2, 1)
        ax3 = ax2.twinx()
        return ax1, ax2, ax3

    def update_plot(self, axes):
        ax1, ax2, ax3 = axes
        self.plot_geo(ax1)
        self.plot_MQ(ax2, ax3)
        ax1.axis('equal')
        ax1.autoscale(tight=True)