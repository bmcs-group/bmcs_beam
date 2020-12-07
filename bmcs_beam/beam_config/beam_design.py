import traits.api as tr
from bmcs_beam.beam_config.boundary_conditions import BoundaryConditions
from bmcs_cross_section.cs_design import CrossSectionDesign
from bmcs_utils.api import InteractiveModel, \
    Int, Item, View, Float, Range, Button, ButtonEditor, mpl_align_yaxis, FloatRangeEditor
from sympy.physics.continuum_mechanics.beam import Beam
import sympy as sp
import numpy as np
from numbers import Number

# Quick fix for [bmcs_shear_zone]
class BeamDesign(CrossSectionDesign):

# class BeamDesign(InteractiveModel):
    name = 'Beam Design'

    # cs_design = tr.Instance(CrossSectionDesign, ())

    n_x = Int(100)
    L = Float(5000)
    F = Float(-1000)
    # TODO [HS]: make the conf a dropdown menu when the implementation is done by [RC]
    beam_conf_name = BoundaryConditions.CONFIG_NAMES[0]
    beam_conf_name_slider = Range(0.)

    beam = tr.Instance(Beam)

    def _beam_default(self):
        return BoundaryConditions.get_configured_beam(self.L, self.F, self.beam_conf_name)

    x = tr.Property(depends_on='L, n_x')

    @tr.cached_property
    def _get_x(self):
        return np.linspace(0, self.L, self.n_x)

    @tr.observe("beam_conf_name_slider")
    def notify_slider_change(self, event):
        self.beam_conf_name = BoundaryConditions.CONFIG_NAMES[int(event.new)]
        self.beam = BoundaryConditions.get_configured_beam(self.L, self.F, self.beam_conf_name)

    @tr.observe("L")
    def notify_L_change(self, event):
        self.beam = BoundaryConditions.get_configured_beam(event.new, self.F, self.beam_conf_name)

    @tr.observe("F")
    def notify_F_change(self, event):
        self.beam = BoundaryConditions.get_configured_beam(self.L, event.new, self.beam_conf_name)

    # add_force_btn = Button()
    #
    # def _add_force_btn_fired(self):
    #     self.L += 100

    ipw_view = View(
        Item('L', latex='L \mathrm{[mm]}'),
        Item('F', latex='F \mathrm{[N]}'),
        Item('n_x', latex='n_x'),
        Item('beam_conf_name_slider', editor=FloatRangeEditor(label='Beam config', n_steps=5, low=0, high=5)),
        # Item('add_force_btn', editor=ButtonEditor(label='Force', icon='plus'))
    )

    def solve_beam_for_reactions(self):
        reactions_names = []

        for load in self.beam.applied_loads:
            load_value = load[0]
            if not isinstance(load_value, Number):
                reactions_names.append(load_value)

        self.beam.solve_for_reaction_loads(*reactions_names)

    def get_M_x(self, solve_beam_first=True):
        if solve_beam_first:
            self.solve_beam_for_reactions()
        x = sp.symbols('x')
        M_ = self.beam.bending_moment().rewrite(sp.Piecewise)
        get_M = sp.lambdify(x, M_, 'numpy')
        return get_M(self.x)

    def get_Q_x(self, solve_beam_first=True):
        if solve_beam_first:
            self.solve_beam_for_reactions()
        x = sp.symbols('x')
        Q_ = self.beam.shear_force().rewrite(sp.Piecewise)
        get_Q = sp.lambdify(x, Q_, 'numpy')
        Q_x = get_Q(self.x)
        return Q_x

    def plot_MQ(self, ax2, ax3):
        x = self.x
        M_scale = 1e+6
        Q_scale = 1000

        M_x = self.get_M_x(solve_beam_first=True) / M_scale
        Q_x = self.get_Q_x(solve_beam_first=False) / Q_scale

        ax2.plot(x, -M_x, color='red', label='moment [kNm]')
        ax2.fill(x, -M_x, color='red', alpha=0.1)
        ax2.set_ylabel('M [kNm]')
        ax2.legend()

        ax3.plot(x, Q_x, lw=0.1, color='green', label='shear [kN]')
        ax3.fill_between(x, Q_x, 0, color='green', alpha=0.1)
        ax3.set_ylabel('Q [kN]')
        ax3.set_xlabel('x [mm]')
        ax3.legend()
        mpl_align_yaxis(ax2, 0, ax3, 0)

    def subplots(self, fig):
        ax1, ax2 = fig.subplots(2, 1)
        ax3 = ax2.twinx()
        return ax1, ax2, ax3

    def update_plot(self, axes):
        ax1, ax2, ax3 = axes
        BoundaryConditions.plot(ax1, self.beam)
        self.plot_MQ(ax2, ax3)
        ax1.axis('equal')
        ax1.autoscale(tight=True)

    # Quick fix: needed for [bmcs_shear_zone]
    def plot_reinforcement(self, ax):
        L = self.L
        for z in self.cross_section_layout.reinforcement.z_j:
            ax.plot([0,L],[z,z], lw=2, color='brown')