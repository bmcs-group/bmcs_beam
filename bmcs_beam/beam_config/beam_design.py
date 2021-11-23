from numbers import Number

import traits.api as tr
from bmcs_beam.beam_config.system.cantilever_system import CantileverDistLoadSystem
from bmcs_beam.beam_config.system.four_pb_system import FourPBBeamSystem
from bmcs_beam.beam_config.system.simple_dist_load_system import SimpleDistLoadBeamSystem
from bmcs_beam.beam_config.system.three_pb_system import ThreePBBeamSystem
from bmcs_cross_section.cs_design import CrossSectionDesign
from bmcs_utils.api import Item, View, EitherType, Str, Bool


class BeamDesign(CrossSectionDesign):

    name = 'BeamSystem Design'

    state_changed = tr.Event
    beam = EitherType(options=[('three_pb', ThreePBBeamSystem),
                                ('four_pb', FourPBBeamSystem),
                               ('simple_beam_dist_load', SimpleDistLoadBeamSystem),
                               ('cantilever_dist_load', CantileverDistLoadSystem),
                               # ('fixed_support_dist_load', CarbonReinfMatMod),
                               # ('fixed_and_roller_support_dist_load', CarbonReinfMatMod),
                               # ('three_span_dist_load', CarbonReinfMatMod),
                               # ('three_pb_fixed_support', CarbonReinfMatMod),
                               # ('single_moment', CarbonReinfMatMod),
                               ])

    tree = []
    L = tr.DelegatesTo('beam_')
    F = tr.DelegatesTo('beam_')
    x = tr.DelegatesTo('beam_')

    ipw_view = View(
        Item('beam'),
    )

    @tr.observe('beam')
    def _update_state(self, event):
        print('beam changed!!')
        self.state_changed = True

    def solve_beam_for_reactions(self):
        # TODO: fix this
        #  The following is a quick fix in order for 4pb to work, without it the self.beam is somehow corrupted after
        #  using the function get_Fw() in DeflectionProfile so it's not able to solve the beam for 4pb!
        # self.beam = BoundaryConditions.get_configured_beam(self.L, self.F, self.beam_conf_name)

        reactions_names = []

        for load in self.beam.applied_loads:
            load_value = load[0]
            if not isinstance(load_value, Number):
                reactions_names.append(load_value)

        self.beam.solve_for_reaction_loads(*reactions_names)

    def get_M_x(self):
        return self.beam_.get_moment()

    def get_Q_x(self):
        return self.beam_.get_shear()

    def subplots(self, fig):
        return self.beam_.subplots(fig)

    # @todo: [HS] this limits the beam design object to just a MQ profile
    # it is, however meant to be only useful for the deflection calculation.
    # For the shear zone - the interface should be generalized and the
    # beam design should provide more services, namely the access to
    # all the components of the beam, draw the plan from the side view,
    # cross sectional view - indeed the whole design including the supports
    # and loading.
    # **[HS]: I would take this second variant**
    # An alternative is to declare this package just to a beam_statical_system
    # which is purely 1D and not to call it design. Then a BeamDesign would
    # contain both - beam_bc (a statical beam) and cs_design. It would provide
    # functionality to change the parameters of the design - including the
    # length, load configuration and supports (BC in one word).
    #
    def update_plot(self, axes):
        print('plot updated!')
        self.beam_.update_plot(axes)

    # Quick fix: needed for [bmcs_shear_zone]
    def plot_reinforcement(self, ax):
        L = self.L
        for z in self.cross_section_layout.reinforcement.z_j:
            ax.plot([0,L],[z,z], lw=2, color='brown')