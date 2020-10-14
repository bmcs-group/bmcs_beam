# from .beam_design import BeamDesign
from beam_design.beam_design import BeamDesign

from .cross_section_shape import CrossSectionShapeBase
from bmcs_utils.api import InteractiveModel, Item, View
import traits.api as tr
import numpy as np


class Reinforcement(InteractiveModel):
    name = 'Reinforcement'

    # TODO->Saeed: prepare the varibles for InteractiveModel (ipw_view and so on...)

    z_j = tr.Array(np.float_, value=[50])
    """z positions of reinforcement layers"""

    A_j = tr.Array(np.float_, value=[np.pi * (16 / 2.) ** 2])
    """cross section area of reinforcement layers"""

    E_j = tr.Array(np.float_, value=[210000])
    """E modulus of reinforcement layers"""

    eps_sy_j = tr.Array(np.float_, value=[500. / 210000.])
    """Steel yield strain"""


class Fabric(Reinforcement):
    E_carbon = tr.Int(200000)
    width = tr.Float(8)
    thickness = tr.Float(1)
    spacing = tr.Float(1)
    n_layers = tr.Int(1)
    A_roving = tr.Float(1)
    f_h = tr.Int(5)

    ipw_view = View(
        Item('E_carbon', param=True, latex='E_r \mathrm{[MPa]}', minmax=(200000, 300000)),
        Item('width', param=True, latex='rov_w \mathrm{[mm]}', minmax=(8, 450)),
        Item('thickness', param=True, latex='rov_t \mathrm{[mm]}', minmax=(1, 100)),
        Item('spacing', param=True, latex='ro_s \mathrm{[mm]}', minmax=(1, 100)),
        Item('n_layers', param=True, latex='n_l \mathrm{[-]}', minmax=(1, 100)),
        Item('A_roving', param=True, latex='A_r \mathrm{[mm^2]}', minmax=(1, 100)),
        Item('f_h', param=True, latex='f_h \mathrm{[mm]}', minmax=(5, 500))
    )


class Bar(Reinforcement):
    """" pass """


class Matrix(InteractiveModel):
    name = 'Matrix'

    ipw_view = View(
        Item('E_ct', minmax=(10, 50000), latex='E_{ct} [N/mm^2]'),
        Item('E_cc', minmax=(10, 50000), latex='E_{cc} [N/mm^2]')
        # TODO->Saeed: complete these
    )

    E_ct = tr.Float(24000)
    """E modulus of matrix on tension"""

    E_cc = tr.Float(25000)
    """E modulus of matrix on compression"""

    eps_cr = tr.Float(0.001)
    """Matrix cracking strain"""

    eps_cy = tr.Float(-0.003)
    """Matrix compressive yield strain"""

    eps_cu = tr.Float(-0.01)
    """Ultimate matrix compressive strain"""

    eps_tu = tr.Float(0.003)
    """Ultimate matrix tensile strain"""

    mu = tr.Float(0.33)
    """Post crack tensile strength ratio (represents how much strength is left after the crack because of short steel 
    fibers in the mixture)"""

    # TODO: what if I don't want to plot anything, just change params? is this a good approach?
    def update_plot(self, axes):
        pass


class CrossSectionLayout(InteractiveModel):
    name = 'CrossSectionLayout'

    matrix = tr.Instance(Matrix, ())
    reinforcement = tr.Instance(Reinforcement, ())

    beam_design = tr.WeakRef
    cross_section_shape = tr.DelegatesTo('beam_design')

    # def _cross_section_shape_default(self):
    #     return self.beam_design.cross_section_shape

    H = tr.DelegatesTo('cross_section_shape')

    # TODO: what params need to show here? is there any?
    ipw_view = View(
    )

    def get_comp_E(self):
        # todo: check it with the bmcs example
        A_composite = self.b * self.H
        n_rovings = self.width / self.spacing  # width or B??
        A_layer = n_rovings * self.A_roving
        A_carbon = self.n_layers * A_layer
        A_concrete = A_composite - A_carbon
        E_comp = (self.E_carbon * A_carbon + self.E_con * A_concrete) / (A_composite)
        return E_comp

    def subplots(self, fig):
        return fig.subplots(1, 1)

    # TODO: what should be plotted here? cross section with steel positions?
    def update_plot(self, ax):
        self.cross_section_shape.update_plot(ax)

        # Quick fix with constant b
        # b_ = 100
        # ax.axis([0, 100, 0, self.H])
        # ax.axis('equal')
        # ax.fill([0, b_, b_, 0, 0], [0, 0, self.H, self.H, 0], color='gray')
        # ax.plot([0, b_, b_, 0, 0], [0, 0, self.H, self.H, 0], color='black')

        # Original from notebook
        # ax.axis([0, self.b, 0, self.H])
        # ax.axis('equal')
        # ax.fill([0, self.b, self.b, 0, 0], [0, 0, self.H, self.H, 0], color='gray')
        # ax.plot([0, self.b, self.b, 0, 0], [0, 0, self.H, self.H, 0], color='black')
        # ax.plot([self.b / 2 - self.width / 2, self.b / 2 + self.width / 2], [self.f_h, self.f_h], color='Blue',
        #         linewidth=self.n_layers * self.thickness)
        # ax.annotate('E_composite = {} GPa'.format(np.round(self.get_comp_E() / 1000), 0),
        #             xy=(self.b / 10, self.f_h * 1.1), color='white')