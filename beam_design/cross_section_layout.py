from beam_design.cross_section_shape import CrossSectionShape
from bmcs_utils.api import InteractiveModel, Item, View
import traits.api as tr
import numpy as np


class Reinforcement(InteractiveModel):
    z_j = tr.Array(np.float_, value=[50])
    """z positions of reinforcement layers"""

    A_j = tr.Array(np.float_, value=[np.pi * (16 / 2.) ** 2])
    """cross section area of reinforcement layers"""

    E_j = tr.Array(np.float_, value=[210000])
    """E modulus of reinforcement layers"""

    eps_sy_j = tr.Array(np.float_, value=[500. / 210000.])
    """Steel yield strain"""


class Matrix(InteractiveModel):
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


class CrossSectionLayout(InteractiveModel):
    name = 'CrossSectionLayout'

    cross_section_shape = tr.Instance(CrossSectionShape, ())

    matrix = tr.Instance(Matrix, ())
    reinforcement = tr.Instance(Reinforcement, ())

    # Reinforcement
    E_carbon = tr.Int(200000)
    width = tr.Float(8)
    thickness = tr.Float(1)
    spacing = tr.Float(1)
    n_layers = tr.Int(1)
    A_roving = tr.Float(1)
    f_h = tr.Int(5)

    # Concerte cross section
    #     L = tr.Int(5000, param=True, latex='L \mathrm{mm}', minmax=(10,10000))
    h = tr.DelegatesTo('cross_section_shape')
    b = tr.DelegatesTo('cross_section_shape')
    E_con = tr.Int(14000)
    F = tr.Float(10)
    n_x = tr.Int(100)

    ipw_view = View(
        Item('E_carbon', param=True, latex='E_r \mathrm{[MPa]}', minmax=(200000, 300000)),
        Item('width', param=True, latex='rov_w \mathrm{[mm]}', minmax=(8, 450)),
        Item('thickness', param=True, latex='rov_t \mathrm{[mm]}', minmax=(1, 100)),
        Item('spacing', param=True, latex='ro_s \mathrm{[mm]}', minmax=(1, 100)),
        Item('n_layers', param=True, latex='n_l \mathrm{[-]}', minmax=(1, 100)),
        Item('A_roving', param=True, latex='A_r \mathrm{[mm^2]}', minmax=(1, 100)),
        Item('f_h', param=True, latex='f_h \mathrm{[mm]}', minmax=(5, 500)),
        Item('E_con', param=True, latex='E \mathrm{[MPa]}', minmax=(14000, 41000)),
        Item('F', param=True, latex='F \mathrm{[N]}', minmax=(10, 1000)),
        Item('n_x', param=True, latex='n_x \mathrm{[-]}', minmax=(1, 1000))
    )

    def get_comp_E(self):
        '''todo: check it with the bmcs example'''
        A_composite = self.b * self.h
        n_rovings = self.width / self.spacing  # width or B??
        A_layer = n_rovings * self.A_roving
        A_carbon = self.n_layers * A_layer
        A_concrete = A_composite - A_carbon
        E_comp = (self.E_carbon * A_carbon + self.E_con * A_concrete) / (A_composite)
        return E_comp

    def subplots(self, fig):
        return fig.subplots(1, 1)

    def update_plot(self, ax):
        ax.axis([0, self.b, 0, self.h])
        ax.axis('equal')
        ax.fill([0, self.b, self.b, 0, 0], [0, 0, self.h, self.h, 0], color='gray')
        ax.plot([0, self.b, self.b, 0, 0], [0, 0, self.h, self.h, 0], color='black')
        ax.plot([self.b / 2 - self.width / 2, self.b / 2 + self.width / 2], [self.f_h, self.f_h], color='Blue',
                linewidth=self.n_layers * self.thickness)
        ax.annotate('E_composite = {} GPa'.format(np.round(self.get_comp_E() / 1000), 0),
                    xy=(self.b / 10, self.f_h * 1.1), color='white')