from .cross_section_layout import CrossSectionLayout
from .cross_section_shape import CrossSectionShape
from bmcs_utils.api import InteractiveModel, Item, View
import traits.api as tr
import numpy as np


class BeamDesign(InteractiveModel):
    cross_section_layout = tr.Instance(CrossSectionLayout, ())
    cross_section_shape = tr.Instance(CrossSectionShape, ())

    name = 'BeamDesign'

    L = tr.Float(5000)
    H = tr.DelegatesTo('cross_section_shape')

    ipw_view = View(
        Item('L', latex='L \mathrm{[mm]}', minmax=(1000, 10000))
    )

    def subplots(self, fig):
        return fig.subplots(1, 1)

    def update_plot(self, ax):
        beam_points_x_coords = [0, self.L, self.L, 0, 0]
        beam_points_y_coords = [0, 0, self.H, self.H, 0]
        ax.fill(beam_points_x_coords, beam_points_y_coords, color='gray')
        ax.plot(beam_points_x_coords, beam_points_y_coords, color='black')
        ax.annotate('L = {} mm'.format(np.round(self.L), 0), xy=(self.L / 2, self.H * 1.1), color='black')
        ax.axis('equal')
