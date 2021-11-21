from bmcs_utils.api import Model, Int, Item, View, Float, mpl_align_yaxis_to_zero, \
    mpl_show_one_legend_for_twin_axes
from anastruct import SystemElements
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import traits.api as tr
from bmcs_beam.beam_config.system.anastruct_custom_plotter import CustomPlotter

class System(Model):

    name = 'System'

    L = Float(3000)
    n_x = Int(20, desc='Number of discretization points along the beam.')
    struct = SystemElements()

    tree = []

    ipw_view = View(
        Item('L', latex='L \mathrm{[mm]}'),
        Item('n_x', latex='n_x'),
    )

    x = tr.Property(depends_on='L, n_x')

    @tr.cached_property
    def _get_x(self):
        return np.linspace(0, self.L, self.n_x)

    def __init__(self, *args, **kw):
        super().__init__(*args, **kw)
        self._update_struct()

    def _get_L(self):
        return (self.struct.element_map[self.n_x].vertex_2.x - self.struct.element_map[1].vertex_1.x)

    def get_new_struct(self):
        struct = SystemElements()
        struct.plotter = CustomPlotter(struct, mesh=50)
        return struct

    def _update_struct(self):
        return NotImplementedError

    def subplots(self, fig):
        axes = fig.subplots(3, 1)
        return axes

    def update_plot(self, axes):
        self._update_struct()
        self.plot_struct(axes[0])
        self.plot_bending_moment(axes[1])
        self.plot_shear(axes[2])

    def plot_struct(self, ax):
        ax.get_yaxis().set_visible(False)
        self.struct.plotter.plot_structure(figsize=(8, 2), verbosity=1, loads=True, ax=ax)
        ax.set_xlabel('x [mm]')

    def plot_shear(self, ax):
        self.struct.plotter.shear_force(figsize=(8, 2), factor=1, adjust_ax_height=False,
                                           verbosity=1, ax=ax, show=False)
        ax.set_ylabel('Q [N]')

    def plot_bending_moment(self, ax):
        self.struct.plotter.bending_moment(figsize=(8, 2), factor=1, adjust_ax_height=False,
                                           verbosity=1, ax=ax, show=False)
        ax.set_ylabel('M [Nmm]')
        # M_x = np.array(self.struct.get_element_result_range(unit = "moment")) / M_scale