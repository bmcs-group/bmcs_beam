from bmcs_utils.api import Model, Int, Item, View, Float
from bmcs_beam.beam_config.system.system import System

class ThreePBSystem(System):

    name = 'ThreePBSystem'

    F = Float(-1000)

    tree = []

    ipw_view = View(
        *System.ipw_view.content,
        Item('F', latex='F \mathrm{[N]}'),
    )

    @staticmethod
    def subplots(fig):
        axes = fig.subplots(2, 1)
        return axes

    def update_plot(self, axes):
        self._update_struct()
        self.plot_struct(axes[0])
        # self._plot_internal_forces(axes[1])

    def _update_struct(self):
        self.struct = self.get_new_struct()

        n_x_to_F = round(0.5 * self.n_x)
        n_x_remains = self.n_x - n_x_to_F
        self.struct.add_multiple_elements([[0, 0], [self.L / 2, 0]], n_x_to_F)
        self.struct.add_multiple_elements([[self.L / 2, 0], [self.L, 0]], n_x_remains)

        self.struct.add_support_hinged(1)
        self.struct.add_support_roll(self.n_x + 1)

        self.struct.point_load(n_x_to_F + 1, Fy=self.F)
        self.struct.solve()