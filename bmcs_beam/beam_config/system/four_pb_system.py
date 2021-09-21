from bmcs_utils.api import Model, Int, Item, View, Float
from bmcs_beam.beam_config.system.system import System

class FourPBSystem(System):

    name = 'FourPBSystem'

    F = Float(-1000)
    L_F = Float(1000, desc='length from support to the force')

    tree = []

    ipw_view = View(
        *System.ipw_view.content,
        Item('F', latex='F \mathrm{[N]}'),
        Item('L_F', latex='L_F \mathrm{[N]}'),
    )

    @staticmethod
    def subplots(fig):
        axes = fig.subplots(4, 2)
        plt.subplots_adjust(hspace=1)
        fig.set_size_inches(7, 8)
        return axes

    def update_plot(self, axes):
        self._update_struct()
        self._plot_struct(axes[0, 0])
        # self._plot_internal_forces(axes[1:, 0])

    def _update_struct(self):
        self.struct = self.get_new_struct()

        n_x_to_F = round((self.L_F / self.L) * self.n_x)
        n_x_between_forces = self.n_x - 2 * n_x_to_F
        self.struct.add_multiple_elements([[0, 0], [self.L_F, 0]], n_x_to_F)
        self.struct.add_multiple_elements([[self.L_F, 0], [self.L - self.L_F, 0]], n_x_between_forces)
        self.struct.add_multiple_elements([[self.L - self.L_F, 0], [self.L, 0]], n_x_to_F)

        self.struct.add_support_hinged(1)
        self.struct.add_support_roll(self.n_x + 1)

        self.struct.point_load(n_x_to_F + 1, Fy=self.F)
        self.struct.point_load(n_x_to_F + n_x_between_forces + 1, Fy=self.F)
        self.struct.solve()
