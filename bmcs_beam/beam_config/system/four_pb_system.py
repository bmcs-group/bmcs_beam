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

    def get_force_from_max_moment(self, m):
        # This should be done numerically, but here for efficiency and
        # because this is a known case it's given directly
        return m / self.L_F
