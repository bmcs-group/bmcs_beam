from bmcs_utils.api import Model, Int, Item, View, Float
from anastruct import SystemElements
import numpy as np
import matplotlib.pyplot as plt

class System(Model):

    name = 'System'

    L = Float(3000)
    n_x = Int(20, desc='Number of discretization points along the beam.')
    struct = SystemElements()
    scale =200

    tree = []

    ipw_view = View(
        Item('L', latex='L \mathrm{[mm]}'),
        Item('n_x', latex='n_x'),
    )

    def __init__(self, *args, **kw):
        super().__init__(*args, **kw)
        self._update_struct()

    def _get_L(self):
        return (self.struct.element_map[self.n_x].vertex_2.x - self.struct.element_map[1].vertex_1.x)

    def get_new_struct(self):
        return SystemElements()

    def _update_struct(self):
        return NotImplementedError

    def _plot_struct(self, axes):
        axes.set_frame_on(False)
        axes.set_aspect('equal', adjustable='box')
        axes.get_yaxis().set_ticks([])

        self._plot_elements(axes)
        self._plot_supports(axes)
        self._plot_loads(axes)

    def _plot_elements(self, axes):
        for i, element in self.struct.element_map.items():
            v1 = element.vertex_1
            v2 = element.vertex_2
            axes.plot([v1.x, v2.x], [v1.y, v2.y], color='k')

    def _plot_supports(self, axes):
        self._plot_hinged_supports(axes)
        self._plot_roll_supports(axes)
        self._plot_fixed_supports(axes)

    def _plot_hinged_supports(self, axes):
        for hinged_support in self.struct.supports_hinged:
            x = hinged_support.vertex.x
            y = hinged_support.vertex.y
            axes.plot((x + self.scale*np.array([0.18, 1.5, -1.5, -0.18])), (y + self.scale*np.array([-0.24, -2, -2, -0.24])),
                      color='k')
            circle = plt.Circle((x, y), self.scale * 0.3, color='k', fill=False)
            axes.add_patch(circle)
            hatch_ys = self.scale * (y + np.array([-2.05, -2.5]))
            hatch_xs = self.scale * (x + np.array([1.4, 1.0]))
            for i in range(5):
                axes.plot(hatch_xs, hatch_ys, color='k')
                hatch_xs -= self.scale * 0.6

    def _plot_roll_supports(self, axes):
        for roll_support in self.struct.supports_roll:
            x = roll_support.vertex.x
            y = roll_support.vertex.y
            axes.plot((x + self.scale*np.array([0.18, 1.5, -1.5, -0.18])), (y + self.scale*np.array([-0.24, -2, -2, -0.24])), color='k')
            axes.plot((x + self.scale*np.array([1.5, -1.5,])), (y + self.scale*np.array([-2.4, -2.4,])), color='k')
            circle = plt.Circle((x, y), self.scale * 0.3, color='k', fill=False)
            axes.add_patch(circle)

    def _plot_fixed_supports(self, axes):
        # TODO: do the implementation to depend on fixed_support positions like in _plot_hinged_supports()
        for fixed_support in self.struct.supports_fixed:
            x = fixed_support.vertex.x
            y = fixed_support.vertex.y
            axes.plot((x + self.scale*np.array([-0.1, -0.1])), (y + self.scale*np.array([-1.5, 1.5])), color='k')
            hatch_ys = self.scale * (y + np.array([1.2, 0.9]))
            hatch_xs = self.scale * (x + np.array([-0.2, -0.6]))
            for i in range(5):
                axes.plot(hatch_xs, hatch_ys, color='k')
                hatch_ys -= self.scale * 0.6


    def _plot_hatch(self,axes):
        # TODO: für Rotation und vereinfachung des Codes
        pass

    def _plot_loads(self, axes):
        self._plot_point_loads(axes)
        self._plot_distributed_loads(axes)

    def _plot_point_loads(self, axes):
        L = self._get_L()
        # for node_i, force in self.struct.loads_point.items():
        #     node = self.struct.node_map[node_i]
        #     if(force[0] != 0 and force[1] == 0):
        #         self._plot_x_arrow(axes,node.vertex.x, node.vertex.y)
        #     elif(force[0] == 0 and force[1] != 0):
        #         self._plot_y_arrow(axes,node.vertex.x, node.vertex.y)
        #     elif(force[0] != 0 and force[1] != 0):
        #         self._plot_x_arrow(axes,node.vertex.x, node.vertex.y)
        #         self._plot_y_arrow(axes,node.vertex.x, node.vertex.y)

        # TODO: This implementation is not adaptive (hard-coded), please change the implementation to something similar
        #  to _plot_hinged_supports(), by taking the forces locations from self.struct and relate all plotting
        #  values to their values
        axes.title.set_text('System')
        axes.axis('equal')
        axes.axis('off')

    def _plot_x_arrow(self, axes, x, y):
        L = self._get_L()
        scale_x = 15
        length_x = 600
        if (self.Fx != 0):
            if (self.Fx > 0):
                arrowFx = mpatches.FancyArrowPatch(
                    (x / L * self.scale * self.n_x - length_x, y -0.5 * self.scale),
                    (x / L * self.scale * self.n_x, y -0.5 * self.scale),
                    mutation_scale=scale_x, facecolor='black',
                    edgecolor='black', arrowstyle='-|>')
                if (self.Fy > 0):
                    if (x < 800):
                        axes.text(-2.75, 0, 'Fx', fontsize=12)
                    elif (x == L):
                        axes.text(self.n_x - 3, -2, 'Fx', fontsize=12)
                    else:
                        axes.text(x / L * self.n_x - length_x, -3, 'Fx', fontsize=12)
                else:
                    axes.text(x / L * self.n_x - length_x, 1, 'Fx', fontsize=12)
            else:
                arrowFx = mpatches.FancyArrowPatch((x / L * self.n_x - length_x, -0.5 * self.scale),
                                                   (x / L * self.n_x, y -0.5 * self.scale),
                                                   mutation_scale=scale_x, facecolor='black',
                                                   edgecolor='black', arrowstyle='-|>')
                if (self.Fy > 0):
                    if (x > 3200 and x < 3800):
                        axes.text(x / L*self.scale * self.n_x - length_x - 2, y-4, 'Fx', fontsize=12)
                    elif (x >= L):
                        axes.text(x / L*self.scale * self.n_x - length_x, y-1, 'Fx', fontsize=12)
                    elif (x == 0):
                        axes.text(x / L*self.scale * self.n_x - length_x - 2, y-4, 'Fx', fontsize=12)
                    else:
                        axes.text(x / L*self.scale * self.n_x - length_x - 2, y-3, 'Fx', fontsize=12)
                else:
                    axes.text(x / L*self.scale * self.n_x - length_x - 2, y + 1, 'Fx', fontsize=12)
            axes.add_patch(arrowFx)

    def _plot_y_arrow(self,axes, x,y):
        L = self.struct.element_map[self.n_x].vertex_2.x - self.struct.element_map[1].vertex_1.x
        scale_y = 15
        length_y = 600
        if (self.Fy != 0):
            if (self.Fy > 0):
                arrowFy = mpatches.FancyArrowPatch((x / L * self.scale * self.n_x, y + length_y),
                                                   (x / L * self.scale * self.n_x, y),
                                                   mutation_scale=scale_y, facecolor='black',
                                                   edgecolor='black', arrowstyle='-|>')
                if (self.Fx > 0):
                    axes.text(x / L * self.scale * self.n_x, y + length_y - 1.25, 'Fy', fontsize=12)
                else:
                    axes.text(x / L * self.scale * self.n_x - 2.5, y + length_y - 1.25, 'Fy', fontsize=12)
            else:
                arrowFy = mpatches.FancyArrowPatch((x / L * self.scale * self.n_x, y - length_y),
                                                   (x / L * self.scale * self.n_x, y),
                                                   mutation_scale=scale_y, facecolor='black',
                                                   edgecolor='black', arrowstyle='-|>')
                if (self.Fx > 0):
                    if (y > 3200):
                        axes.text(x / L * self.scale * self.n_x - 0.5, y-4, 'Fy', fontsize=12)
                    elif (y == 0):
                        axes.text(x / L * self.scale * self.n_x - 0.5, y-4, 'Fy', fontsize=12)
                    else:
                        axes.text(x / L * self.scale * self.n_x + 0.5, y + length_y, 'Fy', fontsize=12)
                else:
                    if (y > 3800):
                        axes.text(x / L * self.scale * self.n_x - 0.5, y-4, 'Fy', fontsize=12)
                    elif (y < 800):
                        axes.text(x / L * self.scale * self.n_x - 0.5, y-4, 'Fy', fontsize=12)
                    else:
                        axes.text(x / L * self.scale * self.n_x - 2, y + length_y, 'Fy', fontsize=12)
            axes.add_patch(arrowFy)

    def _plot_distributed_loads(self, axes):
        return
        # for element_i, q in self.struct.loads_q.items():

    def _plot_internal_forces(self, axes):
        L = self.struct.element_map[self.n_x].vertex_2.x - self.struct.element_map[1].vertex_1.x

        # N,Q,M zeichnen

        # N-Plot
        [nx, ny] = self.struct.show_axial_force(values_only=True, factor=1)
        if (self.x_force_x != 0):
            axes[0].plot(nx, ny, 'g')
            axes[0].fill(nx, ny, 'g', alpha=0.3)
        axes[0].plot([0, L], [0, 0], 'k')

        # Q-Plot
        [qx, qy] = self.struct.show_shear_force(values_only=True, factor=1)
        for node_i, force in self.struct.loads_point.items():
            node = self.struct.node_map[node_i]
            # Force position
            x = node.vertex.x
            #y = node.vertex.y
            if (x != 0 and x != L):
                axes[1].plot(qx, -qy, 'r')
                axes[1].fill(qx, -qy, 'r', alpha=0.3, )
        axes[1].plot([0, L], [0, 0], 'k')

        # M-Plot
        [mx, my] = self.struct.show_bending_moment(values_only=True, factor=1)
        my = [i / -1000 for i in my]
        axes[2].plot(mx, my, 'b')
        axes[2].plot([0, L], [0, 0], 'k')
        axes[2].fill(mx, my, 'b', alpha=0.3)

        ## Optik der Plots

        # N-Plot hübsch
        axes[0].title.set_text('N [KN]')
        axes[0].axis('off')
        n_max = self.get_max_N()
        if (self.Fx > 0):
            axes[0].set_ylim([-n_max - 10, n_max + 10])
            axes[0].plot([0, L], [(n_max + 10) / 10, (n_max + 10) / 10], 'k--')
            axes[0].invert_yaxis()
        else:
            axes[0].set_ylim([n_max + 10, -n_max - 10])
            axes[0].plot([0, L], [(n_max + 10) / 10, (n_max + 10) / 10], 'k--')

        if(self.Fx>0):
            axes[0].text(self.x_force_x/2 -0.7*self.scale, self.Fx + 0.3*self.Fx, self.Fx, fontsize=12)

        if (self.Fx < 0):
            axes[0].text(self.x_force_x / 2 - 0.7 * self.scale, self.Fx + 0.1*self.Fx, self.Fx, fontsize=12)

        # Q-Plot hübsch
        axes[1].title.set_text('Q [KN]')
        axes[1].axis('off')
        q_max = self.get_max_Q()
        axes[1].set_ylim([-q_max - 10, q_max + 10])
        axes[1].plot([0, L], [(q_max + 10) / 9, (q_max + 10) / 9], 'k--')
        axes[1].invert_yaxis()

        # M-Plot hübsch
        axes[2].title.set_text('M [KNm]')
        axes[2].axis('off')
        m_max = self.get_max_M()
        if (self.Fy > 0):
            axes[2].set_ylim([m_max + 15, -m_max - 15])
            axes[2].plot([0, L], [(m_max + 10) / 9, (m_max + 10) / 9], 'k--')
        else:
            axes[2].set_ylim([m_max + 15, -m_max - 15])
            axes[2].plot([0, L], [(m_max + 10) / 9, (m_max + 10) / 9], 'k--')