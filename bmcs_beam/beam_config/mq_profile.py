import matplotlib.patches as mpatches
import numpy as np
import sympy as sp
import traits.api as tr
from bmcs_cross_section.mkappa import MKappa
from bmcs_beam.beam_config.boundary_conditions import BoundaryConditions
from bmcs_utils.api import InteractiveModel, \
    Item, View, Float, Int, FloatEditor, FloatRangeEditor, mpl_align_yaxis
from matplotlib.patches import PathPatch
from matplotlib.path import Path
from sympy.physics.continuum_mechanics.beam import Beam


class MQPProfile(InteractiveModel):
    name = 'M-Q profile'

    mc = tr.Instance(MKappa, ())

    # @todo [SR]: this does not belong here - executable code at the
    #        the logic of boundary conditions belongs into BeamBC
    #        One possibility if to define the classes for elementary
    #        configurations, e.g.:
    #          - Beam3PSimplySupported
    #          - Beam4PSimplySupported
    #          - BeamCustomBC
    #          - Beam3PClamped
    #          - ...
    #
    bc = tr.Instance(BoundaryConditions,())
    bd = tr.DelegatesTo('bc', 'beam_design')
    L = tr.DelegatesTo('bd')
    H = tr.DelegatesTo('bd')

    _GEO = tr.Event
    @tr.on_trait_change('bc, +GEO, bc._GEO')
    def _reset_GEO(self):
        self._GEO = True

    supports_loc = tr.Property(depends_on = '_GEO')

    @tr.cached_property
    def _get_supports_loc(self):
        return  self.bc.get_supports_loc()

    # TODO [SD]: Deprecated style - don't put function calls to the class level
    #            This should be included in a function.
    x, E, I, F_ = sp.symbols('x E I F')
    l = sp.symbols('l', positive=True)  # the l sign
    b3p = Beam(l, E, I)
    R1, R2 = sp.symbols('R1 R2')
    b3p.apply_load(R1, 0, -1)
    b3p.apply_load(R2, l, -1)
    b3p.apply_load(-F_, l / 2, -1)
    b3p.bc_deflection = [(0, 0), (l, 0)]
    b3p.solve_for_reaction_loads(R1, R2)

    conf_name = b3p  # beam configuration name

    F_max = Float(1, desc='maximum load value', BC=True)
    theta_F = Float(1.0, desc='load factor', BC=True)

    F = tr.Property(depends_on='+BC')
    '''Current value of load
    '''
    @tr.cached_property
    def _get_F(self):
        return self.F_max * self.theta_F

    n_x = Int(100)
    G_adj = Float(0.015)

    ipw_view = View(
        Item('theta_F', latex=r'\theta [-]',
             editor=FloatRangeEditor(low=0, high=1)),
        Item('F_max', latex='F_\mathrm{max} [\mathrm{N}]',
             readonly=True),
        Item('G_adj', param=True, latex='G_{adj} [\mathrm{-}]', minmax=(1e-3, 1e-1)),
        Item('n_x', param=True, latex='n_x [\mathrm{-}]'),
    )

    x = tr.Property(depends_on='_GEO')

    @tr.cached_property
    def _get_x(self):
        return np.linspace(0, self.L, self.n_x)

    def get_M_x(self):
        x, F, l = sp.symbols('x F l')
        M_ = self.conf_name.bending_moment().rewrite(sp.Piecewise)
        get_M = sp.lambdify((x, F, l), M_, 'numpy')
        M_x = get_M(self.x, self.F, self.L)
        return M_x

    def get_Q_x(self):
        x, F, l = sp.symbols('x F l')
        Q_ = self.conf_name.shear_force().rewrite(sp.Piecewise)
        get_Q = sp.lambdify((x, F, l), Q_, 'numpy')
        Q_x = get_Q(self.x, self.F, self.L)
        return Q_x

    def plot_geo(self, ax1):
        # beam
        ax1.fill([0, self.L, self.L, 0, 0], [0, 0, self.H, self.H, 0], color='gray')
        #         ax.plot([0,self.L,self.L,0,0], [0,0,self.H,self.H,0],color='black')
        # supports
        vertices = []
        codes = []

        # @TODO [SR]: Check this part of implementation - this should be possible
        #        to transform to a more transparent implementation based on
        #        bc-configurations objects.
        for i in range(0, len(self.supports_loc[0])):

            codes += [Path.MOVETO] + [Path.LINETO] * 2 + [Path.CLOSEPOLY]
            vertices += [(self.supports_loc[0][i] - self.L * (self.G_adj), -self.L * self.G_adj),
                         (self.supports_loc[0][i], 0),
                         (self.supports_loc[0][i] + self.L * (self.G_adj), -self.L * self.G_adj),
                         (self.supports_loc[0][i] - self.L * (self.G_adj), -self.L * self.G_adj)]

            if self.supports_loc[0][i] != 0:
                codes += [Path.MOVETO] + [Path.LINETO] + [Path.CLOSEPOLY]
                vertices += [(self.supports_loc[0][i] + (self.L * self.G_adj), -self.L * self.G_adj * 1.2),
                             (self.supports_loc[0][i] - (self.L * self.G_adj), -self.L * self.G_adj * 1.2),
                             (self.supports_loc[0][i] - (self.L * self.G_adj), -self.L * self.G_adj * 1.2)]

        for i in range(0, len((self.supports_loc[1]))):
            codes += [Path.MOVETO] + [Path.LINETO] * 3 + [Path.CLOSEPOLY]
            vertices += [(self.supports_loc[1][i] - self.L * (self.G_adj), -self.L * self.G_adj * 1.2),
                         (self.supports_loc[1][i] - self.L * (self.G_adj), self.H + self.L * self.G_adj * 1.2),
                         (self.supports_loc[1][i] + self.L * (self.G_adj), self.H + self.L * self.G_adj * 1.2),
                         (self.supports_loc[1][i] + self.L * (self.G_adj), -self.L * self.G_adj * 1.2),
                         (self.supports_loc[1][i] - self.L * (self.G_adj), -self.L * self.G_adj * 1.2)]

        vertices = np.array(vertices, float)
        path = Path(vertices, codes)
        pathpatch = PathPatch(path, facecolor='black', edgecolor='black')
        ax1.add_patch(pathpatch)

        # loads
        load = self.conf_name.applied_loads
        value = [0]
        pos = [0]
        type = [0]

        # reading the value and the positoin of external loads
        for i in range(0, len(load)):
            if load[i][0] == -self.F or load[i][0] == self.F:
                value[0] = (load[i][0])
                pos[0] = (load[i][1])
                load_dic = dict(zip(value, pos))
                load_ = sp.lambdify((self.F, self.l), load_dic)
                Load_x = load_(self.f, self.L)

                # load arrow parameters
                x_tail = (list(Load_x.values())[0])
                x_head = (list(Load_x.values())[0])

                # moment
                if load[i][2] == -2:

                    if load[i][0] == -self.F:

                        ax1.plot([(list(Load_x.values())[0])], [(self.H / 2)],
                                 marker=r'$\circlearrowleft$', ms=self.H / 5)

                    else:

                        ax1.plot([(list(Load_x.values())[0])], [(self.H / 2)],
                                 marker=r'$\circlearrowright$', ms=self.H / 5)

                    ax1.annotate('{} KN.mm'.format(np.round(self.F / 1000), 0),
                                 xy=(list(Load_x.values())[0], self.H * 1.4), color='blue')

                # point load
                elif load[i][2] == -1:

                    if load[i][0] == -self.F:

                        y_tail = self.L / 20 + self.H
                        y_head = 0 + self.H
                        dy = y_head + x_head / 10

                        arrow = mpatches.FancyArrowPatch((x_tail, y_tail), (x_head, y_head),
                                                         color='blue', mutation_scale=self.L / 500)
                        ax1.annotate('{} KN'.format(np.round(self.f / 1000), 0),
                                     xy=(x_tail, y_tail), color='black')
                        ax1.add_patch(arrow)

                    else:

                        y_tail = 0 + self.H
                        y_head = self.L / 20 + self.H
                        dy = y_head + x_head / 10

                        arrow = mpatches.FancyArrowPatch((x_tail, y_tail), (x_head, y_head),
                                                         color='blue', mutation_scale=self.L / 500)
                        ax1.annotate('{} KN'.format(np.round(self.f / 1000), 0),
                                     xy=(x_head, y_head), color='black')
                        ax1.add_patch(arrow)

                else:

                    if load[i][0] == -self.F:
                        y_tail = self.L / 20 + self.H
                        y_head = 0 + self.H
                        dy = y_head + x_head / 10

                        # distributed load
                        l_step = 0
                        while l_step <= self.L:
                            x_tail = (list(Load_x.values())[0]) + l_step
                            y_tail = self.L / 20 + self.H
                            x_head = (list(Load_x.values())[0]) + l_step
                            y_head = 0 + self.H
                            dy = y_head + x_head / 10
                            l_step += self.L / 10
                            arrow = mpatches.FancyArrowPatch((x_tail, y_tail), (x_head, y_head),
                                                             color='blue', mutation_scale=self.L / 500)
                            ax1.add_patch(arrow)

                        ax1.annotate('{} KN'.format(np.round(self.F / 1000), 0),
                                     xy=(self.L / 2, y_tail * 1.1), color='black')
                        ax1.plot([0, self.L], [y_tail, y_tail], color='blue')

                    else:
                        y_tail = 0 + self.H
                        y_head = self.L / 20 + self.H
                        dy = y_head + x_head / 10

                        # distributed load
                        l_step = 0
                        while l_step <= self.L:
                            x_tail = (list(Load_x.values())[0]) + l_step
                            y_tail = 0 + self.H
                            x_head = (list(Load_x.values())[0]) + l_step
                            y_head = self.L / 20 + self.H
                            dy = y_head + x_head / 10
                            l_step += self.L / 10
                            arrow = mpatches.FancyArrowPatch((x_tail, y_tail), (x_head, y_head),
                                                             color='blue', mutation_scale=self.L / 500)
                            ax1.add_patch(arrow)

                        ax1.annotate('{} KN'.format(np.round(self.f / 1000), 0),
                                     xy=(self.L / 2, y_head * 1.1), color='black')
                        ax1.plot([0, self.L], [y_head, y_head], color='blue')

    def plot_MQ(self, ax2, ax3):
        x = self.x

        M_scale = self.mc.M_scale
        M_x = self.get_M_x() / M_scale
        ax2.plot(x, -M_x, color='red', label='moment [kNm]')
        ax2.fill(x, -M_x, color='red', alpha=0.1)
        ax2.set_ylabel('M [kNm]')
        ax2.legend();

        Q_x = self.get_Q_x() / 1000
        ax3.plot(x, Q_x, lw=0.1, color='green', label='shear [kN]')
        ax3.fill_between(x, Q_x, 0, color='green', alpha=0.1)
        ax3.set_ylabel('Q [kN]')
        ax3.set_xlabel('x [mm]')
        ax3.legend();
        mpl_align_yaxis(ax2, 0, ax3, 0)

    def subplots(self, fig):
        ax1, ax2 = fig.subplots(2, 1)
        ax3 = ax2.twinx()
        return ax1, ax2, ax3

    def update_plot(self, axes):
        ax1, ax2, ax3 = axes
        self.plot_geo(ax1)
        self.plot_MQ(ax2, ax3)
        ax1.axis('equal')
        ax1.autoscale(tight=True)
