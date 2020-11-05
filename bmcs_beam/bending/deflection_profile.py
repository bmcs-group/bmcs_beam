import matplotlib.patches as mpatches
import numpy as np
import sympy as sp
from bmcs_beam.beam_bc.mq_profile import MQPProfile
from matplotlib.patches import PathPatch
from matplotlib.path import Path
from scipy.integrate import cumtrapz

class DeflectionProfile(MQPProfile):
    name = 'DeflectionProfile'


    '''Temporary beam definition'''

    def get_kappa_x(self):
        M = self.get_M_x()
        return self.mc.get_kappa_M(M)

    def get_phi_x(self):
        kappa_x = self.get_kappa_x()
        phi_x = cumtrapz(kappa_x, self.x, initial=0)
        phi_L2 = np.interp(self.L / 2, self.x, phi_x)
        phi_x -= phi_L2
        return phi_x

    def get_w_x(self):
        phi_x = self.get_phi_x()
        w_x = cumtrapz(phi_x, self.x, initial=0)
        w_x += w_x[0]
        return w_x

    def set_F_max(self):
        ''''Identify the ultimate limit state based on the maximum moment capacity
        of the cross section.
        '''
        # specific to 3pt bending  - equation should be provided by the
        # BoundaryCondition class - corresponding to the load configuration.
        #
        M_I, kappa_I = self.mc.inv_M_kappa
        self.F_max = 4 * M_I[-1] / self.L

    def subplots(self, fig):
        ax1, ax2, ax3 = fig.subplots(3, 1)
        ax22 = ax2.twinx()
        ax33 = ax3.twinx()
        return ax1, ax2, ax22, ax3, ax33

    def update_plot(self, axes):
        ax1, ax2, ax22, ax3, ax33 = axes

        self.plot_geo(ax1)
        self.plot_MQ(ax2, ax22)

        x = self.x
        kappa_x = self.get_kappa_x()    #self.mc.get_kappa(M)
        ax3.plot(x, kappa_x, color='black', label='$kappa$ [-]')
        ax3.set_ylabel(r'$\kappa [\mathrm{mm}^{-1}]$')
        ax3.set_xlabel(r'$x$')
        ax3.legend();

        w_x = self.get_w_x()
        ax33.plot(x, w_x, color='blue', label='$w$ [mm]')
        ax33.set_ylabel(r'$w [\mathrm{mm}]$')
        ax33.legend();