
import traits.api as tr
import numpy as np
from bmcs_utils.api import Float, mpl_align_yaxis
from bmcs_beam.beam_bc.mq_profile import MQPProfile
from scipy.integrate import cumtrapz
import matplotlib.gridspec as gridspec

class DeflectionProfile(MQPProfile):
    '''
    Deflection model of a BMCS beam
    '''

    name = 'Deflection profile'

    def get_kappa_x(self):
        '''
        Profile of curvature along the beam
        '''
        M = self.get_M_x()
        return self.mc.get_kappa_M(M)

    def get_phi_x(self):
        '''
        Calculate the cross sectional rotation by integrating the curvature
        '''
        kappa_x = self.get_kappa_x()
        phi_x = cumtrapz(kappa_x, self.x, initial=0)
        # resolve the integration constant by requiring zero curvature
        # at the midspan of the beam
        # TODO [SD] this is specific to 3 point bending - generalize
        #           for other loading conditions.
        phi_L2 = np.interp(self.L / 2, self.x, phi_x)
        phi_x -= phi_L2
        return phi_x

    def get_w_x(self):
        '''
        Profile of deflection along the beam
        '''
        phi_x = self.get_phi_x()
        w_x = cumtrapz(phi_x, self.x, initial=0)
        # resolve the integration constant by requiring zero deflection
        # at the left support - the right one comes automatically
        # TODO [SD] this is specific to 3 point bending - generalize
        #           for other loading conditions.
        w_x += w_x[0]
        return w_x

    F_max = tr.Property(Float)
    ''''
    Identify the ultimate limit state based on the maximum moment capacity
    of the cross section.
    '''
    def _get_F_max(self):
        # specific to 3pt bending  - equation should be provided by the
        # BoundaryCondition class - corresponding to the load configuration.
        #
        M_I, kappa_I = self.mc.inv_M_kappa
        return 4 * M_I[-1] / self.L

    def get_Fw(self):
        theta_arr = np.linspace(0,1,41)
        F_arr = self.F_max * theta_arr
        w_list = []
        # @todo [SR,RC]: separate the slider theta_F from the calculation
        #                of the datapoints load deflection curve.
        #                use broadcasting in the functions
        #                get_M_x(x[:,np.newaxis], F[np.newaxis,:] and
        #                in get_Q_x, get_kappa_x, get_w_x, get_phi_x, get_w_x
        #                then, the browsing through the history is done within
        #                the two dimensional array of and now loop over theta is
        #                neeeded then. Theta works just as a slider - as originally
        #                introduced.
        for theta in theta_arr:
            self.theta_F = theta
            w_list.append(np.fabs(np.min(self.get_w_x())))
        return F_arr, np.array(w_list)

    def subplots(self, fig):

        gs = gridspec.GridSpec(2, 2, figure=fig, width_ratios=[0.7,0.3])

        ax_M = fig.add_subplot(gs[0, 0])
        ax_Q = ax_M.twinx()
        # ax2, ax3 = fig.subplots(2, 1)
        ax_w = fig.add_subplot(gs[1, 0])
        ax_k = ax_w.twinx()
        ax_Fw = fig.add_subplot(gs[:,1])
        return ax_M, ax_Q, ax_w, ax_k, ax_Fw

    def update_plot(self, axes):
        ax_M, ax_Q, ax_w, ax_k, ax_Fw = axes

        self.plot_MQ(ax_M, ax_Q)

        x = self.x
        kappa_x = self.get_kappa_x()  # self.mc.get_kappa(M)
        ax_k.plot(x, -kappa_x, color='black', label='$kappa$ [-]')
        ax_k.fill(x, -kappa_x, color='gray', alpha=0.1)
        ax_k.set_ylabel(r'$\kappa [\mathrm{mm}^{-1}]$')
        ax_k.set_xlabel(r'$x$')
        ax_k.legend();

        w_x = self.get_w_x()
        ax_w.plot(x, w_x, color='blue', label='$w$ [mm]')
        ax_w.fill(x, w_x, color='blue', alpha=0.1)
        ax_w.set_ylabel(r'$w [\mathrm{mm}]$')
        ax_w.legend();
        mpl_align_yaxis(ax_w, 0, ax_k, 0)

        ax_Fw.set_xlabel(r'$w_\mathrm{max}$ [mm]')
        ax_Fw.set_ylabel(r'$F$ [kN]')
        F, w = self.get_Fw()
        ax_Fw.plot(w, F / 1000, color='blue', lw=2)
