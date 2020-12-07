
import traits.api as tr
import numpy as np
from bmcs_utils.api import InteractiveModel, View, Item, Button, ButtonEditor, Float, Int, mpl_align_yaxis
from bmcs_beam.beam_config.beam_design import BeamDesign
from bmcs_cross_section.mkappa import MKappa
from scipy.integrate import cumtrapz
import matplotlib.gridspec as gridspec


class DeflectionProfile(InteractiveModel):
    '''
    Deflection model of a BMCS beam
    '''

    name = 'Deflection Profile'

    beam_design = tr.Instance(BeamDesign, ())
    mc = tr.Instance(MKappa, ())
    n_load_steps = Int(31)

    ipw_view = View(
        Item('n_load_steps', latex='n_{\mathrm{load~steps}}')
    )

    def get_kappa_x(self):
        '''
        Profile of curvature along the beam
        '''
        M = self.beam_design.get_M_x()
        return self.mc.get_kappa_M(M)

    def get_phi_x(self):
        '''
        Calculate the cross sectional rotation by integrating the curvature
        '''
        kappa_x = self.get_kappa_x()
        phi_x = cumtrapz(kappa_x, self.beam_design.x, initial=0)
        # resolve the integration constant by requiring zero curvature
        # at the midspan of the beam
        # TODO [SD] this is specific to 3 point bending - generalize
        #           for other loading conditions.
        phi_L2 = np.interp(self.beam_design.L / 2, self.beam_design.x, phi_x)

        if self.beam_design.beam_conf_name == '3pb':
            phi_x -= phi_L2
        elif self.beam_design.beam_conf_name == '4pb':
            phi_x -= phi_L2
        return phi_x

    def get_w_x(self):
        '''
        Profile of deflection along the beam
        '''
        phi_x = self.get_phi_x()
        w_x = cumtrapz(phi_x, self.beam_design.x, initial=0)
        # resolve the integration constant by requiring zero deflection
        # at the left support - the right one comes automatically
        # TODO [SR, HS] this is specific to 3 point bending - generalize
        #           for other loading conditions.
        if self.beam_design.beam_conf_name == '3pb':
            w_x += w_x[0]
        elif self.beam_design.beam_conf_name == '4pb':
            w_x += w_x[0]
        return w_x

    theta_max = tr.Float(1)

    F_max = tr.Property(Float)
    ''''
    Identify the ultimate limit state based on the maximum moment capacity
    of the cross section.
    '''
    def _get_F_max(self):
        # specific to 3pt bending  - equation should be provided by the
        # BoundaryCondition class - corresponding to the load configuration.
        # TODO [SR, HS] this is specific to 3 point bending - generalize
        #           for other loading conditions.
        M_I, kappa_I = self.mc.inv_M_kappa

        if self.beam_design.beam_conf_name == '3pb':
            F_Max_ = 4 * M_I[-1] / self.beam_design.L
    
        elif self.beam_design.beam_conf_name == '4pb':
            F_Max_ = 4 * M_I[-1] / self.beam_design.L

        return F_Max_

    # def run(self):
    #     F_arr = np.linspace(0, self.F_max, self.n_load_steps)
    #     w_list = []
    #     original_F = self.beam_design.F
    #     for F in F_arr:
    #         if F == 0:
    #             w_list.append(0)
    #         else:
    #             self.beam_design.F = -F
    #             # Append the maximum deflection value that corresponds to the new load (F)
    #             w_list.append(np.fabs(np.min(self.get_w_x())))
    #     self.beam_design.F = original_F
    #     return F_arr, np.array(w_list)
    #
    # def reset(self):
    #     self.theta_F = 0

    F_max_old = Float

    def get_Fw(self):
        F_max = self.F_max
        F_arr = np.linspace(0, F_max, self.n_load_steps)
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
        original_F = self.beam_design.F
        for F in F_arr:
            if F == 0:
                w_list.append(0)
            else:
                self.beam_design.F = -F
                # Append the maximum deflection value that corresponds to the new load (F)
                w_list.append(np.fabs(np.min(self.get_w_x())))
        if self.F_max_old == F_max:
            self.beam_design.F = original_F
        self.F_max_old = F_max
        return F_arr, np.array(w_list)

    def subplots(self, fig):
        gs = gridspec.GridSpec(1, 2, figure=fig, width_ratios=[0.7, 0.3])

        # ax2, ax3 = fig.subplots(1, 2)
        ax_w = fig.add_subplot(gs[0, 0])
        ax_k = ax_w.twinx()
        ax_Fw = fig.add_subplot(gs[0, 1])
        return ax_w, ax_k, ax_Fw

    def update_plot(self, axes):
        ax_w, ax_k, ax_Fw = axes

        x = self.beam_design.x

        F_scale = 1000

        # TODO: expensive calculations for all displacements are running with each plot update to produce new
        #  load-displacement curve, this shouldn't be done for example when only the force has changed
        ax_Fw.set_xlabel(r'$w_\mathrm{max}$ [mm]')
        ax_Fw.set_ylabel(r'$F$ [kN]')
        F, w = self.get_Fw()
        ax_Fw.plot(w, F / F_scale, color='blue', lw=2)
        current_F = round(abs(self.beam_design.F / F_scale), 2)
        ax_Fw.axhline(y=current_F, color='r')
        ax_Fw.annotate('F = {} kN'.format(current_F), xy=(0, current_F + 3), color='r')

        kappa_x = self.get_kappa_x()  # self.mc.get_kappa(M)
        ax_k.plot(x, -kappa_x, color='black', label='$kappa$ [-]')
        ax_k.fill(x, -kappa_x, color='gray', alpha=0.1)
        ax_k.set_ylabel(r'$\kappa [\mathrm{mm}^{-1}]$')
        ax_k.set_xlabel(r'$x$')
        ax_k.legend()

        w_x = self.get_w_x()
        ax_w.plot(x, w_x, color='blue', label='$w$ [mm]')
        ax_w.fill(x, w_x, color='blue', alpha=0.1)
        ax_w.set_ylabel(r'$w [\mathrm{mm}]$')
        ax_w.legend(loc='lower right')
        mpl_align_yaxis(ax_w, 0, ax_k, 0)
