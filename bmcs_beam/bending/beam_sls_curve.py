import bmcs_utils.api as bu
from .deflection_profile import DeflectionProfile
import numpy as np
import matplotlib.pyplot as plt
import pickle
import os
import json
from matplotlib.ticker import PercentFormatter
import traits.api as tr

class BeamSLSCurve(bu.ParametricStudy):
    name = 'BeamSystem SLS Curve'

    '''
    - link to data cache identifying the directory where to store the
    interim results
    - dictionary based storage of the individual runs of the study
      which makes it introspectable.
    '''

    sls_to_uls_ratio = bu.Float(0.59)
    l_over_d_min = bu.Int(10)
    l_over_d_max = bu.Int(35)
    rho_min = bu.Float(0.0025)
    rho_max = bu.Float(0.025)

    dp = bu.Instance(DeflectionProfile, ())
    tree = ['dp']

    ipw_view = bu.View(
        bu.Item('l_over_d_min', latex='{l/d}_\mathrm{min}'),
        bu.Item('l_over_d_max', latex='{l/d}_\mathrm{max}'),
        bu.Item('rho_min', latex=r'\rho_\mathrm{min}'),
        bu.Item('rho_max', latex=r'\rho_\mathrm{max}'),
        bu.Item('n_i', latex='n_i'),
        # bu.Item('slenderness_range', latex='l/d', editor=bu.IntRangeEditor(value=(10, 35), low_name='slenderness_low', high_name='slenderness_high', n_steps_name='')),
        bu.Item('sls_to_uls_ratio', latex=r'F_\mathrm{SLS}/F_\mathrm{ULS}'),
        time_editor=bu.ProgressEditor(
            run_method='run',
            reset_method='reset',
            interrupt_var='interrupt',
            time_var='seg',
            time_max='n_seg'
        )
    )

    n_i = bu.Int(3)
    dense_quarter = bu.Bool

    rho_range = tr.Property(depends_on='rho_max, rho_max, n_i')

    @tr.cached_property
    def _get_rho_range(self):
        if self.dense_quarter:
            range = np.linspace(self.rho_min, self.rho_max/4, self.n_i/2)
            np.concatenate(range, np.linspace(self.rho_max/4, self.rho_max, self.n_i/2))
            return np.linspace(self.rho_min, self.rho_max, self.n_i)
        else:
            return np.linspace(self.rho_min, self.rho_max, self.n_i)

    l_over_d_range = tr.Property(depends_on='l_over_d_min, l_over_d_max, n_i')

    @tr.cached_property
    def _get_l_over_d_range(self):
        return np.linspace(self.l_over_d_min, self.l_over_d_max, self.n_i)

    n_seg = bu.Int(5, TIME=True)
    seg = bu.Int(0)
    interrupt = bu.Bool(False)

    def run(self, update_progress=lambda t: t):
        F_u_grid, F_s_grid, rho_grid, sl_grid = self.get_Fu_and_Fs()
        # self.plot_with_ec2_curves(F_u_grid, F_s_grid, rho_grid, sl_grid)

        # while self.seg <= self.n_seg:
        #     if self.interrupt:
        #         break
        #     try:
        #         self.make_incr()
        #     except StopIteration:
        #         print('stopped after', self.seg)
        #         break
        #     if self.seg < self.n_seg:
        #         self.sz_cp.add_x_tip_an(self.sz_cp.sz_ctr.x_tip_ak[:, 0])
        #     self.seg += 1

    def reset(self):
        # reset everything to default values
        pass

    def subplots(self, fig):
        self.ax1 = fig.subplots(1, 1)
        return self.ax1

    def update_plot(self, axes):
        ax = axes
        self.plot_with_ec2_curves(F_u_grid, F_s_grid, rho_grid, sl_grid, ax=ax)
        # if self.should_plot:
        #     ax.plot(self.w_data, self.F_data / 1000)

    def save(self):
        out_dir = os.path.join(bu.data_cache.dir, self.__class__.__name__)
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        with open(os.path.join(out_dir, 'class'), 'wb') as out_file:
            """ Giving TypeError: cannot pickle 'weakref' object"""
            pickle.dump(self, out_file)

    def save_data_vars_in_json(self, f_ck, dp, path):
        mc = dp.mc
        rein = mc.cross_section_layout.items
        output_data = {'mc.n_m': mc.n_m,
                       'mc.n_kappa': mc.n_kappa,
                       'mc.low_kappa': mc.low_kappa,
                       'mc.high_kappa': mc.high_kappa,
                       'mc.E_cc': mc.cs_design.matrix_.E_cc,
                       'mc.E_ct': mc.cs_design.matrix_.E_ct,
                       'mc.eps_tu': mc.cs_design.matrix_.eps_tu,
                       'mc.eps_cr': mc.cs_design.matrix_.eps_cr,
                       'mc.eps_cy': mc.cs_design.matrix_.eps_cy,
                       'mc.eps_cu': mc.cs_design.matrix_.eps_cu,
                       'mc.mu': mc.cs_design.matrix_.mu,
                       'f_ck': f_ck,
                       'rein[0].E': rein[0].matmod_.E,
                       'rein[0].z': rein[0].z,
                       'rein[0].A': rein[0].A,
                       'dp.beam_design.L': dp.beam_design.L, }
        with open(path, 'w') as outfile:
            json.dump(output_data, outfile, sort_keys=True, indent=4)

    def save_resulting_files(self, f_ck, dp, F_u_grid, F_s_grid, rho_grid, sl_grid):
        np.save('F_u_grid_steel_EC2_eq2_tension_c' + str(f_ck) + '.npy', F_u_grid)
        np.save('F_s_grid_steel_EC2_eq2_tension_c' + str(f_ck) + '.npy', F_s_grid)
        np.save('rho_grid_steel_EC2_eq2_tension_c' + str(f_ck) + '.npy', rho_grid)
        np.save('sl_grid_steel_EC2_eq2_tension_c' + str(f_ck) + '.npy', sl_grid)
        self.save_data_vars_in_json(f_ck, dp, 'data_steel_EC2_eq2_tension_c' + str(f_ck) + '.json')

    should_plot = bu.Bool(False)
    # F_data_array = tr.List()
    # w_data_array = tr.List()

    def get_Fu_and_Fs(self, upper_reinforcement=False):
        # self.F_data_array = tr.List()
        # self.w_data_array = tr.List()
        dp = self.dp

        # def get_Fu_and_Fs(dp, rho_range=np.linspace(0.0025, 0.025), l_over_d_range=np.linspace(10, 35), upper_reinforcement = False):
        if upper_reinforcement:
            d = dp.mc.cross_section_layout.items[0].z
        else:
            d = dp.mc.cross_section_shape_.H - dp.mc.cross_section_layout.items[0].z
        b = dp.mc.cross_section_shape_.B
        area_g = b * d

        rho_grid, sl_grid = np.meshgrid(self.rho_range, self.l_over_d_range)
        F_u_grid = np.zeros_like(rho_grid)
        F_s_grid = np.zeros_like(rho_grid)

        _, ax = plt.subplots()
        ax.set_xlabel(r'$w$ [mm]')
        ax.set_ylabel(r'$F$ [KN]')

        for sl_idx in range(0, len(self.l_over_d_range)):
            for rho_idx in range(0, len(self.rho_range)):
                rho = rho_grid[rho_idx, sl_idx]
                sl = sl_grid[rho_idx, sl_idx]

                print('parameter combination', rho, sl)

                # assigning the grid area (area_g) to the reinforcement area variable
                A_j_g = rho * area_g
                dp.mc.cross_section_layout.items[0].A = A_j_g

                # assigning the grid length (L_g) to the beam length variable
                L_g = sl * d
                dp.beam_design.L = L_g

                dp.mc.state_changed = True

                # running the deflection analysis
                F_data, w_data = dp.get_Fw()
                self.F_data_array.append(F_data)
                self.w_data_array.append(w_data)
                # self.should_plot = True

                # plotting, post-processing & saving the data
                ax.plot(w_data, F_data / 1000, label="rho={}%-sl={} ".format(rho * 100, sl))

                w_s = dp.beam_design.L / 250
                F_u = max(F_data)
                F_s = np.interp(w_s, w_data, F_data, right=F_u * 2)

                F_u_grid[rho_idx, sl_idx] = F_u
                F_s_grid[rho_idx, sl_idx] = F_s

        return F_u_grid, F_s_grid, rho_grid, sl_grid

    def plot_with_ec2_curves(self, F_u_grid, F_s_grid, rho_grid, sl_grid, ax=None):
        if not ax:
            fig, ax = plt.subplots()

        z = F_u_grid / F_s_grid - 1. / 0.59
        CS = ax.contour(rho_grid, sl_grid, z)
        ax.clabel(CS, inline=1, fontsize=10)

        # Draw EC2 curve
        self.plot_steel_sls_curves(ax, f_cks=[70], axes_start_from_zero=True)

    def plot_steel_sls_curves(self, ax=None, rho_range=None, f_cks=None, rho_p=0, K=1, axes_start_from_zero=False):
        if not ax:
            fig, ax = plt.subplots()
        if not rho_range:
            rho_range = np.linspace(0.0025, 0.025, 1000)
        if not f_cks:
            f_cks = np.arange(20, 110, 10)

        slenderness = []
        for f_ck in f_cks:
            for rho in rho_range:
                slenderness.append(BeamSLSCurve.get_l_over_d_limit(rho, f_ck, rho_p, K))
            ax.plot(rho_range, slenderness, label=r'$f_{ck} = $' + str(f_ck) + ' MPa')
            slenderness = []

        if axes_start_from_zero:
            ax.set_ylim(0, 35)
            ax.set_xlim(0, rho_range[-1])
        else:
            ax.set_ylim(10, 35)
            ax.set_xlim(rho_range[0], rho_range[-1])

        ax.xaxis.set_major_formatter(PercentFormatter(xmax=1))
        ax.set_ylabel(r'$l/(K \cdot d$)')
        ax.set_xlabel(r'Tensile reinforcement ratio $\rho$ [%]')
        ax.set_title(r'$l/d$ limits according to EC2 eqs. 7.16a & 7.16b')
        ax.grid(color='#e6e6e6', linewidth=0.7)
        ax.legend()

    @staticmethod
    def get_l_over_d_limit(rho, f_ck, rho_p=0, K=1):
        # see EC2, see eqs 7.16
        rho_0 = 0.001 * np.sqrt(f_ck)
        if rho <= rho_0:
            return K * (11 + 1.5 * ((f_ck) ** 0.5) * (rho_0 / rho) + 3.2 * ((f_ck) ** 0.5) * (
                        (rho_0 / rho - 1) ** (3 / 2)))
        else:
            return K * (11 + 1.5 * ((f_ck) ** 0.5) * (rho_0 / (rho - rho_p)) + (1 / 12) * (f_ck ** 0.5) * (
                        (rho_p / rho_0) ** 0.5))

