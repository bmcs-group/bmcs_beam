import pickle
import json
import os
import pickle

import bmcs_utils.api as bu
import matplotlib.pyplot as plt
import numpy as np
import traits.api as tr
from bmcs_beam.api import BoundaryConfig, DeflectionProfile
from bmcs_cross_section.api import MKappa, EC2, ReinfLayer
from matplotlib.ticker import PercentFormatter


class BeamSLSCurve(bu.Model):
    name = 'Beam SLS Curve'

    '''
    - link to data cache identifying the directory where to store the
    interim results
    - dictionary based storage of the individual runs of the study
      which makes it introspectable.
    '''

    sls_to_uls_ratio = bu.Float(0.59)
    slenderness_min = bu.Int(3)
    slenderness_max = bu.Int(50)
    rho_min = bu.Float(0.0002)
    rho_max = bu.Float(0.025)

    dp = bu.Instance(DeflectionProfile)

    def _dp_default(self):
        return self._get_dp()

    tree = ['dp']

    n_i = bu.Int(3)
    dense_quarter = bu.Bool

    rho_range = tr.Property(depends_on='rho_max, rho_max, n_i')

    @tr.cached_property
    def _get_rho_range(self):
        if self.dense_quarter:
            range = np.linspace(self.rho_min, self.rho_max / 4, self.n_i / 2)
            np.concatenate(range, np.linspace(self.rho_max / 4, self.rho_max, self.n_i / 2))
            return np.linspace(self.rho_min, self.rho_max, self.n_i)
        else:
            return np.linspace(self.rho_min, self.rho_max, self.n_i)

    slenderness_range = tr.Property(depends_on='slenderness_min, slenderness_max, n_i')

    @tr.cached_property
    def _get_slenderness_range(self):
        return np.linspace(self.slenderness_min, self.slenderness_max, self.n_i)

    # pseudo time vars
    interrupt = bu.Bool(False)
    n_seg = bu.Int(5, TIME=True)
    seg = bu.Int(0)

    # dp properties:
    apply_material_factors = bu.Bool(False)
    system_type = bu.Str('dist') # system_type can be '4pb' or '3pb' or 'dist'
    concrete_law = bu.Str('EC2') # or 'piecewise linear' or 'EC2 with plateau'
    f_ck = bu.Float(50)
    rein_type = bu.Str('steel')
    use_f_ctm_fl = bu.Bool(True)

    @tr.observe('apply_material_factors')
    @tr.observe('system_type')
    @tr.observe('concrete_law')
    @tr.observe('f_ck')
    @tr.observe('rein_type')
    @tr.observe('use_f_ctm_fl')
    def _update_dp(self, event):
        self.dp = self._get_dp()

    ipw_view = bu.View(
        bu.Item('slenderness_min', latex='{l/d}_\mathrm{min}'),
        bu.Item('slenderness_max', latex='{l/d}_\mathrm{max}'),
        bu.Item('rho_min', latex=r'\rho_\mathrm{min}'),
        bu.Item('rho_max', latex=r'\rho_\mathrm{max}'),
        bu.Item('n_i', latex='n_i'),
        bu.Item('apply_material_factors', latex='\mathrm{Apply}~\gamma'),
        bu.Item('system_type', latex='\mathrm{System}'),
        bu.Item('concrete_law', latex='\mathrm{concrete~law}'),
        bu.Item('f_ck', latex='f_{ck}'),
        bu.Item('rein_type', latex='\mathrm{Rein~type}'),
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

    def _get_dp(self):
        b = 1000
        h = 300
        d = 0.9 * h
        f_ck = self.f_ck
        rein_type = self.rein_type

        E = EC2.get_E_cm(f_ck)
        f_ctm_fl = EC2.get_f_ctm_fl(f_ck, h)
        f_ctm = EC2.get_f_ctm(f_ck)

        eps_cr = f_ctm_fl / E if self.use_f_ctm_fl else f_ctm / E

        # Info: eps_cr =  0.000170 with default concrete law gives good slenderness curve
        mc = MKappa(low_kappa=0, high_kappa=0.00007, n_kappa=100)
        mc.cs_design.matrix = self.concrete_law
        mc.cs_design.matrix_.trait_set(
                                       factor=0.85 / 1.5 if self.apply_material_factors else 1,
                                       eps_cr=eps_cr,
                                       eps_tu=eps_cr,
        )

        if self.concrete_law == 'EC2 with plateau' or self.concrete_law == 'EC2':
            mc.cs_design.matrix_.trait_set(f_cm=EC2.get_f_cm(f_ck))
        elif self.concrete_law == 'piecewise linear':
            mc.cs_design.matrix_.trait_set(
                E_cc = E,
                E_ct = E,
                eps_cy=EC2.get_eps_c3(f_ck),
                eps_cu=EC2.get_eps_cu3(f_ck),
            )

        # The default uses f_ctm. Here, I will use f_ctm_fl (in EC2, they tested with both)
        # mc.cs_design.matrix_.eps_cr = (EC2.get_f_ctm_fl(f_ck, h) * 1.5 / 0.85) / mc.cs_design.matrix_.E_ct
        # mc.cs_design.matrix_.eps_tu = (EC2.get_f_ctm_fl(f_ck, h) * 1.5 / 0.85) / mc.cs_design.matrix_.E_ct

        mc.cross_section_shape_.B = b
        mc.cross_section_shape_.H = h

        # T-section
        # mc.cross_section_shape = 'I-shape'
        # mc.cross_section_shape_.H = 200
        # mc.cross_section_shape_.B_w = 50
        # mc.cross_section_shape_.B_f_bot = 50
        # mc.cross_section_shape_.B_f_top = 150
        # mc.cross_section_shape_.H_f_bot = 50
        # mc.cross_section_shape_.H_f_top = 50

        rho = 0.01
        A_s = rho * b * d

        bl1 = ReinfLayer(name=rein_type, z=h - d, A=A_s, matmod=rein_type)
        if rein_type == 'steel':
            bl1.matmod_.trait_set(E_s=200000, f_sy=500, factor=1 / 1.15 if self.apply_material_factors else 1)
        elif rein_type == 'carbon':
            # carbon material factors :
            # alpha_ft * alpha_f_eff / gamma_frp (see El-Ghadioui2020_PhD P. 122)
            bl1.matmod_.trait_set(E=230000, f_t=3300, factor=0.85 * 0.9 / 1.3 if self.apply_material_factors else 1)
        mc.cross_section_layout.add_layer(bl1)

        dp = DeflectionProfile(mc=mc)

        if self.system_type == '4pb':
            dp.beam_design.beam_conf_name = BoundaryConfig.FOUR_PB
            dp.beam_design.beam_conf_name.first_load_distance = L / 3
        elif self.system_type == '3pb':
            dp.beam_design.beam_conf_name = BoundaryConfig.THREE_PB
        elif self.system_type == 'dist':
            dp.beam_design.beam_conf_name = BoundaryConfig.SIMPLE_BEAM_DIST_LOAD
            dp.F_scale = 1

        dp.beam_design.L = 8 * d
        return dp

    def run(self, update_progress=lambda t: t):
        print('run started...')
        F_u_grid, F_s_grid, rho_grid, sl_grid = self.get_Fu_and_Fs()
        self._plot_with_ec2_curves(F_u_grid, F_s_grid, rho_grid, sl_grid, self.ax1)
        print('run finished...')

    def reset(self):
        # reset everything to default values
        pass

    def subplots(self, fig):
        self.ax1 = fig.subplots(1, 1)
        return self.ax1

    def update_plot(self, axes):
        print('update_plot')

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
    F_data_array = tr.List()
    w_data_array = tr.List()

    def get_Fu_and_Fs(self, upper_reinforcement=False, plot=False):
        self.F_data_array = []
        self.w_data_array = []

        slenderness_range = self.slenderness_range
        rho_range = self.rho_range

        dp = self.dp

        if upper_reinforcement:
            d = dp.mc.cross_section_layout.items[0].z
        else:
            d = dp.mc.cross_section_shape_.H - dp.mc.cross_section_layout.items[0].z

        area_g = dp.mc.cross_section_shape_.get_cs_area()

        rho_grid, sl_grid = np.meshgrid(rho_range, slenderness_range)
        F_u_grid = np.zeros_like(rho_grid)
        F_s_grid = np.zeros_like(rho_grid)

        if plot:
            _, ax = plt.subplots()
            ax.set_xlabel(r'$w$ [mm]')
            ax.set_ylabel(r'$F$ [KN]')

        for sl_idx in range(len(slenderness_range)):
            for rho_idx in range(len(rho_range)):
                if self.interrupt:
                    return

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

                # plotting, post-processing & saving the data
                if plot:
                    ax.plot(w_data, F_data / 1000, label="rho={}%-sl={} ".format(rho * 100, sl))

                w_s = dp.beam_design.L / 250
                F_u = max(F_data)
                F_s = np.interp(w_s, w_data, F_data, right=F_u * 2)

                F_u_grid[rho_idx, sl_idx] = F_u
                F_s_grid[rho_idx, sl_idx] = F_s

        return F_u_grid, F_s_grid, rho_grid, sl_grid

    def _plot_with_ec2_curves(self, F_u_grid, F_s_grid, rho_grid, sl_grid, ax=None, label=None):
        if not ax:
            fig, ax = plt.subplots()

        color = np.random.rand(3,)

        z = F_u_grid / F_s_grid - 1. / self.sls_to_uls_ratio
        CS = ax.contour(rho_grid, sl_grid, z, levels=[0], colors=[color])
        if label:
            CS.collections[0].set_label(label)

        # Draw EC2 curve
        self.plot_steel_sls_curves(ax, f_cks=[self.f_ck],
                                   axes_start_from_zero=True,
                                   color=color)

    def plot_steel_sls_curves(self, ax=None,
                              rho_range=np.linspace(0.0025, 0.025, 1000),
                              f_cks=np.arange(20, 110, 10),
                              rho_p=0,
                              K=1,
                              axes_start_from_zero=False,
                              color=None):
        if not ax:
            fig, ax = plt.subplots()

        slenderness = []
        for f_ck in f_cks:
            for rho in rho_range:
                slenderness.append(BeamSLSCurve.get_slenderness_limit(rho, f_ck, rho_p, K))
            ax.plot(rho_range, slenderness, label=r'$f_{ck} = $' + str(f_ck) + ' MPa (Steel - EC2)', ls='--', c=color)

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
    def get_slenderness_limit(rho, f_ck, rho_p=0, K=1):
        # see EC2, see eqs 7.16
        rho_0 = 0.001 * np.sqrt(f_ck)
        if rho <= rho_0:
            return K * (11 + 1.5 * ((f_ck) ** 0.5) * (rho_0 / rho) + 3.2 * ((f_ck) ** 0.5) * (
                    (rho_0 / rho - 1) ** (3 / 2)))
        else:
            return K * (11 + 1.5 * ((f_ck) ** 0.5) * (rho_0 / (rho - rho_p)) + (1 / 12) * (f_ck ** 0.5) * (
                    (rho_p / rho_0) ** 0.5))


class SLSParamStudy(bu.ParametricStudy):
    def __init__(self, b_sls):
        self.b_sls = b_sls

    def plot(self, ax, plot_title, param_name_value):
        ax.set_xlabel(r'$w_\mathrm{max}$ [mm]')
        ax.set_ylabel(r'$F$ [kN]')

        self.b_sls.dp.mc.state_changed = True

        F_u_grid, F_s_grid, rho_grid, sl_grid = self.b_sls.get_Fu_and_Fs()
        self.b_sls._plot_with_ec2_curves(F_u_grid, F_s_grid, rho_grid, sl_grid, ax, label=param_name_value)
        ax.set_title(plot_title + ' effect')
        ax.legend()

        out_dir = self.get_output_dir()
        np.save(os.path.join(out_dir, 'F_u_grid' + '_' + self.b_sls.rein_type + '__' + param_name_value + '.npy'), F_u_grid)
        np.save(os.path.join(out_dir, 'F_s_grid' + '_' + self.b_sls.rein_type + '__' + param_name_value + '.npy'), F_s_grid)
        np.save(os.path.join(out_dir, 'rho_grid' + '_' + self.b_sls.rein_type + '__' + param_name_value + '.npy'), rho_grid)
        np.save(os.path.join(out_dir, 'sl_grid' + '_' + self.b_sls.rein_type + '__' + param_name_value + '.npy'), sl_grid)

    def plot_all_curves(self):
        f_cks = [50]
        F_u_grids = []
        F_s_grids = []
        rho_grids = []
        sl_grids = []
        reinforcement = 'carbon'
        for f_ck in f_cks:
            f_ck = str(f_ck)
            F_u_grids.append(
                np.load('exported_data/F_u_grid_carbon_EC2_eq2_tension_E230_ft_3000_c' + str(f_ck) + '.npy'))
            F_s_grids.append(
                np.load('exported_data/F_s_grid_carbon_EC2_eq2_tension_E230_ft_3000_c' + str(f_ck) + '.npy'))
            rho_grids.append(
                np.load('exported_data/rho_grid_carbon_EC2_eq2_tension_E230_ft_3000_c' + str(f_ck) + '.npy'))
            sl_grids.append(np.load('exported_data/sl_grid_carbon_EC2_eq2_tension_E230_ft_3000_c' + str(f_ck) + '.npy'))

        _, ax = plt.subplots(1, 1)

        ax.set_ylabel('L/d')
        ax.set_xlabel(r'$\rho$ %')
        ax.set_ylim(0, 35)
        ax.set_xlim(0.0, 0.025)

        for f_ck, F_u_grid, F_s_grid, rho_grid, sl_grid in zip(f_cks, F_u_grids, F_s_grids, rho_grids, sl_grids):
            z = 0.5 * F_u_grid / F_s_grid - 1. / 0.59
            CS = ax.contour(rho_grid, sl_grid, z, colors=[np.random.rand(3, )], levels=[0])
            CS.collections[0].set_label('C' + str(f_ck))
            #     ax.clabel(CS, inline=1, fontsize=10)

            BeamSLSCurve().plot_steel_sls_curves(ax=ax, f_cks=[50], axes_start_from_zero=True)

        ax.legend()

