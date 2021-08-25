import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
import os
import imageio
import natsort
import tempfile
import shutil
from scipy.stats import poisson as poisson
from scipy import interpolate


class SR_system:
    Gamma_0 = 1
    spin = 1 / 2
    m_discretization = 1

    def __init__(self, s, m0, N = 20):
        self.N = N
        self.m0 = m0
        self.s = s
        self.m_values = np.arange(-s, s + .01, self.m_discretization)

        if self.s > self.N * self.spin:
            raise ValueError(f's value must be smaller then N * spin = {self.N * self.spin}')
        if self.m0 not in self.m_values:
            raise ValueError(f'm0 = {self.m0} is not a valid quantum m number')

        self.Pm_0 = {m: 1 if m == m0 else 0 for m in self.m_values}
        self.Pm = pd.DataFrame({'t': [0], **self.Pm_0})  # ToDo make 't' key index

        self.natural_time_unit = 1 / (self.N * self.Gamma_0)  # classical emmition rate
        self.dt = .1 * (1 / (self.N * self.Gamma(0)))

    def Gamma(self, m):
        return self.Gamma_0 * (self.s + m) * (self.s - m + 1)

    def calc_next_Pm(self):
        curr_Pm = self.Pm.iloc[-1][self.m_values]
        # dPm = np.asarray([self.calc_dPm(m, curr_Pm) for m in self.Pm_0.keys()])
        dPm = curr_Pm @ self.dPm_matrix
        next_Pm = curr_Pm + dPm * self.dt
        return next_Pm

    # def calc_dPm(self, m, curr_Pm):
    #     if m + self.m_discretization <= self.s:
    #         dPm = -self.Gamma(m) * curr_Pm[m] + self.Gamma(m + self.m_discretization) * curr_Pm[
    #             m + self.m_discretization]
    #     else:
    #         dPm = -self.Gamma(m) * curr_Pm[m]
    #
    #     return dPm

    def calc_dPm_matrix(self):
        self_decay_gamma = self.Gamma(np.diag(self.m_values)) * np.eye(len(self.m_values))
        pump_gamma = np.append(self_decay_gamma[:, 1:], np.zeros((len(self.m_values), 1)), axis=1)
        return pump_gamma - self_decay_gamma

    def creat_time_evolution(self, t_sim):
        """
        continue time evolution for t_sim period.
        calculated evolution is added to object data even if the function was interrupted.
        """
        self.natural_time_unit = 1 / (self.N * self.Gamma_0)  # classical emmition rate
        self.dt = .1 * (1 / (self.N * self.Gamma(0)))
        self.dPm_matrix = self.calc_dPm_matrix()
        steps = round(t_sim / self.dt)

        for i, _ in enumerate(tqdm(range(steps))):

            if float(1 - self.Pm[self.m_values.min()].iloc[-1]) < 1e-1:
                # The system is at its lowest state
                print('no farther computation needed')
                break
            if i % 100 == 0:
                self.dt = .1 * (1 / (self.N * self.Gamma(self.get_m_expectation_value().values[-1])))

            self.Pm = self.Pm.append({'t': self.Pm['t'].iloc[-1] + self.dt,
                                      **self.calc_next_Pm()},
                                     ignore_index=True)

        self._simulating_sanity_check()

    def _simulating_sanity_check(self):
        P_sum = self.Pm[self.m_values].sum(axis=1)
        Pmin, Pmax = P_sum.max(), P_sum.min()

        if np.abs(Pmin - 1) < 1e-3 and np.abs(Pmax - 1) < 1e-1:
            print('sanity check for numerical approximation passed')
        else:
            raise ValueError('sanity check for numerical approximation didnt pass')

    def get_f_expectation_value(self, f):
        """ given functoin f(m) return its expectation value over time"""
        tmp_df = self.Pm.copy()
        for m in self.m_values:
            tmp_df[m] = f(m) * tmp_df[m]
        return tmp_df[self.m_values].sum(axis=1)

    def get_m_expectation_value(self):
        # tmp_df = self.Pm[np.asarray(list(self.Pm_0.keys()))].copy()
        # for m in np.asarray(list(self.Pm_0.keys())):
        #     tmp_df[m] = m * tmp_df[m]
        # return tmp_df.sum(axis=1)
        return self.get_f_expectation_value(f=lambda m: m)

    def get_n_from_m(self, m):
        return (self.m0 - m) / self.m_discretization

    def get_n_expectation_value(self):
        # tmp_df = self.Pm[np.asarray(list(self.Pm_0.keys()))].copy()
        # for m in np.asarray(list(self.Pm_0.keys())):
        #     tmp_df[m] = self.get_n_from_m(m) * tmp_df[m]
        # return tmp_df.sum(axis=1)
        return self.get_f_expectation_value(f=self.get_n_from_m)

    def get_n_square_expectation_value(self):
        # tmp_df = self.Pm[np.asarray(list(self.Pm_0.keys()))].copy()
        # for m in np.asarray(list(self.Pm_0.keys())):
        #     tmp_df[m] = self.get_n_from_m(m) ** 2 * tmp_df[m]
        # return tmp_df.sum(axis=1)
        return self.get_f_expectation_value(f=lambda m: self.get_n_from_m(m)**2)

    def get_mandel_Q_param(self):
        return ((self.get_n_square_expectation_value() - self.get_n_expectation_value() ** 2) / self.get_n_expectation_value()) - 1

    def get_Gamma_expectation_value(self):
        # tmp_df = self.Pm[np.asarray(list(self.Pm_0.keys()))].copy()
        # for m in np.asarray(list(self.Pm_0.keys())):
        #     tmp_df[m] = self.Gamma(m) * tmp_df[m]
        # return tmp_df.sum(axis=1)
        return self.get_f_expectation_value(f=self.Gamma)

    #     def get_emmistion_rate(self):
    #         emmited_photon = self.m0 - self.get_m_expectation_value()
    #         emmistion_rate = np.diff(emmited_photon, append=[0]) / self.dt
    #         return emmistion_rate

    def plot_m_expectation_value(self, normlize_m_range=False, *args, **kwargs):
        m_expectation_value = self.get_m_expectation_value()
        if normlize_m_range:
            m_expectation_value = m_expectation_value - m_expectation_value.min()
            m_expectation_value = m_expectation_value / m_expectation_value.max()
        plt.plot(self.Pm['t'] / self.natural_time_unit, m_expectation_value, *args, **kwargs)
        plt.xlabel('t [$(N\cdot\Gamma_0)^{-1}$]')
        plt.ylabel('$<m>$')

    def plot_scatter_p_evolution(self, *args, **kwargs):
        pixels = 240
        t, m, p = self.Pm['t'] / self.natural_time_unit, self.m_values, self.Pm[self.m_values].values.transpose()
        pp = interpolate.interp2d(t, m, p)

        im_t, im_m = np.linspace(t.min(), t.max(), pixels), m
        tt, mm = np.meshgrid(im_t, im_m)
        # plt.imshow(p, interpolation='bilinear', origin='lower', extent=(t.min(), t.max(), m.min(), m.max()), aspect="auto")
        plt.scatter(tt, mm, c=pp(im_t, im_m), *args, **kwargs)
        plt.colorbar()
        plt.xlim([t.min(), t.max()])
        self.plot_m_expectation_value()


    def plot_Gamma(self, as_functoin_of_m=False, *args, **kwargs):
        if as_functoin_of_m:
            plt.plot(self.get_m_expectation_value(), self.get_Gamma_expectation_value(), *args, **kwargs)
            plt.xlabel('<m>')
        else:
            plt.plot(self.Pm['t'] / self.natural_time_unit, self.get_Gamma_expectation_value(), *args, **kwargs)
            plt.xlabel('t [$(N\cdot\Gamma_0)^{-1}$]')
        plt.ylabel('$\Gamma_{eff} [ph/t]$')

    def plot_Q(self, as_functoin_of_m=False, *args, **kwargs):
        if as_functoin_of_m:
            plt.plot(self.get_m_expectation_value(), self.get_mandel_Q_param(), *args, **kwargs)
            plt.xlabel('<m>')
        else:
            plt.plot(self.Pm['t'] / self.natural_time_unit, self.get_mandel_Q_param(), *args, **kwargs)
            plt.xlabel('t [$(N\cdot\Gamma_0)^{-1}$]')

        plt.ylabel('$Q$')

    def stem_Pm(self, t, *args, **kwargs):
        t_idx = self.Pm.t.searchsorted(t)

        fig, ax = plt.subplots()
        ax.stem(self.Pm.iloc[t_idx][1:].keys(), self.Pm.iloc[t_idx][1:],
                label='t = ' + "{:.3f}".format(t) + '[$(N\cdot\Gamma_0)^{-1}$]', use_line_collection=True)
        ax.legend(loc='upper left')
        ax.set_xlabel('m')
        ax.set_ylabel('P(m)')
        return fig, ax

    def plot_Pn(self, t, *args, **kwargs):
        t_idx = self.Pm.t.searchsorted(t)

        plt.plot(self.get_n_from_m(self.m_values),
                 self.Pm.iloc[t_idx][self.m_values],
                label='t = ' + "{:.3f}".format(t) + '[$(N\cdot\Gamma_0)^{-1}$]')
        plt.legend(loc='upper left')
        plt.xlabel('n')
        plt.ylabel('P(n)')

    def plot_Pm(self, t, *args, **kwargs):
        t_idx = self.Pm.t.searchsorted(t)

        plt.plot(self.m_values,
                 self.Pm.iloc[t_idx][self.m_values],
                label='t = ' + "{:.3f}".format(t / self.natural_time_unit) + '[$(N\cdot\Gamma_0)^{-1}$]')
        plt.legend(loc='upper left')
        plt.xlabel('m')
        plt.ylabel('P(m)')

    def calc_equivalent_n_pois_dist(self, t, normelize_dist=True):
        """ assuming pois distribution will have the same m expectation value """

        t_idx = self.Pm.t.searchsorted(t)

        expect_m = self.get_m_expectation_value()[t_idx]
        m_vals = self.Pm.iloc[t_idx][1:].keys()

        expect_n = self.get_n_from_m(expect_m)
        n_vals = self.get_n_from_m(m_vals)

        dist = poisson.pmf(k=n_vals, mu=expect_n)
        if normelize_dist:
            dist = dist / dist.sum()

        return dist

    def stem_Pn_vs_poiss_equivalent(self, t, *args, **kwargs):
        t_idx = self.Pm.t.searchsorted(t)

        m_val, Pm = self.Pm.iloc[t_idx][1:].keys(), self.Pm.iloc[t_idx][1:]
        plt.stem(self.get_n_from_m(m_val), Pm, label='SR distribution', use_line_collection=True)
        plt.plot(np.concatenate((self.get_n_from_m(m_val), [-1e-10])),
                 np.concatenate((self.calc_equivalent_n_pois_dist(t), [.0])),
                 '--', label='poisson distribution')
        plt.legend()
        plt.xlabel('n')
        plt.ylabel('P(n)')

    def _diluted_times(self, tot_final_times):
        diluted_times = self.Pm.t
        if tot_final_times >= len(diluted_times):
            return diluted_times

        skip_size = int(np.floor(len(diluted_times) / tot_final_times))
        diluted_times = np.array([diluted_times[skip_size * i] for i in range(tot_final_times)])
        return diluted_times

    def _diluted_idx(self, tot_final_idx):
        diluted_idx = self.Pm.index
        if tot_final_idx >= len(diluted_idx):
            return diluted_idx

        skip_size = int(np.floor(len(diluted_idx) / tot_final_idx))
        diluted_idx = np.array([diluted_idx[skip_size * i] for i in range(tot_final_idx)])
        return diluted_idx

    def animate_Pm(self, dest, sime_times=None, duration=3.5):
        dirpath = tempfile.mkdtemp()

        # if sime_times is None:
        #     sime_times = self.Pm.t
        #     images = int(60 * duration)
        #     skip_group = int(np.floor(len(sime_times) / images))
        #     sime_times = [sime_times[skip_group * i] for i in range(images)]
        sime_times = self._diluted_times(tot_final_times=int(60 * duration))

        for frame, t in enumerate(tqdm(sime_times)):
            fig, ax = self.stem_Pm(t)
            if frame == 0:
                ylim = ax.get_ylim()
            ax.set_ylim(ylim)
            fig.savefig(os.path.join(dirpath, f'fig{frame}.jpg'))
            plt.close(fig)

        convert_images_from_dirc_to_gif(dirpath, duration=duration, dest=dest)

        shutil.rmtree(dirpath)


class CL_system(SR_system):
    def Gamma(self, m):
        return self.Gamma_0 * (m - (- self.s))


def convert_images_from_dirc_to_gif(im_dirc, duration, dest, format='gif', delete_after=True):
    fn_list = [os.path.join(im_dirc, f) for f in natsort.natsorted(os.listdir(im_dirc))]
    frames = [imageio.imread(f) for f in fn_list]
    imageio.mimsave(dest, frames, duration=duration / len(frames))
    if delete_after:
        [os.remove(f) for f in fn_list if not f.startswith('psi_evolution.')]


if __name__ == '__main__':
    N = 20

    SR_m_max = SR_system(N * SR_system.spin, m0=N * SR_system.spin, N=N)
    # SR_m_0 = SR_system(N * SR_system.spin, m0=0, N=N)
    # CL_m_max = CL_system(N * SR_system.spin, m0=N * SR_system.spin, N=N)
    # CL_m_0 = CL_system(N * SR_system.spin, m0=0, N=N)

    ntu = 30
    SR_m_max.creat_time_evolution(ntu * SR_m_max.natural_time_unit)
    # SR_m_0.creat_time_evolution(ntu * SR_m_0.natural_time_unit)
    # CL_m_max.creat_time_evolution(20 * ntu * SR_m_0.natural_time_unit)
    # CL_m_0.creat_time_evolution(20 * ntu * SR_m_0.natural_time_unit)

    SR_m_max.plot_m_expectation_value()
    plt.show()
