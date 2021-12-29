import argparse
import sys
import os
import numpy as np
import xlsxwriter
import openpyxl
import spd_analyze
import pandas as pd
import module_fwhm
#import generate_pdf
from matplotlib import cm
from scipy.fft import fft, ifft, fftfreq

import matplotlib.pyplot as plt
plt.rcParams["figure.figsize"] = (20,3)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams.update({'font.size': 22})
from mpl_toolkits.axes_grid1 import make_axes_locatable
#plt.rcParams["figure.figsize"] = (1, 1)

import warnings
warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning)


class graphics_generation(object):
    def __init__(self, params):
        self.params = params
        self.dir_name = params['dir_name']
        self.dof_name = params['name']
        self.fold = 'report'
        self.params['fold'] = self.fold

        #self.criterion_list = ['SNR', 'QE']
        self.criterion_list = [self.params["criterion"]]

        try:
            os.mkdir(os.path.join(self.dir_name, self.fold))
        except:
            print(os.path.join(self.dir_name, self.fold) + ' folder currently exists')

        for el in self.criterion_list:
            try:
                os.mkdir(os.path.join(self.dir_name, self.fold, el))
            except:
                print(os.path.join(self.dir_name, self.fold, el) + ' folder currently exists')


        #self.SPD = spd_analyze.SPD(self.dir_name)
        #self.SPD = self.SPD.get_spd()
        self.SPD = self.params["SPD"]

        self.temp_list = self.SPD.temp_list
        #Interesting way to determine the qe_sat_list set
        qe_max = 0
        for i in range(len(self.SPD.qe)):
            if np.max(self.SPD.qe[i]) > qe_max:
                qe_max = np.max(self.SPD.qe[i])

        #Settings for common detectors
        #self.qe_sat_list = [10, 15, 20, 25]
        #Settings for gun detector
        #self.qe_sat_list = [20, 25, 30, 35]

        #print('qe_max = ', qe_max)
        # if 18 <= qe_max < 70:
        #     self.qe_sat_list = [10, 15, 20, 25]
        # elif 10 <= qe_max < 18:
        # #    self.qe_sat_list = [round((i + 1) * qe_max / 4, 1) for i in range(4)]
        #     self.qe_sat_list = [5, 7, 10, 15]
        # elif qe_max > 70 :
        #     self.qe_sat_list = [20, 25, 30, 35]

        self.qe_sat_list = [10, 20, 30]

        self.colors = ['b', 'g', 'r', 'magenta', 'k', 'cyan']
        self.params['optimal_params'] = {'head': ['Criterion', 'PDE, %', 'DCR, Hz', 'SNR', 'AP, %', 'DT, mus', 'TR, ps']}
        self.params['settings'] = {'head' : ['Criteria', 'Vg, code', 'Vb, V', 'T, C']}
        self.params['temp'] = self.temp_list
        self.params['type'] = self.SPD.type
        self.params["SPD"] = self.SPD

        self.optimum_snr_tr_hist = {}
        self.optimum_qe_10_tr_hist = {}
        self.optimum_qe_20_tr_hist = {}
        self.optimum_qe_30_tr_hist = {}
        self.optimum_snr_tr_hist_bif = {}
        self.optimum_qe_10_tr_hist_bif = {}
        self.optimum_qe_20_tr_hist_bif = {}
        self.optimum_qe_30_tr_hist_bif = {}
        self.vbias_trhist = {}
        self.vbias_trhist_10 = {}
        self.vbias_trhist_20 = {}
        self.vbias_trhist_30 = {}
        self.vbias_trhist_bif = {}
        self.vbias_trhist_10_bif = {}
        self.vbias_trhist_20_bif = {}
        self.vbias_trhist_30_bif = {}

        #self.criterion = params['criterion']

        if self.SPD.type == 'gated':
            self.params['optimal_params'] = {'head': ['Criterion', 'PDE, %', 'DCR, Hz', 'SNR', 'AP, %', 'DT, mus', 'TR, ps', 'BIF, %']}
            self.stb_list = self.SPD.stb_list
            #if type(self.gate_voltage) != int:
            #    self.params['settings'] = {'head': ['Criterion', 'Vg, V', 'Vb, V', 'T, C']}

    def get_sub_slice(self, qe_sat):
        # First part about PDE
        dcr_sat, ap_sat, vb_sat = [], [], []
        if self.SPD.type == 'gated':
            for i in range(len(self.temp_list)):
                dcr_sat.append([])
                ap_sat.append([])
                vb_sat.append([])
                for j in range(len(self.SPD.stb_list[i])):
                    qe_list = np.array(self.SPD.qe[i][j,:])
                    qe_list[qe_list < 0] = 0
                    ind_vb2 = np.argmin(abs(qe_list * (qe_list >= qe_sat) - qe_sat))
                    ind_vb1 = np.argmin(abs(qe_list * (qe_list < qe_sat) - qe_sat))
                    qe1 = qe_list[ind_vb1]
                    qe2 = qe_list[ind_vb2]

                    if ind_vb1 != ind_vb2:
                        dcr1 = self.SPD.dcr[i][j,ind_vb1]
                        dcr2 = self.SPD.dcr[i][j,ind_vb2]
                        dcr_sat[i].append((dcr2 - dcr1) / (qe2 - qe1) * (qe_sat - qe1) + dcr1)

                        ap1 = self.SPD.ap_list[i][j,ind_vb1]
                        ap2 = self.SPD.ap_list[i][j,ind_vb2]
                        ap_new_val = (ap2 - ap1) / (qe2 - qe1) * (qe_sat - qe1) + ap1
                        if ap_new_val > 50:
                            ap_new_val = 50
                        if ap_new_val < 0:
                            ap_new_val = 0
                        ap_sat[i].append(ap_new_val)

                        vb1 = self.SPD.grid[i][j,ind_vb1]
                        vb2 = self.SPD.grid[i][j,ind_vb2]
                        vb_sat[i].append((vb2 - vb1) / (qe2 - qe1) * (qe_sat - qe1) + vb1)
                    else:
                        qe1 = 0
                        dcr1 = 0
                        dcr2 = self.SPD.dcr[i][j, ind_vb2]
                        dcr_sat[i].append((dcr2 - dcr1) / (qe2 - qe1) * (qe_sat - qe1) + dcr1)

                        ap1 = 0
                        ap2 = self.SPD.ap_list[i][j, ind_vb2]
                        ap_sat[i].append((ap2 - ap1) / (qe2 - qe1) * (qe_sat - qe1) + ap1)

                        vb1 = self.SPD.grid[i][j, ind_vb1] - 1.5
                        vb2 = self.SPD.grid[i][j, ind_vb2]
                        vb_sat[i].append((vb2 - vb1) / (qe2 - qe1) * (qe_sat - qe1) + vb1)

            return np.array(dcr_sat),  np.array(ap_sat), np.array(vb_sat)
        elif self.SPD.type == 'freerun':
            for i in range(len(self.temp_list)):
                qe_list = np.array(self.SPD.qe[i, :])
                ind_vb2 = np.argmin(abs(qe_list * (qe_list >= qe_sat) - qe_sat))
                ind_vb1 = np.argmin(abs(qe_list * (qe_list < qe_sat) - qe_sat))
                qe1 = qe_list[ind_vb1]
                qe2 = qe_list[ind_vb2]

                dcr1 = self.SPD.dcr[i, ind_vb1]
                dcr2 = self.SPD.dcr[i, ind_vb2]
                dcr_sat.append((dcr2 - dcr1) / (qe2 - qe1) * (qe_sat - qe1) + dcr1)

                ap1 = self.SPD.ap_list[i, ind_vb1]
                ap2 = self.SPD.ap_list[i, ind_vb2]
                ap_sat.append((ap2 - ap1) / (qe2 - qe1) * (qe_sat - qe1) + ap1)

                vb1 = self.SPD.grid[i, ind_vb1]
                vb2 = self.SPD.grid[i, ind_vb2]
                vb_sat.append((vb2 - vb1) / (qe2 - qe1) * (qe_sat - qe1) + vb1)
            return np.array(dcr_sat), np.array(ap_sat), np.array(vb_sat)

    def get_slices(self):
        dcr_sat, ap_sat, vb_sat  = [], [], []
        for i in range(len(self.qe_sat_list)):
            dcr_sat_hlp, ap_sat_hlp, vb_sat_hlp  = self.get_sub_slice(self.qe_sat_list[i])
            dcr_sat.append(dcr_sat_hlp)
            ap_sat.append(ap_sat_hlp)
            vb_sat.append(vb_sat_hlp)

        self.dcr_sat = np.array(dcr_sat)
        self.ap_sat = np.array(ap_sat)
        self.vb_sat = np.array(vb_sat)

    def plot_dcr_qe_fr(self, str_spd=""):
        for j in range(len(self.temp_list)):
            fig, ax1 = plt.subplots()
            fig.set_size_inches(8, 6)
            ax1.plot(self.SPD.grid[j], self.SPD.dcr[j], color=self.colors[0], marker='o', label='DCR')
            ax2 = ax1.twinx()
            ax2.plot(self.SPD.grid[j], self.SPD.qe[j], color=self.colors[1], marker='o', label='PDE')

            ax1.set_xlabel(r'$V_b$, V')
            ax1.set_ylabel(r'$DCR$, Hz')
            ax2.set_ylabel(r'$PDE$, \%')
            plt.title(r'DCR and PDE dependence on $V_b$. T = ' + str(self.temp_list[j]) + ' C')
            ax1.legend(fontsize=12, loc='upper left')
            ax2.legend(fontsize=12, loc='upper right')
            ax1.grid(True)
            plt.savefig(os.path.join(self.params['dir_name'], self.fold, self.criterion) + '//DCR_T' + str(self.temp_list[j]) + str_spd + '.png',
                        bbox_inches='tight', transparent=True, dpi=100)
            plt.cla()

    def plot_dcr(self, str_spd=""):
        for j in range(len(self.temp_list)):
            fig, ax1 = plt.subplots()
            fig.set_size_inches(8, 6)

            if type(self.gate_voltage) == int:
                for i in range(len(self.qe_sat_list)):
                    ax1.plot(self.stb_list[j], self.dcr_sat[i][j], color=self.colors[i], marker='o',
                             label='PDE = ' + str(self.qe_sat_list[i]) + '$\%$')
                ax1.set_xlabel(r'$V_g$, code')
                filename = os.path.join(self.params['dir_name'], self.fold, self.criterion) + '//DCR_T' + str(self.temp_list[j]) + '_CD' + str_spd + '.png'
            else:
                for i in range(len(self.qe_sat_list)):
                    ax1.plot(self.gate_voltage.vg.to_numpy()[:len(self.stb_list[j])] / 2, self.dcr_sat[i][j], color=self.colors[i], marker='o',
                             label='PDE = ' + str(self.qe_sat_list[i]) + '$\%$')
                ax1.set_xlabel(r'$V_g$, V')
                filename = os.path.join(self.params['dir_name'], self.fold, self.criterion) + '//DCR_T' + str(
                    self.temp_list[j]) + '_VG' + str_spd + '.png'

            ax1.set_ylabel('DCR, Hz')
            plt.title(r'DCR dependence on $V_g$. T = ' + str(self.temp_list[j]) + ' C')
            ax1.legend(fontsize=14, loc='upper right')
            plt.grid(True)
            plt.savefig(filename, bbox_inches='tight', transparent=True, dpi=100)
            plt.cla()

    def plot_ap_fr(self, str_spd=""):
        for j in range(len(self.temp_list)):
            fig, ax1 = plt.subplots()
            fig.set_size_inches(8, 6)
            ax1.plot(self.SPD.grid[j], self.SPD.ap_list[j], color=self.colors[0], marker='o', label='AP')

            ax1.set_xlabel(r'$V_b$')
            ax1.set_ylabel('AP, \%')
            plt.title(r'AP dependence on $V_b$. T = ' + str(self.temp_list[j]) + ' C')
            ax1.legend(fontsize=14, loc='upper right')
            plt.grid(True)
            plt.savefig(os.path.join(self.params['dir_name'], self.fold, self.criterion) + '//AP_T' + str(self.temp_list[j]) + str_spd + '.png',
                        bbox_inches='tight', transparent=True, dpi=100)
            plt.cla()
            # plt.show()

    def plot_ap(self, str_spd=""):
        for j in range(len(self.temp_list)):
            fig, ax1 = plt.subplots()
            fig.set_size_inches(8, 6)
            if type(self.gate_voltage) == int:
                for i in range(len(self.qe_sat_list)):
                    ax1.plot(self.stb_list[j], self.ap_sat[i][j], color=self.colors[i], marker='o',
                             label='PDE = ' + str(self.qe_sat_list[i]) + '$\%$')
                    ax1.set_xlabel(r'$V_g$, code')
                    filename = os.path.join(self.params['dir_name'], self.fold, self.criterion) + '//AP_T' + str(self.temp_list[j]) + '_CD' + str_spd + '.png'
            else:
                for i in range(len(self.qe_sat_list)):
                    ax1.plot(self.gate_voltage.vg.to_numpy()[:len(self.stb_list[j])] / 2, self.ap_sat[i][j], color=self.colors[i], marker='o',
                             label='PDE = ' + str(self.qe_sat_list[i]) + '$\%$')
                    ax1.set_xlabel(r'$V_g$, V')
                    filename = os.path.join(self.params['dir_name'], self.fold, self.criterion) + '//AP_T' + str(
                        self.temp_list[j]) + '_VG' + str_spd + '.png'

            ax1.set_ylabel('AP, \%')
            plt.title(r'AP dependence on $V_g$. T = ' + str(self.temp_list[j]) + ' C')
            ax1.legend(fontsize=14, loc='upper right')
            plt.grid(True)
            plt.savefig(filename, bbox_inches='tight', transparent=True, dpi=100)
            plt.cla()
            # plt.show()

    def plot_snr_fr(self, str_spd=""):
        for j in range(len(self.temp_list)):
            fig, ax1 = plt.subplots()
            fig.set_size_inches(8, 6)
            ax1.plot(self.SPD.grid[j], self.SPD.snr[j], color=self.colors[0], marker='o', label='AP')

            ax1.set_xlabel(r'$V_b$')
            ax1.set_ylabel('SNR')
            plt.title(r'SNR dependence on $V_b$. T = ' + str(self.temp_list[j]) + ' C')
            ax1.legend(fontsize=14, loc='upper right')
            plt.grid(True)
            plt.savefig(os.path.join(self.params['dir_name'], self.fold, self.criterion) + '//SNR_T' + str(self.temp_list[j]) + str_spd + '.png',
                        bbox_inches='tight', transparent=True, dpi=100)
            plt.cla()
            # plt.show()

    def choose_optimum_fr(self):
        if self.criterion == 'QE':
            for i in range(len(self.temp_list)):
                qe_sat_10 = self.qe_sat_list[0]
                qe_sat_20 = self.qe_sat_list[1]
                qe_sat_30 = self.qe_sat_list[2]

                qe_list = self.SPD.qe[i]
                ind_vb2_10 = np.argmin(abs(qe_list * (qe_list >= qe_sat_10) - qe_sat_10))
                ind_vb1_10 = np.argmin(abs(qe_list * (qe_list < qe_sat_10) - qe_sat_10))
                ind_vb2_20 = np.argmin(abs(qe_list * (qe_list >= qe_sat_20) - qe_sat_20))
                ind_vb1_20 = np.argmin(abs(qe_list * (qe_list < qe_sat_20) - qe_sat_20))
                ind_vb2_30 = np.argmin(abs(qe_list * (qe_list >= qe_sat_30) - qe_sat_30))
                ind_vb1_30 = np.argmin(abs(qe_list * (qe_list < qe_sat_30) - qe_sat_30))

                qe1_10 = qe_list[ind_vb1_10]
                qe2_10 = qe_list[ind_vb2_10]
                qe1_20 = qe_list[ind_vb1_20]
                qe2_20 = qe_list[ind_vb2_20]
                qe1_30 = qe_list[ind_vb1_30]
                qe2_30 = qe_list[ind_vb2_30]

                #DCR block
                dcr1_10 = self.SPD.dcr[i, ind_vb1_10]
                dcr2_10 = self.SPD.dcr[i, ind_vb2_10]
                dcr1_20 = self.SPD.dcr[i, ind_vb1_20]
                dcr2_20 = self.SPD.dcr[i, ind_vb2_20]
                dcr1_30 = self.SPD.dcr[i, ind_vb1_30]
                dcr2_30 = self.SPD.dcr[i, ind_vb2_30]

                dcr_sat_10 = (dcr2_10 - dcr1_10) / (qe2_10 - qe1_10) * (qe_sat_10 - qe1_10) + dcr1_10
                dcr_sat_20 = (dcr2_20 - dcr1_20) / (qe2_20 - qe1_20) * (qe_sat_20 - qe1_20) + dcr1_20
                dcr_sat_30 = (dcr2_30 - dcr1_30) / (qe2_30 - qe1_30) * (qe_sat_30 - qe1_30) + dcr1_30

                #AP block
                ap1_10 = self.SPD.ap_list[i, ind_vb1_10]
                ap2_10 = self.SPD.ap_list[i, ind_vb2_10]
                ap1_20 = self.SPD.ap_list[i, ind_vb1_20]
                ap2_20 = self.SPD.ap_list[i, ind_vb2_20]
                ap1_30 = self.SPD.ap_list[i, ind_vb1_30]
                ap2_30 = self.SPD.ap_list[i, ind_vb2_30]

                ap_sat_10 = (ap2_10 - ap1_10) / (qe2_10 - qe1_10) * (qe_sat_10 - qe1_10) + ap1_10
                ap_sat_20 = (ap2_20 - ap1_20) / (qe2_20 - qe1_20) * (qe_sat_20 - qe1_20) + ap1_20
                ap_sat_30 = (ap2_30 - ap1_30) / (qe2_30 - qe1_30) * (qe_sat_30 - qe1_30) + ap1_30

                #VB block
                vb1_10 = self.SPD.grid[i, ind_vb1_10]
                vb2_10 = self.SPD.grid[i, ind_vb2_10]
                vb1_20 = self.SPD.grid[i, ind_vb1_20]
                vb2_20 = self.SPD.grid[i, ind_vb2_20]
                vb1_30 = self.SPD.grid[i, ind_vb1_30]
                vb2_30 = self.SPD.grid[i, ind_vb2_30]

                vb_sat_10 = (vb2_10 - vb1_10) / (qe2_10 - qe1_10) * (qe_sat_10 - qe1_10) + vb1_10
                vb_sat_20 = (vb2_20 - vb1_20) / (qe2_20 - qe1_20) * (qe_sat_20 - qe1_20) + vb1_20
                vb_sat_30 = (vb2_30 - vb1_30) / (qe2_30 - qe1_30) * (qe_sat_30 - qe1_30) + vb1_30

                #SNR block
                snr_sat_10 = qe_sat_10 / dcr_sat_10
                snr_sat_20 = qe_sat_20 / dcr_sat_20
                snr_sat_30 = qe_sat_30 / dcr_sat_30

                #DT block
                dt_sat_10 = self.SPD.dt_list[i, 0]
                dt_sat_20 = self.SPD.dt_list[i, 0]
                dt_sat_30 = self.SPD.dt_list[i, 0]

                if abs(qe1_10 - qe_sat_10) < abs(qe2_10 - qe_sat_10):
                    index_optimum_qe_10 = ind_vb1_10
                else:
                    index_optimum_qe_10 = ind_vb2_10

                if abs(qe1_20 - qe_sat_20) < abs(qe2_20 - qe_sat_20):
                    index_optimum_qe_20 = ind_vb1_20
                else:
                    index_optimum_qe_20 = ind_vb2_20

                if abs(qe1_30 - qe_sat_30) < abs(qe2_30 - qe_sat_30):
                    index_optimum_qe_30 = ind_vb1_30
                else:
                    index_optimum_qe_30 = ind_vb2_30

                try:
                    self.optimum_qe_10_tr_hist[self.temp_list[i]] = np.array(self.SPD.timeres_dict[i][0][self.SPD.grid[i, index_optimum_qe_10]])
                    self.vbias_trhist_10[self.temp_list[i]] = self.SPD.grid[i, index_optimum_qe_10]
                except:
                    self.optimum_qe_10_tr_hist[self.temp_list[i]] = np.array(self.SPD.timeres_dict[i][0][list(self.SPD.timeres_dict[i][0].keys())[0]])
                    self.vbias_trhist_10[self.temp_list[i]] = list(self.SPD.timeres_dict[i][0].keys())[0]

                try:
                    self.optimum_qe_20_tr_hist[self.temp_list[i]] = np.array(self.SPD.timeres_dict[i][0][self.SPD.grid[i, index_optimum_qe_20]])
                    self.vbias_trhist_20[self.temp_list[i]] = self.SPD.grid[i, index_optimum_qe_20]
                except:
                    self.optimum_qe_20_tr_hist[self.temp_list[i]] = np.array(self.SPD.timeres_dict[i][0][list(self.SPD.timeres_dict[i][0].keys())[0]])
                    self.vbias_trhist_20[self.temp_list[i]] = list(self.SPD.timeres_dict[i][0].keys())[0]

                try:
                    self.optimum_qe_30_tr_hist[self.temp_list[i]] = np.array(self.SPD.timeres_dict[i][0][self.SPD.grid[i, index_optimum_qe_30]])
                    self.vbias_trhist_30[self.temp_list[i]] = self.SPD.grid[i, index_optimum_qe_30]
                except:
                    self.optimum_qe_30_tr_hist[self.temp_list[i]] = np.array(self.SPD.timeres_dict[i][0][list(self.SPD.timeres_dict[i][0].keys())[0]])
                    self.vbias_trhist_30[self.temp_list[i]] = list(self.SPD.timeres_dict[i][0].keys())[0]

                #print(self.optimum_snr_tr_hist[self.temp_list[i]])
                #Use filtration
                cnt = self.optimum_qe_10_tr_hist[self.temp_list[i]]
                T = 10e-12
                N = int(len(cnt))
                freq = fftfreq(N, T)[:N // 2]
                fft_cnt = fft(cnt)

                ind_last = np.where(freq < 1.0e10)[-1][-1]
                # print(ind_last)
                fft_cnt[ind_last::] = np.zeros(N - ind_last)
                self.optimum_qe_10_tr_hist[self.temp_list[i]] = np.real(ifft(fft_cnt))
                self.optimum_qe_10_tr = 10 * module_fwhm.fwtm(self.optimum_qe_10_tr_hist[self.temp_list[i]])

                # Use filtration
                cnt = self.optimum_qe_20_tr_hist[self.temp_list[i]]
                T = 10e-12
                N = int(len(cnt))
                freq = fftfreq(N, T)[:N // 2]
                fft_cnt = fft(cnt)

                ind_last = np.where(freq < 1.0e10)[-1][-1]
                # print(ind_last)
                fft_cnt[ind_last::] = np.zeros(N - ind_last)
                self.optimum_qe_20_tr_hist[self.temp_list[i]] = np.real(ifft(fft_cnt))
                self.optimum_qe_20_tr = 10 * module_fwhm.fwtm(self.optimum_qe_20_tr_hist[self.temp_list[i]])

                # Use filtration
                cnt = self.optimum_qe_30_tr_hist[self.temp_list[i]]
                T = 10e-12
                N = int(len(cnt))
                freq = fftfreq(N, T)[:N // 2]
                fft_cnt = fft(cnt)

                ind_last = np.where(freq < 1.0e10)[-1][-1]
                # print(ind_last)
                fft_cnt[ind_last::] = np.zeros(N - ind_last)
                self.optimum_qe_30_tr_hist[self.temp_list[i]] = np.real(ifft(fft_cnt))
                self.optimum_qe_30_tr = 10 * module_fwhm.fwtm(self.optimum_qe_30_tr_hist[self.temp_list[i]])



                self.params['optimal_params'][self.temp_list[i]] = [['QE = ' + str(round(qe_sat_10, 2)), str(round(qe_sat_10, 2)),
                                                      str(round(dcr_sat_10)), str(round(snr_sat_10, 4)),
                                                      str(round(ap_sat_10, 2)), str(round(dt_sat_10, 2)), str(round(self.optimum_qe_10_tr, 0))],
                                                                    ['QE = ' + str(round(qe_sat_20, 2)),
                                                                     str(round(qe_sat_20, 2)),
                                                                     str(round(dcr_sat_20)), str(round(snr_sat_20, 4)),
                                                                     str(round(ap_sat_20, 2)), str(round(dt_sat_20, 2)),
                                                                     str(round(self.optimum_qe_20_tr, 0))],
                                                                    ['QE = ' + str(round(qe_sat_30, 2)),
                                                                     str(round(qe_sat_30, 2)),
                                                                     str(round(dcr_sat_30)), str(round(snr_sat_30, 4)),
                                                                     str(round(ap_sat_30, 2)), str(round(dt_sat_30, 2)),
                                                                     str(round(self.optimum_qe_30_tr, 0))]
                                                                    ]

                self.params['settings'][self.temp_list[i]] = [['QE = ' + str(round(qe_sat_10, 2)),
                                                              '0', str(round(vb_sat_10, 2)),
                                                              str(round(self.temp_list[i], 1))],
                                                              ['QE = ' + str(round(qe_sat_20, 2)),
                                                               '0', str(round(vb_sat_20, 2)),
                                                               str(round(self.temp_list[i], 1))],
                                                              ['QE = ' + str(round(qe_sat_30, 2)),
                                                               '0', str(round(vb_sat_30, 2)),
                                                               str(round(self.temp_list[i], 1))]
                                                              ]


        if self.criterion == 'SNR':
            for i in range(len(self.temp_list)):
                self.optimum_snr_snr = np.max(self.SPD.snr[i])
                index_optimum_snr = np.where(self.SPD.snr[i] == self.optimum_snr_snr)
                self.index_optimum_snr = [i, index_optimum_snr[0][0]]
                self.optimum_snr_pde = self.SPD.qe[self.index_optimum_snr[0], self.index_optimum_snr[1]]
                self.optimum_snr_dcr = self.SPD.dcr[self.index_optimum_snr[0], self.index_optimum_snr[1]]
                self.optimum_snr_ap = self.SPD.ap_list[self.index_optimum_snr[0], self.index_optimum_snr[1]]
                self.optimum_snr_dt = self.SPD.dt_list[self.index_optimum_snr[0], self.index_optimum_snr[1]]

                self.optimum_snr_vb = self.SPD.grid[self.index_optimum_snr[0], self.index_optimum_snr[1]]
                self.optimum_snr_t = self.SPD.temp_list[self.index_optimum_snr[0]]

                try:
                    self.optimum_snr_tr_hist[self.temp_list[i]] = np.array(self.SPD.timeres_dict[self.index_optimum_snr[0]][0][self.optimum_snr_vb])
                    self.vbias_trhist[self.temp_list[i]] = self.optimum_snr_vb
                except:
                    self.optimum_snr_tr_hist[self.temp_list[i]] = np.array(self.SPD.timeres_dict[self.index_optimum_snr[0]][0][list(self.SPD.timeres_dict[i][0].keys())[0]])
                    self.vbias_trhist[self.temp_list[i]] = list(self.SPD.timeres_dict[i][0].keys())[0]
                #print(self.optimum_snr_tr_hist[self.temp_list[i]])
                #Use filtration
                cnt = self.optimum_snr_tr_hist[self.temp_list[i]]
                T = 10e-12
                N = int(len(cnt))
                freq = fftfreq(N, T)[:N // 2]
                fft_cnt = fft(cnt)

                ind_last = np.where(freq < 1.0e10)[-1][-1]
                # print(ind_last)
                fft_cnt[ind_last::] = np.zeros(N - ind_last)
                self.optimum_snr_tr_hist[self.temp_list[i]] = np.real(ifft(fft_cnt))

                # plt.plot(freq, 2.0 / N * np.abs(fft_cnt[0:N // 2]))
                # plt.plot(cnt)
                # plt.plot(cnt_new)
                # plt.show()

                self.optimum_snr_tr = 10 * module_fwhm.fwtm(self.optimum_snr_tr_hist[self.temp_list[i]])

                self.params['optimal_params'][self.temp_list[i]] = ['SNR', str(round(self.optimum_snr_pde, 2)),
                                                      str(round(self.optimum_snr_dcr)), str(round(self.optimum_snr_snr, 4)),
                                                      str(round(self.optimum_snr_ap, 2)), str(round(self.optimum_snr_dt, 2)), str(round(self.optimum_snr_tr, 0))]

                self.params['settings'][self.temp_list[i]] = ['SNR', '0', str(round(self.optimum_snr_vb, 2)), str(round(self.optimum_snr_t, 1))]

        #print(self.params['optimal_params'])
        #print(self.params['settings'])

    def TR_filtration(self, cnt):
        # Use filtration
        T = 10e-12
        N = int(len(cnt))
        freq = fftfreq(N, T)[:N // 2]
        fft_cnt = fft(cnt)

        ind_last = np.where(freq < 1.0e10)[-1][-1]
        # print(ind_last)
        fft_cnt[ind_last::] = np.zeros(N - ind_last)
        return np.real(ifft(fft_cnt))

    def choose_optimum(self):
        if self.criterion == 'QE':
            for i in range(len(self.temp_list)):
                min_10 = 1e6
                min_20 = 1e6
                min_30 = 1e6

                index_optimum_qe_10 = np.argmin(self.dcr_sat[0][i])
                index_optimum_qe_20 = np.argmin(self.dcr_sat[1][i])
                index_optimum_qe_30 = np.argmin(self.dcr_sat[2][i])

                self.index_optimum_qe_10 = [0, i, index_optimum_qe_10]
                self.index_optimum_qe_20 = [1, i, index_optimum_qe_20]
                self.index_optimum_qe_30 = [2, i, index_optimum_qe_30]

                # (self.index_optimum_snr)
                self.optimum_qe_10_dcr = self.dcr_sat[
                    self.index_optimum_qe_10[0]][self.index_optimum_qe_10[1]][self.index_optimum_qe_10[2]]
                self.optimum_qe_20_dcr = self.dcr_sat[
                    self.index_optimum_qe_20[0]][self.index_optimum_qe_20[1]][ self.index_optimum_qe_20[2]]
                self.optimum_qe_30_dcr = self.dcr_sat[
                    self.index_optimum_qe_30[0]][self.index_optimum_qe_30[1]][self.index_optimum_qe_30[2]]

                self.optimum_qe_10_snr = self.qe_sat_list[0] / self.optimum_qe_10_dcr
                self.optimum_qe_20_snr = self.qe_sat_list[1] / self.optimum_qe_20_dcr
                self.optimum_qe_30_snr = self.qe_sat_list[2] / self.optimum_qe_30_dcr

                self.optimum_qe_10_ap = self.ap_sat[
                    self.index_optimum_qe_10[0]][ self.index_optimum_qe_10[1]][ self.index_optimum_qe_10[2]]
                self.optimum_qe_20_ap = self.ap_sat[
                    self.index_optimum_qe_20[0]][ self.index_optimum_qe_20[1]][ self.index_optimum_qe_20[2]]
                self.optimum_qe_30_ap = self.ap_sat[
                    self.index_optimum_qe_30[0]][self.index_optimum_qe_30[1]][self.index_optimum_qe_30[2]]

                self.optimum_qe_10_vb = self.vb_sat[
                    self.index_optimum_qe_10[0]][ self.index_optimum_qe_10[1]][ self.index_optimum_qe_10[2]]
                self.optimum_qe_20_vb = self.vb_sat[
                    self.index_optimum_qe_20[0]][ self.index_optimum_qe_20[1]][ self.index_optimum_qe_20[2]]
                self.optimum_qe_30_vb = self.vb_sat[
                    self.index_optimum_qe_30[0]][self.index_optimum_qe_30[1]][self.index_optimum_qe_30[2]]


                gate_index_10 = index_optimum_qe_10
                gate_index_20 = index_optimum_qe_20
                gate_index_30 = index_optimum_qe_30

                vb_index_10 = np.argmin(np.abs(np.array(self.SPD.grid[i][gate_index_10]) - self.optimum_qe_10_vb))
                vb_index_20 = np.argmin(np.abs(np.array(self.SPD.grid[i][gate_index_20]) - self.optimum_qe_20_vb))
                vb_index_30 = np.argmin(np.abs(np.array(self.SPD.grid[i][gate_index_30]) - self.optimum_qe_30_vb))

                optimum_grid_10_vb = self.SPD.grid[i][ gate_index_10][ vb_index_10]
                optimum_grid_20_vb = self.SPD.grid[i][ gate_index_20][ vb_index_20]
                optimum_grid_30_vb = self.SPD.grid[i][gate_index_30][vb_index_30]

                self.optimum_qe_10_dt = self.SPD.dt_list[i][ gate_index_10][ vb_index_10]
                self.optimum_qe_20_dt = self.SPD.dt_list[i][ gate_index_20][ vb_index_20]
                self.optimum_qe_30_dt = self.SPD.dt_list[i][gate_index_30][vb_index_30]

                self.optimum_qe_10_vg = self.SPD.stb_list[i][gate_index_10]
                self.optimum_qe_20_vg = self.SPD.stb_list[i][gate_index_20]
                self.optimum_qe_30_vg = self.SPD.stb_list[i][gate_index_30]

                self.optimum_qe_10_t = self.SPD.temp_list[i]
                self.optimum_qe_20_t = self.SPD.temp_list[i]
                self.optimum_qe_30_t = self.SPD.temp_list[i]

                gate_10 = int(self.stb_list[i][gate_index_10])
                gate_20 = int(self.stb_list[i][gate_index_20])
                gate_30 = int(self.stb_list[i][gate_index_30])

                tmp = self.temp_list[i]

                try:
                    self.optimum_qe_10_tr_hist[tmp] = np.array(
                        self.SPD.timeres_dict[i][gate_10][optimum_grid_10_vb])
                    self.vbias_trhist_10[tmp] = optimum_grid_10_vb
                except:
                    self.optimum_qe_10_tr_hist[tmp] = np.array(
                        self.SPD.timeres_dict[i][gate_10][
                            list(self.SPD.timeres_dict[i][gate_10].keys())[0]])
                    self.vbias_trhist_10[tmp] = list(self.SPD.timeres_dict[i][gate_10].keys())[0]


                try:
                    self.optimum_qe_10_tr_hist_bif[tmp] = np.array(
                        self.SPD.timeres_dict_bif[i][gate_10][optimum_grid_10_vb])
                    self.vbias_trhist_10_bif[tmp] = optimum_grid_10_vb
                except:
                    self.optimum_qe_10_tr_hist_bif[tmp] = np.array(
                        self.SPD.timeres_dict_bif[i][gate_10][
                            list(self.SPD.timeres_dict_bif[i][gate_10].keys())[0]])
                    self.vbias_trhist_10_bif[tmp] = list(self.SPD.timeres_dict_bif[i][gate_10].keys())[0]

                try:
                    self.optimum_qe_20_tr_hist[tmp] = np.array(
                        self.SPD.timeres_dict[i][gate_20][optimum_grid_20_vb])
                    self.vbias_trhist_20[tmp] = optimum_grid_20_vb
                except:
                    self.optimum_qe_20_tr_hist[tmp] = np.array(
                        self.SPD.timeres_dict[i][gate_20][
                            list(self.SPD.timeres_dict[i][gate_20].keys())[0]])
                    self.vbias_trhist_20[tmp] = list(self.SPD.timeres_dict[i][gate_20].keys())[0]

                try:
                    self.optimum_qe_20_tr_hist_bif[tmp] = np.array(
                        self.SPD.timeres_dict_bif[i][gate_20][optimum_grid_20_vb])
                    self.vbias_trhist_20_bif[tmp] = optimum_grid_20_vb
                except:
                    self.optimum_qe_20_tr_hist_bif[tmp] = np.array(
                        self.SPD.timeres_dict_bif[i][gate_20][
                            list(self.SPD.timeres_dict_bif[i][gate_20].keys())[0]])
                    self.vbias_trhist_20_bif[tmp] = list(self.SPD.timeres_dict_bif[i][gate_20].keys())[0]

                try:
                    self.optimum_qe_30_tr_hist[tmp] = np.array(
                        self.SPD.timeres_dict[i][gate_30][optimum_grid_30_vb])
                    self.vbias_trhist_30[tmp] = optimum_grid_30_vb
                except:
                    self.optimum_qe_30_tr_hist[tmp] = np.array(
                        self.SPD.timeres_dict[i][gate_30][
                            list(self.SPD.timeres_dict[i][gate_30].keys())[0]])
                    self.vbias_trhist_30[tmp] = list(self.SPD.timeres_dict[i][gate_30].keys())[0]

                try:
                    self.optimum_qe_30_tr_hist_bif[tmp] = np.array(
                        self.SPD.timeres_dict_bif[i][gate_30][optimum_grid_30_vb])
                    self.vbias_trhist_30_bif[tmp] = optimum_grid_30_vb
                except:
                    self.optimum_qe_30_tr_hist_bif[tmp] = np.array(
                        self.SPD.timeres_dict_bif[i][gate_30][
                            list(self.SPD.timeres_dict_bif[i][gate_30].keys())[0]])
                    self.vbias_trhist_30_bif[tmp] = list(self.SPD.timeres_dict_bif[i][gate_30].keys())[0]


                # Use filtration
                self.optimum_qe_10_tr_hist[tmp] = self.TR_filtration(self.optimum_qe_10_tr_hist[tmp])
                self.optimum_qe_10_tr = 10 * module_fwhm.fwtm(self.optimum_qe_10_tr_hist[tmp])

                # Use filtration
                self.optimum_qe_20_tr_hist[tmp] = self.TR_filtration(self.optimum_qe_20_tr_hist[tmp])
                self.optimum_qe_20_tr = 10 * module_fwhm.fwtm(self.optimum_qe_20_tr_hist[tmp])

                # Use filtration
                self.optimum_qe_30_tr_hist[tmp] = self.TR_filtration(self.optimum_qe_30_tr_hist[tmp])
                self.optimum_qe_30_tr = 10 * module_fwhm.fwtm(self.optimum_qe_30_tr_hist[tmp])

                #Slide
                self.bif_10 = self.find_bif(self.optimum_qe_10_tr_hist_bif[tmp])

                #Slide
                self.bif_20 = self.find_bif(self.optimum_qe_20_tr_hist_bif[tmp])

                # Slide
                self.bif_30 = self.find_bif(self.optimum_qe_30_tr_hist_bif[tmp])


                self.params['optimal_params'][tmp] = [['QE = ' + str(round(self.qe_sat_list[0], 2)), str(round(self.qe_sat_list[0], 2)),
                                                                    str(round(self.optimum_qe_10_dcr)),
                                                                    str(round(self.optimum_qe_10_snr, 4)),
                                                                    str(round(self.optimum_qe_10_ap, 2)),
                                                                    str(round(self.optimum_qe_10_dt, 2)),
                                                                    str(round(self.optimum_qe_10_tr, 0)),
                                                                    str(round(self.bif_10, 1))],
                                                      ['QE = ' + str(round(self.qe_sat_list[1], 2)), str(round(self.qe_sat_list[1], 2)),
                                                                   str(round(self.optimum_qe_20_dcr)),
                                                                   str(round(self.optimum_qe_20_snr, 4)),
                                                                   str(round(self.optimum_qe_20_ap, 2)),
                                                                   str(round(self.optimum_qe_20_dt, 2)),
                                                                   str(round(self.optimum_qe_20_tr, 0)),
                                                                    str(round(self.bif_20, 1))],
                                                      ['QE = ' + str(round(self.qe_sat_list[2], 2)),
                                                       str(round(self.qe_sat_list[2], 2)),
                                                       str(round(self.optimum_qe_30_dcr)),
                                                       str(round(self.optimum_qe_30_snr, 4)),
                                                       str(round(self.optimum_qe_30_ap, 2)),
                                                       str(round(self.optimum_qe_30_dt, 2)),
                                                       str(round(self.optimum_qe_30_tr, 0)),
                                                       str(round(self.bif_30, 1))]
                                                      ]

                if type(self.gate_voltage) == int:
                    self.optimum_snr_vg_v = 0
                    self.params['settings'][self.temp_list[i]] = [['QE = ' + str(round(self.qe_sat_list[0], 2)), str(round(self.optimum_qe_10_vg)),
                                                                  str(round(self.optimum_qe_10_vb, 2)),
                                                                  str(round(self.optimum_qe_10_t, 1))],
                                                                  ['QE = ' + str(round(self.qe_sat_list[1], 2)), str(round(self.optimum_qe_20_vg)),
                                                                   str(round(self.optimum_qe_20_vb, 2)),
                                                                   str(round(self.optimum_qe_20_t, 1))],
                                                                  ['QE = ' + str(round(self.qe_sat_list[2], 2)),
                                                                   str(round(self.optimum_qe_30_vg)),
                                                                   str(round(self.optimum_qe_30_vb, 2)),
                                                                   str(round(self.optimum_qe_30_t, 1))]
                                                                  ]
                else:
                    self.optimum_qe_10_vg_v = self.gate_voltage.vg.to_numpy()[gate_index_10] / 2
                    self.optimum_qe_20_vg_v = self.gate_voltage.vg.to_numpy()[gate_index_20] / 2
                    self.optimum_qe_30_vg_v = self.gate_voltage.vg.to_numpy()[gate_index_30] / 2

                    self.params['settings'][self.temp_list[i]] = [['QE = ' + str(round(self.qe_sat_list[0], 2)), str(self.optimum_qe_10_vg_v),
                                                                  str(round(self.optimum_qe_10_vb, 2)),
                                                                  str(round(self.optimum_qe_10_t, 1))],
                                                                  ['QE = ' + str(round(self.qe_sat_list[1], 2)), str(self.optimum_qe_20_vg_v),
                                                                   str(round(self.optimum_qe_20_vb, 2)),
                                                                   str(round(self.optimum_qe_20_t, 1))],
                                                                  ['QE = ' + str(round(self.qe_sat_list[2], 2)),
                                                                   str(self.optimum_qe_30_vg_v),
                                                                   str(round(self.optimum_qe_30_vb, 2)),
                                                                   str(round(self.optimum_qe_30_t, 1))]
                                                                  ]


        if self.criterion == 'SNR':
            for i in range(len(self.temp_list)):
                tmp = self.temp_list[i]
                self.optimum_snr_snr = np.max(self.SPD.snr[i])
                index_optimum_snr = np.where(self.SPD.snr[i] == self.optimum_snr_snr)
                self.index_optimum_snr = [i, index_optimum_snr[0][0], index_optimum_snr[1][0]]
                #(self.index_optimum_snr)
                self.optimum_snr_pde = self.SPD.qe[self.index_optimum_snr[0]][ self.index_optimum_snr[1]][ self.index_optimum_snr[2]]
                self.optimum_snr_dcr = self.SPD.dcr[self.index_optimum_snr[0]][ self.index_optimum_snr[1]][ self.index_optimum_snr[2]]
                self.optimum_snr_ap = self.SPD.ap_list[self.index_optimum_snr[0]][ self.index_optimum_snr[1]][ self.index_optimum_snr[2]]
                self.optimum_snr_dt = self.SPD.dt_list[self.index_optimum_snr[0]][ self.index_optimum_snr[1]][ self.index_optimum_snr[2]]

                self.optimum_snr_vg = self.SPD.stb_list[i][self.index_optimum_snr[1]]
                self.optimum_snr_vb = self.SPD.grid[self.index_optimum_snr[0]][ self.index_optimum_snr[1]][ self.index_optimum_snr[2]]
                self.optimum_snr_t = self.SPD.temp_list[self.index_optimum_snr[0]]

                gate = self.stb_list[i][self.index_optimum_snr[1]]
                #For the simple time resolution
                try:
                    self.optimum_snr_tr_hist[self.temp_list[i]] = np.array(self.SPD.timeres_dict[self.index_optimum_snr[0]][gate][self.optimum_snr_vb])
                    self.vbias_trhist[self.temp_list[i]] = self.optimum_snr_vb
                except:
                    self.optimum_snr_tr_hist[self.temp_list[i]] = np.array(self.SPD.timeres_dict[self.index_optimum_snr[0]][gate][list(self.SPD.timeres_dict[i][gate].keys())[0]])
                    self.vbias_trhist[self.temp_list[i]] = list(self.SPD.timeres_dict[i][gate].keys())[0]
                #print(self.optimum_snr_tr_hist[self.temp_list[i]])

                #For the bif time resolution
                try:
                    self.optimum_snr_tr_hist_bif[tmp] = np.array(
                        self.SPD.timeres_dict_bif[i][gate][self.optimum_snr_vb])
                    self.vbias_trhist_bif[tmp] = self.optimum_snr_vb
                except:
                    self.optimum_snr_tr_hist_bif[tmp] = np.array(
                        self.SPD.timeres_dict_bif[i][gate][
                            list(self.SPD.timeres_dict_bif[i][gate].keys())[0]])
                    self.vbias_trhist_bif[tmp] = list(self.SPD.timeres_dict_bif[i][gate].keys())[0]

                #Use filtration
                self.optimum_snr_tr_hist[tmp] = self.TR_filtration(self.optimum_snr_tr_hist[tmp])
                self.optimum_snr_tr = 10 * module_fwhm.fwtm(self.optimum_snr_tr_hist[tmp])

                # Slide
                self.bif_snr = self.find_bif(self.optimum_snr_tr_hist_bif[tmp])

                self.params['optimal_params'][self.temp_list[i]] = ['SNR',
                                                                str(round(self.optimum_snr_pde, 2)),
                                                                str(round(self.optimum_snr_dcr)),
                                                                str(round(self.optimum_snr_snr, 4)),
                                                                str(round(self.optimum_snr_ap, 2)),
                                                                str(round(self.optimum_snr_dt, 2)),
                                                                str(round(self.optimum_snr_tr, 0)),
                                                                str(round(self.bif_snr, 1))]
                if type(self.gate_voltage) == int:
                    self.optimum_snr_vg_v = 0
                    self.params['settings'][self.temp_list[i]] = ['SNR',
                                                                  str(round(self.optimum_snr_vg)),
                                                                  str(round(self.optimum_snr_vb, 2)),
                                                                  str(round(self.optimum_snr_t, 1))]
                else:
                    self.optimum_snr_vg_v = self.gate_voltage.vg.to_numpy()[self.index_optimum_snr[1]] / 2
                    self.params['settings'][self.temp_list[i]] = ['SNR',
                                                                  str(self.optimum_snr_vg_v),
                                                                  str(round(self.optimum_snr_vb, 2)),
                                                                  str(round(self.optimum_snr_t, 1))]

        #print(self.params['optimal_params'])
        #print(self.params['settings'])

    def find_bif(self, hist):
        slide_ind_10_bif = 0
        min_bif_10 = 1e9
        for j in range(2000, 3000 - 320):
            sm = (np.sum(hist[j: j + 320]) /
                  np.sum(hist[j + 320: j + 640]))
            if sm < min_bif_10:
                min_bif_10 = sm
                slide_ind_10_bif = j
        return 100 * min_bif_10

    def plot_TR(self, str_spd=""):
        #self.optimum_snr_tr_hist
        if self.criterion == 'QE':
            for el in self.temp_list:
                fig, ax1 = plt.subplots()
                fig.set_size_inches(8, 6)
                if type(self.optimum_qe_10_tr_hist[el]) != int:
                    tr_for_plot = self.optimum_qe_10_tr_hist[el][::-1]
                    max_ind = np.argmax(tr_for_plot)
                    tr_for_plot = tr_for_plot[max_ind - 150 : max_ind + 150]

                    ax1.plot(0.01 * np.arange(0, len(tr_for_plot)), tr_for_plot, 'b')
                ax1.set_xlabel(r'$t$, ns')
                ax1.set_ylabel('Counts')
                if self.SPD.type == 'gated' and type(self.gate_voltage) != int:
                    plt.title(r'TR for T = ' + str(el) + ' C. Vb = ' + str(round(self.vbias_trhist_10[el], 2)) + r' V. $V_g$ = ' + str(self.optimum_qe_10_vg_v) + ' V' )
                    filename = os.path.join(self.params['dir_name'], self.fold, self.criterion) + '//TR_T' + str(el) + '_VG' + str_spd + '.png'
                elif self.SPD.type == 'gated':
                    plt.title(r'TR for T = ' + str(el) + ' C. Vb = ' + str(round(self.vbias_trhist_10[el], 2)) + r' V.')
                    filename = os.path.join(self.params['dir_name'], self.fold, self.criterion) + '//TR_T' + str(
                        el) + '_CD' + str_spd + '.png'
                elif self.SPD.type == 'freerun':
                    plt.title(r'TR for T = ' + str(el) + ' C. Vb = ' + str(round(self.vbias_trhist_10[el], 2)) + r' V.')
                    filename = os.path.join(self.params['dir_name'], self.fold, self.criterion) + '//TR_T' + str(
                        el) + str_spd + '.png'
                #plt.legend(fontsize=14, loc='upper right')
                plt.grid(True)
                plt.savefig(filename, bbox_inches='tight', transparent=True, dpi=100)
                plt.cla()
                #plt.show()

        if self.criterion == 'SNR':
            for el in self.temp_list:
                fig, ax1 = plt.subplots()
                fig.set_size_inches(8, 6)
                if type(self.optimum_snr_tr_hist[el]) != int:
                    tr_for_plot = self.optimum_snr_tr_hist[el][::-1]
                    max_ind = np.argmax(tr_for_plot)
                    tr_for_plot = tr_for_plot[max_ind - 150: max_ind + 150]

                    ax1.plot(0.01 * np.arange(0, len(tr_for_plot)), tr_for_plot, 'b')
                ax1.set_xlabel(r'$t$, ns')
                ax1.set_ylabel('Counts')
                if self.SPD.type == 'gated' and type(self.gate_voltage) != int:
                    plt.title(r'TR for T = ' + str(el) + ' C. Vb = ' + str(round(self.vbias_trhist[el], 2)) + r' V. $V_g$ = ' + str(self.optimum_snr_vg_v) + ' V' )
                    filename = os.path.join(self.params['dir_name'], self.fold, self.criterion) + '//TR_T' + str(el) + '_CD' + str_spd + '.png'
                elif self.SPD.type == 'gated':
                    plt.title(r'TR for T = ' + str(el) + ' C. Vb = ' + str(round(self.vbias_trhist[el], 2)) + r' V.')
                    filename = os.path.join(self.params['dir_name'], self.fold, self.criterion) + '//TR_T' + str(
                        el) + '_VG' + str_spd + '.png'
                elif self.SPD.type == 'freerun':
                    plt.title(r'TR for T = ' + str(el) + ' C. Vb = ' + str(round(self.vbias_trhist[el], 2)) + r' V.')
                    filename = os.path.join(self.params['dir_name'], self.fold, self.criterion) + '//TR_T' + str(
                        el) + str_spd +'.png'
                #plt.legend(fontsize=14, loc='upper right')
                plt.grid(True)
                plt.savefig(filename, bbox_inches='tight', transparent=True, dpi=100)
                plt.cla()
                #plt.show()

    def plot_heatmap_bif(self, str_spd=""):
        self.bif_hm = {}
        for i in range(len(self.temp_list)):
            mesh_strobe, grid, bif_hm = [], [], []
            for key in self.SPD.timeres_dict_bif[i].keys():
                mesh_strobe.append([])
                grid.append([])
                bif_hm.append([])
                for key_vb in self.SPD.timeres_dict_bif[i][key].keys():
                    mesh_strobe[-1].append(key)
                    grid[-1].append(key_vb)
                    bif_hm[-1].append(self.find_bif(self.SPD.timeres_dict_bif[i][key][key_vb]))

            self.bif_hm[i] = np.array(bif_hm)
            mesh_strobe = np.array(mesh_strobe)

            fig, ax1 = plt.subplots()
            fig.set_size_inches(8, 6)
            divider = make_axes_locatable(ax1)
            cax = divider.append_axes('right', size='5%', pad=0.05)
            if type(self.gate_voltage) == int:
                im = ax1.pcolormesh(mesh_strobe, grid, self.bif_hm[i], cmap=cm.coolwarm)
                ax1.set_xlabel(r'$V_g$, code')
                filename = os.path.join(self.params['dir_name'], self.fold, self.criterion) + '//HM_BIF_T' + str(
                    self.temp_list[i]) + '_CD' + str_spd +'.png'
            else:
                mesh_strobe_v = []
                for j in range(mesh_strobe.shape[0]):
                    mesh_strobe_v.append([self.gate_voltage.vg.to_numpy()[j] / 2 for _ in range(mesh_strobe.shape[1])])
                mesh_strobe_v = np.array(mesh_strobe_v)
                im = ax1.pcolormesh(mesh_strobe_v, grid, self.bif_hm[i], cmap=cm.coolwarm)
                ax1.set_xlabel(r'$V_g$, V')
                filename = os.path.join(self.params['dir_name'], self.fold, self.criterion) + '//HM_BIF_T' + str(
                    self.temp_list[i]) + '_VG' + str_spd + '.png'

            cb = fig.colorbar(im, cax=cax, orientation='vertical')
            ax1.set_ylabel(r'$V_b$, V')
            ax1.grid(True)
            ax1.set_title('BIF heatmap for T = ' + str(self.temp_list[i]) + ' C')
            plt.savefig(filename, bbox_inches='tight', transparent=True, dpi=100)
            # plt.show()
            cb.remove()
            plt.cla()

    def plot_heatmap(self, str_spd=""):
        for i in range(len(self.temp_list)):
            mesh_strobe = self.SPD.mesh_strobe[i]
            grid = self.SPD.grid[i]

            fig, ax1 = plt.subplots()
            fig.set_size_inches(8, 6)
            divider = make_axes_locatable(ax1)
            cax = divider.append_axes('right', size='5%', pad=0.05)
            if type(self.gate_voltage) == int:
                im = ax1.pcolormesh(mesh_strobe, grid, self.SPD.snr[i], cmap=cm.coolwarm)
                ax1.set_xlabel(r'$V_g$, code')
                filename = os.path.join(self.params['dir_name'], self.fold, self.criterion) + '//HM_SNR_T' + str(self.temp_list[i]) + '_CD' + str_spd + '.png'
            else:
                mesh_strobe_v = [[] for _ in range(len(self.stb_list[i]))]
                for j in range(len(self.stb_list[i])):
                    mesh_strobe_v[j] = [self.gate_voltage.vg.to_numpy()[j] / 2 for _ in range(len(mesh_strobe[j]))]
                mesh_strobe_v = np.array(mesh_strobe_v)
                #print(mesh_strobe)
                #print(mesh_strobe_v)
                im = ax1.pcolormesh(mesh_strobe_v, grid, self.SPD.snr[i], cmap=cm.coolwarm)
                ax1.set_xlabel(r'$V_g$, V')
                filename = os.path.join(self.params['dir_name'], self.fold, self.criterion) + '//HM_SNR_T' + str(
                    self.temp_list[i]) + '_VG' + '.png'

            cb = fig.colorbar(im, cax=cax, orientation='vertical')
            ax1.set_ylabel(r'$V_b$, V')
            ax1.grid(True)
            ax1.set_title('SNR heatmap for T = ' + str(self.temp_list[i]) + ' C')
            plt.savefig(filename, bbox_inches='tight', transparent=True, dpi=100)
            #plt.show()
            cb.remove()
            plt.cla()

    def generate_pdf_gated(self, str_spd=""):
        # Block with gate voltage
        if type(self.SPD.gate_voltage) != int:
            self.params['isvg'] = True
            self.params['settings'] = {'head': ['Criterion', 'Vg, V', 'Vb, V', 'T, C']}
            self.gate_voltage = self.SPD.gate_voltage
            self.get_slices()
            self.plot_dcr(str_spd)
            self.plot_ap(str_spd)
            self.choose_optimum()
            self.plot_TR(str_spd)
            self.plot_heatmap(str_spd)
            self.plot_heatmap_bif(str_spd)
            #generate_pdf.create_pdf(self.params)

        # Block without gate voltage
        self.params['isvg'] = False
        self.params['settings'] = {'head': ['Criterion', 'Vg, code', 'Vb, V', 'T, C']}
        self.gate_voltage = 0
        self.get_slices()
        self.plot_dcr(str_spd)
        self.plot_ap(str_spd)
        self.choose_optimum()
        self.plot_TR(str_spd)
        self.plot_heatmap(str_spd)
        self.plot_heatmap_bif(str_spd)

    def generate_pdf_freerun(self, str_spd=""):
        self.plot_dcr_qe_fr(str_spd)
        self.plot_ap_fr(str_spd)
        self.plot_snr_fr(str_spd)
        self.choose_optimum_fr()
        self.plot_TR(str_spd)

    def generate_pdf(self):
        for el in self.criterion_list:
            self.criterion = el
            self.params['criterion'] = el

            if self.SPD.type == 'gated':
                self.generate_pdf_gated(self.params["str_spd"])

            elif self.SPD.type == 'freerun':
                self.generate_pdf_freerun(self.params["str_spd"])
