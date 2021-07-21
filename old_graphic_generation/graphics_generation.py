import argparse
import sys
import os
import numpy as np
import xlsxwriter
import openpyxl
import spd_analyze
import pandas as pd
import module_fwhm
import generate_pdf
from matplotlib import cm

import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams.update({'font.size': 22})

#Vereya
#dir_name = '!SPDStandResults//DEFNAME//21_04_27__11_36'
#dof_name = 'Vereya'

class graphics_generation(object):
    def __init__(self, params):
        #Uliana
        #dir_name = '!SPDStandResults//Uliana//21_05_13__18_16'
        self.params = params
        self.dir_name = params['dir_name']
        self.dof_name = params['name']
        self.fold = 'report'
        self.params['fold'] = self.fold

        try:
            os.mkdir(os.path.join(self.dir_name, self.fold))
        except:
            print(os.path.join(self.dir_name, self.fold) + ' folder currently exists')

        self.SPD = spd_analyze.SPD(self.dir_name)
        self.temp_list = self.SPD.temp_list
        self.qe_sat_list = [10, 15, 20, 25]
        self.stb_list = self.SPD.stb_list
        self.colors = ['b', 'g', 'r', 'magenta', 'k', 'cyan']
        self.params['optimal_params'] = [['Criteria', 'PDE, %', 'DCR, Hz', 'SNR', 'AP, %', 'DT, mus', 'TR, ps']]
        self.params['settings'] = [['Criteria', 'Vg', 'Vb', 'T']]

    def get_sub_slice(self, qe_sat):
        # First part about PDE

        dcr_sat, ap_sat, vb_sat = [], [], []
        for i in range(len(self.temp_list)):
            dcr_sat.append([])
            ap_sat.append([])
            vb_sat.append([])
            for j in range(len(self.SPD.stb_list)):
                qe_list = np.array(self.SPD.qe[i,j,:])
                ind_vb2 = np.argmin(abs(qe_list * (qe_list >= qe_sat) - qe_sat))
                ind_vb1 = np.argmin(abs(qe_list * (qe_list < qe_sat) - qe_sat))
                qe1 = qe_list[ind_vb1]
                qe2 = qe_list[ind_vb2]

                dcr1 = self.SPD.dcr[i,j,ind_vb1]
                dcr2 = self.SPD.dcr[i,j,ind_vb2]
                dcr_sat[i].append((dcr2 - dcr1) / (qe2 - qe1) * (qe_sat - qe1) + dcr1)

                ap1 = self.SPD.ap_list[i,j,ind_vb1]
                ap2 = self.SPD.ap_list[i,j,ind_vb2]
                ap_sat[i].append((ap2 - ap1) / (qe2 - qe1) * (qe_sat - qe1) + ap1)

                vb1 = self.SPD.grid[i,j,ind_vb1]
                vb2 = self.SPD.grid[i,j,ind_vb2]
                vb_sat[i].append((vb2 - vb1) / (qe2 - qe1) * (qe_sat - qe1) + vb1)
        return np.array(dcr_sat),  np.array(ap_sat), np.array(vb_sat)

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

    def plot_dcr(self):
        for i in range(len(self.qe_sat_list)):
            figure = plt.gcf()
            figure.set_size_inches(8, 6)
            for j in range(len(self.temp_list)):
                plt.plot(self.stb_list, self.dcr_sat[i, j, :], color=self.colors[j], marker='o',
                         label='T = ' + str(self.temp_list[j]) + r'${}^\circ C$')

            plt.xlabel(r'$V_g$')
            plt.ylabel('DCR, Hz')
            plt.title(r'DCR dependence on $V_g$. PDE = ' + str(self.qe_sat_list[i]) + '\%')
            plt.legend(fontsize=14, loc='upper right')
            plt.grid(True)
            plt.savefig(os.path.join(self.params['dir_name'], self.fold) + '//DCR_PDE' + str(self.qe_sat_list[i]) + '.png',
                        bbox_inches='tight', transparent=True, dpi=100)
            plt.cla()

    def plot_ap(self):
        for i in range(len(self.qe_sat_list)):
            figure = plt.gcf()
            figure.set_size_inches(8, 6)
            for j in range(len(self.temp_list)):
                plt.plot(self.stb_list, self.ap_sat[i, j, :], color=self.colors[j], marker='o',
                         label='T = ' + str(self.temp_list[j]) + r'${}^\circ C$')

            plt.xlabel(r'$V_g$')
            plt.ylabel('AP, \%')
            plt.title(r'AP dependence on $V_g$. PDE = ' + str(self.qe_sat_list[i]) + '\%')
            plt.legend(fontsize=14, loc='upper right')
            plt.grid(True)
            plt.savefig(os.path.join(self.params['dir_name'], self.fold) + '//AP_PDE' + str(self.qe_sat_list[i]) + '.png',
                        bbox_inches='tight', transparent=True, dpi=100)
            plt.cla()
            # plt.show()

    def choose_optimum(self, criteria = 'SNR'):
        if criteria == 'SNR':
            self.optimum_snr_snr = np.max(self.SPD.snr)
            index_optimum_snr = np.where(self.SPD.snr == self.optimum_snr_snr)
            self.index_optimum_snr = [index_optimum_snr[0][0], index_optimum_snr[1][0], index_optimum_snr[2][0]]
            self.optimum_snr_pde = self.SPD.qe[self.index_optimum_snr[0], self.index_optimum_snr[1], self.index_optimum_snr[2]]
            self.optimum_snr_dcr = self.SPD.dcr[self.index_optimum_snr[0], self.index_optimum_snr[1], self.index_optimum_snr[2]]
            self.optimum_snr_ap = self.SPD.ap_list[self.index_optimum_snr[0], self.index_optimum_snr[1], self.index_optimum_snr[2]]
            self.optimum_snr_dt = self.SPD.dt_list[self.index_optimum_snr[0], self.index_optimum_snr[1], self.index_optimum_snr[2]]

            self.optimum_snr_vg = self.SPD.stb_list[self.index_optimum_snr[1]]
            self.optimum_snr_vb = self.SPD.grid[self.index_optimum_snr[0], self.index_optimum_snr[1], self.index_optimum_snr[2]]
            self.optimum_snr_t = self.SPD.temp_list[self.index_optimum_snr[0]]

            gate = self.stb_list[self.index_optimum_snr[1]]
            #print(self.SPD.timeres_dict[0][gate][self.optimum_snr_vb])
            try:
                self.optimum_snr_tr_hist = np.array(self.SPD.timeres_dict[self.index_optimum_snr[0]][gate][self.optimum_snr_vb])
            except:
                self.optimum_snr_tr_hist = 0
            self.optimum_snr_tr = 10 * module_fwhm.fwhm(self.optimum_snr_tr_hist)

            self.params['optimal_params'].append(['SNR', str(round(self.optimum_snr_pde, 2)),
                                                  str(round(self.optimum_snr_dcr)), str(round(self.optimum_snr_snr, 4)),
                                                  str(round(self.optimum_snr_ap, 2)), str(round(self.optimum_snr_dt, 2)), str(round(self.optimum_snr_tr, 0))])
            self.params['settings'].append(['SNR', str(round(self.optimum_snr_vg)), str(round(self.optimum_snr_vb, 2)), str(round(self.optimum_snr_t, 1))])

            #print(self.params)
            #print(self.SPD.snr[self.index_optimum_snr])

    def plot_TR(self):
        #self.optimum_snr_tr_hist
        figure = plt.gcf()
        figure.set_size_inches(8, 6)
        if type(self.optimum_snr_tr_hist) != int:
            plt.plot(0.01 * np.arange(0, len(self.optimum_snr_tr_hist)), self.optimum_snr_tr_hist[::-1], 'b')

        plt.xlabel(r'$t$, ns')
        plt.ylabel('Counts')
        plt.title(r'TR for optimal SNR point')
        plt.legend(fontsize=14, loc='upper right')
        plt.grid(True)
        plt.savefig(
            os.path.join(self.params['dir_name'], self.fold) + '//TR.png', bbox_inches='tight', transparent=True, dpi=100)
        plt.cla()
        #plt.show()

    def plot_heatmap(self):
        for i in range(len(self.temp_list)):
            mesh_strobe = self.SPD.mesh_strobe
            grid = self.SPD.grid[i]
            figure = plt.gcf()
            figure.set_size_inches(8, 6)
            plt.pcolormesh(mesh_strobe, grid, self.SPD.snr[i], cmap=cm.coolwarm)

            #if i == self.index_optimum_snr[0]:
            #   plt.scatter(self.optimum_snr_vg, self.optimum_snr_vb, 1 + self.optimum_snr_snr, color = 'magenta')

            cb = plt.colorbar()
            plt.ylabel(r'$V_b$, V')
            plt.xlabel(r'$V_g$')
            plt.grid(True)
            plt.title('SNR heatmap for T = ' + str(self.temp_list[i]))
            plt.savefig(os.path.join(self.params['dir_name'], self.fold) + '//HM_SNR_T' + str(round(self.temp_list[i])) + '.png', bbox_inches='tight', transparent=True, dpi=100)
            #plt.show()
            cb.remove()
            plt.cla()

    def generate_pdf(self):
        generate_pdf.create_pdf(self.params)