import numpy as np
import os
from matplotlib import cm
import pandas as pd
from matplotlib.ticker import LinearLocator
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams.update({'font.size': 22})
import afterpulse as aft
from time_resolution import time_resolution
import module_pde

import warnings
warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning)


class SPD(object):
    def __init__(self, dir_name):
        #Dowload temperature files
        self.dir_name = dir_name
        #self.type = 'gated'
        temp_name = 'Grid_Temperature.csv'
        temp_list = []
        with open(os.path.join(dir_name, temp_name)) as fn:
            for line in fn:
                line = line.strip().split(';')
                for el in line:
                    temp_list.append(float(el))
        self.temp_list = temp_list

        #Download grid for gate:
        #And derive the SPD type
        self.type = 'freerun'
        self.stb_name = 'Grid_Stb_Code.csv'
        for fn in os.listdir(dir_name):
            if fn == self.stb_name:
                self.type = 'gated'

        if self.type == 'freerun':
            self.load_freerun()
        elif self.type == 'gated':
            self.load_gated()

    def func_exp(self, x, a, b):
        return a * np.exp(b * x)

    def load_freerun(self):
        n = len(self.temp_list)
        temp_fold = ['T_' + str(el) for el in self.temp_list]

        # Download the afterpulse csv
        timeres_dict = [[] for i in range(n)]
        hv_timeres_dict = [[] for i in range(n)]
        ap_dict = [[] for i in range(n)]
        dt_dict = [[] for i in range(n)]
        hv_dict = [[] for i in range(n)]
        ap_list = [[] for i in range(n)]
        dt_list = [[] for i in range(n)]
        dcr = [[] for i in range(n)]
        sigma_dcr = [[] for i in range(n)]
        R = [[] for i in range(n)]
        QE = [[] for i in range(n)]
        grid = [[] for i in range(n)]


        # Download grid
        file_grid = 'Grid_HV_bias.csv'
        for i in range(n):
            with open(os.path.join(self.dir_name, temp_fold[i], file_grid)) as fn:
                for line in fn:
                    line = line.strip().split(';')
                    grid[i] = [round(float(line[k]), 4) for k in range(len(line))]
        grid = np.array(grid)
        #print(grid)

        # Download Histograms for AP and time resolution
        for i in range(n):
            ap_gates = '0x000'
            conv_ap_gates = int(ap_gates.replace('0x', ''), 16)

            timeres_dict_hlp = {}
            ap_dict_hlp = {}
            dt_dict_hlp = {}
            hv_dict_hlp = {}

            # Read histograms for AP/DT and for time resolution
            path = os.path.join(self.dir_name, temp_fold[i], ap_gates)
            AP_list, DT_list, HV_list = [0], [], [0]
            HV_list_tr = []
            tr_dict_hlp2 = {}
            for file_name in os.listdir(path):
                if file_name.startswith('AP_'):
                    path2 = os.path.join(path, file_name)
                    DT, AP = aft.Ap_analyze(path2, self.type)
                    AP_list.append(AP)
                    DT_list.append(DT)
                    cur_hv = float(file_name.replace('AP_', '').replace('.csv', ''))
                    HV_list.append(cur_hv - grid[i][0])
                if file_name.startswith('TR_'):
                    path2 = os.path.join(path, file_name)
                    trhist = time_resolution(path2)
                    cur_hv = float(file_name.replace('TR_', '').replace('.csv', ''))
                    HV_list_tr.append(grid[i][np.argmin(np.abs(grid[i][:] - cur_hv))])
                    tr_dict_hlp2[HV_list_tr[-1]] = trhist

            red_gt = conv_ap_gates
            timeres_dict_hlp[red_gt] = tr_dict_hlp2
            ap_dict_hlp[red_gt] = AP_list
            dt_dict_hlp[red_gt] = DT_list
            hv_dict_hlp[red_gt] = HV_list

            timeres_dict[i] = timeres_dict_hlp
            ap_dict[i] = ap_dict_hlp
            dt_dict[i] = dt_dict_hlp
            hv_dict[i] = hv_dict_hlp

            '''
            print(timeres_dict)
            print(ap_dict)
            print(dt_dict)
            print(hv_dict)
            '''

        for i in range(n):
            vbias_small = hv_dict[i][0]
            AP_small = ap_dict[i][0]
            popt, pcov = curve_fit(self.func_exp, vbias_small, AP_small, method='lm', maxfev=1800, p0=[1, 1])
            req_hvlist = np.array(grid[i][:])
            req_hvlist = req_hvlist - req_hvlist[0]
            ap_list[i] = self.func_exp(req_hvlist, *popt)
            DT_aver = np.mean(dt_dict[i][0])
            dt_list[i] = [DT_aver for i in range(len(req_hvlist))]

        # Download DCR, sigmaDCR and QE:
        for i in range(n):
            file_dcr = 'DCR.csv'
            file_sigma_dcr = 'DCRSigma.csv'
            file_R = 'QE.csv'

            with open(os.path.join(self.dir_name, temp_fold[i], file_dcr)) as fn:
                for line in fn:
                    line = line.strip().split(';')
                    dcr[i] = [float(line[k]) for k in range(len(line))]
                    for k in range(len(dcr[i])):
                        if dcr[i][k] == 0:
                            dcr[i][k] = 30

            with open(os.path.join(self.dir_name, temp_fold[i], file_sigma_dcr)) as fn:
                for line in fn:
                    line = line.strip().split(';')
                    sigma_dcr[i] = [float(line[k]) for k in range(len(line))]

            with open(os.path.join(self.dir_name, temp_fold[i], file_R)) as fn:
                for line in fn:
                    line = line.strip().split(';')
                    R[i]= [float(line[k]) for k in range(len(line))]
                    QE[i] = [module_pde.pde(ap_list[i][k], dt_list[i][k], R[i][k], dcr[i][k])
                                  for k in range(len(line))]

        dcr = np.array(dcr)
        sigma_dcr = np.array(sigma_dcr)
        R = np.array(R)
        qe = np.array(QE)
        grid = np.array(grid)
        ap_list = np.array(ap_list)
        dt_list = np.array(dt_list)

        # Save to the class:
        self.dcr = dcr
        self.R = R
        self.qe = 100 * qe
        self.sigma_dcr = sigma_dcr
        self.grid = grid
        self.ap_list = 100 * ap_list
        self.dt_list = 1e6 * dt_list
        self.timeres_dict = timeres_dict
        self.hv_timeres_dict = hv_timeres_dict
        self.snr = self.qe / self.dcr
        #print(self.timeres_dict)
        #print(self.dt_list)

    def load_gated(self):
        n = len(self.temp_list)

        try:
            self.gate_voltage = pd.read_excel(os.path.join(self.dir_name, 'gate_voltage.xlsx'))
        except:
            self.gate_voltage = 0
            print('there is no gate voltage tables')
        temp_fold = ['T_' + str(el) for el in self.temp_list]
        stb_list = []
        with open(os.path.join(self.dir_name, self.stb_name)) as fn:
            for line in fn:
                line = line.strip().split(';')
                for el in line:
                    stb_list.append(float(el))
        N = len(stb_list)

        # Download the afterpulse csv
        timeres_dict = [[] for i in range(n)]
        timeres_dict_bif = [[] for i in range(n)]
        hv_timeres_dict = [[] for i in range(n)]
        hv_timeres_dict_bif = [[] for i in range(n)]
        ap_dict = [[] for i in range(n)]
        dt_dict = [[] for i in range(n)]
        hv_dict = [[] for i in range(n)]
        ap_list = [[] for i in range(n)]
        dt_list = [[] for i in range(n)]
        dcr = [[] for i in range(n)]
        sigma_dcr = [[] for i in range(n)]
        R = [[] for i in range(n)]
        QE = [[] for i in range(n)]
        grid = [[] for i in range(n)]
        n_apgates = []

        '''
        for i in range(n):
            for j in range(N):
                dcr[i] = [[] for j in range(N)]
                sigma_dcr[i] = [[] for j in range(N)]
                R[i] = [[] for j in range(N)]
                QE[i] = [[] for j in range(N)]
                grid[i] = [[] for j in range(N)]
        '''

        #Download grid
        file_grid = 'Grid_HV_bias.csv'
        for i in range(n):
            with open(os.path.join(self.dir_name, temp_fold[i], file_grid)) as fn:
                ind = 0
                for line in fn:
                    line = line.strip().split(';')
                    if line[0] != 'NaN':
                        grid[i].append([round(float(line[k]), 4) for k in range(len(line))])
                        ind += 1
        grid = np.array(grid)
        #print(grid)

        #Download Histograms for AP and time resolution
        for i in range(n):
            ap_gates = []
            for file_name in os.listdir(os.path.join(self.dir_name, temp_fold[i])):
                path = os.path.join(self.dir_name, temp_fold[i], file_name)
                if os.path.isdir(path):
                    ap_gates.append(file_name)

            conv_ap_gates = [int(el.replace('0x', ''), 16) for el in ap_gates]
            n_apgates.append(len(ap_gates))
            #print(n_apgates)

            timeres_dict_hlp = {}
            timeres_dict_hlp_bif = {}
            ap_dict_hlp = {}
            dt_dict_hlp = {}
            hv_dict_hlp = {}

            #Read histograms for AP/DT and for time resolution
            for j in range(n_apgates[i]):
                path = os.path.join(self.dir_name, temp_fold[i], ap_gates[j])
                AP_list, DT_list, HV_list = [0], [], [0]
                HV_list_tr, HV_list_tr_bif = [], []
                tr_dict_hlp2, tr_dict_hlp2_bif = {}, {}
                for file_name in os.listdir(path):
                    if file_name.startswith('AP_'):
                        path2 = os.path.join(path, file_name)
                        DT, AP = aft.Ap_analyze(path2, self.type)
                        AP_list.append(AP)
                        DT_list.append(DT)
                        cur_hv = float(file_name.replace('AP_', '').replace('.csv', ''))
                        HV_list.append(cur_hv - grid[i][j][0])
                    if file_name.startswith('TR_'):
                        path2 = os.path.join(path, file_name)
                        trhist = time_resolution(path2)
                        if file_name.startswith('TR_Bif_'):
                            cur_hv = float(file_name.replace('TR_Bif_', '').replace('.csv', ''))
                            #print(cur_hv)
                            #print(np.array(grid[i][j]) - cur_hv)
                            HV_list_tr_bif.append(grid[i][j][np.argmin(np.abs(np.array(grid[i][j][:]) - cur_hv))])
                            tr_dict_hlp2_bif[HV_list_tr_bif[-1]] = trhist
                        else:
                            cur_hv = float(file_name.replace('TR_', '').replace('.csv', ''))
                            #print(cur_hv)
                            #print(np.array(grid[i][j]) - cur_hv)
                            HV_list_tr.append(grid[i][j][np.argmin(np.abs(np.array(grid[i][j][:]) - cur_hv))])
                            tr_dict_hlp2[HV_list_tr[-1]] = trhist

                red_gt = int(conv_ap_gates[j])
                timeres_dict_hlp[red_gt] = tr_dict_hlp2
                timeres_dict_hlp_bif[red_gt] = tr_dict_hlp2_bif
                ap_dict_hlp[red_gt] = AP_list
                dt_dict_hlp[red_gt] = DT_list
                hv_dict_hlp[red_gt] = HV_list
            timeres_dict[i] = timeres_dict_hlp
            timeres_dict_bif[i] = timeres_dict_hlp_bif
            ap_dict[i] = ap_dict_hlp
            dt_dict[i] = dt_dict_hlp
            hv_dict[i] = hv_dict_hlp

        for i in range(n):
            #print('***')
            #print(n_apgates[i])
            #print(self.temp_list[i])
            for j in range(n_apgates[i]):
                try:
                    vbias_small = hv_dict[i][stb_list[j]]
                    AP_small = ap_dict[i][stb_list[j]]
                    popt, pcov = curve_fit(self.func_exp, vbias_small, AP_small, method='lm', maxfev=1800, p0=[1, 1])
                    req_hvlist = np.array(grid[i][j][:])
                    req_hvlist = req_hvlist - req_hvlist[0]
                    #print(req_hvlist)
                    ap_list[i].append(self.func_exp(req_hvlist, *popt))
                    DT_aver = np.mean(dt_dict[i][stb_list[j]])
                    dt_list[i].append([DT_aver for i in range(len(req_hvlist))])
                except:
                    print('There is no required AP histogram for T = ' + str(self.temp_list[i]) + ' and gate = ' + str(hex(int(stb_list[j]))))
                    n_apgates[i] = j
                    del grid[i][j]

        #Download DCR, sigmaDCR and QE:
        for i in range(n):
            file_dcr = 'DCR.csv'
            file_sigma_dcr = 'DCRSigma.csv'
            file_R = 'QE.csv'

            with open(os.path.join(self.dir_name, temp_fold[i], file_dcr)) as fn:
                ind = 0
                for line in fn:
                    line = line.strip().split(';')
                    if line[0] != 'NaN' and ind < n_apgates[i]:
                        dcr[i].append([float(line[k]) for k in range(len(line))])
                        for k in range(len(dcr[i][ind])):
                            if dcr[i][ind][k] == 0:
                                dcr[i][ind][k] = 30
                        ind += 1

            #print(dcr)

            with open(os.path.join(self.dir_name, temp_fold[i], file_sigma_dcr)) as fn:
                ind = 0
                for line in fn:
                    line = line.strip().split(';')
                    if line[0] != 'NaN' and ind < n_apgates[i]:
                        sigma_dcr[i].append([float(line[k]) for k in range(len(line))])
                        ind += 1

            #print(sigma_dcr)



            with open(os.path.join(self.dir_name, temp_fold[i], file_R)) as fn:
                ind = 0
                for line in fn:
                    line = line.strip().split(';')
                    if line[0] != 'NaN' and ind < n_apgates[i]:
                        R[i].append([float(line[k]) for k in range(len(line))])
                        QE[i].append([module_pde.pde(ap_list[i][ind][k], dt_list[i][ind][k], R[i][ind][k], dcr[i][ind][k]) for k in range(len(line))])
                        ind += 1

        #Preparing for plot the Surface:

        mesh_strobe = []
        for i in range(n):
            NN = len(dcr[0][0][:])
            mesh_strobe.append([stb_list[:n_apgates[i]] for j in range(NN)])
            mesh_strobe[i] = np.array(mesh_strobe[i])
            mesh_strobe[i] = mesh_strobe[i].transpose()
        #print(mesh_strobe)
        self.mesh_strobe = mesh_strobe

        self.dt_list = []
        self.dcr = []
        self.R = []
        self.qe = []
        self.sigma_dcr = []
        self.grid = []
        self.stb_list = []
        self.ap_list = []
        self.snr = []
        self.n_apgates = n_apgates

        for i in range(n):
            self.dt_list.append(1e6 * np.array(dt_list[i]))
            self.dcr.append(np.array(dcr[i]))
            self.R.append(np.array(R[i]))
            self.qe.append(100 * np.array(QE[i]))
            self.sigma_dcr.append(np.array(sigma_dcr[i]))
            self.grid.append(np.array(grid[i]))
            self.stb_list.append(np.array(stb_list[:n_apgates[i]]))
            self.ap_list.append(100 * np.array(ap_list[i]))
            self.snr.append(self.qe[i] / self.dcr[i])
        #Save to the class:

        '''
        self.dcr = dcr
        self.R = R
        self.qe = 100 * qe
        self.sigma_dcr = sigma_dcr
        self.grid = grid
        self.stb_list = stb_list
        self.ap_list = 100 * ap_list
        '''

        self.timeres_dict = timeres_dict
        self.hv_timeres_dict = hv_timeres_dict
        self.timeres_dict_bif = timeres_dict_bif
        self.hv_timeres_dict_bif = hv_timeres_dict_bif
        #self.snr = self.qe / self.dcr

        #print(self.gate_voltage.code.iloc[5])
        #print(self.stb_list[5])
        #print(self.gate_voltage.code.iloc[5] == self.stb_list[5])

    def plot_section(self, choosed_T, choosed_graph):
        fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

        #choosed_T = 1
        # ablut graphs: 0 == dcr; 1 == qe; 2 == qe / dcr
        #choosed_graph = 1

        # Plot the DCR
        if choosed_graph == 0:
            surf = ax.plot_surface(self.mesh_strobe, self.grid[choosed_T], self.dcr[choosed_T], cmap=cm.coolwarm, linewidth=0,
                                   antialiased=False)
            plt.title('DCR dependence for T = ' + str(self.temp_list[choosed_T]))
            ax.set_zlabel('DCR, Hz')
        # Plot the QE
        if choosed_graph == 1:
            surf = ax.plot_surface(self.mesh_strobe, self.grid[choosed_T], self.qe[choosed_T], cmap=cm.coolwarm, linewidth=0,
                                   antialiased=False)
            plt.title('QE dependence for T = ' + str(self.temp_list[choosed_T]))
            ax.set_zlabel('QE, \%')
        # Plot the QE/DCR
        if choosed_graph == 2:
            surf = ax.plot_surface(self.mesh_strobe, self.grid[choosed_T], self.qe[choosed_T] / self.dcr[choosed_T], cmap=cm.coolwarm,
                                   linewidth=0, antialiased=False)
            plt.title('QE/DCR dependence for T = ' + str(self.temp_list[choosed_T]))
            ax.set_zlabel('QE/DCR')
        ax.set_xlabel('Gate number')
        ax.set_ylabel('Bias')

        # fig.colorbar(surf, shrink=0.5, aspect=5)
        plt.show()
    def slice_return(self, plot_index = 0, tempind = 0, ind = 0):
        if self.type == 'gated':
            if plot_index == 0:
                ret = 1 * self.dcr[tempind][ind]
            elif plot_index == 1:
                ret = 1 * self.qe[tempind][ind]
            elif plot_index == 2:
                ret = 1 * np.array(self.qe[tempind][ind]) / np.array(self.dcr[tempind][ind])
            elif plot_index == 3:
                ret = 1 * self.ap_list[tempind][ind]
            elif plot_index == 4:
                ret = 1 * self.dt_list[tempind][ind]
            n = len(ret)
            return [self.stb_list[tempind][ind] for i in range(n)], 1 * self.grid[tempind][ind], ret
        elif self.type == 'freerun':
            if plot_index == 0:
                ret = 1 * self.dcr[tempind, :]
            elif plot_index == 1:
                ret = 1 * self.qe[tempind, :]
            elif plot_index == 2:
                ret = 1 * self.qe[tempind, :] / self.dcr[tempind, :]
            elif plot_index == 3:
                ret = 1 * self.ap_list[tempind, :]
            elif plot_index == 4:
                ret = 1 * self.dt_list[tempind, :]
            n = len(ret)
            return  self.grid[tempind, :], ret

    def timeres_return(self, tempind = 0, gateind = 0, Vb = 0):
        #hv = round(self.grid[tempind, gateind, hvind], 4)
        if self.type == 'gated':
            gate = self.stb_list[tempind][gateind]
            exitcode = 0
            if self.timeres_dict[tempind].get(gate) != None:
                tr_hist = self.timeres_dict[tempind][gate][Vb]
            else:
                exitcode = 1
        elif self.type == 'freerun':
            gate = 0
            exitcode = 0
            if self.timeres_dict[tempind].get(gate) != None:
                tr_hist = self.timeres_dict[tempind][gate][Vb]
            else:
                exitcode = 1

        if exitcode == 0:
            return tr_hist
        elif exitcode == 1:
            print('There is not such TR histogram!')
            return 0
