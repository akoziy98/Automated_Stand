import numpy as np
import os
from matplotlib import cm
import pandas as pd
from scipy.optimize import curve_fit
import warnings
warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning)
from abc import ABC, abstractmethod
import shutil
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams.update({'font.size': 22})

import afterpulse as aft
from time_resolution import time_resolution
import module_pde


DEFAULT_GATED_CHECK = "Grid_Stb_Code.csv"
DEFAULT_BATTERFLY_CHECK = "Grid_DT.csv"
DEFAULT_TEMP_NAME = 'Grid_Temperature.csv'
DEFAULT_GRID_DT = "Grid_DT.csv"
DEFAULT_GRID_INSENS = "Grid_Insens.csv"
DEFAULT_DEAD_TIME_FOLDER = "DT"
DEFAULT_DEAD_TIME_NAME = "DT.csv"

DEFAULT_DCR_SIGMA_NAME = "DCRSigma.csv"
DEFAULT_DCR_SIGMA_BUTTERFLY_NAME = "DCR_Sigma.csv"
DEFAULT_PDE_SIGMA_NAME = "QESigma.csv"
DEFAULT_PDE_SIGMA_BUTTERFLY_NAME = "QE_Sigma.csv"

class SPD():
    def __init__(self, dir_name):
        self.dir_name = dir_name

        self.is_freerun = True
        self.is_gated = False
        self.is_batterfly = False

        for fn in os.listdir(dir_name):
            if fn == DEFAULT_GATED_CHECK:
                self.is_freerun = False
                self.is_gated = True
            if fn == DEFAULT_BATTERFLY_CHECK:
                self.is_batterfly = True

        if self.is_batterfly:
            if self.is_gated:
                self.spd = SPDButterfly(dir_name, "gated")
            elif self.is_freerun:
                self.spd = SPDButterfly(dir_name, "freerun")
        elif self.is_gated:
            self.spd = SPDGated(dir_name)
        elif self.is_freerun:
            self.spd = SPDFreerun(dir_name)

    def get_spd(self):
        return self.spd

    def get_type(self):
        return self.is_batterfly, self.is_gated, self.is_freerun


def read_single_row_tables(filename):
    array = []
    with open(filename) as fn:
        for line in fn:
            line = line.strip().split(';')
            for el in line:
                array.append(float(el))

            break

    return array


class SPDCommon(ABC):
    def __init__(self, dir_name):
        self.dir_name = dir_name
        temp_filename = os.path.join(dir_name, DEFAULT_TEMP_NAME)
        self.temp_list = read_single_row_tables(temp_filename)

    def func_exp(self, x, a, b):
        return a * np.exp(b * x)

    @abstractmethod
    def load_spd(self):
        pass

    @abstractmethod
    def slice_return(self):
        pass

    @abstractmethod
    def timeres_return(self):
        pass


class SPDButterfly(SPDCommon):
    def __init__(self, dir_name, type):
        super().__init__(dir_name)
        self.type = "butterfly"
        self.type_of_spd = type
        #copy files to DT folder
        self.copy_required_files()
        self.load_spd()

    def copy_required_files(self):
        root_path = os.getcwd()
        main_data_folder_path = os.path.join(root_path, self.dir_name)
        temp_main_folder_filename = os.path.join(main_data_folder_path, DEFAULT_TEMP_NAME)
        if self.type_of_spd == "gated":
            gates_main_folder_filename = os.path.join(main_data_folder_path, DEFAULT_GATED_CHECK)

        for file_name in os.listdir(self.dir_name):
            if file_name.startswith('DT_'):
                path = os.path.join(main_data_folder_path, file_name, DEFAULT_TEMP_NAME)
                shutil.copy(temp_main_folder_filename, path)
                if self.type_of_spd == "gated":
                    path = os.path.join(main_data_folder_path, file_name, DEFAULT_GATED_CHECK)
                    shutil.copy(gates_main_folder_filename, path)


    def load_spd(self):
        self.spd_dict = {}
        self.path_dict = {}

        grid_insens_filename = os.path.join(self.dir_name, DEFAULT_GRID_INSENS)
        grid_dt_filename = os.path.join(self.dir_name, DEFAULT_GRID_DT)
        dead_time_filename = os.path.join(self.dir_name, DEFAULT_DEAD_TIME_FOLDER, DEFAULT_DEAD_TIME_NAME)

        self.grid_insens = read_single_row_tables(grid_insens_filename)
        self.grid_dt = read_single_row_tables(grid_dt_filename)
        self.grid_dt = [round(el) for el in self.grid_dt]
        self.dead_time_list = read_single_row_tables(dead_time_filename)
        self.dead_time_list = 1e6 * np.array(self.dead_time_list)

        print(self.grid_dt)

        for ind, dead_time in enumerate(self.dead_time_list):
            pass

        for file_name in os.listdir(self.dir_name):
            #if file_name.startswith('DT_'):
            dt_code = file_name.lower().replace("dt_0x", "")
            try:
                dt_code = int(dt_code, 16)
            except:
                dt_code = None

            if dt_code is not None:
                ind_of_dt_code = self.grid_dt.index(dt_code)
                dead_time = self.dead_time_list[ind_of_dt_code]
                path = os.path.join(self.dir_name, file_name)
                print(f"load butterfly path {path}")
                self.path_dict[dead_time] = file_name

                if self.type_of_spd == "gated":
                    self.spd_dict[dead_time] = SPDGated(path)
                if self.type_of_spd == "freerun":
                    self.spd_dict[dead_time] = SPDFreerun(path)

        self.spd_deadtimes = self.spd_dict.keys()

    def slice_return(self, plot_index=0, tempind=0, ind=0, spd_dead_time = 0):
        spd = self.get_specific_spd(spd_dead_time)
        if spd is not None:
            return spd.slice_return(plot_index, tempind, ind)

    def timeres_return(self, tempind=0, gateind=0, Vb=0, spd_dead_time = 0):
        spd = self.get_specific_spd(spd_dead_time)
        if spd is not None:
            return spd.timeres_return(tempind, gateind, Vb)

    def get_specific_spd(self, spd_dead_time = 0):
        spd = self.spd_dict.get(spd_dead_time)
        if spd is None:
            print("there is no such spd dead time in butterfly spd")
            return None
        else:
            return spd

class SPDFreerun(SPDCommon):
    def __init__(self, dir_name):
        super().__init__(dir_name)
        self.type = "freerun"
        self.load_spd()

    def load_spd(self):
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

            try:
                with open(os.path.join(self.dir_name, temp_fold[i], file_sigma_dcr)) as fn:
                    for line in fn:
                        line = line.strip().split(';')
                        sigma_dcr[i] = [float(line[k]) for k in range(len(line))]
            except:
                file_sigma_dcr = DEFAULT_DCR_SIGMA_BUTTERFLY_NAME
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

    def slice_return(self, plot_index = 0, tempind = 0, ind = 0):
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


class SPDGated(SPDCommon):
    def __init__(self, dir_name):
        super().__init__(dir_name)
        stb_list_filename = os.path.join(self.dir_name, DEFAULT_GATED_CHECK)
        self.stb_list = read_single_row_tables(stb_list_filename)
        self.type = "gated"
        self.load_spd()

    def load_spd(self):
        n = len(self.temp_list)

        try:
            self.gate_voltage = pd.read_excel(os.path.join(self.dir_name, 'gate_voltage.xlsx'))
        except:
            self.gate_voltage = 0
            print('there is no gate voltage tables')
        temp_fold = ['T_' + str(el) for el in self.temp_list]

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
                    vbias_small = hv_dict[i][self.stb_list[j]]
                    AP_small = ap_dict[i][self.stb_list[j]]
                    popt, pcov = curve_fit(self.func_exp, vbias_small, AP_small, method='lm', maxfev=1800, p0=[1, 1])
                    req_hvlist = np.array(grid[i][j][:])
                    req_hvlist = req_hvlist - req_hvlist[0]
                    #print(req_hvlist)
                    ap_list[i].append(self.func_exp(req_hvlist, *popt))
                    DT_aver = np.mean(dt_dict[i][self.stb_list[j]])
                    dt_list[i].append([DT_aver for i in range(len(req_hvlist))])
                except:
                    print('There is no required AP histogram for T = ' + str(self.temp_list[i]) + ' and gate = ' + str(hex(int(self.stb_list[j]))))
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
            try:
                with open(os.path.join(self.dir_name, temp_fold[i], file_sigma_dcr)) as fn:
                    ind = 0
                    for line in fn:
                        line = line.strip().split(';')
                        if line[0] != 'NaN' and ind < n_apgates[i]:
                            sigma_dcr[i].append([float(line[k]) for k in range(len(line))])
                            ind += 1
            except:
                file_sigma_dcr = DEFAULT_DCR_SIGMA_BUTTERFLY_NAME
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
            mesh_strobe.append([self.stb_list[:n_apgates[i]] for j in range(NN)])
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
        self.stb_list_new = []
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
            self.stb_list_new.append(np.array(self.stb_list[:n_apgates[i]]))
            self.ap_list.append(100 * np.array(ap_list[i]))
            self.snr.append(self.qe[i] / self.dcr[i])
        #Save to the class:

        self.stb_list = self.stb_list_new

        self.timeres_dict = timeres_dict
        self.hv_timeres_dict = hv_timeres_dict
        self.timeres_dict_bif = timeres_dict_bif
        self.hv_timeres_dict_bif = hv_timeres_dict_bif

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

    def timeres_return(self, tempind = 0, gateind = 0, Vb = 0):
        gate = self.stb_list[tempind][gateind]
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
