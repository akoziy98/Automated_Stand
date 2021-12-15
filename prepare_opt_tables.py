import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

from graphics_generation import graphics_generation


class prepare_table(graphics_generation):
    def __init__(self, params):
        super().__init__(params)

        self.qe_sat_list = np.arange(5, 16)
        self.DEFAULT_OPT_TABLE_NAME = "optimization_table.csv"
        self.process_bif_hist()

    def func_exp(self, x, a, b):
        return a * np.exp(b * x)

    def get_sub_slice(self, qe_sat):
        # First part about PDE
        dcr_sat, ap_sat, vb_sat, bif_sat = {}, {}, {}, {}
        if self.SPD.type == 'gated':
            for i, temp in enumerate(self.temp_list):
                dcr_sat[temp] = []
                ap_sat[temp] = []
                vb_sat[temp] = []
                bif_sat[temp] = []

                for j, vg in enumerate(self.SPD.stb_list[i]):
                    qe_list = np.array(self.SPD.qe[i][j,:])
                    qe_list[qe_list < 0] = 0
                    ind_vb2 = np.argmin(abs(qe_list * (qe_list >= qe_sat) - qe_sat))
                    ind_vb1 = np.argmin(abs(qe_list * (qe_list < qe_sat) - qe_sat))
                    qe1 = qe_list[ind_vb1]
                    qe2 = qe_list[ind_vb2]

                    if ind_vb1 != ind_vb2:
                        dcr1 = self.SPD.dcr[i][j,ind_vb1]
                        dcr2 = self.SPD.dcr[i][j,ind_vb2]
                        dcr_sat[temp].append((dcr2 - dcr1) / (qe2 - qe1) * (qe_sat - qe1) + dcr1)

                        ap1 = self.SPD.ap_list[i][j,ind_vb1]
                        ap2 = self.SPD.ap_list[i][j,ind_vb2]
                        ap_new_val = (ap2 - ap1) / (qe2 - qe1) * (qe_sat - qe1) + ap1
                        if ap_new_val > 50:
                            ap_new_val = 50
                        if ap_new_val < 0:
                            ap_new_val = 0
                        ap_sat[temp].append(ap_new_val)

                        vb1 = self.SPD.grid[i][j,ind_vb1]
                        vb2 = self.SPD.grid[i][j,ind_vb2]
                        vb_sat[temp].append((vb2 - vb1) / (qe2 - qe1) * (qe_sat - qe1) + vb1)
                    else:
                        qe1 = 0
                        dcr1 = 0
                        dcr2 = self.SPD.dcr[i][j, ind_vb2]
                        dcr_sat[temp].append((dcr2 - dcr1) / (qe2 - qe1) * (qe_sat - qe1) + dcr1)

                        ap1 = 0
                        ap2 = self.SPD.ap_list[i][j, ind_vb2]
                        ap_sat[temp].append((ap2 - ap1) / (qe2 - qe1) * (qe_sat - qe1) + ap1)

                        vb1 = self.SPD.grid[i][j, ind_vb1] - 1.5
                        vb2 = self.SPD.grid[i][j, ind_vb2]
                        vb_sat[temp].append((vb2 - vb1) / (qe2 - qe1) * (qe_sat - qe1) + vb1)

                    # analyze bif
                    x = []
                    y = []
                    for k, v in self.bif_postproc[temp][vg].items():
                        x.append(k)
                        y.append(v)
                    x_min = min(self.SPD.grid[i][j, :])
                    popt, pcov = curve_fit(self.func_exp, np.array(x) - x_min, y, method='lm', maxfev=1800, p0=[1, 1])
                    bif_sat[temp].append(self.func_exp(vb_sat[temp][-1] - x_min, *popt))
                    # x_list = np.linspace(min(self.SPD.grid[i][j, :]), max(self.SPD.grid[i][j, :]), 100)
                    # y_list = self.func_exp(x_list - x_min, *popt)
                    # plt.plot(x, y)
                    # plt.plot(x_list, y_list, "--")
                    #plt.scatter(vb_sat[temp][-1], bif_sat[temp][-1], marker="o", color="magenta")
                    # plt.show()

                # plt.plot(vb_sat[temp], bif_sat[temp], linestyle="-", marker="o")
                # plt.show()
                # end analyze bif

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

        return dcr_sat, ap_sat, vb_sat, bif_sat

    def get_slices(self):
        self.dcr_sat = {}
        self.ap_sat = {}
        self.vb_sat = {}
        self.bif_sat = {}

        for i, qe in enumerate(self.qe_sat_list):
            dcr_sat_hlp, ap_sat_hlp, vb_sat_hlp, bif_sat_hlp  = self.get_sub_slice(qe)
            self.dcr_sat[qe] = dcr_sat_hlp
            self.ap_sat[qe] = ap_sat_hlp
            self.vb_sat[qe] = vb_sat_hlp
            self.bif_sat[qe] = bif_sat_hlp

    def find_bif(self, hist):
        min_bif_10 = 1e9
        for j in range(2000, 3000 - 320):
            sm = (np.sum(hist[j: j + 320]) /
                  np.sum(hist[j + 320: j + 640]))
            if sm < min_bif_10:
                min_bif_10 = sm
        return 100 * min_bif_10

    def process_bif_hist(self):
        self.bif_postproc = {}
        for ind_t, t in enumerate(self.temp_list):
            self.bif_postproc[t] = {}
            for key_vg in self.SPD.timeres_dict_bif[ind_t].keys():
                self.bif_postproc[t][key_vg] = {}
                for key_vb in self.SPD.timeres_dict_bif[ind_t][key_vg].keys():
                    hist = self.SPD.timeres_dict_bif[ind_t][key_vg][key_vb]
                    self.bif_postproc[t][key_vg][key_vb] = self.find_bif(hist)

    def create_opt_table(self):
        self.get_slices()
        self.df_opt = pd.DataFrame(columns=["t", "vb", "vg", "qe", "dcr", "dt", "ap", "bif"])

        for ind_qe, qe in enumerate(self.qe_sat_list):
            for ind_t, t in enumerate(self.temp_list):
                for ind_vg, vg in enumerate(self.SPD.stb_list[ind_t]):
                    dict_to_append = {"t" : t,
                                      "vb" : self.vb_sat[qe][t][ind_vg],
                                      "vg" : vg,
                                      "qe" : qe / 100,
                                      "dcr" : self.dcr_sat[qe][t][ind_vg],
                                      "dt" : np.mean(self.SPD.dt_list[ind_t][ind_vg]) / 1e6,
                                      "ap" : self.ap_sat[qe][t][ind_vg],
                                      "bif" : self.bif_sat[qe][t][ind_vg] / 100}
                    self.df_opt = self.df_opt.append(dict_to_append, ignore_index=True)


        #print(self.df_opt)
        self.df_opt.to_csv(f"{os.path.join(self.dir_name, self.fold, self.DEFAULT_OPT_TABLE_NAME)}")






#Params for Rodion
# params = {}
# params['name'] = 'Rodion'
# params['date'] = '21.07.2021'
# params['dir_name'] = '!SPDStandResults//Rodion//21_07_21__18_14'
# params['type'] = 'gated'
#
# tst_opt_table = prepare_table(params)
# tst_opt_table.create_opt_table()