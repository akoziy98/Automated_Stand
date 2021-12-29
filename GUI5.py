import argparse
import sys
import os

import numpy as np

from matplotlib import cm
from matplotlib.pyplot import figure
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from mpl_toolkits.mplot3d import Axes3D

from PyQt5.QtCore import Qt, pyqtSlot
from PyQt5 import QtWidgets, QtCore
from PyQt5.QtWidgets import QApplication as Application, QWidget, QPushButton as Button, QFileDialog, QComboBox
from PyQt5.QtWidgets import QLabel as Label, QGridLayout, QDesktopWidget, QSlider

import spd_analyze
import design5
import generate_pdf
import prepare_opt_tables

DEFAULT_GATED_NAME = "gated"
DEFAULT_FREERUN_NAME = "freerun"
DEFAULT_BUTTERFLY_NAME = "butterfly"



class ProgramGUI(QtWidgets.QMainWindow, design5.Ui_MainWindow):

    def __init__(self):
        super().__init__()
        self.setupUi(self)
        # GUI window specific values
        self.title = 'SPD view programm'
        self.setWindowTitle(self.title)
        self.left = 100
        self.top = 100
        self.width = 1600
        self.height = 1000
        # Other object specific values
        self.plot_status = u'a'
        self.X_plot_val = None
        self.Y_plot_val = None
        self.Z_plot_val = None

        self.directory = '!SPDStandResults\\'
        # initialize UI
        self.init_ui()
        #self.horizontalLayout.setStretchFactor(self.horizontalLayout.itemAt(0), 1)
        #self.horizontalLayout.setStretchFactor(self.horizontalLayout.itemAt(1), 2)

    def init_ui(self):
        #
        # Setup Window Title and geometry
        #
        '''
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)
        #self.center()
        '''
        self.combobox_plotnumber.addItem('DCR')
        self.combobox_plotnumber.addItem('QE')
        self.combobox_plotnumber.addItem('QE / DCR')
        self.combobox_plotnumber.addItem('AP')
        self.combobox_plotnumber.addItem('DT')

        self.combobox_2dplot.addItem(r'Vg = const')

        self.directory_view()

        self.combobox_name.activated.connect(self.date_view)
        self.combobox_date.activated.connect(self.download_spd)
        self.combobox_temp.activated.connect(self.shift_slider2_newplot)
        self.combobox_plotnumber.activated.connect(self.shift_slider2_newplot)
        self.combobox_deadtime.activated.connect(self.choose_deadtime_setting)
        self.slider1.valueChanged.connect(self.shift_slider2_newplot)
        self.slider2.valueChanged.connect(self.shift_slider2)
        self.report_button.clicked.connect(self.create_report)


    def create_report(self):
        report = generate_pdf.PdfDocument(self.params)
        report.generate_pdf()

        try:
            opt_table = prepare_opt_tables.prepare_table(self.params)
            opt_table.create_opt_table()
        except:
            pass


    def center(self):
        """centers the window on the screen"""
        screen = QDesktopWidget().screenGeometry()
        size = self.geometry()

    def directory_view(self):
        for file_name in os.listdir(self.directory):  # для каждого файла в директории
            self.combobox_name.addItem(file_name)
        self.date_view(0)

    def date_view(self, index):
        if self.combobox_date.count() > 0:
            self.combobox_date.clear()
            self.combobox_temp.clear()
        for file_name in os.listdir(self.directory + self.combobox_name.itemText(index) + '/'):
            self.combobox_date.addItem(file_name)
        self.download_spd()
        #self.slice = self.shift_slider(0)
        #self.plot3d_click()


    def download_spd(self):
        #self.directory = QFileDialog.getExistingDirectory(self, "Выберите папку")
        # открыть диалог выбора директории и установить значение переменной
        # равной пути к выбранной директории
        self.params = {}
        dir_name = os.getcwd() + '\\' + self.directory + self.combobox_name.currentText() + '\\' + self.combobox_date.currentText()

        self.params['dir_name'] = dir_name
        self.params['name'] = self.combobox_name.currentText()
        date = self.combobox_date.currentText()
        date = date[:date.index('__')]
        date_day = date[date.rindex('_') + 1:]
        date_year = '20' + date[: date.index('_')]
        date_month = date[date.index('_') + 1:date.rindex('_')]
        date = date_day + '.' + date_month + '.' + date_year
        self.params['date'] = date

        #print(os.getcwd())
        print(dir_name)
        self.download_files(dir_name)
        self.shift_slider2_newplot()

    def download_files(self, dir_name):
        try:
            del self.SPD
            print('Old SPD is deleted!')
            self.combobox_temp.clear()
            self.combobox_deadtime.clear()
        except:
            print('SPD is loaded!')

        self.SPD = spd_analyze.SPD(dir_name)
        self.SPD = self.SPD.get_spd()
        self.params["type"] = self.SPD.type

        #print(self.SPD.temp_list)
        if self.SPD.type == DEFAULT_GATED_NAME or self.SPD.type == DEFAULT_FREERUN_NAME:
            self.SPD_butterfly = None
        elif self.SPD.type == DEFAULT_BUTTERFLY_NAME:
            self.SPD_butterfly = self.SPD
            deadtime_key = list(self.SPD.spd_dict.keys())[0]
            self.SPD = self.SPD_butterfly.spd_dict[deadtime_key]
            for ind, el in enumerate(self.SPD_butterfly.spd_deadtimes):
                self.combobox_deadtime.addItem("DT = " + str(el))

        self.set_start_setup_after_loading()

    def choose_deadtime_setting(self):
        deadtime_ind = self.combobox_deadtime.currentIndex()
        list_of_keys = list(self.SPD_butterfly.spd_dict.keys())
        cur_key = list_of_keys[deadtime_ind]
        self.SPD = self.SPD_butterfly.spd_dict[cur_key]
        self.combobox_temp.clear()

        self.set_start_setup_after_loading()
        self.shift_slider2_newplot()

    def set_start_setup_after_loading(self):
        self.ntemps = len(self.SPD.temp_list)
        for i in range(len(self.SPD.temp_list)):
            self.combobox_temp.addItem('T = ' + str(self.SPD.temp_list[i]))

        self.label_type.setText('SPD type: ' + self.SPD.type)
        self.slider1.setMinimum(0)
        temp_index = self.combobox_temp.currentIndex()
        if self.SPD.type == 'gated':
            self.slider1.setMaximum(len(self.SPD.stb_list[temp_index]) - 1)
        elif self.SPD.type == 'freerun':
            self.slider1.setMaximum(0)

        self.slider1.setSingleStep(1)
        self.maxdcr = np.max([np.max(self.SPD.dcr[i]) for i in range(len(self.SPD.temp_list))])
        self.maxqe = np.max([np.max(self.SPD.qe[i]) for i in range(len(self.SPD.temp_list))])
        self.maxap = np.max([np.max(self.SPD.ap_list[i]) for i in range(len(self.SPD.temp_list))])
        self.maxdt = np.max([np.max(self.SPD.dt_list[i]) for i in range(len(self.SPD.temp_list))])
        self.maxqedcr = -100
        for i in range(self.ntemps):
            qedcr = self.SPD.snr[i]
            maxqedcr = np.max(qedcr)
            if maxqedcr > self.maxqedcr:
                self.maxqedcr = maxqedcr


    def shift_slider2_newplot(self):
        self.shift_slider(is_newplot=1)

    def shift_slider2(self, index=0):
        self.shift_slider(is_newplot=0)
        self.fill_labels()

        plot_index = self.combobox_plotnumber.currentIndex()
        temp_index = self.combobox_temp.currentIndex()
        slice_index = self.slider1.value()
        slice2_index = self.slider2.value()


        Vb = self.grid_slider2[slice2_index]

        if self.SPD.type == 'gated':
            gate = self.SPD.stb_list[temp_index][slice_index]
            Vb_ind = np.where(self.SPD.grid[temp_index][slice_index] == Vb)

            self.label_slider2.setText('Vb = ' + str(Vb))
            if Vb in self.keys_tr:
                tr_hist = self.SPD.timeres_return(temp_index, slice_index, Vb)
            else:
                tr_hist = 0

            if type(tr_hist) != int:
                self.plot_container.plot_point([gate, Vb, self.slice[2][Vb_ind]])
                self.twodplot_container.plot_slice(Vb)
                self.trplot_container.draw_graph(
                    [np.arange(0, len(tr_hist)), 0.01 * np.arange(0, len(tr_hist)), tr_hist],
                    title=r'TR. $V_g = \ $' + str(gate) + r'. $V_b$ = \ ' + str(Vb),
                    style='-', labels=['$t$, ns', '$Counts$'], limit=None)
            else:
                self.plot_container.plot_point([gate, Vb, self.slice[2][Vb_ind]])
                self.twodplot_container.plot_slice(Vb)
                self.trplot_container.draw_graph(0, title=r'TR. $V_g = \ $' + str(gate) + r'. $V_b$ = \ ' + str(Vb),
                                                 style='-', labels=['$t$, ns', '$Counts$'], limit=None)

        elif self.SPD.type == 'freerun':
            gate = 0
            Vb_ind = np.where(self.SPD.grid[temp_index, :] == Vb)

            self.label_slider2.setText('Vb = ' + str(Vb))
            if Vb in self.keys_tr:
                tr_hist = self.SPD.timeres_return(temp_index, 0, Vb)
            else:
                tr_hist = 0

            if type(tr_hist) != int:
                self.twodplot_container.plot_slice(Vb)
                self.trplot_container.draw_graph(
                    [np.arange(0, len(tr_hist)), 0.01 * np.arange(0, len(tr_hist)), tr_hist],
                    title=r'TR. $V_g = \ $' + str(gate) + r'. $V_b$ = \ ' + str(Vb),
                    style='-', labels=['$t$, ns', '$Counts$'], limit=None)
            else:
                self.twodplot_container.plot_slice(Vb)
                self.trplot_container.draw_graph(0, title=r'TR. $V_g = \ $' + str(gate) + r'. $V_b$ = \ ' + str(Vb),
                                                 style='-', labels=['$t$, ns', '$Counts$'], limit=None)


        width_cur = self.centralwidget.width()
        height_cur = self.centralwidget.height()
        if width_cur != 100 and height_cur != 30:
            self.centralwidget.setFixedWidth(width_cur - 1)
            self.centralwidget.setFixedWidth(width_cur)
            self.centralwidget.setFixedHeight(height_cur - 1)
            self.centralwidget.setFixedHeight(height_cur)
            self.center()

    def shift_slider(self, index = 0, is_newplot = 0):
        if self.SPD.type == 'gated':
            #set slider 1 parameters
            plot_index = self.combobox_plotnumber.currentIndex()
            temp_index = self.combobox_temp.currentIndex()
            slice_index = self.slider1.value()
            gate = self.SPD.stb_list[temp_index][slice_index]
            self.label_slider1.setText('Vg = ' + str(gate))
            self.slice = self.SPD.slice_return(tempind=temp_index, ind=slice_index, plot_index=plot_index)

            #set slider2 parameters
            self.slider2.setMinimum(0)
            if (is_newplot == 1) and (self.SPD.timeres_dict[temp_index].get(gate) != None):
                self.keys_tr = list(self.SPD.timeres_dict[temp_index][gate].keys())
                self.grid_slider2 = self.slice[1]
                #print(self.keys_tr)
                self.slider2.setMaximum(len(self.grid_slider2) - 1)
                self.slider2.setSingleStep(1)
                self.label_slider2.setText('Vb = ' + str(self.grid_slider2))
                #self.slider2.setValue(self.keys_tr[0])
                return self.shift_slider2()

            self.plot_container.axes.clear()
            self.plot3d_click()
            self.plot_container.plot_slice(self.slice)
            self.plot2d_click()

        elif self.SPD.type == 'freerun':
            # set slider 1 parameters
            plot_index = self.combobox_plotnumber.currentIndex()
            temp_index = self.combobox_temp.currentIndex()
            self.label_slider1.setText('Vg = 0')
            self.slice = self.SPD.slice_return(tempind=temp_index, plot_index=plot_index)

            # set slider2 parameters
            self.slider2.setMinimum(0)
            if (is_newplot == 1):
                self.keys_tr = list(self.SPD.timeres_dict[temp_index][0].keys())
                self.grid_slider2 = self.slice[0]
                # print(self.keys_tr)
                self.slider2.setMaximum(len(self.grid_slider2) - 1)
                self.slider2.setSingleStep(1)
                self.label_slider2.setText('Vb = ' + str(self.grid_slider2))
                # self.slider2.setValue(self.keys_tr[0])
                return self.shift_slider2()

            #self.plot_container.axes.clear()
            self.plot_container.clear_graphics()
            self.plot2d_click()

    def fill_labels(self):
        temp_index = self.combobox_temp.currentIndex()
        slice_index = self.slider1.value()
        slice2_index = self.slider2.value()
        Vb = self.grid_slider2[slice2_index]

        if self.SPD.type == 'gated':
            gate = self.SPD.stb_list[temp_index][slice_index]
            Vb_ind = np.where(self.SPD.grid[temp_index][ slice_index] == Vb)
            # Efficiency parameters
            pde_cur = self.SPD.qe[temp_index][slice_index][Vb_ind]
            pde_cur = round(pde_cur[0], 2)
            self.QE_label.setText('PDE = ' + str(pde_cur) + ' %')
            dcr_cur = self.SPD.dcr[temp_index][slice_index][Vb_ind]
            dcr_cur = round(dcr_cur[0])
            self.DCR_label.setText('DCR = ' + str(dcr_cur) + ' Hz')
            pde_dcr_cur = round(pde_cur / dcr_cur, 3)
            self.QE_DCR_label.setText('QE / DCR = ' + str(pde_dcr_cur))
            ap_cur = self.SPD.ap_list[temp_index][slice_index][Vb_ind]
            ap_cur = round(ap_cur[0], 2)
            self.AP_label.setText('AP = ' + str(ap_cur) + ' %')
            dt_cur = self.SPD.dt_list[temp_index][slice_index][Vb_ind]
            dt_cur = round(dt_cur[0], 2)
            self.DT_label.setText('DT = ' + str(dt_cur) + ' mus')
        elif self.SPD.type == 'freerun':
            gate = 0
            Vb_ind = np.where(self.SPD.grid[temp_index] == Vb)
            # Efficiency parameters
            pde_cur = self.SPD.qe[temp_index][Vb_ind]
            pde_cur = round(pde_cur[0], 2)
            self.QE_label.setText('PDE = ' + str(pde_cur) + ' %')
            dcr_cur = self.SPD.dcr[temp_index][Vb_ind]
            dcr_cur = round(dcr_cur[0])
            self.DCR_label.setText('DCR = ' + str(dcr_cur) + ' Hz')
            pde_dcr_cur = round(pde_cur / dcr_cur, 3)
            self.QE_DCR_label.setText('QE / DCR = ' + str(pde_dcr_cur))
            ap_cur = self.SPD.ap_list[temp_index][Vb_ind]
            ap_cur = round(ap_cur[0], 2)
            self.AP_label.setText('AP = ' + str(ap_cur) + ' %')
            dt_cur = self.SPD.dt_list[temp_index][Vb_ind]
            dt_cur = round(dt_cur[0], 2)
            self.DT_label.setText('DT = ' + str(dt_cur) + ' mus')

        #Technical parameters
        self.T_label.setText('T = ' + str(self.SPD.temp_list[temp_index]) + ' C')
        self.VB_label.setText('Vb = ' + str(round(Vb, 2)) + ' V')
        self.Gate_label.setText('Vg = ' +  'dec: ' + str(int(gate)) + '|| hex:' + str(hex(int(gate)).split('x')[-1]))

    def plot2d_click(self):
        gate_index = self.slider1.value()
        plot2d_index = self.combobox_2dplot.currentIndex()
        plot_index = self.combobox_plotnumber.currentIndex()
        temp_index = self.combobox_temp.currentIndex()

        if self.SPD.type == 'gated':
            if plot2d_index == 0:
                if plot_index == 0:
                    limit = self.maxdcr
                    labels = [r'$V_b$', r'$DCR$, Hz']
                    title = r'$DCR$. $V_g =$ ' + str(self.SPD.stb_list[temp_index][gate_index])
                elif plot_index == 1:
                    limit = self.maxqe
                    labels = [r'$V_b$', r'$QE, \%$']
                    title = r'$QE$ \%. $V_g =$ ' + str(self.SPD.stb_list[temp_index][gate_index])
                elif plot_index == 2:
                    limit = self.maxqedcr
                    labels = [r'$V_b$', r'$QE / DCR$']
                    title = r'$QE / DCR$. $V_g =$ ' + str(self.SPD.stb_list[temp_index][gate_index])
                elif plot_index == 3:
                    limit = self.maxap
                    labels = [r'$V_b$', r'$AP, \ \%$']
                    title = r'$AP$ \%. $V_g =$ ' + str(self.SPD.stb_list[temp_index][gate_index])
                elif plot_index == 4:
                    limit = self.maxdt
                    labels = [r'$V_b$', r'$DT, \ \mu s$']
                    title = r'$DT, \ \mu s$. $V_g =$ ' + str(self.SPD.stb_list[temp_index][gate_index])

                self.twodplot_container.draw_graph(self.slice, limit=limit, labels=labels, title = title)

        elif self.SPD.type == 'freerun':
            if plot2d_index == 0:
                if plot_index == 0:
                    limit = self.maxdcr
                    labels = [r'$V_b$', r'$DCR$, Hz']
                    title = r'$DCR$. $V_g = 0$'
                elif plot_index == 1:
                    limit = self.maxqe
                    labels = [r'$V_b$', r'$QE, \%$']
                    title = r'$QE$ \%. $V_g = 0$'
                elif plot_index == 2:
                    limit = self.maxqedcr
                    labels = [r'$V_b$', r'$QE / DCR$']
                    title = r'$QE / DCR$. $V_g = 0$'
                elif plot_index == 3:
                    limit = self.maxap
                    labels = [r'$V_b$', r'$AP, \ \%$']
                    title = r'$AP$ \%. $V_g = 0$'
                elif plot_index == 4:
                    limit = self.maxdt
                    labels = [r'$V_b$', r'$DT, \ \mu s$']
                    title = r'$DT, \ \mu s$. $V_g = 0$'

                self.twodplot_container.draw_graph(self.slice, limit=limit, labels=labels, title=title)
    @pyqtSlot()
    def plot3d_click(self):
        print('try to plot!')
        temp_index = self.combobox_temp.currentIndex()
        plot_index = self.combobox_plotnumber.currentIndex()
        #print(temp_index, plot_index)
        mesh_strobe = self.SPD.mesh_strobe[temp_index]
        grid = self.SPD.grid[temp_index]

        if plot_index == 0:
            dcr = self.SPD.dcr[temp_index]
            #print(mesh_strobe, grid, dcr)
            title = 'DCR dependence for T = ' + str(self.SPD.temp_list[temp_index])
            labels = ['Gate number', 'Bias', 'DCR, Hz']
            self.plot_container.draw_graph(mesh_strobe, grid, dcr, title, labels)
        elif plot_index == 1:
            qe =  self.SPD.qe[temp_index]
            title = 'QE dependence for T = ' + str(self.SPD.temp_list[temp_index])
            labels = ['Gate number', 'Bias', 'QE, \%']
            self.plot_container.draw_graph(mesh_strobe, grid, qe, title, labels)
        elif plot_index == 2:
            qe =  self.SPD.qe[temp_index]
            dcr = self.SPD.dcr[temp_index]
            title = 'QE/DCR dependence for T = ' + str(self.SPD.temp_list[temp_index])
            labels = ['Gate number', 'Bias', 'QE/DCR']
            self.plot_container.draw_graph(mesh_strobe, grid, qe / dcr, title, labels)
        elif plot_index == 3:
            AP =  self.SPD.ap_list[temp_index]
            title = 'AP dependence for T = ' + str(self.SPD.temp_list[temp_index])
            labels = ['Gate number', 'Bias', 'AP, \%']
            self.plot_container.draw_graph(mesh_strobe, grid, AP, title, labels)
        elif plot_index == 4:
            DT =  self.SPD.dt_list[temp_index]
            title = 'DT dependence for T = ' + str(self.SPD.temp_list[temp_index])
            labels = ['Gate number', 'Bias', 'DT, mus']
            self.plot_container.draw_graph(mesh_strobe, grid, DT, title, labels)

if __name__ == '__main__':
    app = Application(sys.argv)
    gui = ProgramGUI()

    qr = gui.frameGeometry()
    cp = QDesktopWidget().availableGeometry().center()
    qr.moveCenter(cp)
    gui.move(qr.topLeft())
    app.processEvents()

    #gui.show()
    gui.showMaximized()

    exit_val = app.exec_()

    # behaviour to trigger on exit
    sys.exit(exit_val)
