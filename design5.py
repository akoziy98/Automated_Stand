# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'design5.ui'
#
# Created by: PyQt5 UI code generator 5.9.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(988, 1022)
        MainWindow.setStyleSheet("background-color: rgb(250, 250, 255);")
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout(self.centralwidget)
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.verticalLayout_2 = QtWidgets.QVBoxLayout()
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.label_5 = QtWidgets.QLabel(self.centralwidget)
        self.label_5.setMinimumSize(QtCore.QSize(160, 0))
        self.label_5.setMaximumSize(QtCore.QSize(160, 16777215))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(11)
        self.label_5.setFont(font)
        self.label_5.setObjectName("label_5")
        self.verticalLayout_2.addWidget(self.label_5)
        self.line_2 = QtWidgets.QFrame(self.centralwidget)
        self.line_2.setMinimumSize(QtCore.QSize(160, 0))
        self.line_2.setMaximumSize(QtCore.QSize(140, 16777215))
        self.line_2.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_2.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_2.setObjectName("line_2")
        self.verticalLayout_2.addWidget(self.line_2)
        self.label = QtWidgets.QLabel(self.centralwidget)
        self.label.setMinimumSize(QtCore.QSize(160, 0))
        self.label.setMaximumSize(QtCore.QSize(160, 16777215))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(11)
        self.label.setFont(font)
        self.label.setObjectName("label")
        self.verticalLayout_2.addWidget(self.label)
        self.combobox_name = QtWidgets.QComboBox(self.centralwidget)
        self.combobox_name.setMinimumSize(QtCore.QSize(160, 0))
        self.combobox_name.setMaximumSize(QtCore.QSize(160, 16777215))
        self.combobox_name.setObjectName("combobox_name")
        self.verticalLayout_2.addWidget(self.combobox_name)
        self.label_2 = QtWidgets.QLabel(self.centralwidget)
        self.label_2.setMinimumSize(QtCore.QSize(160, 0))
        self.label_2.setMaximumSize(QtCore.QSize(160, 16777215))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(11)
        self.label_2.setFont(font)
        self.label_2.setObjectName("label_2")
        self.verticalLayout_2.addWidget(self.label_2)
        self.combobox_date = QtWidgets.QComboBox(self.centralwidget)
        self.combobox_date.setMinimumSize(QtCore.QSize(160, 0))
        self.combobox_date.setMaximumSize(QtCore.QSize(160, 16777215))
        self.combobox_date.setObjectName("combobox_date")
        self.verticalLayout_2.addWidget(self.combobox_date)
        self.label_13 = QtWidgets.QLabel(self.centralwidget)
        self.label_13.setMinimumSize(QtCore.QSize(160, 0))
        self.label_13.setMaximumSize(QtCore.QSize(160, 16777215))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(11)
        self.label_13.setFont(font)
        self.label_13.setObjectName("label_13")
        self.verticalLayout_2.addWidget(self.label_13)
        self.combobox_deadtime = QtWidgets.QComboBox(self.centralwidget)
        self.combobox_deadtime.setMinimumSize(QtCore.QSize(160, 0))
        self.combobox_deadtime.setMaximumSize(QtCore.QSize(160, 16777215))
        self.combobox_deadtime.setObjectName("combobox_deadtime")
        self.verticalLayout_2.addWidget(self.combobox_deadtime)
        self.label_3 = QtWidgets.QLabel(self.centralwidget)
        self.label_3.setMinimumSize(QtCore.QSize(160, 0))
        self.label_3.setMaximumSize(QtCore.QSize(160, 16777215))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(11)
        self.label_3.setFont(font)
        self.label_3.setObjectName("label_3")
        self.verticalLayout_2.addWidget(self.label_3)
        self.combobox_temp = QtWidgets.QComboBox(self.centralwidget)
        self.combobox_temp.setMinimumSize(QtCore.QSize(160, 0))
        self.combobox_temp.setMaximumSize(QtCore.QSize(160, 16777215))
        self.combobox_temp.setObjectName("combobox_temp")
        self.verticalLayout_2.addWidget(self.combobox_temp)
        self.label_4 = QtWidgets.QLabel(self.centralwidget)
        self.label_4.setMinimumSize(QtCore.QSize(160, 0))
        self.label_4.setMaximumSize(QtCore.QSize(160, 16777215))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(11)
        self.label_4.setFont(font)
        self.label_4.setObjectName("label_4")
        self.verticalLayout_2.addWidget(self.label_4)
        self.combobox_plotnumber = QtWidgets.QComboBox(self.centralwidget)
        self.combobox_plotnumber.setMinimumSize(QtCore.QSize(160, 0))
        self.combobox_plotnumber.setMaximumSize(QtCore.QSize(160, 16777215))
        self.combobox_plotnumber.setObjectName("combobox_plotnumber")
        self.verticalLayout_2.addWidget(self.combobox_plotnumber)
        self.line = QtWidgets.QFrame(self.centralwidget)
        self.line.setMinimumSize(QtCore.QSize(160, 0))
        self.line.setMaximumSize(QtCore.QSize(140, 16777215))
        self.line.setFrameShape(QtWidgets.QFrame.HLine)
        self.line.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line.setObjectName("line")
        self.verticalLayout_2.addWidget(self.line)
        self.label_6 = QtWidgets.QLabel(self.centralwidget)
        self.label_6.setMinimumSize(QtCore.QSize(160, 0))
        self.label_6.setMaximumSize(QtCore.QSize(160, 16777215))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(11)
        self.label_6.setFont(font)
        self.label_6.setObjectName("label_6")
        self.verticalLayout_2.addWidget(self.label_6)
        self.line_3 = QtWidgets.QFrame(self.centralwidget)
        self.line_3.setMinimumSize(QtCore.QSize(160, 0))
        self.line_3.setMaximumSize(QtCore.QSize(140, 16777215))
        self.line_3.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_3.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_3.setObjectName("line_3")
        self.verticalLayout_2.addWidget(self.line_3)
        self.label_7 = QtWidgets.QLabel(self.centralwidget)
        self.label_7.setMinimumSize(QtCore.QSize(160, 0))
        self.label_7.setMaximumSize(QtCore.QSize(160, 16777215))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(11)
        self.label_7.setFont(font)
        self.label_7.setObjectName("label_7")
        self.verticalLayout_2.addWidget(self.label_7)
        self.combobox_2dplot = QtWidgets.QComboBox(self.centralwidget)
        self.combobox_2dplot.setMinimumSize(QtCore.QSize(160, 0))
        self.combobox_2dplot.setMaximumSize(QtCore.QSize(160, 16777215))
        self.combobox_2dplot.setObjectName("combobox_2dplot")
        self.verticalLayout_2.addWidget(self.combobox_2dplot)
        self.label_slider1 = QtWidgets.QLabel(self.centralwidget)
        self.label_slider1.setMinimumSize(QtCore.QSize(160, 0))
        self.label_slider1.setMaximumSize(QtCore.QSize(160, 16777215))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(11)
        self.label_slider1.setFont(font)
        self.label_slider1.setText("")
        self.label_slider1.setObjectName("label_slider1")
        self.verticalLayout_2.addWidget(self.label_slider1)
        self.slider1 = QtWidgets.QScrollBar(self.centralwidget)
        self.slider1.setMinimumSize(QtCore.QSize(160, 0))
        self.slider1.setMaximumSize(QtCore.QSize(140, 100))
        self.slider1.setOrientation(QtCore.Qt.Horizontal)
        self.slider1.setObjectName("slider1")
        self.verticalLayout_2.addWidget(self.slider1)
        self.line_4 = QtWidgets.QFrame(self.centralwidget)
        self.line_4.setMinimumSize(QtCore.QSize(160, 0))
        self.line_4.setMaximumSize(QtCore.QSize(140, 16777215))
        self.line_4.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_4.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_4.setObjectName("line_4")
        self.verticalLayout_2.addWidget(self.line_4)
        self.label_10 = QtWidgets.QLabel(self.centralwidget)
        self.label_10.setMinimumSize(QtCore.QSize(160, 0))
        self.label_10.setMaximumSize(QtCore.QSize(160, 16777215))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(11)
        self.label_10.setFont(font)
        self.label_10.setObjectName("label_10")
        self.verticalLayout_2.addWidget(self.label_10)
        self.label_slider2 = QtWidgets.QLabel(self.centralwidget)
        self.label_slider2.setMinimumSize(QtCore.QSize(160, 0))
        self.label_slider2.setMaximumSize(QtCore.QSize(160, 16777215))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(11)
        self.label_slider2.setFont(font)
        self.label_slider2.setText("")
        self.label_slider2.setObjectName("label_slider2")
        self.verticalLayout_2.addWidget(self.label_slider2)
        self.slider2 = QtWidgets.QScrollBar(self.centralwidget)
        self.slider2.setMinimumSize(QtCore.QSize(160, 0))
        self.slider2.setMaximumSize(QtCore.QSize(140, 16777215))
        self.slider2.setOrientation(QtCore.Qt.Horizontal)
        self.slider2.setObjectName("slider2")
        self.verticalLayout_2.addWidget(self.slider2)
        self.line_8 = QtWidgets.QFrame(self.centralwidget)
        self.line_8.setMinimumSize(QtCore.QSize(160, 0))
        self.line_8.setMaximumSize(QtCore.QSize(160, 16777215))
        self.line_8.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_8.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_8.setObjectName("line_8")
        self.verticalLayout_2.addWidget(self.line_8)
        self.label_type = QtWidgets.QLabel(self.centralwidget)
        self.label_type.setMinimumSize(QtCore.QSize(160, 0))
        self.label_type.setMaximumSize(QtCore.QSize(160, 16777215))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(11)
        self.label_type.setFont(font)
        self.label_type.setObjectName("label_type")
        self.verticalLayout_2.addWidget(self.label_type)
        self.line_5 = QtWidgets.QFrame(self.centralwidget)
        self.line_5.setMinimumSize(QtCore.QSize(160, 0))
        self.line_5.setMaximumSize(QtCore.QSize(140, 16777215))
        self.line_5.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_5.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_5.setObjectName("line_5")
        self.verticalLayout_2.addWidget(self.line_5)
        self.label_11 = QtWidgets.QLabel(self.centralwidget)
        self.label_11.setMinimumSize(QtCore.QSize(160, 0))
        self.label_11.setMaximumSize(QtCore.QSize(160, 16777215))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(11)
        self.label_11.setFont(font)
        self.label_11.setObjectName("label_11")
        self.verticalLayout_2.addWidget(self.label_11)
        self.QE_label = QtWidgets.QLabel(self.centralwidget)
        self.QE_label.setMinimumSize(QtCore.QSize(160, 0))
        self.QE_label.setMaximumSize(QtCore.QSize(160, 16777215))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(11)
        self.QE_label.setFont(font)
        self.QE_label.setObjectName("QE_label")
        self.verticalLayout_2.addWidget(self.QE_label)
        self.DCR_label = QtWidgets.QLabel(self.centralwidget)
        self.DCR_label.setMinimumSize(QtCore.QSize(160, 0))
        self.DCR_label.setMaximumSize(QtCore.QSize(160, 16777215))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(11)
        self.DCR_label.setFont(font)
        self.DCR_label.setObjectName("DCR_label")
        self.verticalLayout_2.addWidget(self.DCR_label)
        self.QE_DCR_label = QtWidgets.QLabel(self.centralwidget)
        self.QE_DCR_label.setMinimumSize(QtCore.QSize(160, 0))
        self.QE_DCR_label.setMaximumSize(QtCore.QSize(160, 16777215))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(11)
        self.QE_DCR_label.setFont(font)
        self.QE_DCR_label.setObjectName("QE_DCR_label")
        self.verticalLayout_2.addWidget(self.QE_DCR_label)
        self.AP_label = QtWidgets.QLabel(self.centralwidget)
        self.AP_label.setMinimumSize(QtCore.QSize(160, 0))
        self.AP_label.setMaximumSize(QtCore.QSize(160, 16777215))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(11)
        self.AP_label.setFont(font)
        self.AP_label.setObjectName("AP_label")
        self.verticalLayout_2.addWidget(self.AP_label)
        self.DT_label = QtWidgets.QLabel(self.centralwidget)
        self.DT_label.setMinimumSize(QtCore.QSize(160, 0))
        self.DT_label.setMaximumSize(QtCore.QSize(160, 16777215))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(11)
        self.DT_label.setFont(font)
        self.DT_label.setObjectName("DT_label")
        self.verticalLayout_2.addWidget(self.DT_label)
        self.line_7 = QtWidgets.QFrame(self.centralwidget)
        self.line_7.setMinimumSize(QtCore.QSize(160, 0))
        self.line_7.setMaximumSize(QtCore.QSize(140, 16777215))
        self.line_7.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_7.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_7.setObjectName("line_7")
        self.verticalLayout_2.addWidget(self.line_7)
        self.label_12 = QtWidgets.QLabel(self.centralwidget)
        self.label_12.setMinimumSize(QtCore.QSize(160, 0))
        self.label_12.setMaximumSize(QtCore.QSize(160, 16777215))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(11)
        self.label_12.setFont(font)
        self.label_12.setObjectName("label_12")
        self.verticalLayout_2.addWidget(self.label_12)
        self.T_label = QtWidgets.QLabel(self.centralwidget)
        self.T_label.setMinimumSize(QtCore.QSize(160, 0))
        self.T_label.setMaximumSize(QtCore.QSize(160, 16777215))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(11)
        self.T_label.setFont(font)
        self.T_label.setObjectName("T_label")
        self.verticalLayout_2.addWidget(self.T_label)
        self.VB_label = QtWidgets.QLabel(self.centralwidget)
        self.VB_label.setMinimumSize(QtCore.QSize(160, 0))
        self.VB_label.setMaximumSize(QtCore.QSize(160, 16777215))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(11)
        self.VB_label.setFont(font)
        self.VB_label.setObjectName("VB_label")
        self.verticalLayout_2.addWidget(self.VB_label)
        self.Gate_label = QtWidgets.QLabel(self.centralwidget)
        self.Gate_label.setMinimumSize(QtCore.QSize(160, 0))
        self.Gate_label.setMaximumSize(QtCore.QSize(160, 16777215))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(11)
        self.Gate_label.setFont(font)
        self.Gate_label.setObjectName("Gate_label")
        self.verticalLayout_2.addWidget(self.Gate_label)
        self.line_6 = QtWidgets.QFrame(self.centralwidget)
        self.line_6.setMinimumSize(QtCore.QSize(160, 0))
        self.line_6.setMaximumSize(QtCore.QSize(140, 16777215))
        self.line_6.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_6.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_6.setObjectName("line_6")
        self.verticalLayout_2.addWidget(self.line_6)
        self.report_button = QtWidgets.QPushButton(self.centralwidget)
        self.report_button.setMinimumSize(QtCore.QSize(160, 0))
        self.report_button.setMaximumSize(QtCore.QSize(160, 16777215))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(11)
        self.report_button.setFont(font)
        self.report_button.setObjectName("report_button")
        self.verticalLayout_2.addWidget(self.report_button)
        spacerItem = QtWidgets.QSpacerItem(20, 297, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.verticalLayout_2.addItem(spacerItem)
        self.label_9 = QtWidgets.QLabel(self.centralwidget)
        self.label_9.setMinimumSize(QtCore.QSize(160, 0))
        self.label_9.setMaximumSize(QtCore.QSize(140, 92))
        self.label_9.setText("")
        self.label_9.setPixmap(QtGui.QPixmap("pic_labels/koziy.PNG"))
        self.label_9.setScaledContents(True)
        self.label_9.setObjectName("label_9")
        self.verticalLayout_2.addWidget(self.label_9)
        self.label_8 = QtWidgets.QLabel(self.centralwidget)
        self.label_8.setMinimumSize(QtCore.QSize(160, 0))
        self.label_8.setMaximumSize(QtCore.QSize(140, 42))
        self.label_8.setText("")
        self.label_8.setPixmap(QtGui.QPixmap("pic_labels/logo.png"))
        self.label_8.setScaledContents(True)
        self.label_8.setWordWrap(False)
        self.label_8.setOpenExternalLinks(False)
        self.label_8.setObjectName("label_8")
        self.verticalLayout_2.addWidget(self.label_8)
        self.horizontalLayout.addLayout(self.verticalLayout_2)
        self.plot_container = ThreeDSurfaceGraphWindow(self.centralwidget)
        self.plot_container.setMinimumSize(QtCore.QSize(0, 0))
        self.plot_container.setMaximumSize(QtCore.QSize(1111111, 16777215))
        self.plot_container.setObjectName("plot_container")
        self.horizontalLayout.addWidget(self.plot_container)
        self.verticalLayout = QtWidgets.QVBoxLayout()
        self.verticalLayout.setObjectName("verticalLayout")
        self.twodplot_container = TwoDGraphicWindow(self.centralwidget)
        self.twodplot_container.setMinimumSize(QtCore.QSize(0, 0))
        self.twodplot_container.setObjectName("twodplot_container")
        self.verticalLayout.addWidget(self.twodplot_container)
        self.trplot_container = TwoDGraphicWindow(self.centralwidget)
        self.trplot_container.setMinimumSize(QtCore.QSize(0, 0))
        self.trplot_container.setMaximumSize(QtCore.QSize(16777215, 1111111))
        self.trplot_container.setObjectName("trplot_container")
        self.verticalLayout.addWidget(self.trplot_container)
        self.horizontalLayout.addLayout(self.verticalLayout)
        self.horizontalLayout_2.addLayout(self.horizontalLayout)
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 988, 22))
        self.menubar.setObjectName("menubar")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow"))
        self.label_5.setText(_translate("MainWindow", "SPD selection"))
        self.label.setText(_translate("MainWindow", "SPD name"))
        self.label_2.setText(_translate("MainWindow", "Measurement date"))
        self.label_13.setText(_translate("MainWindow", "Dead time setup"))
        self.label_3.setText(_translate("MainWindow", "Temperature"))
        self.label_4.setText(_translate("MainWindow", "Plot setup"))
        self.label_6.setText(_translate("MainWindow", "Slice setup"))
        self.label_7.setText(_translate("MainWindow", "2d plot"))
        self.label_10.setText(_translate("MainWindow", "Time resolution"))
        self.label_type.setText(_translate("MainWindow", "SPD type:"))
        self.label_11.setText(_translate("MainWindow", "Efficiency parameters"))
        self.QE_label.setText(_translate("MainWindow", "TextLabel"))
        self.DCR_label.setText(_translate("MainWindow", "TextLabel"))
        self.QE_DCR_label.setText(_translate("MainWindow", "TextLabel"))
        self.AP_label.setText(_translate("MainWindow", "TextLabel"))
        self.DT_label.setText(_translate("MainWindow", "TextLabel"))
        self.label_12.setText(_translate("MainWindow", "Technical parameters"))
        self.T_label.setText(_translate("MainWindow", "TextLabel"))
        self.VB_label.setText(_translate("MainWindow", "TextLabel"))
        self.Gate_label.setText(_translate("MainWindow", "TextLabel"))
        self.report_button.setText(_translate("MainWindow", "Create report"))

from threedsurfacegraphwindow import ThreeDSurfaceGraphWindow
from twodgraphicwindow import TwoDGraphicWindow
