from matplotlib import cm
from matplotlib.pyplot import figure
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from mpl_toolkits.mplot3d import Axes3D
from PyQt5.QtWidgets import QWidget

class TwoDGraphicWindow(FigureCanvas):  # Class for 3D window
    def __init__(self, parent = 0):
        #super().__init__()
        self.plot_colorbar = None
        self.plot_figure = figure(figsize=(7, 7))
        FigureCanvas.__init__(self, self.plot_figure)  # creating FigureCanvas
        self.axes = self.plot_figure.add_subplot(111)
        #self.axes = self.plot_figure.gca(projection='3d')  # generates 3D Axes object
        #self.toolbar = NavigationToolbar(self.plot_figure, self)
        #self.fig, self.axes  = self.plot_figure.subplots()
        self.setWindowTitle("figure")  # sets Window title

    def draw_graph(self, slice, title = 'title', labels = ['x', 'y'], style = 'o-', limit = None):  # Function for graph plotting
        self.axes.clear()

        if type(slice) == int:
            self.axes.set_title(title)
            self.axes.set_xlabel(labels[0])
            self.axes.set_ylabel(labels[1])
            self.axes.grid(True)
            self.draw()
        elif len(slice) == 3:
            # plot
            x = slice[1]
            y = slice[2]
            plot_stuff = self.axes.plot(x, y, style)

            self.axes.set_title(title)
            self.axes.set_xlabel(labels[0])
            self.axes.set_ylabel(labels[1])
            self.axes.grid(True)
            if limit != None:
                self.axes.set_ylim((0, 1.05 * limit))
            self.draw()
        elif len(slice) == 2:
            # plot
            x = slice[0]
            y = slice[1]
            plot_stuff = self.axes.plot(x, y, style)

            self.axes.set_title(title)
            self.axes.set_xlabel(labels[0])
            self.axes.set_ylabel(labels[1])
            self.axes.grid(True)
            if limit != None:
                self.axes.set_ylim((0, 1.05 * limit))
            self.draw()

    def plot_slice(self, slice):
        self.axes.axvline(x=slice, color='m', linestyle='--')