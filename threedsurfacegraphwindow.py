from matplotlib import cm
from matplotlib.pyplot import figure
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from mpl_toolkits.mplot3d import Axes3D
from PyQt5.QtWidgets import QWidget

class ThreeDSurfaceGraphWindow(FigureCanvas):  # Class for 3D window
    def __init__(self, parent = 0):
        #super().__init__()
        self.plot_colorbar = None
        self.plot_figure = figure(figsize=(7, 7))
        FigureCanvas.__init__(self, self.plot_figure)  # creating FigureCanvas
        self.axes = self.plot_figure.gca(projection='3d')  # generates 3D Axes object
        self.setWindowTitle("figure")  # sets Window title

    def draw_graph(self, x, y, z, title = 'title', labels = ['x', 'y', 'z']):  # Function for graph plotting
        self.axes.clear()
        if self.plot_colorbar is not None:  # avoids adding one more colorbar at each draw operation
            self.plot_colorbar.remove()
        # plots the 3D surface plot
        plot_stuff = self.axes.plot_surface(x, y, z,
                                            cmap=cm.coolwarm, linewidth=0, alpha = 0.7, antialiased=False)
        self.axes.zaxis.set_major_locator(LinearLocator(10))
        self.axes.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

        self.axes.set_title(title)
        self.axes.set_xlabel(labels[0])
        self.axes.set_ylabel(labels[1])
        #self.axes.set_zlabel(labels[2])
        # Add a color bar which maps values to colors.
        #self.plot_colorbar = self.plot_figure.colorbar(plot_stuff, shrink=0.5, aspect=5)
        # draw plot
        self.draw()

    def clear_graphics(self):
        self.axes.clear()
        if self.plot_colorbar is not None:  # avoids adding one more colorbar at each draw operation
            self.plot_colorbar.remove()

    def plot_slice(self, slice):
        #slice = [x, y, z]
        '''
        try:
            line.remove()
        except NameError:
            print('there is not such variable line')
        '''
        self.line = self.axes.plot(slice[0],slice[1],slice[2], color = 'k', linewidth=5)
        self.draw()

    def plot_point(self, point):
        self.line2 = self.axes.scatter(point[0], point[1], point[2], color='m', linewidth=10)