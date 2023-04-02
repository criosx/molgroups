# This Python file uses the following encoding: utf-8
import sys

from PySide6.QtWidgets import QApplication, QWidget

# Important:
# You need to run the following command to generate the ui_form.py file
#     pyside6-uic form.ui -o ui_form.py, or
#     pyside2-uic form.ui -o ui_form.py
from ui_form import Ui_ML_widget

import matplotlib
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar)
matplotlib.use("Qt5Agg")

class Widget(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.ui = Ui_ML_widget()
        self.ui.setupUi(self)

        # Initialize the figure in our window
        self.fig = None
        self.ax = None
        self.canvas = None
        self.toolbar = None
        figure = Figure()  # Prep empty figure
        axis = figure.add_subplot(111)  # Prep empty plot
        self.initialize_figure(figure, axis)  # Initialize!

    def initialize_figure(self, fig, ax):
        """ Initializes a matplotlib figure inside a GUI container.
            Only call this once when initializing.
        """
        # Figure creation (self.fig and self.ax)
        self.fig = fig
        self.ax = ax
        # Canvas creation
        self.canvas = FigureCanvas(self.fig)
        self.ui.plotLayout.addWidget(self.canvas)
        self.canvas.draw()
        # Toolbar creation
        self.toolbar = NavigationToolbar(self.canvas, self.ui.w_PlotField,
                                         coordinates=True)
        self.ui.plotLayout.addWidget(self.toolbar)


if __name__ == "__main__":
    app = QApplication(sys.argv)
    widget = Widget()
    widget.show()
    sys.exit(app.exec())
