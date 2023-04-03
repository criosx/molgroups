# This Python file uses the following encoding: utf-8

import numpy
import os
import pickle
import sys
import tensorflow as tf

from sasmodels.data import load_data
from PySide6.QtWidgets import QApplication, QWidget, QFileDialog
from PySide6.QtCore import QTimer
from PySide6.QtCore import SIGNAL

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


def fnLoadObject(sFileName):
    with open(sFileName, 'rb') as file:
        load_object = pickle.load(file)
    return load_object


class Widget(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.ui = Ui_ML_widget()
        self.ui.setupUi(self)

        # initialize file and folder memory
        self.old_sans_folder = '.'
        self.old_ml_folder = '.'

        # Initialize data
        self.sans_file_name = None
        self.Iq = None
        self.Q = None
        self.dI = None
        self.dQ = None

        self.background = 0.0
        self.solvent_sld = 6.4
        self.ui.lE_Background.setText(str(self.background))
        self.ui.lE_SolventSLD.setText(str(self.solvent_sld))

        # Initialize ML model
        dirname = 'ml_model'
        self.sans_models = None
        self.par_names = None
        self.ml_model = None
        self.ui.lE_MLFileName.setText(dirname)
        self.load_ml_model()

        # Initialize the figure in our window
        self.fig = None
        self.ax = None
        self.canvas = None
        self.toolbar = None
        figure = Figure()  # Prep empty figure
        axis = figure.add_subplot(111)  # Prep empty plot
        self.initialize_figure(figure, axis)  # Initialize!

        # input timer
        self.timer = QTimer(None)
        self.timer.setSingleShot(True)
        self.timer.timeout.connect(self.process_window)

        # functionalities
        self.ui.tB_OpenSANSFile.clicked.connect(self.open_sans_file_dialog)
        self.ui.tB_OpenMLFile.clicked.connect(self.open_ml_folder_dialog)

        self.ui.pB_AutoBackground.clicked.connect(self.auto_background)

        self.ui.lE_MLFileName.returnPressed.connect(self.load_ml_model)

        self.ui.lE_SANSFileName.textChanged.connect(self.entries_changed)
        self.ui.lE_Background.textChanged.connect(self.entries_changed)
        self.ui.lE_SolventSLD.textChanged.connect(self.entries_changed)

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

    def auto_background(self):
        if self.Iq is None:
            return

        # average of last n points that is within the error bar of the n-1 th point
        background = None
        for i in range(-5, -30, -1):
            background_new = numpy.average(self.Iq[i:])
            if numpy.abs(self.Iq[i-1] - background_new) > 2 * self.dI[i-1]:
                break
            background = background_new

        if background is not None:
            self.ui.lE_Background.setText(str(background))


    def open_sans_file_dialog(self):
        file_name, _ = QFileDialog.getOpenFileName(parent=self, caption='OPEN SANS file', dir=self.old_sans_folder)
        self.ui.lE_SANSFileName.setText(os.path.relpath(file_name))
        self.old_sans_folder = os.path.dirname(file_name)

    def open_ml_folder_dialog(self):
        folder_name = QFileDialog.getExistingDirectory(parent=self, caption='OPEN ML folder', dir=self.old_ml_folder)
        self.ui.lE_MLFileName.setText(os.path.relpath(folder_name))
        self.old_ml_folder = folder_name
        self.load_ml_model()

    def entries_changed(self):
        self.timer.start(500)

    def process_window(self):
        sans_file_name = self.ui.lE_SANSFileName.text()

        if sans_file_name != self.sans_file_name:
            self.sans_file_name = sans_file_name
            if not os.path.isfile(sans_file_name):
                self.ui.tB_Prediction.setText("Select valid SANS Data file.")
                return
            try:
                ds = load_data(sans_file_name)
                self.Q = ds.x
                self.Iq = ds.y
                self.dI = ds.dy
                self.dQ = ds.dx
            except IOError:
                self.ui.tB_Prediction.setText("Failed to load SANS Data file.")
                return

        solvent_sld = self.ui.lE_SolventSLD.text()
        background = self.ui.lE_Background.text()

        try:
            self.solvent_sld = float(solvent_sld)
        except ValueError:
            self.ui.tB_Prediction.setText("Solvent SLD is not a number.")
            return
        try:
            self.background = float(background)
        except ValueError:
            self.ui.tB_Prediction.setText("Background is not a number.")
            return
        self.ui.tB_Prediction.clear()

        self.update_plot()

        Iq = self.Iq - self.background
        Iq = numpy.log10(numpy.abs(Iq))

        # interpolation of SANS data to hardcoded grid
        qmin = 0.01
        qmax = 0.8
        numpoints = int((numpy.log10(qmax) - numpy.log10(qmin)) * 60)
        qvec = numpy.logspace(numpy.log10(qmin), numpy.log10(qmax), num=numpoints, endpoint=True)
        qvec = qvec[:105]

        intensity = numpy.interp(qvec, self.Q, Iq)
        intensity = intensity[numpy.newaxis, :]
        intensity = tf.convert_to_tensor(intensity, dtype=tf.float32)

        sup = numpy.array([self.background, self.solvent_sld]).astype('float32')
        sup = sup[numpy.newaxis, :]

        y_pred = self.ml_model.predict([intensity, sup])

        self.ui.tB_Prediction.clear()

        self.ui.tB_Prediction.append("---Classification---")
        for i, model in enumerate(self.sans_models):
            pstr = f'{y_pred[-1][0][i]:.2f}' + ' ' + model
            self.ui.tB_Prediction.append(pstr)
        self.ui.tB_Prediction.append("")

        self.ui.tB_Prediction.append("---Regression---")
        for i in range(len(y_pred) - 1):
            pstr = 'Model: ' + self.sans_models[i]
            self.ui.tB_Prediction.append(pstr)
            for j in range(len(self.par_names[i])):
                parname = self.par_names[i][j]
                if 'sld' in parname:
                    correction = 0.1
                else:
                    correction = 1
                pstr = parname + ' ' + f'{y_pred[i][0][j] * correction:.4f}'
                self.ui.tB_Prediction.append(pstr)
            self.ui.tB_Prediction.append("")

        self.ui.tB_Prediction.verticalScrollBar().setValue(0)

    def load_ml_model(self):
        dirname = self.ui.lE_MLFileName.text()
        try:
            self.sans_models = fnLoadObject(os.path.join(dirname, 'sans_models.dat'))
            self.par_names = fnLoadObject(os.path.join(dirname, 'par_names.dat'))
            # Load model for prediction. Compile = False avoids supplying the custom loss function.
            self.ml_model = tf.keras.models.load_model(dirname, compile=False)
        except IOError:
            self.ui.tB_Prediction.setText('Could not load ML model')
            self.ml_model = None
            self.sans_models = None
            self.par_names = None

        self.ui.tB_Prediction.clear()
        self.process_window()

    def update_plot(self):

        # Clear whatever was in the plot before
        self.ax.clear()
        # Plot data, add labels, change colors, ...
        self.ax.errorbar(self.Q, self.Iq, self.dI, ls='none', color='deepskyblue')
        self.ax.scatter(self.Q, self.Iq, s=30, marker='o', facecolors='none', edgecolors='deepskyblue', label='orignal')
        self.ax.errorbar(self.Q, self.Iq-self.background, self.dI, ls='none', color='darkred')
        self.ax.scatter(self.Q, self.Iq-self.background, s=30, marker='o', facecolors='none', edgecolors='darkred',
                        label='-background')

        self.ax.legend(fontsize=16)
        self.ax.set_ylabel("$Iq$ (cm$^{-1}$)", fontsize=16)
        self.ax.set_yscale('log')
        self.ax.set_xscale('log')
        self.ax.minorticks_on()
        self.ax.tick_params(which="both", direction="in", labelsize=16)
        self.ax.tick_params(bottom=True, top=True, left=True, right=True, which="both")
        self.ax.set_xlabel("$q$ (Å$^{-1}$)", fontsize=16)
        # self.ax.figure.set_size_inches(8, 5)
        # Make sure everything fits inside the canvas
        self.fig.tight_layout()
        # Show the new figure in the interface
        self.canvas.draw()


if __name__ == "__main__":
    app = QApplication(sys.argv)
    widget = Widget()
    widget.show()
    sys.exit(app.exec())
