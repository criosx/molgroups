# -*- coding: utf-8 -*-

################################################################################
## Form generated from reading UI file 'form.ui'
##
## Created by: Qt User Interface Compiler version 6.4.3
##
## WARNING! All changes made in this file will be lost when recompiling UI file!
################################################################################

from PySide6.QtCore import (QCoreApplication, QDate, QDateTime, QLocale,
    QMetaObject, QObject, QPoint, QRect,
    QSize, QTime, QUrl, Qt)
from PySide6.QtGui import (QAction, QBrush, QColor, QConicalGradient,
    QCursor, QFont, QFontDatabase, QGradient,
    QIcon, QImage, QKeySequence, QLinearGradient,
    QPainter, QPalette, QPixmap, QRadialGradient,
    QTransform)
from PySide6.QtWidgets import (QApplication, QGridLayout, QGroupBox, QLabel,
    QLineEdit, QSizePolicy, QTextBrowser, QToolButton,
    QVBoxLayout, QWidget)

class Ui_ML_widget(object):
    def setupUi(self, ML_widget):
        if not ML_widget.objectName():
            ML_widget.setObjectName(u"ML_widget")
        ML_widget.setWindowModality(Qt.ApplicationModal)
        ML_widget.resize(650, 656)
        self.actionShow = QAction(ML_widget)
        self.actionShow.setObjectName(u"actionShow")
        self.gridLayout = QGridLayout(ML_widget)
        self.gridLayout.setObjectName(u"gridLayout")
        self.tB_Prediction = QTextBrowser(ML_widget)
        self.tB_Prediction.setObjectName(u"tB_Prediction")

        self.gridLayout.addWidget(self.tB_Prediction, 5, 0, 1, 2)

        self.w_PlotField = QWidget(ML_widget)
        self.w_PlotField.setObjectName(u"w_PlotField")
        self.w_PlotField.setMinimumSize(QSize(20, 200))
        self.w_PlotField.setBaseSize(QSize(20, 20))
        self.plotLayout = QVBoxLayout(self.w_PlotField)
        self.plotLayout.setObjectName(u"plotLayout")

        self.gridLayout.addWidget(self.w_PlotField, 3, 0, 1, 2)

        self.label_4 = QLabel(ML_widget)
        self.label_4.setObjectName(u"label_4")

        self.gridLayout.addWidget(self.label_4, 4, 0, 1, 2)

        self.groupBox = QGroupBox(ML_widget)
        self.groupBox.setObjectName(u"groupBox")
        self.groupBox.setAlignment(Qt.AlignCenter)
        self.gridLayout_2 = QGridLayout(self.groupBox)
        self.gridLayout_2.setObjectName(u"gridLayout_2")
        self.lineEdit = QLineEdit(self.groupBox)
        self.lineEdit.setObjectName(u"lineEdit")
        self.lineEdit.setMinimumSize(QSize(515, 0))

        self.gridLayout_2.addWidget(self.lineEdit, 0, 1, 1, 1)

        self.label_2 = QLabel(self.groupBox)
        self.label_2.setObjectName(u"label_2")

        self.gridLayout_2.addWidget(self.label_2, 2, 0, 1, 1)

        self.label = QLabel(self.groupBox)
        self.label.setObjectName(u"label")

        self.gridLayout_2.addWidget(self.label, 0, 0, 1, 1)

        self.label_3 = QLabel(self.groupBox)
        self.label_3.setObjectName(u"label_3")

        self.gridLayout_2.addWidget(self.label_3, 3, 0, 1, 1)

        self.toolButton = QToolButton(self.groupBox)
        self.toolButton.setObjectName(u"toolButton")

        self.gridLayout_2.addWidget(self.toolButton, 0, 2, 1, 1)

        self.lineEdit_2 = QLineEdit(self.groupBox)
        self.lineEdit_2.setObjectName(u"lineEdit_2")

        self.gridLayout_2.addWidget(self.lineEdit_2, 2, 1, 1, 2)

        self.lineEdit_3 = QLineEdit(self.groupBox)
        self.lineEdit_3.setObjectName(u"lineEdit_3")

        self.gridLayout_2.addWidget(self.lineEdit_3, 3, 1, 1, 2)


        self.gridLayout.addWidget(self.groupBox, 0, 0, 2, 2)

        self.label_5 = QLabel(ML_widget)
        self.label_5.setObjectName(u"label_5")

        self.gridLayout.addWidget(self.label_5, 2, 0, 1, 2)


        self.retranslateUi(ML_widget)

        QMetaObject.connectSlotsByName(ML_widget)
    # setupUi

    def retranslateUi(self, ML_widget):
        ML_widget.setWindowTitle(QCoreApplication.translate("ML_widget", u"Widget", None))
        self.actionShow.setText(QCoreApplication.translate("ML_widget", u"Show", None))
        self.label_4.setText(QCoreApplication.translate("ML_widget", u"Prediction:", None))
        self.groupBox.setTitle(QCoreApplication.translate("ML_widget", u"SANS ML Prediction Widget", None))
        self.label_2.setText(QCoreApplication.translate("ML_widget", u"Solvent SLD", None))
        self.label.setText(QCoreApplication.translate("ML_widget", u"File name", None))
        self.label_3.setText(QCoreApplication.translate("ML_widget", u"Background", None))
        self.toolButton.setText(QCoreApplication.translate("ML_widget", u"...", None))
        self.label_5.setText(QCoreApplication.translate("ML_widget", u"SANS data with and without background subtraction", None))
    # retranslateUi

