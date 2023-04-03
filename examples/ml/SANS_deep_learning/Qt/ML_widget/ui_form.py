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
from PySide6.QtWidgets import (QApplication, QFrame, QGridLayout, QGroupBox,
    QHBoxLayout, QLabel, QLineEdit, QPushButton,
    QSizePolicy, QSpacerItem, QTextBrowser, QToolButton,
    QVBoxLayout, QWidget)

class Ui_ML_widget(object):
    def setupUi(self, ML_widget):
        if not ML_widget.objectName():
            ML_widget.setObjectName(u"ML_widget")
        ML_widget.setWindowModality(Qt.ApplicationModal)
        ML_widget.resize(1055, 692)
        self.actionShow = QAction(ML_widget)
        self.actionShow.setObjectName(u"actionShow")
        self.horizontalLayout = QHBoxLayout(ML_widget)
        self.horizontalLayout.setObjectName(u"horizontalLayout")
        self.groupBox = QGroupBox(ML_widget)
        self.groupBox.setObjectName(u"groupBox")
        font = QFont()
        font.setPointSize(14)
        font.setUnderline(False)
        self.groupBox.setFont(font)
        self.groupBox.setAlignment(Qt.AlignCenter)
        self.gridLayout_2 = QGridLayout(self.groupBox)
        self.gridLayout_2.setObjectName(u"gridLayout_2")
        self.w_PlotField = QWidget(self.groupBox)
        self.w_PlotField.setObjectName(u"w_PlotField")
        self.w_PlotField.setMinimumSize(QSize(600, 200))
        self.w_PlotField.setMaximumSize(QSize(16777215, 1000))
        self.w_PlotField.setBaseSize(QSize(20, 20))
        self.plotLayout = QVBoxLayout(self.w_PlotField)
        self.plotLayout.setObjectName(u"plotLayout")

        self.gridLayout_2.addWidget(self.w_PlotField, 9, 0, 1, 4)

        self.tB_OpenMLFile = QToolButton(self.groupBox)
        self.tB_OpenMLFile.setObjectName(u"tB_OpenMLFile")

        self.gridLayout_2.addWidget(self.tB_OpenMLFile, 13, 3, 1, 1)

        self.label_2 = QLabel(self.groupBox)
        self.label_2.setObjectName(u"label_2")
        sizePolicy = QSizePolicy(QSizePolicy.Fixed, QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_2.sizePolicy().hasHeightForWidth())
        self.label_2.setSizePolicy(sizePolicy)

        self.gridLayout_2.addWidget(self.label_2, 6, 0, 1, 1)

        self.tB_OpenSANSFile = QToolButton(self.groupBox)
        self.tB_OpenSANSFile.setObjectName(u"tB_OpenSANSFile")

        self.gridLayout_2.addWidget(self.tB_OpenSANSFile, 4, 3, 1, 1)

        self.lE_MLFileName = QLineEdit(self.groupBox)
        self.lE_MLFileName.setObjectName(u"lE_MLFileName")

        self.gridLayout_2.addWidget(self.lE_MLFileName, 13, 1, 1, 1)

        self.label = QLabel(self.groupBox)
        self.label.setObjectName(u"label")
        sizePolicy.setHeightForWidth(self.label.sizePolicy().hasHeightForWidth())
        self.label.setSizePolicy(sizePolicy)

        self.gridLayout_2.addWidget(self.label, 4, 0, 1, 1)

        self.lE_SANSFileName = QLineEdit(self.groupBox)
        self.lE_SANSFileName.setObjectName(u"lE_SANSFileName")
        self.lE_SANSFileName.setMinimumSize(QSize(0, 0))

        self.gridLayout_2.addWidget(self.lE_SANSFileName, 4, 1, 1, 2)

        self.lE_Background = QLineEdit(self.groupBox)
        self.lE_Background.setObjectName(u"lE_Background")

        self.gridLayout_2.addWidget(self.lE_Background, 7, 1, 1, 2)

        self.label_3 = QLabel(self.groupBox)
        self.label_3.setObjectName(u"label_3")
        sizePolicy.setHeightForWidth(self.label_3.sizePolicy().hasHeightForWidth())
        self.label_3.setSizePolicy(sizePolicy)

        self.gridLayout_2.addWidget(self.label_3, 7, 0, 1, 1)

        self.label_4 = QLabel(self.groupBox)
        self.label_4.setObjectName(u"label_4")
        font1 = QFont()
        font1.setPointSize(14)
        font1.setBold(True)
        font1.setUnderline(False)
        self.label_4.setFont(font1)

        self.gridLayout_2.addWidget(self.label_4, 3, 0, 1, 1)

        self.label_5 = QLabel(self.groupBox)
        self.label_5.setObjectName(u"label_5")
        self.label_5.setFont(font1)

        self.gridLayout_2.addWidget(self.label_5, 12, 0, 1, 1)

        self.pB_AutoBackground = QPushButton(self.groupBox)
        self.pB_AutoBackground.setObjectName(u"pB_AutoBackground")

        self.gridLayout_2.addWidget(self.pB_AutoBackground, 7, 3, 1, 1)

        self.lE_SolventSLD = QLineEdit(self.groupBox)
        self.lE_SolventSLD.setObjectName(u"lE_SolventSLD")

        self.gridLayout_2.addWidget(self.lE_SolventSLD, 6, 1, 1, 2)

        self.label_6 = QLabel(self.groupBox)
        self.label_6.setObjectName(u"label_6")

        self.gridLayout_2.addWidget(self.label_6, 13, 0, 1, 1)

        self.line = QFrame(self.groupBox)
        self.line.setObjectName(u"line")
        self.line.setFrameShape(QFrame.HLine)
        self.line.setFrameShadow(QFrame.Sunken)

        self.gridLayout_2.addWidget(self.line, 11, 0, 1, 4)

        self.verticalSpacer = QSpacerItem(20, 40, QSizePolicy.Minimum, QSizePolicy.Expanding)

        self.gridLayout_2.addItem(self.verticalSpacer, 10, 0, 1, 4)


        self.horizontalLayout.addWidget(self.groupBox)

        self.tB_Prediction = QTextBrowser(ML_widget)
        self.tB_Prediction.setObjectName(u"tB_Prediction")

        self.horizontalLayout.addWidget(self.tB_Prediction)


        self.retranslateUi(ML_widget)

        QMetaObject.connectSlotsByName(ML_widget)
    # setupUi

    def retranslateUi(self, ML_widget):
        ML_widget.setWindowTitle(QCoreApplication.translate("ML_widget", u"Widget", None))
        self.actionShow.setText(QCoreApplication.translate("ML_widget", u"Show", None))
        self.groupBox.setTitle(QCoreApplication.translate("ML_widget", u"SANS ML Prediction ", None))
        self.tB_OpenMLFile.setText(QCoreApplication.translate("ML_widget", u"...", None))
        self.label_2.setText(QCoreApplication.translate("ML_widget", u"Solvent SLD", None))
        self.tB_OpenSANSFile.setText(QCoreApplication.translate("ML_widget", u"...", None))
        self.label.setText(QCoreApplication.translate("ML_widget", u"File Name", None))
        self.label_3.setText(QCoreApplication.translate("ML_widget", u"Background", None))
        self.label_4.setText(QCoreApplication.translate("ML_widget", u"SANS Data", None))
        self.label_5.setText(QCoreApplication.translate("ML_widget", u"Machine Learning Model", None))
        self.pB_AutoBackground.setText(QCoreApplication.translate("ML_widget", u"Auto", None))
        self.label_6.setText(QCoreApplication.translate("ML_widget", u"Folder Name", None))
    # retranslateUi

