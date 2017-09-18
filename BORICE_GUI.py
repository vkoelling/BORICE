#!/usr/bin/env python2
"""
BORICE is software for the estimation of population mean outcrossing rates and
inbreeding coefficients using Bayesian methods.
Copyright (C) 2012 Vanessa A. Koelling, Patrick J. Monnahan, John K. Kelly

This program is free software: you can redistribute it and/or modify it under 
the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>
"""
from PyQt4 import QtCore, QtGui
import sys
import os
from MyDoubleValidator import MyDoubleValidator
import BORICE

class MainWindow(QtGui.QMainWindow):
	"""
	This class defines the graphical user interface for the BORICE program,
	allowing for easy user input and data display.
	"""
	def __init__(self,parent,app):

		#initialize superclass and references to the application and parent
		super(MainWindow, self).__init__(parent)
		self.parent = parent
		self.app = app

		#set window size to default and set the icon and window title
		self.resize(700,550)
		self.setWindowTitle("BORICE")
		self.setWindowIcon(QtGui.QIcon(os.path.join('Images','dna.png')))

		#
		# Input/Output Data
		#
		#Outcrossing Rate. Can range from 0-1
		self.outcrossingRate = 0.50
		#Number of steps
		self.numSteps = 100000
		#Number of steps to burn in
		self.numBurnInSteps = 9999
		#Number of steps taken
		self.numStepsTaken = None
		#Number of marker loci
		self.numMarkerLoci = None
		#Marker loci names
		self.markerLociNames = []
		#Are there population names?
		self.popNamesExist = None
		#Are there subgroup names?
		self.subNamesExist = None

		self.dataFileName = ""

		self.outcrossingRateTuningParam = 0.05
		self.alleleFreqTuningParam = 0.1
		
		self.writeOutput2 = True
		self.writeOutput3 = True
		self.writeOutput4 = True

		#build the user interface
		self.setupUI()

	def setupUI(self):
		"""
		Creates the layout and all widgets for the MainWindow.
		"""
		#Central Widget
		centralWidget = QtGui.QWidget(self)
		self.setCentralWidget(centralWidget)

		#Central Layout
		centralLayout = QtGui.QVBoxLayout()
		centralLayout.setAlignment(QtCore.Qt.AlignCenter)
		centralWidget.setLayout(centralLayout)

		#Actions for the Menu Bar
		#Open Data File
		fileOpenAction = QtGui.QAction("&Open Data File", self)
		fileOpenAction.setShortcut(QtGui.QKeySequence.Open)
		helpText = "Open a data file"
		fileOpenAction.setToolTip(helpText)
		fileOpenAction.setStatusTip(helpText)
		fileOpenAction.triggered.connect(self.fileOpen)
		#Quit
		fileQuitAction = QtGui.QAction("&Quit", self)
		fileQuitAction.setShortcut(QtGui.QKeySequence.Quit)
		helpText = "Exit"
		fileQuitAction.setToolTip(helpText)
		fileQuitAction.setStatusTip(helpText)
		fileQuitAction.triggered.connect(self.close)

		#Menu Bar
		self.fileMenu = self.menuBar().addMenu("&File")
		self.helpMenu = self.menuBar().addMenu("&Help")

		#Add actions to Menu Bar
		self.fileMenu.addAction(fileOpenAction)
		self.fileMenu.addAction(fileQuitAction)
		#self.helpMenu.addAction()
		
		#Validator for entries to be from 0-1
		zeroOneValidator = MyDoubleValidator(0,1,2,self)
		
		#Validator for only positive integers(for numbers of steps)
		maxStep = 999999999
		posValidator = QtGui.QIntValidator(0, maxStep, self)
		
		#Max Widths
		maxLineWidth = 250
		
		self.numStepsText = QtGui.QLineEdit("100000")
		self.numStepsText.textChanged.connect(self.setNumSteps)
		self.numStepsText.setMaximumWidth(maxLineWidth)
		self.numStepsText.setAlignment(QtCore.Qt.AlignRight)
		self.numStepsText.setValidator(posValidator)
		#Number of steps to Burn In input
		self.numStepsBurnInText = QtGui.QLineEdit("9999")
		self.numStepsBurnInText.textChanged.connect(self.setBurnInSteps)
		self.numStepsBurnInText.setMaximumWidth(maxLineWidth)
		self.numStepsBurnInText.setAlignment(QtCore.Qt.AlignRight)
		self.numStepsBurnInText.setValidator(posValidator)
		#Outcrossing Rate Tuning Parameter
		self.outcrossingRateTuningParamText = QtGui.QLineEdit("0.05")
		self.outcrossingRateTuningParamText.textChanged.connect(self.setOutcrossingRateTuningParam)
		self.outcrossingRateTuningParamText.setMaximumWidth(maxLineWidth)
		self.outcrossingRateTuningParamText.setAlignment(QtCore.Qt.AlignRight)
		self.outcrossingRateTuningParamText.setValidator(zeroOneValidator)
		#Initial Population Outcrossing Rate
		self.initialPopulationOutcrossingRateText = QtGui.QLineEdit("0.5")
		self.initialPopulationOutcrossingRateText.textChanged.connect(self.setInitialPopulationOutcrossingRate)
		self.initialPopulationOutcrossingRateText.setValidator(zeroOneValidator)
		self.initialPopulationOutcrossingRateText.setMaximumWidth(maxLineWidth)
		self.initialPopulationOutcrossingRateText.setAlignment(QtCore.Qt.AlignRight)
		#Allele Frequency Tuning Parameter
		self.AlleleFreqTuningParamText = QtGui.QLineEdit("0.1")
		self.AlleleFreqTuningParamText.textChanged.connect(self.setAlleleFreqTuningParam)
		self.AlleleFreqTuningParamText.setMaximumWidth(maxLineWidth)
		self.AlleleFreqTuningParamText.setAlignment(QtCore.Qt.AlignRight)
		self.AlleleFreqTuningParamText.setValidator(zeroOneValidator)

		#Run Button
		runBox = QtGui.QWidget()
		self.runLayout = QtGui.QVBoxLayout()
		self.runLayout.setAlignment(QtCore.Qt.AlignCenter)
		runBox.setLayout(self.runLayout)
		self.runButton = QtGui.QPushButton("Run")
		self.runLayout.addWidget(self.runButton)
		
		#Posterior Distributions of Maternal Inbreeding Histories
		self.outputOptionTwoCheckBox = QtGui.QCheckBox()
		self.outputOptionTwoCheckBox.setChecked(True)
		self.outputOptionTwoCheckBox.toggled.connect(self.setWrite2)
		
		#List of t and F values
		self.outputOptionThreeCheckBox = QtGui.QCheckBox()
		self.outputOptionThreeCheckBox.setChecked(True)
		self.outputOptionThreeCheckBox.toggled.connect(self.setWrite3)
		
		#Posterior Distributions for each maternal genotype at each locus in each family
		self.outputOptionFourCheckBox = QtGui.QCheckBox()
		self.outputOptionFourCheckBox.setChecked(True)
		self.outputOptionFourCheckBox.toggled.connect(self.setWrite4)
		
		#Connect run button to spawn a progress dialog when clicked
		self.runButton.clicked.connect(self.progress)
		
		#Add input widgets to the layout
		settingsBox = QtGui.QGroupBox("Settings")
		centralLayout.addWidget(settingsBox)
		settingsLayout = QtGui.QVBoxLayout()
		settingsBox.setLayout(settingsLayout)
		centralLayout.addWidget(runBox)
		
		self.tabWidget = QtGui.QTabWidget(settingsBox)
		settingsLayout.addWidget(self.tabWidget)
		
		#General settings box (Number of settings known a priori)
		genSettingsBox = QtGui.QWidget()
		genSettingsLayout = QtGui.QFormLayout(genSettingsBox)
		genSettingsLayout.setSizeConstraint(QtGui.QLayout.SetMinimumSize)
		genSettingsLayout.setFieldGrowthPolicy(QtGui.QFormLayout.AllNonFixedFieldsGrow)
		genSettingsLayout.setContentsMargins(0, 0, 0, 0)
		genSettingsLayout.setFormAlignment(QtCore.Qt.AlignCenter)
		genSettingsLayout.setLabelAlignment(QtCore.Qt.AlignLeft)
		genSettingsBox.setLayout(genSettingsLayout)
		genSettingsLayout.addRow("Number of Steps: ", self.numStepsText)
		genSettingsLayout.addRow("Number of Burn In Steps: ", self.numStepsBurnInText)
		genSettingsLayout.addRow("Outcrossing Rate Tuning Parameter: ", self.outcrossingRateTuningParamText)
		genSettingsLayout.addRow("Allele Frequency Tuning Parameter: ", self.AlleleFreqTuningParamText)
		genSettingsLayout.addRow("Initial Population Outcrossing Rate: ", self.initialPopulationOutcrossingRateText)
		self.tabWidget.addTab(genSettingsBox, "General Settings")
		
		outputOptionsBox = QtGui.QWidget()
		outputOptionsLayout = QtGui.QFormLayout(outputOptionsBox)
		outputOptionsLayout.setContentsMargins(0, 0, 0, 0)
		outputOptionsLayout.setFormAlignment(QtCore.Qt.AlignCenter)
		outputOptionsLayout.setLabelAlignment(QtCore.Qt.AlignLeft)
		outputOptionsBox.setLayout(outputOptionsLayout)
		outputOptionsLayout.addRow("Posterior Distributions of Maternal Inbreeding Histories: ", self.outputOptionTwoCheckBox)
		outputOptionsLayout.addRow("List of t and F values", self.outputOptionThreeCheckBox)
		outputOptionsLayout.addRow("Posterior Distributions for each maternal genotype at each locus in each family: ", self.outputOptionFourCheckBox)
		self.tabWidget.addTab(outputOptionsBox, "File Output Options")
		
		infoBox = QtGui.QWidget()
		infoLayout = QtGui.QFormLayout()
		infoLayout.setContentsMargins(0, 0, 0, 0)
		infoLayout.setFormAlignment(QtCore.Qt.AlignCenter)
		infoLayout.setLabelAlignment(QtCore.Qt.AlignLeft)
		infoBox.setLayout(infoLayout)
		
		self.numLociLabel = QtGui.QLabel("")
		self.numLociLabel.setSizePolicy(QtGui.QSizePolicy.Maximum,
										QtGui.QSizePolicy.Maximum)
		self.numFamiliesLabel = QtGui.QLabel("")
		self.numFamiliesLabel.setSizePolicy(QtGui.QSizePolicy.Maximum,
										QtGui.QSizePolicy.Maximum)
		self.numIndividualsLabel = QtGui.QLabel("")
		self.numIndividualsLabel.setSizePolicy(QtGui.QSizePolicy.Maximum,
										QtGui.QSizePolicy.Maximum)

		infoLayout.addRow("Number of Marker Loci: ", self.numLociLabel)
		infoLayout.addRow("Number of Families: ", self.numFamiliesLabel)
		infoLayout.addRow("Number of Individuals: ", self.numIndividualsLabel)
		self.tabWidget.addTab(infoBox, "Input Data Summary")


		#Hide the interface until the user opens a file
		self.centralWidget().hide()
	
	def addUI(self):
		"""Adds the Locus Settings Window to the Main Window.
		"""
		#Locus CheckBox
		locusCheckBoxList = []
		for n, locus in enumerate(self.locusModel):
			locusCheckBox = QtGui.QCheckBox(str(n), self)
			palette = locusCheckBox.palette()
			palette.setColor(QtGui.QPalette.WindowText, QtGui.QColor(0,0,0,0))
			locusCheckBox.setPalette(palette)
			locusCheckBox.setChecked(False)
			locusCheckBox.toggled.connect(self.setLocus)
			locusCheckBoxList.append(locusCheckBox)
		self.locusCheckBoxList = locusCheckBoxList

		#Locus settings box (Number of settings determined by input data file)
		locusSettingsBox = QtGui.QWidget()
		locusSettingsLayout = QtGui.QFormLayout(locusSettingsBox)
		locusSettingsLayout.setContentsMargins(0, 0, 0, 0)
		locusSettingsLayout.setFormAlignment(QtCore.Qt.AlignCenter)
		locusSettingsLayout.setLabelAlignment(QtCore.Qt.AlignLeft)
		locusSettingsBox.setLayout(locusSettingsLayout)
		for n, locusCheckBox in enumerate(self.locusCheckBoxList):
			locusSettingsLayout.addRow("Locus %s: " % (n + 1), locusCheckBox)
		self.tabWidget.addTab(locusSettingsBox, "Locus Settings")

	def setBurnInSteps(self, steps):
		if steps:			
			self.numBurnInSteps = int(str(steps))

	def setNumSteps(self, steps):
		if steps:
			self.numSteps = int(str(steps))
			
	def setOutcrossingRateTuningParam(self, tuningParam):
		if tuningParam:
			self.outcrossingRateTuningParam = tuningParam
			
	def setAlleleFreqTuningParam(self, tuningParam):
		if tuningParam:
			self.alleleFreqTuningParam = tuningParam
	
	def setInitialPopulationOutcrossingRate(self, outcrossingRate):
		if outcrossingRate:
			self.outcrossingRate = outcrossingRate
	
	def setLocus(self, null_locus):
		sender = self.sender()
		i = int(sender.text())
		self.locusModel[i] = null_locus
	
	def setWrite2(self, writeOutput2):
		self.writeOutput2 = writeOutput2
	
	def setWrite3(self, writeOutput3):
		self.writeOutput3 = writeOutput3
	
	def setWrite4(self, writeOutput4):
		self.writeOutput4 = writeOutput4
	
	def fileOpen(self):
		self.dataFileName = QtGui.QFileDialog.getOpenFileName(self, "Open Data File", ".", "Data Files (*.csv)")
		#Populate data entry tables for Allele Frequency Selection and Inbreeding History Selection
		if(self.dataFileName):
			dataFile = open(self.dataFileName, "rU")
			try:
				marker_names, families = BORICE.parse_csv(dataFile, ',')	
			except BORICE.CSVFileParseException as x:
				sys.exit(str(x))

			dataFile = open(self.dataFileName, "rU")
			numIndividuals = 0
			for line in dataFile:
				if line.strip():
					numIndividuals += 1

			numIndividuals -= 2
			
			# the default model for each locus is False, meaning null alleles are not considered at that locus
			self.locus_list = marker_names
			locusModelList = []
			for locus in self.locus_list:
				null_locus = False
				locusModelList.append(null_locus)
			self.locusModel = locusModelList
			
			self.numLociLabel.setText(str(len(marker_names)))
			self.numFamiliesLabel.setText(str(len(families)))
			self.numIndividualsLabel.setText(str(numIndividuals))

			self.addUI()
			self.centralWidget().show()
		else:
			pass
				
	def progress(self):
		#Progress Dialog
		if not self.numSteps:
			message = QtGui.QMessageBox(self)
			message.setWindowTitle("Error")
			message.setText("Please enter a value for the number of steps!")
			message.exec_()
			return

		if not self.numBurnInSteps:
			message = QtGui.QMessageBox(self)
			message.setWindowTitle("Error")
			message.setText("Please enter a value for the number of burn in steps!")
			message.exec_()
			return

		numSteps = self.numSteps
		self.progress = QtGui.QProgressDialog("Calculating...", "Stop", 0, numSteps, self)
		self.progress.setWindowModality(QtCore.Qt.WindowModal)
		self.progress.resize(300,125)
		self.progress.setWindowTitle("Calculating...")
		self.progress.show()
		self.runButton.hide()

		thread = BoriceThread(self, self.dataFileName, self.locusModel, self.numSteps, self.numBurnInSteps, self.outcrossingRateTuningParam, self.alleleFreqTuningParam, self.outcrossingRate, self.writeOutput2, self.writeOutput3, self.writeOutput4)

		thread.start()
		canceled = False

		while thread.getStep() != self.numSteps - 1:
			self.progress.setValue(100 * thread.getStep()/self.numSteps)
			QtGui.QApplication.instance().processEvents()			
			if(self.progress.wasCanceled()):
				thread.quit()
				canceled = True				
				break
		self.progress.setValue(self.numSteps)
		self.runButton.show()
	
		if canceled:
			return

		message = QtGui.QMessageBox(self)
		message.setWindowTitle("Complete!")
		message.setText("Done! Please look in your current working directory for output files.")
		message.exec_()
		return
		
class BoriceThread(QtCore.QThread):
	def __init__(self, parent, dataFileName, locusModel, numSteps, numBurnInSteps, outcrossingRateTuningParam, alleleFreqTuningParam, outcrossingRate, writeOutput2, writeOutput3, writeOutput4):
		super(BoriceThread, self).__init__(parent)
		self.dataFileName = dataFileName
		self.locusModel = locusModel
		print(self.locusModel)
		self.numSteps = numSteps
		self.numBurnInSteps = numBurnInSteps
		self.outcrossingRateTuningParam = outcrossingRateTuningParam
		self.alleleFreqTuningParam = alleleFreqTuningParam
		self.outcrossingRate = outcrossingRate
		self.writeOutput2 = writeOutput2
		self.writeOutput3 = writeOutput3
		self.writeOutput4 = writeOutput4
		#print self.writeOutput2, self.writeOutput3, self.writeOutput4
		self.Borice = BORICE.BORICE()

	def run(self):
		self.Borice.main(self.dataFileName, self.locusModel, self.numSteps, self.numBurnInSteps, self.outcrossingRateTuningParam, self.alleleFreqTuningParam, self.outcrossingRate, self.writeOutput2, self.writeOutput3, self.writeOutput4)

	def getStep(self):
		return self.Borice.getStep()
		
		
		
		
