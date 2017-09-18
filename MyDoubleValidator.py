#!/usr/bin/python2
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
#
# Modules
#
from PyQt4 import QtGui

class MyDoubleValidator(QtGui.QDoubleValidator):
	"""
	A stricter form of the double validator to limit user
	input in certain editable widgets
	"""
	def __init__(self, bottom, top, decimals, parent = None):
		QtGui.QDoubleValidator.__init__(self, bottom, top, decimals, parent)

	def validate(self, input, pos):
		state, pos = QtGui.QDoubleValidator.validate(self, input, pos)
		if input.isEmpty() or input == '.':
			return QtGui.QValidator.Intermediate, pos
		if state != QtGui.QValidator.Acceptable:
			return QtGui.QValidator.Invalid, pos
		return QtGui.QValidator.Acceptable, pos
