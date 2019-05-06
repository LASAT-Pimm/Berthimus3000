# SCRIPT BASE SUR LA V19

from __future__ import division
import tools
import collections
from collections import defaultdict
from collections import Counter

import GUIV3 as GUI
import operator
import weakref
import matplotlib.pyplot as plt
import shelve
import numpy as np
import time as zoipejfsdlkjnsd
import sys 
import os
import Tkinter as tk

error = 0
matNb = 0

##################################################

def codePatch(a = [], b = []):
	"""taking care of minor problems"""
	c= defaultdict()
	# for i in range(len(b)):
	c[tuple(b)] = []
	for j in range(len(a)):
		try:
			if set(b).issubset(a[j]):
				if len(b) == len(a[j]):
					pass
				else:
					c[tuple(b)].append(a[j])
		except:
			pass
	c[tuple(b)].sort(key=len)

	return c[tuple(b)]

def pressureInterface(motherId, motherCurve, calculated):
	""" When a single curve crosses at an intersection """
	# id of the interface where the mothercurve lands
	interEnd = motherCurve.interEnd
	intersecInterface = setting.interface[interEnd].position
	direction = motherCurve.direction
	daughterCurve0Id = curveHeritage[motherId][0][0]
	daughterCurve1Id = curveHeritage[motherId][1][0]
	daughter0 = curve[daughterCurve0Id]
	daughter1 = curve[daughterCurve1Id]

	if (direction == 1) :
		motherMat = setting.interface[interEnd].previousMat
		# daughterMat = material next to the mother where one of the shockwave will go
		daughterMat = setting.interface[interEnd].nextMat
	elif (direction == 0):
		motherMat = setting.interface[interEnd].nextMat
		daughterMat = setting.interface[interEnd].previousMat

	motherMatZ = material[motherMat].Z
	daughterMatZ = material[daughterMat].Z
	motherVel = motherCurve.velocity
	motherPres = motherCurve.pressure
	presSign = motherPres/abs(motherPres)

	# pressure of the outgoing curve that has the same direction as the mothercurve
	p1 = [motherVel, motherPres]
	p2 = [motherVel + 1.E9, presSign*(abs(motherPres) - motherMatZ*1.E9)]
	l1 = tools.line(p1,p2)
	# next material hugoniot
	p3 = [0.,0.]
	p4 = [1.E9,presSign*daughterMatZ*1E9]
	l2 = tools.line(p3,p4)
	# intersection of both hugoniot
	newPoint = tools.intersection(l1,l2)
	daughter1.velocity = newPoint[0]
	daughter1.pressure = newPoint[1]
	calculated[daughterCurve1Id] = newPoint
	presSign = daughter1.pressure/abs(daughter1.pressure)
	# pressure of the outgoing curve that has the opposite direction than the mothercurve
	p1 = [daughter1.velocity, daughter1.pressure]
	p2 = [daughter1.velocity + 1.E9, daughter1.pressure + presSign*motherMatZ*1.E9]
	l1 = tools.line(p1,p2)
	p3 = [0., 0.]
	p4 = [1.E9, presSign*(-motherMatZ)*1.E9]
	l2 = tools.line(p3,p4)
	newPoint = tools.intersection(l1,l2)
	calculated[daughterCurve0Id] = newPoint
	daughter0.velocity = newPoint[0]
	daughter0.pressure = newPoint[1]
	
	return

def pressureMultInt(point, inc):# value of pressure and velocity for incomming curves, list inc curves, list outgoing curves
	""" When multiple curve crosses at an intersection """
	pres1 = point[1] + 1E-13
	vel1 = point[0]
	direction = curve[inc[0]].direction

	refCurve = curve[inc[0]]
	interEnd = refCurve.interEnd

	if (direction == 1) :
		mat1 = setting.interface[interEnd].previousMat
		mat2 = setting.interface[interEnd].nextMat
	elif (direction == 0):
		mat1 = setting.interface[interEnd].nextMat
		mat2 = setting.interface[interEnd].previousMat

	mat1Z = material[mat1].Z
	mat2Z = material[mat2].Z
	presSign = pres1/abs(pres1)

	p1 = [vel1, pres1]
	if vel1 >= 0:
		p2 = [vel1 + 1.E9, presSign*(abs(pres1) - mat1Z*1.E9)]
	else:
		p2 = [vel1 + 1.E9, presSign*(abs(pres1) + mat1Z*1.E9)]
	l1 = tools.line(p1,p2)
	p3 = [0.,0.]
	if vel1 >= 0 :
		p4 = [1.E9,presSign*mat2Z*1.E9]
	else:
		p4 = [-1.E9,presSign*mat2Z*1.E9]
	l2 = tools.line(p3,p4)
	temp1 = tools.intersection(l1,l2)

	p1 = [temp1[0], temp1[1]]
	if vel1 >=0 :
		p2 = [temp1[0] + 1.E9, temp1[1] + presSign*mat1Z*1.E9]
	else :
		p2 = [temp1[0] + 1.E9, temp1[1] - presSign*mat1Z*1.E9]
	l1 = tools.line(p1,p2)
	# next material hugoniot
	p3 = [0., 0.]
	if vel1 >= 0 :
		p4 = [1.E9, presSign*(-mat1Z)*1.E9]
	else:
		p4 = [1.E9, presSign*(mat1Z)*1.E9]
	l2 = tools.line(p3,p4)
	temp2 = tools.intersection(l1,l2)

	if (direction == 0) :
		if refCurve.cType == 0 :
			if mat1Z < mat2Z:
				if abs(temp2[1]) < abs(temp1[1]):
					topL = temp1
					topR = temp2
				else:
					topL = temp2
					topR = temp1
			else:
				if temp2[1]/pres1 < 0:
					topL = temp1
					topR = temp2
				else:
					topL = temp2
					topR = temp1
		else :
			if mat1Z < mat2Z:
				if abs(temp2[1]) < abs(temp1[1]):
					topL = temp1
					topR = temp2
				else:
					topL = temp2
					topR = temp1
			else:
				if temp2[1]/pres1 < 0:
					topL = temp1
					topR = temp2
				else:
					topL = temp2
					topR = temp1
	elif (direction == 1):
		if refCurve.cType == 0 :
			if mat1Z < mat2Z:
				if abs(temp2[1]) < abs(temp1[1]):
					topR = temp1
					topL = temp2
				else:
					topR = temp2
					topL = temp1
			else:
				if temp2[1]/pres1 < 0:
					topR = temp1
					topL = temp2
				else:
					topR = temp2
					topL = temp1
		else :
			if mat1Z < mat2Z:
				if abs(temp2[1]) < abs(temp1[1]):
					topR = temp1
					topL = temp2
				else:
					topR = temp2
					topL = temp1
			else:
				if temp2[1]/pres1 < 0:
					topR = temp1
					topL = temp2
				else:
					topR = temp2
					topL = temp1

	return topL, topR 

def pressureIntersec(point1, refCurve1, point2): # refCurve1 is required just to get the direction and material
	pres1 = point1[1] + 1.E-20
	vel1 = point1[0] + 1.E-20
	pres2 = point2[1] + 1.E-20
	vel2 = point2[0] + 1.E-20

	interEnd1 = refCurve1.interEnd
	intersecInterface = setting.interface[interEnd1].position
	direction =refCurve1.direction

	if (direction == 1) :
		mat1 = setting.interface[interEnd1].previousMat
		mat2 = setting.interface[interEnd1].nextMat
	elif (direction == 0):
		mat1 = setting.interface[interEnd1].nextMat
		mat2 = setting.interface[interEnd1].previousMat

	mat1Z = material[mat1].Z
	mat2Z = material[mat2].Z
	presSign = pres1/abs(pres1)

	p1 = [vel1, pres1]
	p2 = [vel1 + 1.E9, abs(pres1) - mat1Z*1.E9]
	p3 = [vel2, pres2]
	p4 = [vel2 + 1.E9, abs(pres2) + mat2Z*1.E9]


	l1 = tools.line(p1,p2)
	l2 = tools.line(p3,p4)

	newPoint = tools.intersection(l1, l2)

	return newPoint

def pressureCross(point1, refCurve1, point2,refCurve2): # refCurve1 is required just to get the material
	pres1 = point1[1] + 1.E-13
	vel1 = point1[0] + 1.E-13
	pres2 = point2[1] + 1.E-13
	vel2 = point2[0] + 1.E-13

	interEnd1 = refCurve1.interEnd
	interEnd2 = refCurve2.interEnd

	mat1 = setting.interface[interEnd1].previousMat
	mat2 = setting.interface[interEnd2].nextMat

	mat1Z = material[mat1].Z
	mat2Z = material[mat2].Z
	presSign = pres1/abs(pres1)
	velSign = vel1/abs(vel1)

	p1 = [vel1, pres1]
	p2 = [vel1 + 1.E9, abs(pres1) - mat1Z*1.E9]
	p3 = [vel2, pres2]
	p4 = [vel2 + 1.E9, abs(pres2) + mat2Z*1.E9]
	
	if mat1Z == mat2Z:
		if pres1 / pres2 > 0:
			if vel1 < vel2:
				p2 = [vel1 + 1.E9, abs(pres1) + presSign*mat1Z*1.E9]
				p4 = [vel2 + 1.E9, abs(pres2) - presSign*mat2Z*1.E9]
			else:
				p2 = [vel1 + 1.E9, abs(pres1) - presSign*mat1Z*1.E9]
				p4 = [vel2 + 1.E9, abs(pres2) + presSign*mat2Z*1.E9]
		elif vel1 / vel2 > 0:
			if pres1 < pres2:
				p2 = [vel1 + 1.E9, abs(pres1) + presSign*mat1Z*1.E9]
				p4 = [vel2 + 1.E9, abs(pres2) - presSign*mat2Z*1.E9]
			else:
				p2 = [vel1 + 1.E9, abs(pres1) - presSign*mat1Z*1.E9]
				p4 = [vel2 + 1.E9, abs(pres2) + presSign*mat2Z*1.E9]
	l1 = tools.line(p1,p2)
	l2 = tools.line(p3,p4)

	newPoint = tools.intersection(l1, l2)



	return newPoint

def pressureParallel(curvelist):
	point = [0.0,0.0]
	for i in curvelist:
		try:
			point[1] += calculated[i][1]
			point[0] += calculated[i][0]
		except:
			pass
	return point


def pressureReflection(motherId, motherCurve, calculated):
	daughterCurveId = curveHeritage[motherId][0][0]
	daughter.pressure = -motherCurve.pressure
	daughter.velocity = motherCurve.velocity
	calculated[daughterCurveId] = [daughter.velocity, daughter.pressure]

	return

def pressureMultRef(point, inc):
	refCurve = curve[inc[0]]
	interEnd = refCurve.interEnd

	if (refCurve.direction == 1) :
		mat = setting.interface[interEnd].previousMat
	else:
		mat = setting.interface[interEnd].nextMat
	matZ = material[mat].Z
	presSign = point[1]/abs(point[1]+1.e-13)

	if (point[0] < 0) :
		p2 = [point[0] + 1., presSign*(abs(point[1]) + matZ)]
	else :
		p2 = [point[0] + 1., presSign*(abs(point[1]) - matZ)]
	p1 = [point[0], point[1]]
	p3 = [0,0]
	p4 = [1,0]
	l1 = tools.line(p1,p2)
	l2 = tools.line(p3,p4)

	newPoint = tools.intersection(l1, l2)

	return newPoint	


def separation(refCurve, curveList, pulseDuration):
	top = []
	bot = []
	for nb in curveList:
		if (curve[refCurve].end - pulseDuration <= curve[nb].origin) and (curve[nb].origin <= curve[refCurve].end + pulseDuration) :
			top.append(nb)
		else :
			bot.append(nb)

	return top, bot

def findMother(curveList):
	mother = []
	for i in range(len(curveList)) :
		mother.append(curve[curveList[i]].previousCurve)
	
	return sorted(mother)

def findDaughter(curveList):
	daughter = []
	for i in range(len(curveList)) :
		daughter.append(curveHeritage[curveList[i]][0][0])
		try :
			daughter.append(curveHeritage[curveList[i]][1][0])
		except:
			pass
	
	return daughter

def dirSeparation(curveList):
	right = []
	left = []
	for i in range(len(curveList)) :
		if curve[curveList[i]].direction == 1 :
			right.append(curveList[i])
		else :
			left.append(curveList[i])
	
	return sorted(left), sorted(right)

def initialParam(setup = "mono", delay = 0., pressure = []):
	curve = defaultdict()
	point = defaultdict()
	curveDB = defaultdict()
	pointMat = defaultdict()
	calculated = defaultdict()
	for ply in range(len(setting.plySeq)):
		pointMat[ply]=[]
	time = []
	curveId = 0
	#parameters for the initial curve coming from the left
	interIniLeft = 0
	tIniLeft = 0.0
	posIniLeft = setting.interface[interIniLeft].position
	interEndLeft = interIniLeft + 1
	posEndLeft = setting.interface[interEndLeft].position
	tEndLeft = posEndLeft/material[setting.plySeq[0]].soundVel
	# initial velocity left
	print setting.interface[interEndLeft].previousMat
	vel0 = pressure[0]/material[setting.interface[interEndLeft].previousMat].Z
	

	# adding left curve 
	pointMat[0].append([curveId,[posIniLeft,tIniLeft],[posEndLeft,tEndLeft]])
	time.append(tEndLeft)
	point[tEndLeft] = [(interEndLeft,tEndLeft)]
	curve[curveId] = Curve(interIniLeft, tIniLeft, interEndLeft, tEndLeft, -2, 0, 1, pressure[0], vel0)
	curveDB[(interEndLeft,tEndLeft)] = [curveId]
	calculated[curveId] = [vel0,pressure[0]]
	curveId += 1

	if (setup == "mono"):

		return pointMat, time, point, curve, curveDB, curveId, calculated

	if (setup == "sym"):
		#parameters for the initial curve coming from the right
		interIniRight = len(setting.plySeq)
		tIniRight = delay
		posIniRight = setting.interface[interIniRight].position
		interEndRight = interIniRight - 1
		posEndRight = setting.interface[interEndRight].position
		tEndRight = abs(setting.xmax-posEndRight)/material[setting.plySeq[-1]].soundVel + delay
		# initial velocity right
		vel1 = -pressure[1]/material[setting.interface[interEndRight].nextMat].Z

		pointMat[len(setting.plySeq) - 1].append([curveId,[posIniRight,tIniRight],[posEndRight,tEndRight]])
		if tEndRight not in time :
			time.append(tEndRight)
			point[tEndRight] = [(interEndRight,tEndRight)]
		else :
			point[tEndRight].append((interEndRight,tEndRight))
		curve[curveId] = Curve(interIniRight, tIniRight, interEndRight, tEndRight, -1, 0, 0, pressure[1], vel1)
		curveDB[(interEndRight,tEndRight)] = [curveId]
		calculated[curveId] = [vel1,pressure[1]]
		curveId += 1

		return pointMat, time, point, curve, curveDB, curveId, calculated

	if (setup == "double"):
		#parameters for the initial curve coming from the right
		interIniLeft = 0
		tIniLeft = delay
		posIniLeft = setting.interface[interIniLeft].position
		interEndLeft = interIniLeft + 1
		posEndLeft = setting.interface[interEndLeft].position
		tEndLeft = posEndLeft/material[setting.plySeq[0]].soundVel + tIniLeft
		# initial velocity left
		vel1 = pressure[1]/material[setting.interface[interEndLeft].previousMat].Z

		pointMat[0].append([curveId,[posIniLeft,tIniLeft],[posEndLeft,tEndLeft]])
		if tEndLeft not in time:
			time.append(tEndLeft)
			point[tEndLeft] = [(interEndLeft,tEndLeft)]
		else:
			point[tEndLeft].append((interEndLeft,tEndLeft))
		curve[curveId] = Curve(interIniLeft, tIniLeft, interEndLeft, tEndLeft, -1, 0, 1, pressure[1], vel1)
		try :
			curveDB[(interEndLeft,tEndLeft)].append(curveId)
		except:
			curveDB[(interEndLeft,tEndLeft)] = [curveId]
		calculated[curveId] = [vel1,pressure[1]]
		curveId += 1

		return pointMat, time, point, curve, curveDB, curveId, calculated

class Material(object):
	"""Def material properties"""
	instances = []

	def __init__(self, matNb, name, young, poisson, density, soundVel, minBound = -1, maxBound = -1):

		self.id = matNb
		self.name = name
		self.young = young
		self.poisson = poisson
		self.density = density
		self.soundVel = soundVel
		self.minBound = minBound
		self.maxBound = maxBound
		self.Z = density * soundVel
		self.__class__.instances.append(weakref.proxy(self))

	def bound(self, minBound, maxBound):
		self.minBound = minBound
		self.maxBound = maxBound

	def getParam(self):
		"""If minBound or maxBound are equal to -1, values where not set"""
		print "Material ID = " + str(self.id)
		print "Name = " + self.name
		print "Young Modulus = " + str(self.young)
		print "Poisson ration = " + str(self.poisson)
		print "Density = " + str(self.density)
		print "Sound speed = " + str(self.soundVel)
		print "Min bound = " + str(self.minBound)
		print "Max bound = " + str(self.maxBound)
		print "Impedance = " + str(self.Z)

class Interface(object):
	"""Create an interface, and allocate properties"""

	def __init__(self, id, position = -1337, behaviour = [-1, -1], previousMat = -1, nextMat = -1):
		self.id = id
		self.position = position
		self.behaviour = behaviour
		self.previousMat = previousMat
		self.nextMat = nextMat

	def pos(self, position):
		self.position = position

	def bhv(self, behaviour):
		# if behaviour = 1, transmission if behaviour = 0 reflection
		self.behaviour = behaviour

	def getParam(self):
		print "ID = " + str(self.id)		
		print "position = " + str(self.position)		
		print "behaviour = " + str(self.behaviour)		

class Intersection(object):
	"""Create an intersection"""

	def __init__(self, matList, interfaceList):
		self.id = id
		self.position = position
		self.time = time
		self.interfaceNb = interfaceNb

	def getParam(self):
		print "ID = " + str(self.id)
		print "Position = " + str(self.position)
		print "Time = " + str(self.time)
		print "interfaceNb = " + str(self.interfaceNb)

class Setting(object):
	"""Setting up the ply sequence, boundaries, ..."""

	def __init__(self):
		super(Setting, self).__init__()
		

	def plySequence(self, plies):
		"""Sets the ply sequence based on material ID"""
		self.plySeq = plies
		# for ply in arg :
		# 	if ply in material.keys() :
		# 		self.plySeq.append(ply)
		# 	else:
		# 		print "Material ID " + str(ply) + " does not exist"
		# 		self.plySeq.append(-1)
		# 		break
		return self.plySeq

	def interfaceProp(self,previousMat,nextMat):
		if (previousMat >= nextMat):
			inter = [0, 1]
		else :
			inter = [1, 0]

		return inter

	def plyThickness(self, tckness):
		global error
		if (len(tckness) < len(self.plySeq)):
			print "Error, not enough thicknesses defined"
			error = 1
			return
		elif (len(tckness) > len(self.plySeq)):
			print "Error, too many thicknesses defined"
			error = 2
			return
		elif -1 in self.plySeq:
			print "Error in the definition of the ply sequence"
			error = 3
			return
		else :
			self.plyThick = []
			position = 0.
			k = 0
			self.pos2int = defaultdict()
			self.xmax = 0.
			self.interface = defaultdict()
			self.interface[0] = Interface(0,position,[1,1],-1,self.plySeq[0])
			self.pos2int[0] = 0.
			for ply in tckness :
				position += ply
				self.plyThick.append(ply)
				if (k == len(self.plySeq)-1):
					bhv = [1, 1]
					self.interface[k+1] = Interface(k+1,position,bhv,material[self.plySeq[(k)]].id)
					self.pos2int[position] = k+1
				else :
					bhv = self.interfaceProp(material[self.plySeq[k]].Z, material[self.plySeq[k+1]].Z)
					self.interface[k+1] = Interface(k+1,position,bhv,material[self.plySeq[(k)]].id,  material[self.plySeq[(k+1)]].id)
					self.pos2int[position] = k+1
					k += 1
			self.xmax = position	
			return self.plyThick, self.interface, self.xmax, self.pos2int

	def getParam(self):
		global error
		if error != 0 :
			print "Whoups error " + str(error)
		else :
			print "Ply sequence : " + str(self.plySeq)
			print "Ply thickness : " + str(self.plyThick)
			for key in self.interface.keys():
				print key
				print "Interface nb: " + str(self.interface[key].id) + ", behaviour = " + str(self.interface[key].behaviour) + " and position = " + str(self.interface[key].position)

class Curve(object):
	"""Creates a curve with its properties"""

	def __init__(self, interIni, origin, interEnd, end, previousCurve, cType, direction, pressure = 1.e30, velocity = 1.e30):
		"""cType = 1 release, = 0 shock
	   		direction = 1 right, direction = 0 left""" 
		super(Curve, self).__init__()
		self.interIni = interIni
		self.origin = origin
		self.interEnd = interEnd
		self.end = end
		self.previousCurve = previousCurve
		self.cType = cType
		self.direction = direction
		self.pressure = pressure
		self.velocity = velocity
############################################################################################################################################################################################
############################################################################################################################################################################################
############################################################################################################################################################################################
############################################################################################################################################################################################
############################################################################################################################################################################################
############################################################################################################################################################################################
############################################################################################################################################################################################
############################################################################################################################################################################################

################### GUI PART #####################

xstart = GUI.GUI()
winSelect = xstart.launch()
if winSelect == 0 :
	simulation = GUI.Simulation()
	FOLDER, OUTPUT, PULSE, P1, P2, DELAY, DURATION, XRES, TRES, MATSEQ, MATERIAL, MATTCK, TYPESIM = simulation.launch()
	path = os.path.normpath(FOLDER + "/" + OUTPUT)

	##################################################
	presIni = 1.0

	# mat def
	global material
	material = defaultdict()
	for key in MATERIAL.keys():
		mat = MATERIAL[key][0]
		material[mat[0]] = Material(mat[0], mat[1],mat[2],mat[3],mat[4],mat[5])


	# material[matNb] = Material("Alu",72.0,0.36,2700.0,5300.0) #1
	# material[matNb] = Material("Cu",112.0,0.36,8930.0,3933.0) #2
	# material[matNb] = Material("Colle",5.2,0.3,1260.,2600.) #3
	# material[matNb] = Material("CarboneL",12.6,0.3,1630.,3000.) #4
	# material[matNb] = Material("Carbone3Plis",10.4,0.31,1710.,2970) #5
	# material[matNb] = Material("Polymere",1.,1.,950.,936.) #6
	# material[matNb] = Material("Ceramique",1.4,0.31,5770.,1927.) #7

	setting = Setting()

	setting.plySequence(MATSEQ)
	setting.plyThickness(MATTCK)
	# setting.getParam() 



	timeLimit = DURATION
	pulseDuration = PULSE
	xRes = XRES	
	tRes = TRES



	# intialisation
	curve = defaultdict()
	point = defaultdict()
	time = []
	treatedCurve = []
	pointMat = defaultdict()
	for ply in range(len(setting.plySeq)):
		pointMat[ply]=[]
	# Associate a "mother curve" to its descendents: curveHeritage[orinCurve] = [descendent1, descendent2, ...] Normally two descendent max.
	# An indice is link to each curve : 1 daughter curve same direction as mother, else 0
	curveHeritage = defaultdict()
	curveDB = defaultdict()
	calculated = defaultdict()
	curveId = 0

	pointMat, time, point, curve, curveDB, curveId, calculated = initialParam(setup = TYPESIM, delay = DELAY, pressure = [P1,P2])

	loop = 1
	while (time[0] < timeLimit) :
		print " TIME = ", time[0]/timeLimit*100., " %          \r",

		for nb in range(len(point[time[0]])):
			interfaceNb = point[time[0]][nb][0] # ID of the current interface
			interCoord = point[time[0]][nb][1] # y coord of the current  inersection point
			interCurrent = setting.interface[interfaceNb].position # position current interface position
			# looping through initial curves
			index = 0
			while curveDB[(interfaceNb,interCoord)][index] in treatedCurve :
				index += 1
			curveIni = curveDB[(interfaceNb,interCoord)][index] # gives id of the curve

			if interCurrent > 0 :
				# parameters for the curve class
				matLeft = setting.interface[interfaceNb].previousMat # material on the left of the current interface
				interLeft = setting.interface[interfaceNb - 1].position # position next interface on the left
				endLeft = tools.interLoc(material[matLeft].soundVel,interCurrent,interLeft,interCoord) # time coordinate of the intersection between the curve and the left interface
				# generation of curve created at the intersection point
				# curve going to the left :
				currentCurveDir = 0
				cType = tools.bounceDet(curve[curveIni].cType,setting.interface[interfaceNb].behaviour,curve[curveIni].direction,currentCurveDir) 
				curve[curveId] = Curve(interfaceNb, interCoord, interfaceNb - 1, endLeft, curveIni, cType, currentCurveDir)
				# adding the curve to the database
				try:
					curveDB[(interfaceNb-1,endLeft)].append(curveId)
				except:
					curveDB[(interfaceNb-1,endLeft)] = [curveId]
				# verifying this time has not already been stored
				if endLeft not in time : 
					# adding the intersection to the time list
					time.append(endLeft)
				# adding point to list of points
				if endLeft in point.keys():
					point[endLeft].append((interfaceNb-1, endLeft))
				else:
					point[endLeft] = [(interfaceNb-1, endLeft)]
				# updating heritage curve list
				try:
					if curve[curveIni].direction == currentCurveDir: 
						curveHeritage[curveIni].append([curveId,1])
					else:
						curveHeritage[curveIni].append([curveId,0])
				except:
					if curve[curveIni].direction == currentCurveDir: 
						curveHeritage[curveIni] = [[curveId,1]]
					else:
						curveHeritage[curveIni] = [[curveId,0]]
				# adding point to the point list
				pointMat[interfaceNb-1].append([curveId,[interCurrent,interCoord],[interLeft,endLeft]])
				curveId += 1

			if interCurrent < setting.xmax :
				matRight = setting.interface[interfaceNb].nextMat # material on the right of the current interface
				interRight = setting.interface[interfaceNb + 1].position # position next interface on the right
				endRight = tools.interLoc(material[matRight].soundVel,interCurrent,interRight,interCoord) # time coordinate of the intersection between the curve and the right interface
				# generation of both curve created at the intersection point
				# curve going to the right :
				currentCurveDir = 1 # direction of the curve, 1 for right and 0 for left
				cType = tools.bounceDet(curve[curveIni].cType,setting.interface[interfaceNb].behaviour,curve[curveIni].direction,currentCurveDir) # type of curve, 0 for shock and 1 for release
				curve[curveId] = Curve(interfaceNb, interCoord, interfaceNb + 1, endRight, curveIni, cType, currentCurveDir) # creation + indexation of the curve
				# adding the curve to the database
				try:
					curveDB[(interfaceNb+1,endRight)].append(curveId)
				except:
					curveDB[(interfaceNb+1,endRight)] = [curveId]
				# verifying this time has not already been stored
				if endRight not in time : 
					# adding the intersection to the time list
					time.append(endRight)
				# adding point to list of points
				if endRight in point.keys():
					point[endRight].append((interfaceNb+1, endRight))
				else:
					point[endRight] = [(interfaceNb+1, endRight)]
				# updating heritage curve list
				try:
					if curve[curveIni].direction == currentCurveDir: 
						curveHeritage[curveIni].append([curveId,1])
					else:
						curveHeritage[curveIni].append([curveId,0])
				except:
					if curve[curveIni].direction == currentCurveDir: 
						curveHeritage[curveIni] = [[curveId,1]]
					else:
						curveHeritage[curveIni] = [[curveId,0]]
				# adding point to the point list
				pointMat[interfaceNb].append([curveId,[interCurrent,interCoord],[interRight,endRight]])
				curveId += 1
			treatedCurve.append(curveIni) 

		# removing the first time of the time list
		time.pop(0)
		# sorting time list to have the next one on top
		time.sort()
		loop += 1

	for key in pointMat.keys():
		pointMat[key].sort(key=operator.itemgetter(1))
	for key in curveHeritage.keys():
		curveHeritage[key].sort(key=operator.itemgetter(1))
	 # computation pressure
	treatedCurve.sort()
	pressure = []

	print " END CURVE LOOP"


	sorted_curveDB = sorted(curveDB.items(), key=operator.itemgetter(0))

	pouet = defaultdict()
	curveNbList = []
	for key in curveDB.keys():
		localDirection = [] # empty list where direction will be stored. Goal: make sur there is only one cruve per direction (avoid curve overlapping)
		temp = [] # list where curve are temporally stored
		for i in range(len(curveDB[key])):
			curveDirection = curve[curveDB[key][i]].direction
			if curveDirection not in localDirection:
				localDirection.append(curveDirection)
				temp.append(curveDB[key][i])
				curveNbList.append(curveDB[key][i])
		curveDB[key] = temp

	sorted_curveDB = sorted(curveDB.items(), key=operator.itemgetter(0))


	#### Ploting curve with number #####
	# for nb in curveNbList:
	# # for nb in [34,35,29,28,31,30,25,14]:
	# 	if curve[nb].cType == 0:
	# 		c = 'r'
	# 	else :############3 -33.0130000114

	# 		c = 'b'
	# 	plt.plot([setting.interface[curve[nb].interIni].position,setting.interface[curve[nb].interEnd].position],[curve[nb].origin,curve[nb].end], color = c)
	# 	plt.text((setting.interface[curve[nb].interIni].position+setting.interface[curve[nb].interEnd].position)/2,(curve[nb].origin+curve[nb].end)/2,str(nb))
	# plt.show()
	# sys.exit()
	###################################

	globGrid = Counter()


	globGrid = tools.grid(xRes, tRes)
	globGrid2 = tools.grid(xRes, tRes)
	i = 0
	start = zoipejfsdlkjnsd.time()
	# print pointMat
	print " START CURVE DISCRETISATION "
	for matId in pointMat.keys():
		# print "matId = " + str(matId)
		print " DISCRETISATION = ", i/len(pointMat.keys())*100, " %           \r",
		j = 0
		for couple in range(len(pointMat[matId])):
			print "------------ ", j/len(pointMat[matId])*100, " %           \r",
			initialPoint = [pointMat[matId][couple][1][0],pointMat[matId][couple][1][1]]
			endPoint = [pointMat[matId][couple][2][0],pointMat[matId][couple][2][1]]
			tools.disc(initialPoint, endPoint, pulseDuration, xRes, tRes, setting.xmax, time[0], pointMat[matId][couple][0], globGrid)
			j += 1
		i += 1

	vmax = globGrid[max(globGrid.iteritems(), key=operator.itemgetter(1))[0]]
	vmin = globGrid[min(globGrid.iteritems(), key=operator.itemgetter(1))[0]]

	print " END CURVE DISCRETISATION "

	stop = zoipejfsdlkjnsd.time() - start

	##### Create an array listing all elements #####

	eltList = []
	for key in globGrid.keys():
		eltList.append(key)

	eltList.sort(key = operator.itemgetter(1, 0)) # ordering the element by increasing x and increasing t




	##### Finding interface snapped coordinates #####

	snappedInter = []
	for inter in setting.interface.keys():
		xPos = setting.interface[inter].position
		snappedInter.append(int(xPos*float(xRes)/setting.xmax))


	##### List all curve combinasion that can be found #####


	i = 0
	j = 0

	uniqueElt = []
	print " FINGIND UNIQUE ELEMENTS "
	for key in eltList:
		if key[0] not in snappedInter:
			if len(globGrid[key]) > 0:
				if globGrid[key] not in uniqueElt:
					uniqueElt.append(globGrid[key])

	##### Separte curve incoming / out going #####
	outgoing = defaultdict()
	for i in range(len(uniqueElt)):
		outgoing[i] = []
	parallel = []
	cross = []
	inter = []

	for group in range(len(uniqueElt)):
		uniqueElt[group].sort()
		if len(uniqueElt[group]) != 1 :
			refCurve = np.amin(uniqueElt[group])
			for curveNb in uniqueElt[group]:
				if curve[curveNb].interIni == curve[refCurve].interEnd:
					outgoing[group].append(curveNb)
			if len(outgoing[group]) == 0:
				outgoing[group] = uniqueElt[group]
			if len(outgoing[group]) == len(uniqueElt[group]): # it means it is not at an interface or an intersection
				parallel.append(uniqueElt[group])
			if uniqueElt[group] not in parallel:
				temp = []
				for curveNb in uniqueElt[group]:
					if curve[curveNb].origin == curve[refCurve].end:
						temp.append(curveNb)
				if len(temp) != 0 :
					inter.append(uniqueElt[group])
				else :
					cross.append(uniqueElt[group])
		else :
			outgoing[group] = uniqueElt[group]

	print " END CURVE SORTING "
	# print "inter ", inter
	# print "cross ", cross
	# print "parallel ", parallel
	# sys.exit()
	print " STARTING PRESSURE CALCULATION "
	for i in range(len(uniqueElt)):
		" POURCENTATGE = ", i/len(uniqueElt)*100., " %           \r",
		group = uniqueElt[i]
		try :
			if group in inter:
				refCurve = group[0]
				out, inc = separation(refCurve, group, pulseDuration/10.)
				bot = findMother(out)
				botR, botL = dirSeparation(bot)
				top = findDaughter(bot)
				topL, topR = dirSeparation(top) # swap to keep R/L as indices to specify one SIDE of the interface (and not a direction)
				if len(out) == len(top) : # describes a reflection on one of the free surface
					if (len(out) == 1) and (len(inc) == 1) :
						calculated[out[0]] = [calculated[inc[0]][0], -calculated[inc[0]][1]]
					elif len(out) == 1 :
						calculated[out[0]] = [calculated[tuple(inc)][0], -calculated[tuple(inc)][1]]
					elif len(inc) == 1 :
						calculated[tuple(out)] = [calculated[inc[0]][0], -calculated[inc[0]][1]]
					else :
						calculated[tuple(out)] = [calculated[tuple(inc)][0], -calculated[tuple(inc)][1]]

					if len(bot) == 1 :
						calculated[tuple(group)] = pressureMultRef(calculated[bot[0]],bot)
					else :
						calculated[tuple(group)] = pressureMultRef(calculated[tuple(bot)],bot)
				elif (len(botR) == 0) or (len(botL) == 0): #generation of 2 waves after incontering an interface
					temp1 = bot + top
					temp1.sort()
					tot = tuple(temp1)
					if len(bot) == 1 :
						calculated[topL[0]], calculated[topR[0]] = pressureMultInt(calculated[bot[0]], bot)
						if curve[bot[0]].direction == 1 :
							calculated[tuple(tot)] = calculated[topR[0]]
							calculated[tuple(group)] = calculated[topR[0]]
						else :
							calculated[tuple(tot)] = calculated[topL[0]]
							calculated[tuple(group)] = calculated[topL[0]]
					else :
						calculated[tuple(topL)], calculated[tuple(topR)] = pressureMultInt(calculated[tuple(bot)], bot)
						if curve[bot[0]].direction == 1 :
							calculated[tuple(tot)] = calculated[tuple(topL)]
							calculated[tuple(group)] = calculated[tuple(topL)]
						else :
							calculated[tuple(tot)] = calculated[tuple(topR)]
							calculated[tuple(group)] = calculated[tuple(topR)]
				else :
					temp1 = tuple(botL + topL)
					temp2 = tuple(botR + topR)
					temp3 = bot + top
					temp3.sort()
					tot = tuple(temp3)
					topL1 = [0.,0.]
					topR1 = [0.,0.]
					topL2 = [0.,0.]
					topR2 = [0.,0.]
					if (len(botL) == 1) and (len(botR) == 1):
						calculated[temp1] = pressureCross(calculated[botL[0]], curve[botL[0]], calculated[botR[0]],curve[botR[0]])
						calculated[temp2] = calculated[temp1]
						calculated[tot] = calculated[temp1]
					elif (len(botL) == 1):
						calculated[temp1] = pressureCross(calculated[botL[0]], curve[botL[0]], calculated[tuple(botR)], curve[botR[0]])
						calculated[temp2] = calculated[temp1]
						calculated[tot] = calculated[temp1]
					elif (len(botR) == 1):
						calculated[temp1] = pressureCross(calculated[tuple(botL)], curve[botL[0]], calculated[botR[0]], curve[botR[0]])
						calculated[temp2] = calculated[temp1]
						calculated[tot] = calculated[temp1]
					else :
						calculated[temp1] = pressureCross(calculated[tuple(botL)], curve[botL[0]], calculated[tuple(botR)], curve[botR[0]])
						calculated[temp2] = calculated[temp1]
						calculated[tot] = calculated[temp1]

					sortedBot = sorted(bot)
					try: 
						if (len(botL) == 1):
							topL1, topR1 = pressureMultInt(calculated[botL[0]], botL)
						else:
							topL1, topR1 = pressureMultInt(calculated[tuple(botL)], botL)
					except:
						pass
					try:
						if (len(botR) == 1):
							topL2, topR2 = pressureMultInt(calculated[botR[0]], botR)
						else:
							topL2, topR2 = pressureMultInt(calculated[tuple(botR)], botR)
					except:
						pass
					if (len(topL) == 1) and (len(topR) == 1) :
						calculated[topL[0]] = [topL1[0]+topL2[0], topL1[1]+topL2[1]]
						calculated[topR[0]] = [topR1[0]+topR2[0], topR1[1]+topR2[1]]
					elif len(topR) == 1 :
						calculated[tuple(topL)] = [topL1[0]+topL2[0], topL1[1]+topL2[1]]
						calculated[topR[0]] = [topR1[0]+topR2[0], topR1[1]+topR2[1]]
					elif len(topL) == 1:
						calculated[topL[0]] = [topL1[0]+topL2[0], topL1[1]+topL2[1]]
						calculated[tuple(topR)] = [topR1[0]+topR2[0], topR1[1]+topR2[1]]
					else :
						calculated[tuple(topL)] = [topL1[0]+topL2[0], topL1[1]+topL2[1]]
						calculated[tuple(topR)] = [topR1[0]+topR2[0], topR1[1]+topR2[1]]

			elif group in cross:
				left = []
				right = []
				left, right = dirSeparation(group)
				# print "GROUP = ", group
				# if group == [24, 38]:
				# 	print "CALCULATED ", calculated
				# 	print "CALCULATED[24]", calculated[24]
				# 	print "CALCULATED[38]", calculated[38]
				if (len(right) == 1) and (len(left) == 1):
					# if group == [4, 10]:
						# print "LEFT = ", left
						# print "RIGHT = ", right
						# print "calculated[right[0]] ", calculated[right[0]]
						# print "calculated[left[0]] ", calculated[left[0]]
					calculated[tuple(group)] = pressureCross(calculated[right[0]], curve[right[0]], calculated[left[0]], curve[left[0]])
				elif (len(right) == 1) :
					calculated[tuple(group)] = pressureCross(calculated[right[0]], curve[right[0]], calculated[tuple(left)], curve[left[0]])
				elif (len(left) == 1):
					calculated[tuple(group)] = pressureCross(calculated[tuple(right)], curve[right[0]], calculated[left[0]], curve[left[0]])
				else :
					calculated[tuple(group)] = pressureCross(calculated[tuple(right)], curve[right[0]], calculated[tuple(left)], curve[left[0]])
			elif group in parallel:
				if tuple(group) not in calculated.keys():
					# print "POEUTPOUETPOUET"
					# print "group ", group
					calculated[tuple(group)] = pressureParallel(group)
			# else:
				# print "GROUPGROUPGROUP ", group
			# try:
			# 	print "GROUP ", group
			# 	print "CALCULaTED[4] ", calculated[4]
			# except:
			# 	pass
			# 	try:
			# 		print "CALCULATED ", calculated[tuple(group)]
			# 	except:
			# 		print "CALCULATED ", calculated[group[0]]
			# except:
			# 	print "###################"
		except:
			print "## GROUP = ", group
			pass

	print " END PRESSURE CALCULATION"
	patch = defaultdict()
	problem = []
	ok = []
	j = 0

	print " BEGINNING GRID UPDATE "
	# while j<100 :
	for key in globGrid.keys():
		if len(globGrid[key]) > 1 :
			sortedkey = globGrid[key]
			sortedkey.sort()
			if sortedkey not in ok:
				ok.append(sortedkey)
			try :
				globGrid2[key] = calculated[tuple(sortedkey)]
			except:
				try:
					if sortedkey not in problem:
						problem.append(sortedkey)
						equivList = codePatch(calculated.keys(), sortedkey)
						patch[tuple(sortedkey)] = equivList[0]
					globGrid2[key] = calculated[tuple(patch[tuple(sortedkey)])]
				except:
					globGrid2[key] = [0.,0.]
		elif len(globGrid[key]) == 1 :
			try:
				globGrid2[key] = calculated[globGrid[key][0]]
			except:
				globGrid2[key] = [0.,0.]
		else :
			globGrid2[key] = [0.,0.]



	vmin, vmax = tools.extremity(calculated)


	d = shelve.open(path + ".HUGO")
	d['globGrid2'] = globGrid2
	d['vmin'] = vmin
	d['vmax'] = vmax
	d['time'] = timeLimit
	d['position'] = setting.xmax
	d['interface'] = setting.interface
	d.close()

	tools.discPlot(globGrid2,xRes, tRes, timeLimit, setting.xmax, vmin, vmax, setting.interface)

else:
	read = GUI.ReadHUGO()
	globGrid2,xRes, tRes, timeLimit, xmax, vmin, vmax, interface = read.launch()
	print "LOADING PLEASE WAIT ......"
	tools.discPlot(globGrid2,xRes, tRes, timeLimit, xmax, vmin, vmax, interface)