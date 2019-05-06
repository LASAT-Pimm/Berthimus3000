
import scipy.misc as smp
import scipy
import matplotlib
from matplotlib import pyplot as plt
from matplotlib import rc
import time as sdfjsdlkfj
import numpy as np
from matplotlib.patches import Rectangle  
from collections import defaultdict
from collections import Counter
from decimal import Decimal


def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
    '''
    Function to offset the "center" of a colormap. Useful for
    data with a negative min and positive max and you want the
    middle of the colormap's dynamic range to be at zero

    Input
    -----
      cmap : The matplotlib colormap to be altered
      start : Offset from lowest point in the colormap's range.
          Defaults to 0.0 (no lower ofset). Should be between
          0.0 and `midpoint`.
      midpoint : The new center of the colormap. Defaults to 
          0.5 (no shift). Should be between 0.0 and 1.0. In
          general, this should be  1 - vmax/(vmax + abs(vmin))
          For example if your data range from -15.0 to +5.0 and
          you want the center of the colormap at 0.0, `midpoint`
          should be set to  1 - 5/(5 + 15)) or 0.75
      stop : Offset from highets point in the colormap's range.
          Defaults to 1.0 (no upper ofset). Should be between
          `midpoint` and 1.0.
    '''
    cdict = {
        'red': [],
        'green': [],
        'blue': [],
        'alpha': []
    }

    # regular index to compute the colors
    reg_index = np.linspace(start, stop, 257)

    # shifted index to match the data
    shift_index = np.hstack([
        np.linspace(0.0, midpoint, 128, endpoint=False), 
        np.linspace(midpoint, 1.0, 129, endpoint=True)
    ])

    for ri, si in zip(reg_index, shift_index):
        r, g, b, a = cmap(ri)

        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))

    newcmap = matplotlib.colors.LinearSegmentedColormap(name, cdict)
    plt.register_cmap(cmap=newcmap)

    return newcmap


def __init__(self, arg):
	super(tools, self).__init__()
	self.arg = arg
	

def bounceDet(cType, interBhv, prevCurvDir, currentCurveDir):
	"""defines the type of curve created at an interface"""
	if prevCurvDir == currentCurveDir :
		output = cType
		# else :
	else :
		if interBhv[prevCurvDir] == 0 :
			output = 1 if cType == 1 else 0
		else :
			output = 1 if cType == 0 else 0


	return output

def interLoc(matVel, posIntStart, posIntEnd, pointIni):
	"""intersection between curve and interface"""
	end = abs(posIntEnd - posIntStart)/matVel + pointIni

	return end

def line(p1, p2):
	"""simple line definition"""
	A = (p1[1] - p2[1])
	B = (p2[0] - p1[0])
	C = (p1[0]*p2[1] - p2[0]*p1[1])

	return A, B, -C

def intersection(L1, L2):
	"""intersectino of two lines"""
	D  = L1[0] * L2[1] - L1[1] * L2[0]
	Dx = L1[2] * L2[1] - L1[1] * L2[2]
	Dy = L1[0] * L2[2] - L1[2] * L2[0]
	if D != 0:
		x = Dx / D
		y = Dy / D
		return x,y
	else:
		return False

def genGrid(xMax, tMax, xRes, tRes):
	"""generating x and t coordinate for the grid"""

	xDiv = int(xMax/xRes)+1
	tDiv = int(tMax/tRes)+1

	X = np.linspace(0, xMax, xDiv)
	T = np.linspace(0, tMax, tDiv)
	xv, tv = np.meshgrid(X,T)

	grid = np.dstack([xv.ravel(),tv.ravel()])[0]

	return grid

def snapGrid(grid,points = []):
	"""snapping a point to a grid"""
	mytree = scipy.spatial.cKDTree(grid)
	points_list = list(points.transpose())
	dist, indexes = mytree.query(points_list)

	return indexes

def snapGrid2(xRes, tRes, xMax, tMax, point =[]):
	indexes = list()
	indexes.append([int(point[0]/xMax*xRes), int(point[1]/tMax*tRes)])

	return indexes

def disc(initialPoint, endPoint, pulseDuration, xRes, tRes, xMax, tMax, value, globGrid):
	""" adds discret line to the global grid"""
	snapped = snapGrid2(xRes, tRes, xMax, tMax, initialPoint)
	x0 = snapped[0][0]
	t0 = snapped[0][1]
	snapped = snapGrid2(xRes, tRes, xMax, tMax, endPoint)
	xf = snapped[0][0]
	tf = snapped[0][1]

	length = abs(xf-x0)
	direction = (xf-x0)/length
	slope = (float(tf)-float(t0))/(float(xf)-float(x0))
	pulseDiv = int(pulseDuration*float(tRes)/tMax)

	error = 0
	temp = []
	start = sdfjsdlkfj.time()
	# TIMEA = sdfjsdlkfj.time()
	for i in range(0,length+1):
		try:
			position = x0+(i)*direction
			time = int(t0+direction*(i)*slope)
			if (position,time) not in temp:
				globGrid[position,time].append(value)
				temp.append((position,time))
				# print "############1", sdfjsdlkfj.time() - TIMEA 
				# TIMEB = sdfjsdlkfj.time()
				
				for j in range(1,pulseDiv):
					if (position,time+j) not in temp:
						try:
							globGrid[position,time + j].append(value)
							temp.append((position,time+j))
						except:
							error += 1
				# print "############2", sdfjsdlkfj.time() - TIMEB 
				# TIMEC = sdfjsdlkfj.time()

		except:
			pass

	# print "############3", sdfjsdlkfj.time() - TIMEC 
	# sys.exit()
	return


def isReflection(curveHeritage, motherCurve):
	if len(curveHeritage[motherCurve]) == 1:
		reflection = True
	else:
		reflection = False

	return reflection

def separate():
	pass

def grid(xRes, tRes):
	grid = Counter()
	for j in range(tRes + 1):
		for i in range(xRes + 1):
			grid[i,j] = []

	return grid

def genGridDic(xRes, tRes, tMax, xMax):
	i = 0. 
	j = 0. # initialising indices for the 2D matrix
	grid = defaultdict() # dictionnary that will act as matrix
	while j <= tMax : # looping as long as the last computed time is not part of the grid
		i = 0 # resetting row variable
		while i <= xMax : # going all the way to the last interface (back face)
			grid[float('{:0.2e}'.format(i)),float('{:0.2e}'.format(j))] = []
			i += xRes # increasing i by one x resolution
		j += tRes # increasing j by one tRes
	return grid

def extremity(dic = defaultdict()):
	valueRange = []
	for key in dic.keys():
		valueRange.append(dic[key][1])
	vmax = max(valueRange)
	vmin = min(valueRange)
	return vmin, vmax

def format_e(n):
    a = '%.2E' % n
    return a.split('E')[0].rstrip('0').rstrip('.') + 'E' + a.split('E')[1]

def discPlot(grid,xResolution, tResolution, tMax, xMax, vmin, vmax, interface):
	"""plot sum of discretized line"""
	discX = []
	discT = []
	z = []

	font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 12}

	matplotlib.rc('font', **font)

	for x, t in grid.keys():
		discX.append(x)
		discT.append(t)
		z.append(grid[(x,t)][1])

		# print z

	dx = [xResolution] * len(discX)
	dt = [tResolution] * len(discT)

	norm = plt.Normalize(vmin=vmin, vmax=vmax)
	midpoint = 1 - vmax/(vmax + abs(vmin))
	orig_cmap = matplotlib.cm.get_cmap('seismic')
	shifted_cmap = shiftedColorMap(orig_cmap, midpoint=midpoint, name='shifted')

	zi, yi, xi = np.histogram2d(discT, discX, bins=(xResolution,tResolution), weights=z, normed=False)
	counts, _, _ = np.histogram2d(discT, discX, bins=(xResolution,tResolution))
	zi = zi / counts
	zi = np.ma.masked_invalid(zi)

	fig, ax = plt.subplots()
	ax.pcolormesh(xi, yi, zi, edgecolors=None, cmap = shifted_cmap)
	scat = ax.scatter(discX, discT, c=z, s=0, cmap = shifted_cmap)
	fig.colorbar(scat)
	ax.margins(0.05)

	#### RECALAGE DES TICKS AVEC LES VRAIES VALEURS ####
	xticks = ["0"]
	tticks = ["0"]
	xlabels = [item.get_text() for item in ax.get_xticklabels()]
	tlabels = [item.get_text() for item in ax.get_yticklabels()]

	plt.locator_params(axis='x', nbins=5)
	for i in range(1, 6):
	    xticks.append(format_e(Decimal(float(i)*xResolution/float(4)*xMax/xResolution)))

	for i in range(1, len(tlabels)):
	    tticks.append(format_e(Decimal(float(i)*tResolution/float(len(tlabels)-1)*(tMax)/tResolution)))

	# print xticks
	# print tticks
	ax.set_xticklabels(xticks)
	ax.set_yticklabels(tticks)


	for key in interface.keys():
		x = interface[key].position * xResolution / xMax
		if (x != 0) or (x != xMax): 
			plt.plot([x, x], [0., tResolution], color='k', linestyle='-', linewidth=1)
	
	cid = fig.canvas.mpl_connect('button_press_event', lambda event: onclick(event, grid, tMax, xResolution, tResolution))

	plt.gcf().subplots_adjust(left=0.15)
	plt.xlabel("Thickness (m)")
	plt.ylabel("Time (s)")


	plt.show()

	return

def onclick(event, grid, tMax, xResolution, tResolution):
	x = int(event.xdata)
	y = int(event.ydata)

	if event.button == 1 :
		print "vitesse  ", format_e(grid[(x,y)][0]), " m/s"
		print "pression ", format_e(grid[(x,y)][1]), " Pa"
	   
	if event.button == 3 :
	    plotVelocity(tMax, tResolution, xResolution, grid, x)

	if event.button == 2 :
	    plotPressure(tMax, tResolution, xResolution, grid, x)

def plotVelocity(tMax, tRes, xRes, grid, position):
	tCoord = []
	vel = []
	sortedKeys = sorted(grid.keys())
	for key in sortedKeys:
		if key[0] == position:
			tCoord.append(key[1]*tMax/tRes)
			vel.append(grid[key][0])

    # plt.figure()
	fig, ax = plt.subplots()
	# ax.set_xticklabels(tticks)
	plt.plot(tCoord, vel)
	# plt.locator_params(axis='x', nbins=5)
	plt.xlabel("Time (s)")
	plt.ylabel("Velocity (m/s)")
	plt.show()

	return

def plotPressure(tMax, tRes, xRes, grid, position):
	tCoord = []
	pres = []
	sortedKeys = sorted(grid.keys())
	for key in sortedKeys:
		if key[0] == position:
			tCoord.append(key[1]*tMax/tRes)
			pres.append(grid[key][1])

    # plt.figure()
	fig, ax = plt.subplots()
	# ax.set_xticklabels(tticks)
	plt.plot(tCoord, pres)
	plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
	# plt.locator_params(axis='x', nbins=5)
	plt.xlabel("Time (s)")
	plt.ylabel("Pressure (GPa)")
	plt.show()

	return

