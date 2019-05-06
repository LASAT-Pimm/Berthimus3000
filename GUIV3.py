# Read DB file generated with pyHugo (or Berthimus 3000)
from collections import defaultdict
import shelve
import tkFileDialog
import tkMessageBox as messagebox
import ttk
import os
import sys
import inspect
import sqlite3
if sys.version_info[0] < 3:
	from Tkinter import *
	import Tkinter as Tk
else:
	from tkinter import *
	import tkinter as Tk

class HoverButton(Tk.Button):
    def __init__(self, master, **kw):
        Tk.Button.__init__(self,master=master,**kw)
        self.defaultBackground = self["background"]
        self.bind("<Enter>", self.on_enter)
        self.bind("<Leave>", self.on_leave)

    def on_enter(self, e):
        self['background'] = self['activebackground']

    def on_leave(self, e):
        self['background'] = self.defaultBackground


class GUI(object):
	"""docstring for GUI"""
	def __init__(self):
		super(GUI, self).__init__()

	def simulation(self):
		self.value = 0
		self.window.destroy()
		return

	def read(self):
		self.value = 1
		self.window.destroy()
		return

	def launch(self):
		self.value = -1
		# window property
		self.window = Tk.Tk()
		self.window.title("Berthimus 3000")
		self.window.geometry('367x50')
		simBtn = HoverButton(self.window, width = 25, height = 3, text="Launch simulation", activebackground='grey', command = self.simulation)
		simBtn.grid(column=0, row=0)
		simBtn = HoverButton(self.window, width = 25, height = 3, text="Open simulation", activebackground='grey', command = self.read)
		simBtn.grid(column=1, row=0)

		self.window.mainloop()

		return self.value

class ReadHUGO(GUI):
	"""docstring for ReadHUGO"""
	def __init__(self):
		super(GUI, self).__init__()
		
	def launch(self):
		pos = 0

		# window property
		self.window = Tk.Tk()
		self.window.title("Berthimus 3000")
		self.window.geometry('373x355')

		pos = 0
		chocLbl = Tk.Label(self.window,text="")# add blank 
		chocLbl.grid(row=pos)
		pos += 1
		chocLbl = Tk.Label(self.window,text="")# add blank 
		chocLbl.grid(row=pos)
		pos += 1
		chocLbl = Tk.Label(self.window,text="")# add blank 
		chocLbl.grid(row=pos)
		pos += 1

		# select DB
		fileLbl = Tk.Label(self.window, text="Select .HUGO file")
		fileLbl.grid(column=0, row=pos)
		self.inputFile = Tk.Entry(self.window)
		self.inputFile.grid(column=1, row=pos)
		fileBtn = Tk.Button(self.window, text="Select", width = 20, command = self.fileLoc)
		fileBtn.grid(column=2, row=pos)
		pos += 1

		chocLbl = Tk.Label(self.window,text="")# add blank 
		chocLbl.grid(row=pos)
		pos += 1


		# sim resoltion
		chocLbl = Tk.Label(self.window,text="") # add blank 
		chocLbl.grid(row=pos)
		pos += 1
		resXOpenLbl = Tk.Label(self.window, text="Res X")
		resXOpenLbl.grid(column=0, row=pos)
		self.inputResXOpen = Tk.Entry(self.window)
		self.inputResXOpen.insert(0,"800")
		self.inputResXOpen.grid(column=1, row=pos)
		pos += 1
		resYOpenLbl = Tk.Label(self.window, text="Res Y")
		resYOpenLbl.grid(column=0, row=pos)
		self.inputResYOpen = Tk.Entry(self.window)
		self.inputResYOpen.insert(0,"800")
		self.inputResYOpen.grid(column=1, row=pos)
		pos += 1
		# plot
		chocLbl = Tk.Label(self.window,text="")# add blank 
		chocLbl.grid(row=pos)
		pos += 1
		plotBtn = Tk.Button(self.window, text="Plot", width = 30,  command = self.plot)
		plotBtn.grid(column=1, row=pos, columnspan=2)

		self.window.mainloop()

		return self.globGrid2,self.xRes, self.tRes, self.timeLimit, self.xmax, self.vmin, self.vmax, self.interface

	def plot(self):

		d = shelve.open(self.inputFile.get())
		self.globGrid2 = d['globGrid2'] 
		self.vmin = d['vmin']
		self.vmax = d['vmax']
		self.timeLimit = d['time']
		self.xmax = d['position']
		self.interface = d['interface']
		self.xRes = int(self.inputResXOpen.get())
		self.tRes = int(self.inputResYOpen.get())
		d.close()

		self.window.destroy()

		return

	def fileLoc(self):
		#current work directory
		filename = inspect.getframeinfo(inspect.currentframe()).filename
		workdir = os.path.dirname(os.path.abspath(filename))
		file = tkFileDialog.askopenfilename(defaultextension='.HUGO', initialdir = os.path.normpath(workdir))
		self.inputFile.insert(10,file)

		return		

class Simulation(GUI):
	"""docstring for FileOperation"""
	def __init__(self):
		super(GUI, self).__init__()

	def file_selection(self):
		root = Tk.Tk()
		root.title("DB reader")
		root.withdraw()

		# Current file dir
		filename = inspect.getframeinfo(inspect.currentframe()).filename
		workdir = os.path.dirname(os.path.abspath(filename))
		# print "filename ", filename
		# print "workdir ", workdir
		self.db = tkFileDialog.askopenfilename(defaultextension='.png', initialdir = os.path.normpath(workdir))
		return self.db
		# root.mainloop()

	def search_inp(self):
		self.matches = []
		for root, dirnames, filenames in os.walk(self.root_folder):
			for filename in fnmatch.filter(filenames, '*.f90'):
				self.matches.append(os.path.normpath(os.path.join(root, filename)))

		return os.path.normpath(self.db), self.matches

# acousticApprox = GUI()

	def addMat(self):
		# add material to the table
		# extracting values
		name = self.inputName.get()
		E = self.inputE.get()
		poisson = self.inputPoisson.get()
		density = self.inputDensity.get()
		velocity = self.inputVelocity.get()

		for i in range(self.nbRow+1):
			self.mat.execute('SELECT * FROM material WHERE nb = ?', (i,))
			if len(self.mat.fetchall()) == 0:
				formatName = (name,)
				# check if entry already exists
				self.mat.execute('SELECT * FROM material WHERE name = ?', formatName)

				if len(self.mat.fetchall()) != 0:
					self.error = 1

					return

				else:
					self.mat.execute('''INSERT INTO material VALUES (?, ?, ?, ?, ?, ?)''', (i, name, float(E), float(poisson), float(density), float(velocity)))
					self.materialDB.commit()
					self.nbRow += 1

				break
		
		self.win.destroy()
		self.printTable()
		return

	def delMat(self):
		delID = int(self.delInput.get())
		self.mat.execute('SELECT * FROM material WHERE nb=?', (int(delID),))
		self.mat.execute('DELETE FROM material WHERE nb = ?', (int(delID),))
		self.materialDB.commit()
		self.win.destroy()
		self.printTable()

	def printTable(self):
		width = 400
		div = 42
		self.win = Tk.Tk()
		self.error = 0

		# setting up scroll bar on the right
		sb = Tk.Scrollbar(self.win, orient=VERTICAL, bd=0, highlightthickness=0)
		sb.grid(column=6, sticky=N+S)
		c = Canvas(self.win,yscrollcommand=sb.set)
		c.grid(row=0, column=0, sticky="news")
		self.win.title("Window")
		self.win.geometry("450x"+str(4*38))
		sb.config(command=c.yview)
		self.win.grid_rowconfigure(0, weight=1)
		self.win.grid_columnconfigure(0, weight=1)
		fr = Frame(c)


	# 	Setting header
		self.mat.execute('PRAGMA table_info(material)')
		header = self.mat.fetchall()

	# 	listing material
		line = 0
		matID = Tk.Label(fr, text=header[0][1], width = width/div, borderwidth = 5)
		matID.grid(row=line, column=0)
		matName = Tk.Label(fr, text=header[1][1], width = width/div)
		matName.grid(row=line, column=1)
		matE = Tk.Label(fr, text=header[2][1], width = width/div)
		matE.grid(row=line, column=2)
		unitE = Tk.Label(fr, text="(GPa)", width = width/div)
		unitE.grid(row=line+1, column=2)
		matPoisson = Tk.Label(fr, text=header[3][1], width = width/div)
		matPoisson.grid(row=line, column=3)
		matDensity = Tk.Label(fr, text=header[4][1], width = width/div)
		matDensity.grid(row=line, column=4)
		unitDensity = Tk.Label(fr, text="(kg/m^3)", width = width/div)
		unitDensity.grid(row=line+1, column=4)
		matVelocity = Tk.Label(fr, text=header[5][1], width = width/div)
		matVelocity.grid(row=line, column=5)
		unitVelocity = Tk.Label(fr, text="(m/s)", width = width/div)
		unitVelocity.grid(row=line+1, column=5)
		line += 2

		for row in self.mat.execute('SELECT * FROM material ORDER BY nb'):
			matID = Tk.Label(fr, text=row[0], width = width/div, borderwidth = 5)
			matID.grid(row=line, column=0)
			matName = Tk.Label(fr, text=row[1], width = width/div)
			matName.grid(row=line, column=1)
			matE = Tk.Label(fr, text=row[2], width = width/div)
			matE.grid(row=line, column=2)
			matPoisson = Tk.Label(fr, text=row[3], width = width/div)
			matPoisson.grid(row=line, column=3)
			matDensity = Tk.Label(fr, text=row[4], width = width/div)
			matDensity.grid(row=line, column=4)
			matVelocity = Tk.Label(fr, text=row[5], width = width/div)
			matVelocity.grid(row=line, column=5)
			line += 1

		# add material
		self.inputName = Tk.Entry(fr, width = width/div)
		self.inputName.grid(column=1, row=line)
		self.inputE = Tk.Entry(fr, width = width/div)
		self.inputE.grid(column=2, row=line)
		self.inputPoisson = Tk.Entry(fr, width = width/div)
		self.inputPoisson.grid(column=3, row=line)
		self.inputDensity = Tk.Entry(fr, width = width/div)
		self.inputDensity.grid(column=4, row=line)
		self.inputVelocity = Tk.Entry(fr, width = width/div)
		self.inputVelocity.grid(column=5, row=line)
		line += 1

		# add button
		add = Tk.Button(fr, text="Add mat", width = width/div, command = self.addMat)
		add.grid(row=line, column=4)

		# close button
		b = Tk.Button(fr, text="Close", command=self.win.destroy, width = width/div)
		b.grid(row=line, column=5)

		# delete section
		delLbl = Tk.Label(fr, text="ID", width = width/div)
		delLbl.grid(row=line, column=0)
		self.delInput = Tk.Entry(fr, width = width/div)
		self.delInput.grid(row=line, column=1)
		delBtn = Tk.Button(fr, text = "Delete", width = width/div, command = self.delMat)
		delBtn.grid(row=line, column=2)

		c.create_window(0, 0,  window=fr)
		fr.update_idletasks()
		c.config(scrollregion=c.bbox("all"))

# 		Resetting view
		c.yview_moveto(0)

		return

	def folderLoc(self):
		#current work directory
		filename = inspect.getframeinfo(inspect.currentframe()).filename
		workdir = os.path.dirname(os.path.abspath(filename))
		root_folder = tkFileDialog.askdirectory(initialdir = os.path.normpath(workdir))
		self.inputFolder.insert(10,root_folder)

		return

	def submit(self):
		self.FOLDER = self.inputFolder.get()
		if self.FOLDER == "" :
			filename = inspect.getframeinfo(inspect.currentframe()).filename
			workdir = os.path.dirname(os.path.abspath(filename))
			self.inputFolder.insert(10,workdir)
			messagebox.showwarning('Folder location', 'folder set to current folder')
			return

		self.OUTPUT = self.inputOutput.get()
		if self.OUTPUT == "":
			messagebox.showerror('OUTPUT missing', 'please add output name')
			return

		try:
			self.PULSE = float(self.inputPulse.get()) * 1.E-9
		except:
			messagebox.showerror('Pulse duration', 'Check pulse duration format')
			return
		try:	
			self.P1 = float(self.inputP1.get()) * 1.E9
		except:
			messagebox.showerror('P1', 'Check P1 format')
			return

		try:
			self.P2 = float(self.inputP2.get()) * 1.E9
		except:
			if self.typeSim.get() != "mono":
				messagebox.showerror('P2', 'Check P2 format')
				return
			else:
				self.P2 = 0.0

		try:
			self.DELAY = float(self.inputDelay.get()) * 1.E-9
		except:
			messagebox.showerror('Delay', 'Check delay format')
			return
		
		try:
			self.DURATION = float(self.inputTime.get()) * 1.E-6
		except:
			messagebox.showerror('Sim duration', 'Check sim duration format')
			return
		
		try:
			self.XRES = int(self.inputResX.get())
		except:
			messagebox.showerror('Xres', 'Requires Res X integer')
			return

		try:
			self.YRES = int(self.inputResY.get())
		except:
			messagebox.showerror('Yres', 'Requires Res Y integer')
			return

		try:
			self.MATSEQ = list(map(int, self.inputMat.get().split(',')))
			self.MATERIAL = defaultdict()
			matList = list(set(sorted(self.MATSEQ)))
			for i in range(len(matList)):
				self.MATERIAL[matList[i]] = self.getMatParam(matList[i])
				if len(self.MATERIAL[matList[i]]) == 0 :
					messagebox.showerror('Mat sequence', 'Error, material ID ' + str(self.MATSEQ[i]) + 'does not exist, please refere to material list')
					return
		except:
			messagebox.showerror('Mat sequence', 'Check mat sequence format, should be: mat#1,mat#2,mat#3....')
			return

		try:
			self.MATTCK = list(map(float, self.inputTck.get().split(',')))
			for i in range(len(self.MATTCK)):
				self.MATTCK[i] = self.MATTCK[i] * 1.E-6
		except:
			messagebox.showerror('Mat thickness', 'Check format, should be: thicknessMat#1,thicknessMat#2,thicknessMat#3....')
			return

		self.TYPESIM = self.typeSim.get()

		self.window.destroy()
		
		return

	def getMatParam(self, matID):

		self.mat.execute('''SELECT * FROM material where nb = ? ''', (matID,))

		return self.mat.fetchall()


	def launch(self):

		self.error=0
		pos = 0

		# window property
		self.window = Tk.Tk()
		self.window.title("Berthimus 3000")
		self.window.geometry('373x325')

		# folder definition
		folderLbl = Tk.Label(self.window, text="Define output folder")
		folderLbl.grid(column=0, row=pos)
		self.inputFolder = Tk.Entry()
		self.inputFolder.grid(column=1, row=pos)
		folderBtn = Tk.Button(self.window, text="Select", command = self.folderLoc)
		folderBtn.grid(column=2, row=pos)
		pos += 1

		# output name
		outputLbl = Tk.Label(self.window, text="Name output")
		outputLbl.grid(column=0, row=pos)
		self.inputOutput = Tk.Entry()
		self.inputOutput.grid(column=1, row=pos)
		pos += 1

		# choc parameters
		chocLbl = Tk.Label(self.window,text="")# add blank 
		chocLbl.grid(row=pos)
		pos += 1
		pulseLbl = Tk.Label(self.window, text="Pulse duration (ns)")
		pulseLbl.grid(column=0, row=pos)
		self.inputPulse = Tk.Entry()
		self.inputPulse.grid(column=0, row=pos+1)

		p1Lbl = Tk.Label(self.window, text="P1 (GPa)")
		p1Lbl.grid(column=1, row=pos)
		self.inputP1 = Tk.Entry()
		self.inputP1.grid(column=1, row=pos+1)

		p2Lbl = Tk.Label(self.window, text="P2 (GPa)")
		p2Lbl.grid(column=2, row=pos)
		self.inputP2 = Tk.Entry()
		self.inputP2.grid(column=2, row=pos+1)
		pos = pos + 2

		# choc parameters
		delayLbl = Tk.Label(self.window, text="Delay (ns)")
		delayLbl.grid(column=0, row=pos)
		self.inputDelay = Tk.Entry()
		self.inputDelay.grid(column=0, row=pos+1)

		timeLbl = Tk.Label(self.window, text="Sim duration (us)")
		timeLbl.grid(column=1, row=pos)
		self.inputTime = Tk.Entry()
		self.inputTime.grid(column=1, row=pos+1)
		
		pos = pos + 2

		chocLbl = Tk.Label(self.window,text="")# add blank 
		chocLbl.grid(row=pos)
		pos += 1

		# material
		self.materialDB = sqlite3.connect('material.db') #connect to material.db
		self.mat = self.materialDB.cursor()
		try:
			# create material DB if does not exists
			self.mat.execute('''CREATE TABLE material(nb, name, young, poisson, density, velocity)''')
			# mat.execute('''INSERT INTO material VALUEs (?, ?, ?, ?, ?, ?)''', ("nb", "name", "E (GPa)", "poisson", "density (kg/cm^3)", "velocity (m/s)"))
		except:
			pass
		# find number of entries in the database
		self.mat.execute('SELECT * FROM material')
		results = self.mat.fetchall()
		self.nbRow = len(results)

		matLbl = Tk.Label(self.window, text="Mat sequence")
		matLbl.grid(column=0, row=pos)
		self.inputMat = Tk.Entry()
		self.inputMat.insert(0,"ex. 1,2,1")
		self.inputMat.grid(column=1, row=pos)
		matBtn = Tk.Button(self.window, text="See material list", command = self.printTable)
		matBtn.grid(column=2, row=pos)
		pos += 1

		# thicknes
		tckLbl = Tk.Label(self.window, text="Thickness in um")
		tckLbl.grid(column=0, row=pos)
		self.inputTck = Tk.Entry()
		self.inputTck.insert(0, "ex. 300,100,300")
		self.inputTck.grid(column=1, row=pos)

		pos += 1
		# radio button type sym selection
		self.typeSim = StringVar()
		self.typeSim.set("mono")
		monoBtn = Tk.Radiobutton(self.window, variable=self.typeSim, text="mono", value="mono")
		monoBtn.grid(column=0, row=pos)
		symBtn = Tk.Radiobutton(self.window, variable=self.typeSim, text="sym", value="sym")
		symBtn.grid(column=1, row=pos)
		doubleBtn = Tk.Radiobutton(self.window, variable=self.typeSim, text="double", value="double")
		doubleBtn.grid(column=2, row=pos)
		pos += 1

		# sim resoltion
		chocLbl = Tk.Label(self.window,text="") # add blank 
		chocLbl.grid(row=pos)
		pos += 1
		resXLbl = Tk.Label(self.window, text="Res X")
		resXLbl.grid(column=0, row=pos)
		self.inputResX = Tk.Entry()
		self.inputResX.insert(0,"800")
		self.inputResX.grid(column=0, row=pos+1)

		resYLbl = Tk.Label(self.window, text="Res Y")
		resYLbl.grid(column=1, row=pos)
		self.inputResY = Tk.Entry()
		self.inputResY.insert(0,"800")
		self.inputResY.grid(column=1, row=pos+1)

		pos = pos + 2

		# submit
		subBtn = Tk.Button(self.window, text="Submit", command = self.submit)
		subBtn.grid(column=2, row=pos)


		self.materialDB.commit()
		self.window.mainloop()

		return self.FOLDER, self.OUTPUT, self.PULSE, self.P1, self.P2, self.DELAY, self.DURATION, self.XRES, self.YRES, self.MATSEQ, self.MATERIAL, self.MATTCK, self.TYPESIM
		



# window = Window()
# window.mainWindow()
# gui = GUI()
# gui.launch()
# FOLDER, OUTPUT, PULSE, P1, P2, DELAY, DURATION, XRES, YRES, MATSEQ, MATERIAL, MATTCK, TYPESIM = window.mainWindow()
# print FOLDER
# print OUTPUT
# print PULSE
# print P1
# print P2
# print DELAY
# print XRES
# print YRES
# print MATSEQ
# print MATERIAL
# print MATTCK
# print TYPESIM
