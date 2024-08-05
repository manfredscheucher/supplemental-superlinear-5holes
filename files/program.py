#!/usr/bin/python
"""
	This program uses functions from the 
	python order type library pyotlib.
	For sakes of shortness, relevant functions were 
	copied from the library into this program 
	and all unused code parts were removed.
	(c) 2013-2016 Manfred Scheucher <mscheuch@ist.tugraz.at>
"""



from itertools import permutations
from datetime import datetime
import sys
import struct



sign = lambda x: (x>0) - (x<0)



# class for point sets 
class PointSet(object):


	# initialize a point set
	def __init__(self,n,points):
		assert(len(points) == n)
		self.n = n
		self.points = points


	# compute the orientation of the triple i,j,k 
	def calcOrientation(self,i,j,k):
		ix,iy = self.points[i]
		jx,jy = self.points[j]
		kx,ky = self.points[k]
		return sign(ix*(jy-ky)+jx*(ky-iy)+kx*(iy-jy))


	# compute the big lambda matrix (i.e., the set of triple-orientations)
	def toBigLambda(self): 
		orientations = [[[0 for i in range(self.n)] for j in range(self.n)] for k in range(self.n)]
		for i in range(self.n):
			for j in range(i,self.n):
				for k in range(j,self.n):
					o = self.calcOrientation(i,j,k) if (i<j and j<k) else 0
					orientations[i][j][k] = orientations[j][k][i] = orientations[k][i][j] = o
					orientations[i][k][j] = orientations[j][i][k] = orientations[k][j][i] = -o
		return BigLambda(self.n, orientations)


	# select a subset of the given point set
	def selectPoints(self,perm):
		k = len(perm)
		return PointSet(k,[self.points[perm[i]] for i in range(k)])



# class for big lambda matrices (i.e., triple-orientations)
class BigLambda(object):


	# initialize
	def __init__(self,n,o):
		self.n = n
		self.o = o


	# compute extremal points from triple-orientations
	def getExtremalPoints(self):
		hull = []
		for i in range(self.n):
			for j in range(self.n):
				if j != i and [k for k in range(self.n) if self.o[i][j][k] > 0] == []:
					hull.append(i)
					break
		return hull


	# compute the rotation system from triple-orientations
	def getRotationSystem(self):	
		n = self.n
		o = self.o
		order = []
		for i in range(n):
			order.append([])
			j = (i+1)%n
			sign = +1
			while True:
				if sign == +1:
					if j in order[i]: 
						break
					else:
						order[i].append(j)

				jleft = {k for k in range(n) if o[i][j][k] > 0}
				if not jleft: sign *= -1 # current point is extremal
				
				jnext = None
				for j2 in range(n):
					if j2 == j or j2 == i: continue

					j2left  = {k for k in range(n) if o[i][j2][k] > 0}
					j2right = {k for k in range(n) if o[i][j2][k] < 0}

					if j2 in jleft and (j2left|{j2}) == jleft:
						jnext = j2
						break

					if j2 not in jleft and j2right == jleft:
						sign *= -1
						jnext = j2
						break

				assert(jnext != None) # there is always such a consecutive point 
				j = jnext

			assert(len(order[i]) == n-1)
		return order
		


# enumerate all k-holes from triple-orientations  
def enumerateKHoles(BL,k):
	# ensure natural labeling, that is, 
	# the first point p is extremal
	# and the others are sorted around p
	for i in range(1,BL.n-1):
		assert(BL.o[0][1][2] == BL.o[0][i][i+1])

	for pt in range(k-1,BL.n):
		for poly in _enumerateKHolesInner(BL,k-1,selection=[pt],PotentialPoints=range(pt)):
			yield poly



# auxiliary function for enumerating all k-holes from triple-orientations  
def _enumerateKHolesInner(BL,k,selection,PotentialPoints):
	if len(PotentialPoints) < k:
		return # not enough points for a k-hole

	if k == 0:
		yield selection
	else:
		b = selection[-1]
		s = selection[0]
		for a in PotentialPoints:
			newPotentialPoints = [c for c in PotentialPoints if c != a and BL.o[s][a][c] == 1 and BL.o[b][a][c] == 1]
			if b != s:
				# at least 2 points must be present, otherwise no triangle exists
				if not isEmptyTriangle(BL,s,b,a,PotentialPoints): continue

			for poly in _enumerateKHolesInner(BL,k-1,selection+[a],newPotentialPoints):
				yield poly



# test whether three points form an empty triangle
def isEmptyTriangle(BL,a,b,c,PotentialInner):
	for p in enumerateTriangleInnerPoints(BL,a,b,c,PotentialInner=PotentialInner):
		return False # inner point found
	return True



# enumerate points that lie inside the triangle a,b,c
def enumerateTriangleInnerPoints(BL,a,b,c,PotentialInner):
	if BL.o[a][b][c] == 1:
		for p in PotentialInner:
			if p!=a and p!=b and p!=c and BL.o[a][b][p] == 1 and BL.o[b][c][p] == 1 and BL.o[c][a][p] == 1:
				yield p
	else:
		assert(BL.o[a][c][b]==1)
		for p in enumerateTriangleInnerPoints(BL,a,c,b,PotentialInner): 
			yield p



# class for reading point sets from a file
class PointSetBinaryReader(object):


	# initialize the file-reader
	def __init__(self,n,bytes,filepath):
		self.n = n
		self.bytes = bytes
		self.filepath = filepath
		self.f = open(filepath)


	# auxiliary function for reading bytes
	def _unpack(self,string):
		# for documentation, see 
		# http://docs.python.org/2/library/struct.html#format-characters
		if self.bytes == 1: return struct.unpack("<B",string)[0]
		if self.bytes == 2: return struct.unpack("<H",string)[0]
		exit("Invalid number of bytes!")


	# enumerate all point sets read from the file 
	def readAll(self):
		while True:
			l = self.readNext()
			if l is None: 
				break 
			yield l


	# read the next point set from the file
	def readNext(self):
		coordinates = []
		for i in range(self.n):
			px = self.f.read(self.bytes)
			py = self.f.read(self.bytes)
			if px == "" or py == "": return None # end of the file reached
			point = (self._unpack(px),self._unpack(py))
			coordinates.append(point)			
		return PointSet(self.n, coordinates)



# class for arguments from the command line 
class ArgumentReader(object):	


	# initialize the argument-reader
	def __init__(self,arguments):
		self.d = dict()
		l = len(arguments)
		for i in range(1,l-1,2):
			self.d[arguments[i]] = arguments[i+1] 


	# get the value of a command line parameter
	def get(self,name,default=None):
		if name in self.d: return self.d[name]
		if default != None: return default		
		exit("No parameter \""+name+"\" found!")



# class for testing the statements
class HasDividedFiveHoleScript(object):


	# return progress text
	def progressText(self): 
		return " ** count: "+str(self.count)+" **"


	# return time stamp
	def timestamp(self): 
		return "["+str(datetime.now())+"]"


	# print usage information
	def printUsageInfo(self):
		print "inputfile specific"
		print "\tn - number of points"
		print "\tbytes - bytes"
		print "\tfp - filepath"
		print "set sizeA to: the number of points in A"
		print "set convA to:"
		print "\t 0 if A can be arbitrary" 
		print "\t+1 if A must be in convex position"
		print "\t-1 if A must not be in convex position"
		print "set testwedges:"
		print "\t0 if only l-divided 5-holes should be tested"
		print "\t1 if also a-wedges should be tested" 
		print "set nonconvexwedgeempty to:"		
		print "\t0 if the non-convex a-wedge can contain points" 
		print "\t1 if the non-convex a-wedge must be empty" 


	# initialize
	def __init__(self):
		if len(sys.argv) == 1 or len(sys.argv)%2 == 0:
			self.printUsageInfo()
			exit(-1)

		self.arguments = ArgumentReader(sys.argv)

		self.n = int(self.arguments.get("n"))
		self.bytes = int(self.arguments.get("bytes"))
		self.filepath = self.arguments.get("fp")

		if self.n > 6:
			self.sizeA = int(self.arguments.get("sizeA"))
			self.sizeB = self.n - self.sizeA
			self.convA = int(self.arguments.get("convA"))
			self.testwedges = int(self.arguments.get("testwedges"))
			if self.testwedges:
				self.nonconvexwedgeempty = int(self.arguments.get("nonconvexwedgeempty"))

		self.reader = PointSetBinaryReader(self.n, self.bytes, self.filepath)
		self.fail = 0
	

	# iterate over all point sets from the given file and test each set
	def action(self):
		self.count = 0
		self.show = 1
		print self.timestamp(),"loop started"
		
		for PS in self.reader.readAll():
			if self.count == self.show:
				if self.show < 1000:
					self.show = min(1000,2*self.show)
				else:
					self.show += 1000
				print self.timestamp()+self.progressText()

			self.actionInner(PS)
			self.count += 1

		print self.timestamp(),"loop done"
		self.actionEnd()


	# test a given point set
	def actionInner(self,PS):
		n = self.n
		BL = PS.toBigLambda()
		RS = BL.getRotationSystem()

		extremal = list(BL.getExtremalPoints())
		c5hs = list(enumerateKHoles(BL,5))
		
		# run the program for order types P with |P|=6 to verify 
		# h_5(P) in {0,1,2,6}
		if n == 6:
			assert(len(c5hs) in [0,1,2,6])
			return

		# try to partition P = A cup B by a line l
		for a,b in permutations(range(n),2):
			A = [i for i in range(n) if BL.o[a][b][i] > 0]+[a]
			B = [i for i in range(n) if BL.o[a][b][i] < 0]+[b]
			if len(A) != self.sizeA: continue

			PSA = PS.selectPoints(A)
			BLA = PSA.toBigLambda()
			
			extremals = BL.getExtremalPoints()
			extremals_A = [A[i] for i in BLA.getExtremalPoints()]
			extremals_A_P = [aa for aa in extremals if aa in A]
			
			# test whether A is (not) in convex position 
			if self.convA == +1 and len(extremals_A) != len(A): continue
			if self.convA == -1 and len(extremals_A) == len(A): continue

			# check whether there is an l-divided 5-hole in P
			found_ldivided_5hole = False
			for c5h in c5hs:
				found_c5h_incident_to_A = False
				found_c5h_incident_to_B = False
				for p in c5h:
					if p in A: found_c5h_incident_to_A = True
					if p in B: found_c5h_incident_to_B = True

				if found_c5h_incident_to_A and found_c5h_incident_to_B: 
					found_ldivided_5hole = True
					break
			if found_ldivided_5hole: continue

			if not self.testwedges:
				self.fail += 1
				print "A counterexample was found. ",self.fail,"of",self.count,"failed."
				return

			else:
				for alpha in A:
					perm = RS[alpha]
					cperm = perm+perm 

					convex_wedge_has_3_points = False
					nonconvex_wedge_is_empty = True

					prev_A = None
					inner = 0
					for pt in cperm:
						if pt in A:
							if prev_A != None: 
								# test whether the non-convex alpha-wedge (if it exist) 
								# is empty of point of P
								if inner >= 1 and BL.o[alpha][prev_A][pt] < 0:
									nonconvex_wedge_is_empty = False
								# test whether a convex alpha-wedge contains more than two points of P
								if inner >= 3 and BL.o[alpha][prev_A][pt] > 0:
									convex_wedge_has_3_points = True

							inner = 0
							prev_A = pt
						else:
							inner += 1

					if convex_wedge_has_3_points and (not self.nonconvexwedgeempty or nonconvex_wedge_is_empty):
						self.fail += 1
						print "A counterexample was found. ",self.fail,"of",self.count,"failed."
						return


	# print status after the computations are done
	def actionEnd(self):
		print "status:",self.fail,"of",self.count,"failed."
		if self.fail == 0:
			print "The statement is verified."
		else:
			print "Counterexamples were found."



# runs the program
HasDividedFiveHoleScript().action()
