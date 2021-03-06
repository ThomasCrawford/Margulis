import subprocess
import csv
import cmath
import numpy as np
from scipy.optimize import fsolve
#import verify


import time
start_time = time.time()


muGuess = 1.2
#if muGuess is too low, try again with muGuess + 
increment = .1
#with too high an increment, snap is more likely to crash
getOrtholineIncrement = .5

#interact with snap, return a list of geodesics of lenth less than cutoff
def get_geodesics(manifold_number, cutoff):
	snap_interaction = "r closed {}\n print geodesics {}\n".format(str(manifold_number), str(cutoff))
	p = subprocess.Popen("snap", stdin=subprocess.PIPE, stdout=subprocess.PIPE, shell=True, universal_newlines=True)  
	output, err = p.communicate(input=snap_interaction)
	output = output.split('\n')
	if output[2] == "Problem computing a Dirichlet domain for this group.": return False
	else:
		geodesics = []
		for item in output[2:-1]:
			geodesics.append(complex(item[item.find("]")+1:item.rfind("*")]+"j"))
		return geodesics

#interact with snap, return a list of ortholines between geodesic1 and geodesic2, where both lenght of the
#ortholine and length of the geodesics are less than cutoff
#note: has to call geodesics first to load list of geodesics
def get_ortholines(manifold_number, cutoff, geodesic1, geodesic2):
	assert len(get_geodesics(manifold_number,cutoff)) >= max(geodesic1,geodesic2)
	snap_interaction = "r closed {}\n print geodesics {}\n print ortholines {}\n".format(str(manifold_number),str(cutoff),str(cutoff) +" "+ str(geodesic1) +" "+ str(geodesic2))
	p = subprocess.Popen("snap", stdin=subprocess.PIPE, stdout=subprocess.PIPE, shell=True, universal_newlines=True)  
	output, err = p.communicate(input=snap_interaction)
	output_list = output.split('\n')[2:-1]
	ortholine_count = int(output_list[-1][:output_list[-1].find(" ")])
	just_ortholines = output_list[-ortholine_count-1:-1]
	return just_ortholines

#interact with snap, returns name of manifold and volume of manifold
def get_name_and_volume(manifold_number):
	snap_interaction = "r closed {}\n print name \n print volume\n".format(str(manifold_number))
	p = subprocess.Popen("snap", stdin=subprocess.PIPE, stdout=subprocess.PIPE, shell=True, universal_newlines=True)  
	output, err = p.communicate(input=snap_interaction)
	output_list = output.split('\n')[2:-1]
	return [output_list[0][7:], output_list[1][11:]]


#returns just the shortest (real length) ortholine between geodesic1 and geodesic2
#has to distinguish between ortholines going to different geodesics from ortholines going back to original geodesic
def get_shortest_ortholine(manifold_number, cutoff, geodesic1, geodesic2):
	ortholines = get_ortholines(manifold_number, cutoff, geodesic1, geodesic2)
	if len(ortholines) == 0: 
		return get_shortest_ortholine(manifold_number, cutoff+getOrtholineIncrement, geodesic1, geodesic2)
	else:
		if geodesic1 == geodesic2:
			return crop_ortholine_text(ortholines[0])
		else:
			for ortholine in ortholines:
				if ortholine_index_different(ortholine):
					return crop_ortholine_text(ortholine)
			return get_shortest_ortholine(manifold_number, cutoff+getOrtholineIncrement, geodesic1, geodesic2)

#outputs the complex length of an ortholine
def crop_ortholine_text(ortholine_string):
	return complex(ortholine_string[1:ortholine_string.find("*")]+"j")

#checks if an ortholine has endpoints on different geodesics, 
#useful because specifying multiple geodesics in snap's print ortholines command
def ortholine_index_different(ortholine_string):
	lines = ortholine_string.split("i")
	return int(lines[1][:lines[1].find(":")]) !=int(lines[2][:lines[2].find(":")])

#Verifies that a number is margulis
def isMargulis(manifoldNumber, number):
	geodesics = get_geodesics(manifoldNumber,number)
	if len(geodesics)==0: return True
	else: 
		tube_radii =[]
		for geo in geodesics:
			tube_radii.append(tubeRadius(geo,number))
		for geo0 in range(len(geodesics)):
			for geo1 in range(geo0,len(geodesics)):
				if tube_radii[geo0]+ tube_radii[geo1] >= get_shortest_ortholine(manifoldNumber,number,geo0,geo1).real:
					return geo0, tube_radii[geo0], geo1, tube_radii[geo1], get_shortest_ortholine(manifoldNumber,number,geo0,geo1).real
		return True

#calculates the radius of a tube around a geodesic, keeping in mind that the tube sometimes might be larger when considering integer multiples of the geodesic 
def tubeRadius(geodesic_length, number):
	r = geodesic_length.real
	im = geodesic_length.imag
	if number - r < 0.00000001: return 0
	else:
		possibleRadius = []
		for n in range(1,int(np.ceil(number/r))):
			possibleRadius.append(np.arccosh(np.sqrt((np.cosh(number)-np.cos(n*im))/(np.cosh(n*r)-np.cos(n*im)))))
		return max(possibleRadius)

def solveForMu(geoLength0, geoLength1, ortholength):
	if geoLength0.real > geoLength1.real:
		r0 = geoLength1.real
		im0 = geoLength1.imag
		r1 = geoLength0.real
		im1 = geoLength0.imag
	else:
		r0 = geoLength0.real
		im0 = geoLength0.imag
		r1 = geoLength1.real
		im1 = geoLength1.imag
	#potentially, a tube around geodesic0, dictated by mu which is less than r1 can intersect geodesic1.  
	#in this case, the cutoff mu is r1  (this messes with the later solve functions)
	#we can assume r1>r0
	if tubeRadius(geoLength0,r1) >= ortholength: 
		return [r1]
	func = lambda mu : -1000 * ortholength if mu<r0 else \
		np.arccosh(np.sqrt((np.cosh(mu)-np.cos(im0))/(np.cosh(r0)-np.cos(im0))))-ortholength if mu<=r1 else \
		np.arccosh(np.sqrt((np.cosh(mu)-np.cos(im0))/(np.cosh(r0)-np.cos(im0))))+\
		np.arccosh(np.sqrt((np.cosh(mu)-np.cos(im1))/(np.cosh(r1)-np.cos(im1)))) - ortholength
	initialAnswer = fsolve(func,max(r0,r1)+.01, factor = .001)
	#now we take into account the possibliity of a multiple of one of the geodsics doing better
	potentialAnswers = []
	for n0 in range(1,int(np.ceil((initialAnswer+.01)/r0))):
		for n1 in range(1,int(np.ceil((initialAnswer+.01)/r1))):
			func = lambda mu : -1000*ortholength if mu<r0 else \
				np.arccosh(np.sqrt((np.cosh(mu)-np.cos(n0*im0))/(np.cosh(n0*r0)-np.cos(n0*im0))))-ortholength if mu<=r1 else \
				np.arccosh(np.sqrt((np.cosh(mu)-np.cos(n0*im0))/(np.cosh(n0*r0)-np.cos(n0*im0))))+\
				np.arccosh(np.sqrt((np.cosh(mu)-np.cos(n1*im1))/(np.cosh(n1*r1)-np.cos(n1*im1)))) - ortholength
			potentialAnswers.append(fsolve(func,max(n0*r0,n1*r1)+.01,factor = .001))
	return min(potentialAnswers)

#Primary function.  Organizes the search for minimum mu
#if margulisGuess is an overestimation, some tubes intersect, one of these intersections is the cutoff  
#Return the number that makes them just barely intersect
def findCutoff(manifoldNumber, margulisGuess):
	geodesics = get_geodesics(manifoldNumber,margulisGuess)
	if not geodesics: return [0,0,0,"No Dirichlet Domain"]
	# print manifoldNumber, len(geodesics)
	cutoffCandidates = []
	assert len(geodesics)!=0
	tube_radii =[]
	for geo in geodesics:
		tube_radii.append(tubeRadius(geo,margulisGuess))
	for geo0 in range(len(geodesics)):
		for geo1 in range(geo0,len(geodesics)):
			ortholine = get_shortest_ortholine(manifoldNumber,margulisGuess,geo0,geo1).real
			if tube_radii[geo0]+ tube_radii[geo1] >= ortholine:
				cutoffCandidates.append([geodesics[geo0],geodesics[geo1],ortholine])
	if len(cutoffCandidates)==0: 
		return findCutoff(manifoldNumber,margulisGuess+increment)
	else:
		candidatesOptimized = []
		for geo0, geo1, ortho in cutoffCandidates:
			candidatesOptimized.append([geo0, geo1, ortho, solveForMu(geo0,geo1,ortho)[0]])
#if margulisGuess is an underestimation, we potentially didn't consider all relevent pairs of geodesics, run again with output as a new guess 
		if min(candidatesOptimized, key=lambda x: x[3])[3] > margulisGuess:
			return findCutoff(manifoldNumber, min(margulisGuess + increment, min(candidatesOptimized, key=lambda x: x[3])[3]))
		return min(candidatesOptimized, key=lambda x: x[3])

#reorders and reformats for .csv output
def organize(manifoldNumber, margulisGuess):
	name, volume = get_name_and_volume(manifoldNumber)
	geoLength0, geoLength1, ortholength, margulis = findCutoff(manifoldNumber, margulisGuess)
	return[manifoldNumber, name, volume, margulis, geoLength0, geoLength1, ortholength]


def main():
	# print solveForMu(1.21928314545+2.75158841936j,1.21928314545+2.75158841936j,0.033831524241448)

	with open('margulis_incorrect.csv','r') as file:
		file_reader = csv.reader(file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
		for row in file_reader:
			try: 
				print row[0],findCutoff(row[0], muGuess)[3] - float(row[9])
			except:
				print row[0]
	print findCutoff(1176, muGuess)

	# for i in range(12440000, 12000):
	# 	# print "------------------", time.time() - start_time, "seconds ---------------"
	# 	csvLine = organize(i,muGuess)
	# 	# print csvLine[3]
	# 	with open('margulis.csv','a') as file:
	# 		file_writer = csv.writer(file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
	# 		file_writer.writerow(csvLine)
	# 		file.close()


if __name__ == "__main__":
	main()

