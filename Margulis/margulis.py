import subprocess
import cmath
import numpy as np

import time
start_time = time.time()



#interact with snap, return a list of geodesics of lenth less than cutoff
def get_geodesics(manifold_number, cutoff):
	snap_interaction = "r closed {}\n print geodesics {}\n".format(str(manifold_number), str(cutoff))
	p = subprocess.Popen("snap", stdin=subprocess.PIPE, stdout=subprocess.PIPE, shell=True, universal_newlines=True)  
	output, err = p.communicate(input=snap_interaction)
	output = output.split('\n')
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

#returns just the shortest (real length) ortholine between geodesic1 and geodesic2
#has to distinguish between ortholines going to different geodesics from ortholines going back to original geodesic
def get_shortest_ortholine(manifold_number, cutoff, geodesic1, geodesic2):
	ortholines = get_ortholines(manifold_number, cutoff, geodesic1, geodesic2)
	if len(ortholines) == 0: return get_shortest_ortholine(manifold_number, cutoff+1, geodesic1, geodesic2)
	else:
		if geodesic1 == geodesic2:
			return crop_ortholine_text(ortholines[0])
		else:
			for ortholine in ortholines:
				if ortholine_index_different(ortholine):
					return crop_ortholine_text(ortholine)
			return get_shortest_ortholine(manifold_number, cutoff+1, geodesic1, geodesic2)

#outputs the complex lengt of an ortholine
def crop_ortholine_text(ortholine_string):
	return complex(ortholine_string[1:ortholine_string.find("*")]+"j")

#checks if an ortholine has endpoints on different geodesics
def ortholine_index_different(ortholine_string):
	return int(ortholine_string[ortholine_string.find(":")-1]) != int(ortholine_string[ortholine_string.rfind(":")-1])







def isMargulis(manifoldNumber, number):
	geodesics = get_geodesics(manifoldNumber,number)
	print len(geodesics)
	if len(geodesics)==0: return True
	else: 
		tube_radii =[]
		for geo in geodesics:
			tube_radii.append(naiveTubeRadius(geo,number))
		for geo0 in range(len(geodesics)):
			for geo1 in range(geo0,len(geodesics)):
				print tube_radii[geo0], tube_radii[geo1], get_shortest_ortholine(manifoldNumber,number,geo0,geo1).real
				if tube_radii[geo0]+ tube_radii[geo1] >= get_shortest_ortholine(manifoldNumber,number,geo0,geo1).real:
					return False
		return True

def naiveTubeRadius(geodesic_lenght, number):
	r = geodesic_lenght.real
	im = geodesic_lenght.imag
	return np.arccosh(np.sqrt((np.cosh(number)-np.cos(im))/(np.cosh(r)-np.cos(im))))

for i in range(1,10):
	print isMargulis(i,.9)





print "------------------", time.time() - start_time, "seconds ---------------"


