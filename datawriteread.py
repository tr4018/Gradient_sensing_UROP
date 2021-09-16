import random
import numpy as np
#help from: https://stackoverflow.com/questions/7618858/how-to-to-read-a-matrix-from-a-given-filewith 

def params_string(pick_dir,
receptornum = 10,
recepsurface_ratio = 100,
particlenum = 20,
sourcenum = 10,
sourceexact = -1,
random_yn  = 0,
diffusionnum = 5,
diffusionexact = -1,
distancenum = 5,
maxdistance = 10,
distanceexact = -1,
radiusnum = 5,
maxradius = 1,
radiusexact = -1,
ratenum = 5,
maxrate = 1,
rateexact = -1, 
receptor_seed = 1,
initial_source_seed = 1,
particle_seed= 1):#deal with source exact later
    return(['#'+pick_dir,'#datacreate() called:','#receptornum= '+str(receptornum),'#recepsurface_ratio= '+str(recepsurface_ratio), \
    '#particlenum= '+str(particlenum),'#sourcenum= '+str(sourcenum), '#sourceexact=','#random_yn= '+str(random_yn),'#diffusionnum= '+str(diffusionnum), \
    '#diffusionexact= '+ str(diffusionexact),'#distancenum= '+ str(distancenum),'#maxdistance= '+ str(maxdistance), '#distanceexact= '+ str(distanceexact),\
    '#radiusnum= '+ str(radiusnum), '#maxradius= '+ str(maxradius), '#radiusexact= '+ str(radiusexact), '#ratenum= '+ str(ratenum), '#maxrate= '+ str(maxrate),\
    '#rateexact= '+ str(rateexact), '#receptor_seed= '+ str(receptor_seed), '#initial_source_seed='+ str(initial_source_seed), '#particle_seed= '+ str(particle_seed)])


def write_datafile(filename, params, data):
    with open(filename +'.txt', 'w') as filehandle:
        filehandle.writelines("%s\n" % p for p in params)
        data = np.matrix(data)
        for line in data:
            np.savetxt(filehandle, line, fmt='%.2f')
    return

def read_datafile(filename):
    with open(filename+'.txt','r') as fin:
        a = []
        for line in fin:
            if line.startswith('#'): pass
            else:
                b = []
                line = line.split()
                for x in line: 
                    b.append(float(x))
                a.append(b)
    return(a)

#test
"""
params_1 = ['#pick_direction(0, 10)', '#datacreate() called:' , '#receptornum=10', '#particlenum=20']
X = np.random.rand(3,3)
write_datafile('X',params_1,X)
X = read_datafile('X.txt')
print(X)
"""
