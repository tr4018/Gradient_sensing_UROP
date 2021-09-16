from os import read
from numpy.core.function_base import linspace
from ReceptorNeuralNetwork import *
from datacreation import *
from IdealDirection import *
from datawriteread import *
from sklearn.model_selection import train_test_split
import numpy as np
import pickle


######## FIRST CHECK ##################
#direction picking for cell: 
#create some kind of print command for user to select the directions
#direction_sphcoords = pick_direction(0, 5)
#print(direction_sphcoords)

#why index 0 more activated than the rest?
#X, Y = datacreate(direction_sphcoords, sourcenum=5 ,receptornum=10,particlenum=10, recepsurface_ratio=10, distancenum=3,radiusnum=3,diffusionnum=3,ratenum=2,receptor_seed=1,initial_source_seed=1)
#params = params_string('pick_direction(0,10)', sourcenum=5 ,receptornum=10,particlenum=10, recepsurface_ratio=10, distancenum=3,radiusnum=3,diffusionnum=3,ratenum=2,receptor_seed=1,initial_source_seed=1)
#write_datafile('X_test=',params, X)
#write_datafile('Y_test=',params, Y)
#check double index in Y
#coords = cart2spherical_point(-4.651031627755137,2.903085818191935,-0.43588706087042833)
#print(direction_sphcoords)
#print(coords)
#ind = ideal_direction(coords[1],coords[2],direction_sphcoords,1)
#print(ind)
#check 

################### compare particle number accuracy  ###########################
####### THIS WORKS!! EXACT BUT DIFFERENT SEEDS###########
#particletest= [5,10,20,30]
#particletest= [40,50,60]


distancetest= [2,3,4,5,6,7,8,9,10,11,12]
direction_sphcoords = pick_direction(0,10) #same as sourcenum
seeds = 1000

Xfinal = np.zeros((10*seeds,11)) #((sourcenum*seeds,receptornum))
Yfinal = np.zeros((10*seeds,len(direction_sphcoords)))
for i in distancetest:
    
    for seed_particle in range(1,seeds+1): #10 sources corresponding each to the 10 directions -> therefore 10 arrays for each seed of fifo
        print(seed_particle, i)
        X, Y = datacreate(direction_sphcoords, sourcenum=10 ,sourceexact=direction_sphcoords,receptornum=10,particlenum=50, recepsurface_ratio=10, distanceexact=i,radiusexact=1,diffusionexact=1,rateexact=1,receptor_seed=1,particle_seed=seed_particle)
        Xfinal[10*(seed_particle-1):(10*(seed_particle-1)+ 9),:] = X
        Yfinal[10*(seed_particle-1):(10*(seed_particle-1)+ 9),:] = Y
    params = params_string('pick_direction(0,10)', sourcenum=10 ,receptornum=10,particlenum=50, recepsurface_ratio=10, distanceexact=i,radiusexact=1,diffusionexact=1,rateexact=1,receptor_seed=1,initial_source_seed=1,particle_seed=seed_particle)
    write_datafile('Xseed_distance='+str(i), params , Xfinal)
    write_datafile('Yseed_distance='+ str(i), params , Yfinal)


accuracy =[]
accuracynn = []
for i in distancetest:
    print(i)
    X = read_datafile('Xseed_distance='+str(i))
    Y = read_datafile('Yseed_distance='+str(i))
    training_x, predict_x,training_y, predict_y = train_test_split(X, Y)
    #print('training data amount percentage:' + str((len(training_x)/len(X))*100))
    mlp = train(training_x, training_y, layers_tuple = (100,50), max_iterations=5000)
    save_neural_network(mlp, distance=i)
    acc, probs, score, accnn = test(mlp, predict_x, predict_y,direction_sphcoords,0.5)
    accuracy.append(acc)
    accuracynn.append(accnn)
print(accuracy, 'accuracy harsh')
print(accuracynn, 'accuracy nearest neighbour')
