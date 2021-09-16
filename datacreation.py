from IdealDirection import *
import numpy as np
from ReceptorMap import *
from ML_Brownian_Interface import *
from random_3d_rotation import random_3d_rotation
from sphericaltransf import *

def datacreate(
direction_sphcoords,
receptornum = 10,
recepsurface_ratio = 10,
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
particle_seed= 1): 
############## Parameters for receptor map ###################
#the variables ...exact are a shortcut to skipping all the for loops and plugging in one value per variable
#random_yn decides if we want to take a random uniform approach for the data (1) or if we want to do equally spaced values between 
#chosen boundaries (0).
    if np.all(sourceexact) == -1: pass #we create a random source
    else: 
        sourcenum = len(sourceexact)
        source_theta = sourceexact[:,0]
        source_phi = sourceexact[:,1]

    if diffusionexact== -1:
        if random_yn==0:
            diffusion_constants  = np.linspace(0.5,1,diffusionnum)
        elif random_yn==1:
            diffusion_constants = np.random.default_rng().uniform(0.5,1, diffusionnum) #too slow if close to 0
        else: 
            raise ValueError("Pick if diffusion constants should be equally spaced (random_yn = 0) or randomly chosen (random_yn = 1)")
    else: diffusion_constants = np.array([diffusionexact])
    if distanceexact== -1:
        if random_yn==0:
            distance_from_source = np.linspace((maxradius+maxdistance)/2,maxdistance,distancenum)
        elif random_yn==1:
            distance_from_source = np.random.default_rng().uniform(0,maxdistance, distancenum)
        else:
            raise ValueError("Pick if distance constants should be equally spaced (random_yn = 0) or randomly chosen (random_yn = 1)")
    else: distance_from_source = np.array([distanceexact])
    if rateexact == -1:
        if random_yn==0:
            rate = np.linspace(0.1,maxrate,ratenum)
        elif random_yn==1:
            rate = np.random.default_rng().uniform(0.5,maxrate, ratenum)
        else:
            raise ValueError("Pick if rate constants should be equally spaced (random_yn = 0) or randomly chosen (random_yn = 1)")
    else: rate=np.array([rateexact])
    if radiusexact== -1:
        if random_yn==0:
            radius_sphere = np.linspace(0.1,maxradius,radiusnum)
        elif random_yn==1:
            radius_sphere = np.random.default_rng().uniform(0.5,maxradius, radiusnum)
        else:
            raise ValueError("Pick if radius constants should be equally spaced (random_yn = 0) or randomly chosen (random_yn = 1)")
    else: radius_sphere = np.array([radiusexact])


########### LOOPs for data #################
    #fix number of receptors for each training data, it's like fixing the number of eyes the cell has... makes sense, I think.
    receptor_sphcoords,receptor_cartcoords, activation_array = init_Receptors(1,receptornum,0,receptor_seed)
    #print('receptornum in regular method:'+ str(len(receptor_sphcoords)))
    loops = sourcenum*len(radius_sphere)*len(distance_from_source)*len(rate)*len(diffusion_constants)
    X = np.zeros((loops,len(receptor_sphcoords)))
    Y = np.zeros((loops,direction_sphcoords.shape[0]))
    loop_count = 0
    np.random.seed(seed=initial_source_seed)
    source_theta_init = np.random.uniform(0,math.pi)
    source_phi_init = np.random.uniform(0,2*math.pi)   
    for s in range(1,sourcenum+1):
        if np.all(sourceexact) == -1:
            source_theta,source_phi = random_3d_rotation(source_phi_init,source_theta_init,s)
            sourcex,sourcey,sourcez = spherical2cart_point(source_theta,source_phi)
        #function to relate source coordinates to action direction -> make Y vector
            move = ideal_direction(source_theta,source_phi,direction_sphcoords, 1)
        else:  #use exact coordinates of source that we can control
            sourcex,sourcey,sourcez = spherical2cart_point(source_theta[s-1],source_phi[s-1])
            move = ideal_direction(source_theta[s-1],source_phi[s-1],direction_sphcoords,1)
        for r in radius_sphere:
            mindistance = r*math.pi/recepsurface_ratio
            #visualize_Receptors(receptor_cartcoords,r,mindistance)            

            for distance in distance_from_source:
                sx = sourcex * distance
                sy = sourcey * distance
                sz = sourcez * distance
                for ra in rate:
                    for dif in diffusion_constants:
                        loop_count+=1
                        activation_array = np.zeros((1,len(receptor_sphcoords)))
                        #needs source position and radius to be included in parameters
                        brownian_pipe,received,source = init_BrownianParticle(sx,sy,sz,rate=ra,radius=r,diffusion=dif, use_seed=particle_seed) 
                            #what is the seed for?
                            #do we fix parameters training,cutoff,events,iterations? 
                        count = 1 #count how many particles in one activation array measure. Starts with 1 particle.
                        while(count <= particlenum):
                            theta_mol = received[0]
                            phi_mol = received[1]
                            ind = activation_Receptors(theta_mol,phi_mol,receptor_sphcoords,r,mindistance)
                            if ind == -1: pass
                            else: activation_array[0,ind[0]] += 1
                            received,source = update_BrownianParticle(brownian_pipe)
                            count+=1
                        stop_BrownianParticle(brownian_pipe)
                        X[loop_count-1,:] = activation_array
                        if (activation_array==0).all():
                            Y[loop_count-1,:] = np.zeros((1,len(direction_sphcoords)))
                        else: 
                            Y[loop_count-1,:] = move
                        #print(str(X[loop_count-1,:])+'\t'+str(Y[loop_count-1,:]))
                            
    return X, Y
