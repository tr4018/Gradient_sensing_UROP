import numpy as np

def analyticalreg(x,y): 
    #method to solve for beta analytically
    #(X^tX)^-1X^ty
    
    return np.matmul(np.matmul(np.linalg.inv(np.matmul(np.transpose(x),x)),np.transpose(x)),y)