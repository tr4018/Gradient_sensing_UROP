from datawriteread import read_datafile
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import MinMaxScaler
from sklearn.neural_network import MLPClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
import numpy as np
from haversine import * 
import pickle

def fit_mlp(X, Y, layers_tuple, max_iterations):
    #function to scale X to between 0 and 1 which is best for neural network, then creates mlp object 
    X = MinMaxScaler().fit_transform(X)
    #mlp = MLPClassifier(hidden_layer_sizes=layers_tuple,random_state=0, max_iter=max_iterations, solver='sgd', learning_rate='constant',\
                     #   momentum=0, learning_rate_init=0.2) 
    mlp = MLPClassifier(hidden_layer_sizes=layers_tuple, solver='adam', max_iter=max_iterations) # lbfgs for small data
    #layers_tuple: Each element in the tuple is the number of nodes at the ith position. 
    #Length of tuple denotes the total number of layers.
    mlp.fit(X, Y)
    
    return mlp

def predict(mlp, X):
    #function to scale X to between 0 and 1 and then use the Neural Network to produce predictions for Y
    X = MinMaxScaler().fit_transform(X)
    y = mlp.predict(X)
    
    return y

def accuracy(true_y, predicted_y):
    score = 0
    for true, predicted in zip(true_y, predicted_y):
        is_there_a_zero = np.linalg.norm(true - predicted)
        if is_there_a_zero == 0:
            score += 1
        
    return score/len(true_y)

def separate_train_set(X,Y):
    #if data is ordered use this to randomise it so every source position is trained on, then use first half of data to train the network
    #use second half of the data to predict if network is working
    num_rows, num_cols = np.shape(X)
    indices = np.random.choice(range(num_rows), num_rows, replace=False) 
    n_train = int(len(indices)/2)
    training_x = np.array(X)[indices[:n_train].astype(int)]
    training_y = np.array(Y)[indices[:n_train].astype(int)]
    predict_x = np.array(X)[indices[n_train:].astype(int)] 
    predict_y = np.array(Y)[indices[n_train:].astype(int)]
    return training_x, training_y, predict_x, predict_y

def train(training_x, training_y, layers_tuple, max_iterations):
    #function to train our neural network with half the data given
    mlp = fit_mlp(training_x, training_y, layers_tuple, max_iterations)
    return mlp

def test(mlp, predict_x, predict_y, direction_sphcoords, frac):
    pred = predict(mlp, predict_x)
    score = mlp.score(predict_x, predict_y)
    acc = accuracy_score(predict_y,pred)
    directprob = direction_probabilities(mlp, predict_x)
    accnn = nearest_neighbours_accuracy(direction_sphcoords,predict_y,pred,frac)
    #print("Accuracy of MLPClassifier : ", acc)
    #print("Probabilities of each direction : ", directprob)
    return acc, directprob, score, accnn
    
def save_neural_network(mlp, particlenum=None,receptornum=None,diffusion=None,rate=None,cutoff=None,events=None,iterations=None, distance=None):
    filename = 'MLPClassifier'
    if particlenum != None:
        filename += ' -p '+str(particlenum)
    if receptornum != None:
        filename += ' -r '+str(receptornum)
    if rate != None:
        filename += ' -R '+str(rate)
    if cutoff != None:
        filename += ' -c '+str(cutoff)
    if events != None:
        filename += ' -E '+str(events)
    if iterations != None:
        filename += ' -N '+str(iterations)
    if diffusion != None:
        filename += ' -d '+str(diffusion)
    if distance != None:
        filename += ' -D '+str(distance)
    pickle.dump(mlp, open(filename, 'wb'))
    return filename
    
    
def load_neural_network(filename):
    restored_mlp = pickle.load(open(filename, 'rb'))
    return restored_mlp
          
def direction_probabilities(mlp, X):
    return mlp.predict_proba(X)

def nearest_neighbours_accuracy(direction_sphcoords, true_y, predicted_y, frac_area, radius=1):
        # uses the fraction of area of sphere you want covered to locate nearest directions, and see whether the predicted value was included in these

    true_y = list(true_y)
    predicted_y = list(predicted_y)

    neighbours = np.array(find_nearest_neighbours(frac_area, radius, direction_sphcoords), dtype=object)
    i = -1
    score = 0
    for trued, predicted in zip(true_y, predicted_y):
        trued = list(trued)
        i+=1
        # we have to deal with move arrays of all 0's
        
        if (sum(trued)==0): 
            if (sum(predicted)==0): 
                score += 1
            else:
                score += 0
        else:    
            idx = trued.index(1) 
            is_there_a_zero = np.linalg.norm(neighbours[idx] - predicted, axis=1)
            if np.all(is_there_a_zero) == 0: #if there is then score augments
                score += 1

    return score/len(true_y)
    
def find_nearest_neighbours(frac_area, radius, direction_sphcoords):
    # uses the fraction of area of sphere you want covered to locate nearest directions. Stores these as a list of lists, where the element i is a list of the nearest neighbours for direction i
    cap_area = frac_area * 4 * np.pi * np.power(radius,2)
    dtheta = np.arccos(1-cap_area/(2 * np.pi * np.power(radius,2)))
    max_distance = haversine(radius,0,0,dtheta,0)
    directionnum=len(direction_sphcoords)
    distances = []
    neighbours = []
    
    for coords in direction_sphcoords:
        distances.append(haversine(radius,coords[0],coords[1],direction_sphcoords[:,0].reshape(directionnum,1), direction_sphcoords[:,1].reshape(directionnum,1)))
    for d in distances:
        idx = []
        j = 0
        
        for elem in d:
            j += 1
            if elem <= max_distance:
                idx.append(j-1)
    
        best_directions = []
        for i in idx:
            best_direction = np.zeros(len(direction_sphcoords))
            best_direction[i] = 1
            best_directions.append(best_direction)

        neighbours.append(best_directions)
            
    return neighbours

