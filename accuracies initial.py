import numpy as np
import matplotlib.pyplot as plt

acc1 = [0.9, 0.752, 0.668, 0.5, 0.424, 0.356, 0.312, 0.372, 0.264, 0.288, 0.296, 0.236, 0.244, 0.244, 0.224, 0.204, 0.22, 0.22, 0.228, 0.196, 0.204, 0.22, 0.26]
acc2 = [0.908, 0.764, 0.696, 0.536, 0.436, 0.364, 0.316, 0.32, 0.264, 0.248, 0.284, 0.224, 0.26, 0.248, 0.26, 0.228, 0.192, 0.216, 0.2, 0.228, 0.22, 0.204, 0.188]
acc3 = [0.932, 0.816, 0.652, 0.512, 0.464, 0.388, 0.308, 0.308, 0.288, 0.224, 0.276, 0.256, 0.248, 0.244, 0.268, 0.256, 0.232, 0.2, 0.216, 0.236, 0.276, 0.184, 0.244]
accnn25 = [0.94, 0.788, 0.712, 0.544, 0.488, 0.38, 0.352, 0.416, 0.348, 0.336, 0.352, 0.276, 0.292, 0.264, 0.256, 0.236, 0.252, 0.264, 0.26, 0.236, 0.248, 0.252, 0.3] 
accnn35 = [0.948, 0.892, 0.82, 0.708, 0.596, 0.564, 0.5, 0.532, 0.512, 0.408, 0.464, 0.412, 0.464, 0.464, 0.456, 0.408, 0.368, 0.408, 0.384, 0.392, 0.432, 0.364, 0.396]
accnn50 = [0.968, 0.904, 0.836, 0.74, 0.7, 0.684, 0.604, 0.536, 0.572, 0.5, 0.484, 0.52, 0.508, 0.492, 0.504, 0.508, 0.508, 0.436, 0.44, 0.512, 0.488, 0.42, 0.508]


distancetest= [2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10,10.5,11,11.5,12,12.5,13]


plt.plot(distancetest, acc1, color = 'black')
plt.plot(distancetest, acc2, color = 'black')
plt.plot(distancetest, acc3, color = 'black')
plt.plot(distancetest, accnn25, color = 'red')
plt.plot(distancetest, accnn35, color = 'green')
plt.plot(distancetest, accnn50, color = 'blue')

plt.show()
print('hello world')
