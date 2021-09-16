import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


# Load Data
path = 'BrownianParticle_RUNRUN01/'
df1 = pd.read_csv(path+'BrownianParticle_RUNRUN01_1_1000_P01_20210705_173431.txt', header = None, sep=" ")
df1.head()

df2 = pd.read_csv(path+'BrownianParticle_RUNRUN01_2_1000_P01_20210705_173431.txt', header = None, sep=" ")
df2.head()

df3 = pd.read_csv(path+'BrownianParticle_RUNRUN01_3_1000_P01_20210705_173431.txt', header = None, sep=" ")
df3.head()

df4 = pd.read_csv(path+'BrownianParticle_RUNRUN01_4_1000_P01_20210705_173431.txt', header = None, sep=" ")
df4.head()

df5 = pd.read_csv(path+'BrownianParticle_RUNRUN01_5_1000_P01_20210705_173431.txt', header = None, sep=" ")
df5.head()

df6 = pd.read_csv(path+'BrownianParticle_RUNRUN01_6_1000_P01_20210705_173431.txt', header = None, sep=" ")
df6.head()

df7 = pd.read_csv(path+'BrownianParticle_RUNRUN01_7_1000_P01_20210705_173431.txt', header = None, sep=" ")
df7.head()

df8 = pd.read_csv(path+'BrownianParticle_RUNRUN01_8_1000_P01_20210705_173431.txt', header = None, sep=" ")
df8.head()

df9 = pd.read_csv(path+'BrownianParticle_RUNRUN01_9_1000_P01_20210705_173431.txt', header = None, sep=" ")
df9.head()


theta = pd.concat([ df1.iloc[:, 2] , df2.iloc[:, 2], df3.iloc[:, 2], df4.iloc[:, 2], df5.iloc[:, 2] , df6.iloc[:, 2] , df7.iloc[:, 2] , df8.iloc[:, 2] , df9.iloc[:, 2]  ], axis=1)
theta.columns = ['0', '1', '2', '3', '4', '5', '6', '7', '8']
phi = pd.concat([df1.iloc[:, 3] , df2.iloc[:, 3], df3.iloc[:, 3], df4.iloc[:, 3], df5.iloc[:, 3] , df6.iloc[:, 3] , df7.iloc[:, 3] , df8.iloc[:, 3] , df9.iloc[:, 3]  ], axis=1)
phi.columns = ['0', '1', '2', '3', '4', '5', '6', '7', '8']
# ======================== Part 1 - Visualise the data ======================

fig, axes = plt.subplots(3,3, sharex=False, sharey=False, figsize=(15,14))
mylegend = [r'$RUN = 1$', r'$RUN = 2$', r'$RUN = 3$', r'$RUN = 4$', r'$RUN = 5$', r'$RUN = 6$', r'$RUN = 7$', r'$RUN = 8$', r'$RUN = 9$']
color = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']

for i, ax in enumerate(axes.flatten()):    
    n, bins = np.histogram(theta.iloc[:,i], bins=100, range=(0,np.pi))
    width = 0.99 * (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2
    ax.set_ylim(0,max(n/np.sin(center)))
    ax.bar(center, n/np.sin(center), align='center', width=width, color=color[i], label=mylegend[i])
    ax.set_xlabel(r'$\Theta$', fontsize='large')
    ax.set_ylabel(r'$Frequency/sin({\Theta})$')
    ax.legend(fontsize = 'large',  loc='upper right')
    plt.tight_layout()
#plt.show()
plt.legend()
plt.savefig('Theta Distribution for Gradient Sensing.png')

plt.clf()
    
fig, axes = plt.subplots(3,3, sharex=False, sharey=False, figsize=(15,14))
mylegend = [r'$RUN = 1$', r'$RUN = 2$', r'$RUN = 3$', r'$RUN = 4$', r'$RUN = 5$', r'$RUN = 6$', r'$RUN = 7$', r'$RUN = 8$', r'$RUN = 9$']
color = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']

for i, ax in enumerate(axes.flatten()):    
    n, bins = np.histogram(phi.iloc[:,i], bins=100, range=(0,np.pi))
    width = 0.99 * (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2
    ax.set_ylim(0,max(n))
    ax.bar(center, n, align='center', width=width, color=color[i], label=mylegend[i])
    ax.set_xlabel(r'$\Phi$', fontsize='large')
    ax.set_ylabel(r'$Frequency$')
    ax.legend(fontsize = 'large',  loc='upper right')
    plt.tight_layout()
#plt.show()

plt.legend()
plt.savefig('Phi Distribution for Gradient Sensing.png')