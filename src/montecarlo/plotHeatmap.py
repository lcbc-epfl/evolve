'''
plotHeatmap.py

@author: Nicholas Browning

'''


from numpy.random import uniform, seed
from matplotlib.mlab import griddata
import matplotlib.pyplot as plt
import numpy as np
# make up data.
# npts = int(raw_input('enter # of random points to plot:'))

from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import LogNorm
import time

general_data2 = np.fromfile('psiphi/rama500-general-nosec.data', sep=' ')
gly_data = np.fromfile('psiphi/rama500-gly-sym-nosec.data', sep=' ')
prepro_data = np.fromfile('psiphi/rama500-prepro.data', sep=' ')
pro_data = np.fromfile('psiphi/rama500-pro.data', sep=' ')

general_data = prepro_data

general_data = general_data.reshape(len(general_data) / 3, 3)
# x = data[:, 0]
# y = data[:, 1]
general_probability = general_data[:, 2]

general_probability = general_probability / max(general_probability)


psi_arr = np.zeros(200)
phi_arr = np.zeros(200)


N = int(len(general_probability) ** .5)


general_probability = general_probability.reshape(N, N)

iterN = 10000

np.random.seed()

plt.imshow(general_probability.T, cmap=plt.get_cmap("rainbow"), extent=[0, 180, 0, 180], origin="lower")
plt.colorbar()

# plt.ion()

for j in xrange (0, len(psi_arr)):
    phi_arr[j] = np.random.random_integers(low=0, high=179)
    psi_arr[j] = np.random.random_integers(low=0, high=179)


start_phi = np.copy(phi_arr)
start_psi = np.copy(psi_arr)

cl = np.arange(len(psi_arr))

print "Starting MonteCarlo: ", iterN, "steps"
for i in xrange (0, iterN):
    for j in xrange (0, len(psi_arr)):
        y_change = np.random.random_integers(low=-10, high=10)
        x_change = np.random.random_integers(low=-10, high=10)

        new_x = phi_arr[j] + x_change
        new_y = psi_arr[j] + y_change

        if (new_y < 0):
            new_y = 180 + new_y
        if (new_y > 179):
            new_y = new_y - 180
        if (new_x < 0):
            new_x = 180 + new_x
        if (new_x > 179):
            new_x = new_x - 180

        if (np.random.rand() < general_probability[new_x][new_y]):
            phi_arr[j] = new_x
            psi_arr[j] = new_y
print "Finished MonteCarlo"


plt.scatter(phi_arr, psi_arr, c=cl, s=15)

# phi_arr = -180 + 2 *phi_arr
# psi_arr =  -180 + 2 *psi_arr

plt.scatter(start_phi, start_psi, c='black', s=15)
# print indices


plt.show()

#
# fig = plt.figure()
# #ax = fig.add_subplot(111, projection='3d')
#
# # define grid.
#
# xi = np.linspace(-180, 180, 90)
# yi = np.linspace(-180, 180, 90)
#
# # grid the data.
# zi = griddata(x, y, z, xi, yi, interp='linear')
#
# print zi
# CS = plt.contourf(xi, yi, zi, 30, cmap=plt.cm.rainbow,
#                   vmax=abs(zi).max(), vmin=abs(zi).min())
# plt.colorbar()  # draw colorbar
# # plot data points.
#
#
# plt.show()
