#variable_to_save1 = variable_to_save2 = []
#variable_to_save1.append([indiv.phi_dihedrals[j], indiv.psi_dihedrals[j]])

'''           
            variable_to_save2.append([indiv.phi_dihedrals[j], indiv.psi_dihedrals[j]])

    # Save an object in a file
    import pickle

    # Write
    with open('dev_utils/MC_Dihedrals', 'wb') as my_file1: # write and binary, file data is created if not existing
        my_pickler1 = pickle.Pickler(my_file1)
        my_pickler1.dump(variable_to_save1)

    with open('dev_utils/MCBasilisk_Dihedrals', 'wb') as my_file2:
        my_pickler2 = pickle.Pickler(my_file2)
        my_pickler2.dump(variable_to_save2)
'''

# libraries
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import kde
import pickle

# Read
with open('dev_utils/MC_Dihedrals', 'rb') as my_file1:
    my_unplicker = pickle.Unpickler(my_file1)
    MC_dihedrals = my_unplicker.load()

with open('dev_utils/MCBasilisk_Dihedrals', 'rb') as my_file2:
    my_unplicker = pickle.Unpickler(my_file2)
    MCBasilisk_dihedrals = my_unplicker.load()
    
# create data
psi_MC, phi_MC  = np.array(map(list, zip(*MC_dihedrals)))
psi_MC_basilisk, phi_MC_basilisk = np.array(map(list, zip(*MCBasilisk_dihedrals)))


x = psi_MC
y = phi_MC
 
# Evaluate a gaussian kde on a regular grid of nbins x nbins over data extents
nbins=300
k = kde.gaussian_kde([x,y])
xi, yi = np.mgrid[x.min():x.max():nbins*1j, y.min():y.max():nbins*1j]
zi = k(np.vstack([xi.flatten(), yi.flatten()]))
 
# Add color bar
#plt.figure(1)
plt.subplot(2,1,1)
plt.pcolormesh(xi, yi, zi.reshape(xi.shape))

plt.xlabel('$\phi$ [degrees]')
plt.ylabel('$\psi$ [degrees]')
plt.title('Distribution of Dihedrals from MC generator')
#plt.text(60, .025, r'$\mu=100,\ \sigma=15$')
plt.axis([-180, 180, -180, 180])
plt.grid(False)
plt.colorbar()

x = psi_MC_basilisk
y = phi_MC_basilisk
 
# Evaluate a gaussian kde on a regular grid of nbins x nbins over data extents
nbins=300
k = kde.gaussian_kde([x,y])
xi, yi = np.mgrid[x.min():x.max():nbins*1j, y.min():y.max():nbins*1j]
zi = k(np.vstack([xi.flatten(), yi.flatten()]))
 
# Add color bar
plt.subplot(2,1,2)
plt.pcolormesh(xi, yi, zi.reshape(xi.shape))

plt.xlabel('$\phi$ [degrees]')
plt.ylabel('$\psi$ [degrees]')
plt.title('Distribution of Dihedrals after MC-Basilisk')
#plt.text(60, .025, r'$\mu=100,\ \sigma=15$')
plt.axis([-180, 180, -180, 180])
plt.grid(False)
plt.colorbar()
plt.show()