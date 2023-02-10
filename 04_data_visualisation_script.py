from matplotlib import pyplot as plt
import numpy as np
import matplotlib.colors as mcolors
from matplotlib import cm
import matplotlib as mpl
from matplotlib import rc
import pylab
rc('mathtext', default='regular')

####
#### Copy the data from the Excel file and make the final plot
####
energy = np.array([500, 90, 95, 200, 62, 300, 400, 800, 70, 800, 225, 15, 127])
intensity = np.array([3e9, 1e13, 1e7, 2e12, 1e13, 1e9, 1e10, 1e10, 1e9, 1.5e9, 1e11, 1e9, 1e9])
dose_rate = np.array([300, 500, 0.03, 40, 300, 1e9, 120, 40, 10, 15, 330, 1.8e9, 9.7e8])# previous values 7e8, 2.6e9])
n_ion_species = np.array([10, 10, 4, 2, 3, 4, 7, 5, 3, 5, 1, 2, 4])
facilities = [r'FAIR', r'KVI-CART', r'GANIL', r'iTHEMBA', r'LNS', r'ELIMAIA', r'CNAO', r'HIMAC', r'LNL-SPES', r'NICA', r'TIFPA', r'LhARA (s1)', r'LhARA'+ "\n" + r'(s2)']
mks = ['o','o','o','o','o','o','o','o','o','o','o','o','o']

# 50 values for later use from 0 to 1
#greens = cm.cividis(np.linspace(-2,3, num=100))
# 50 values red for later use from 1.5 to 2.5 
#red = [(1,0,0,1)]*len(np.linspace(0,0.5))

#colors = np.vstack((greens, red))
# in total we now have 175 colors in the colormap
#mycmap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)

# cmap = mpl.colors.LinearSegmentedColormap.from_list('custom', 
#                                              [(0,    'red'),
#                                               (0.5, 'green'),
#                                               (1,    'Orange')], N=126)
#norm = mpl.colors.Normalize(vmin=2, vmax=10)
#mpl.colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm) 


#cmap = mpl.colors.ListedColormap([(0.687, 0.525, 0.931), (0.829, 0.473, 0.965), (0.753, 0.0, 0.561), (0.933, 0.337, 0.349), (1.0, 0.667, 0.71), (0.906, 1.0, 0.0)])
cmap = plt.get_cmap('plasma', 6)

# cNorm = mcolors.Normalize(vmin=7, vmax=13) #re-wrapping normalization
# scalarMap = cm.ScalarMappable(norm=cNorm, cmap=cmap)

# Scatter the points, using size and color but no label
fig, ax = plt.subplots()
ax.tick_params(direction='in', length=10, width=1, color='grey', labelsize=26)

plt.scatter(dose_rate[0:1], energy[0:1], label=None, c=np.log10(intensity)[0:1], cmap=cmap, norm=mcolors.Normalize(vmin=7, vmax=13), s=300*n_ion_species[0:1], linewidth=1, marker='o', alpha=0.8, edgecolor='black', zorder=10) #planned
plt.scatter(dose_rate[1:2], energy[1:2], label=None, c=np.log10(intensity[1:2]), cmap=cmap, norm=mcolors.Normalize(vmin=7, vmax=13),s=300*n_ion_species[1:2], linewidth=1, marker='^', alpha=0.8, edgecolor='black', zorder=10) #in operation
plt.scatter(dose_rate[2:3], energy[2:3], label=None, c=np.log10(intensity[2:3]), cmap=cmap, norm=mcolors.Normalize(vmin=7, vmax=13),s=300*n_ion_species[2:3], linewidth=1, marker='D', alpha=0.8, edgecolor='black', zorder=10) #comissioning
plt.scatter(dose_rate[3:4], energy[3:4], label=None, c=np.log10(intensity[3:4]), cmap=cmap, norm=mcolors.Normalize(vmin=7, vmax=13),s=300*n_ion_species[3:4], linewidth=1, marker='D', alpha=0.8, edgecolor='black', zorder=10) #re-comissioning
plt.scatter(dose_rate[4:5], energy[4:5], label=None, c=np.log10(intensity[4:5]), cmap=cmap, norm=mcolors.Normalize(vmin=7, vmax=13),s=300*n_ion_species[4:5], linewidth=1, marker='o', alpha=0.8, edgecolor='black', zorder=10) #upgrading
plt.scatter(dose_rate[5:6], energy[5:6], label=None, c=np.log10(intensity[5:6]), cmap=cmap, norm=mcolors.Normalize(vmin=7, vmax=13),s=300*n_ion_species[5:6], linewidth=1, marker='o', alpha=0.8, edgecolor='black', zorder=10) #planned
plt.scatter(dose_rate[6:7], energy[6:7], label=None, c=np.log10(intensity[6:7]), cmap=cmap, norm=mcolors.Normalize(vmin=7, vmax=13),s=300*n_ion_species[6:7], linewidth=1, marker='^', alpha=0.8, edgecolor='black', zorder=10) #in operation
plt.scatter(dose_rate[7:8], energy[7:8], label=None, c=np.log10(intensity[7:8]), cmap=cmap, norm=mcolors.Normalize(vmin=7, vmax=13),s=300*n_ion_species[7:8], linewidth=1, marker='^', alpha=0.8, edgecolor='black', zorder=10) #in operation
plt.scatter(dose_rate[8:9], energy[8:9], label=None, c=np.log10(intensity[8:9]), cmap=cmap, norm=mcolors.Normalize(vmin=7, vmax=13),s=300*n_ion_species[8:9], linewidth=1, marker='D', alpha=0.8, edgecolor='black', zorder=10) #comissioning
plt.scatter(dose_rate[9:10], energy[9:10], label=None, c=np.log10(intensity[9:10]), cmap=cmap, norm=mcolors.Normalize(vmin=7, vmax=13),s=300*n_ion_species[9:10], linewidth=1, marker='o', alpha=0.8, edgecolor='black', zorder=10) #planned
plt.scatter(dose_rate[10:11], energy[10:11], label=None, c=np.log10(intensity[10:11]), cmap=cmap, norm=mcolors.Normalize(vmin=7, vmax=13),s=300*n_ion_species[10:11], linewidth=1, marker='^', alpha=0.8, edgecolor='black', zorder=10) #in operation
plt.scatter(dose_rate[11:12], energy[11:12], label=None, c=np.log10(intensity[11:12]), cmap=cmap, norm=mcolors.Normalize(vmin=7, vmax=13),s=300*n_ion_species[11:12], linewidth=1, marker='o', alpha=0.8, edgecolor='black', zorder=10) #planned
plt.scatter(dose_rate[12:13], energy[12:13], label=None, c=np.log10(intensity[12:13]), cmap=cmap, norm=mcolors.Normalize(vmin=7, vmax=13),s=300*n_ion_species[12:13], linewidth=1, marker='o', alpha=0.8, edgecolor='black', zorder=10) #planned


#plt.scatter(dose_rate, energy, label=None, c=np.log10(intensity), cmap=cmap, s=300*n_ion_species, linewidth=1, marker='o', alpha=0.8, edgecolor='black', zorder=10)

ax.annotate(facilities[0], (dose_rate[0]+300, energy[0]+40), zorder=20, fontsize=28) #FAIR
ax.annotate(facilities[1], (dose_rate[1]+250, energy[1]+40), zorder=20, fontsize=28) #KVI-CART
ax.annotate(facilities[2], (dose_rate[2]-0.027, energy[2]+45), zorder=20, fontsize=28) #GANIL
ax.annotate(facilities[3], (dose_rate[3]-39.5, energy[3]+30), zorder=20, fontsize=28) #iTHEMBA
ax.annotate(facilities[4], (dose_rate[4]-225, energy[4]-70), zorder=20, fontsize=28) #LNS 
ax.annotate(facilities[5], (dose_rate[5]/100, energy[5]+20), zorder=20, fontsize=28) # ELIMAIA
ax.annotate(facilities[6], (dose_rate[6]-117, energy[6]+10), zorder=20, fontsize=28) # CNAO
ax.annotate(facilities[7], (dose_rate[7]+60, energy[7]-50), zorder=20, fontsize=28) #HIMAC
ax.annotate(facilities[8], (dose_rate[8]-9.8, energy[8]+45), zorder=20, fontsize=28) #LNL-SPES
ax.annotate(facilities[9], (dose_rate[9]-14.5, energy[9]-50), zorder=20, fontsize=28) #NICA
ax.annotate(facilities[10], (dose_rate[10], energy[10]+30), zorder=20, fontsize=28) #TIFPA
ax.annotate(facilities[11], (dose_rate[11]/300, energy[11]-10), zorder=20, fontsize=28) #LhARA stage 1
ax.annotate(facilities[12], (dose_rate[12]/20, energy[12]), zorder=20, fontsize=28) #LhARA stage 2

plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.xlim([10e-4, 10e9])

plt.axis(aspect='equal')
plt.xlabel(r'Dose rate (Gy/s)', fontsize=26)
plt.ylabel(r'Energy (Mev/u)',fontsize=26)
plt.semilogx()
cb = plt.colorbar()
plt.grid(True, alpha=0.7, which='both', zorder=-10)
ax.set_facecolor((0.965, 0.878, 0.839))

cb.ax.tick_params(labelsize=26)
#cb.ax.set_yticklabels(['7', '8', '9', '10', '11', '12', '13'], )
cb.set_label(label=r'log$_{10}$(Intensity) (pps)',size=26)
ax.set_facecolor('bisque')
#ax.set_facecolor('aliceblue')
#ax.set_facecolor('lightcyan')

line1 = pylab.Line2D(range(1),range(1),color='darkblue',marker='^',markersize=20, markerfacecolor="snow",alpha=1.0,linewidth=0, markeredgewidth=2)
line2 = pylab.Line2D(range(1),range(1),color='darkblue',marker='D',markersize=20, markerfacecolor="snow",alpha=1.0,linewidth=0, markeredgewidth=2)
line3 = pylab.Line2D(range(1),range(1),color='darkblue',marker='o',markersize=20, markerfacecolor="snow",alpha=1.0,linewidth=0, markeredgewidth=2)

plt.legend((line1,line2,line3),('in operation','comissioning','planned'),numpoints=1,loc=1, fontsize=26)

# lgnd.legendHandles[0]._sizes = [400]
# lgnd.legendHandles[1]._sizes = [400]
# lgnd.legendHandles[2]._sizes = [400]

plt.show()