from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from H_surface import Hsurface
import numpy as np


def plotHsurface(Hsup: Hsurface, cmap, plotHs = [True, True], log = False, pointSize = 5.0, figSize = (10, 8), alfa = 1.0, orientation = [30, 45], project = False):

    H = Hsup.H
    if log: H = np.log10(np.abs(H))

    minH = H.min()
    maxH = H.max()
    H = (H-minH)/(maxH-minH)

    znorm = H - minH
    znorm /= znorm.ptp()
    znorm.min(), znorm.max()
    color = cmap(znorm)
    
    if project:
        if plotHs[0]:
            # *** 2D H surface projection ***
            fig1 = plt.figure(figsize=figSize)
            ax1 = plt.axes()
            # H surface plot:
            ax1.scatter(Hsup.e1, Hsup.deltaW, c=color, s=pointSize, alpha = alfa)

        if plotHs[1]:
            # *** 2D H surface projection ***
            fig2 = plt.figure(figsize=figSize)
            ax2 = plt.axes()
            # H surface plot:
            ax2.scatter(Hsup.e2, Hsup.deltaW, c=color, s=pointSize, alpha = alfa)
    else:
        if plotHs[0]:
            # *** 3D classic matplotlib plot (not interactive) ***
            fig1 = plt.figure(figsize=figSize)
            ax1 = plt.axes(projection='3d')
            # H surface plot:
            ax1.scatter3D(Hsup.e1, Hsup.deltaW, Hsup.sigma, c=color, s=pointSize, alpha = alfa)
            ax1.view_init(elev=orientation[0], azim=orientation[1])

        if plotHs[1]:
            # *** 3D classic matplotlib plot (not interactive) ***
            fig2 = plt.figure(figsize=figSize)
            ax2 = plt.axes(projection='3d')
            # H surface plot:
            ax2.scatter3D(Hsup.e2, Hsup.deltaW, Hsup.sigma, c=color, s=pointSize, alpha = alfa)
            ax2.view_init(elev=orientation[0], azim=orientation[1])

    plt.show()

# ************************

# Customized colormap to plot the model in grayscale and the numerical integration in colors
class Hcolormap:
    def __init__(self, stack = True, resol = 1, Ncm = 40, Nseg2 = 1, cmap_used = 'gray'):
        if stack:
            # Get the colormaps for the top and bottom halves
            top = cm.get_cmap(cmap_used + '_r', resol)
            bottom = cm.get_cmap(cmap_used, resol)

            # Concatenate the colormaps to create newcolors
            newcolors = np.vstack((top(np.linspace(0, 1, resol)), bottom(np.linspace(0, 1, resol))))

            # Calculate the number of times to repeat the colormap
            repetitions = Ncm // (2 * Nseg2 + 1)

            # Concatenate the colormap multiple times
            newcolors = [newcolors] * repetitions

            # Flatten the list of arrays and concatenate them
            newcolors = np.concatenate(newcolors)

            self.newcm = ListedColormap(newcolors, name='Ngrayscale') # newcm is the new colormap used in the next cell

        else:

            self.newcm = cm.get_cmap(cmap_used) # use the already existent matplotlib library cmap

