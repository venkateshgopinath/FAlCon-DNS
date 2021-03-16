# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

def equatContour( data, radius, label=None, levels=65, cm='seismic',
                  normed=True, vmax=None, vmin=None, cbar=True, tit=True,
                  normRad=False):
    """
    Plot the equatorial cut of a given field

    :param data: the input data (an array of size (nphi,nr))
    :type data: numpy.ndarray
    :param radius: the input radius
    :type radius: numpy.ndarray
    :param label: the name of the input physical quantity you want to
                  display
    :type label: str
    :param normRad: when set to True, the contour levels are normalised
                    radius by radius (default is False)
    :type normRad: bool
    :param levels: the number of levels in the contour
    :type levels: int
    :param cm: name of the colormap ('jet', 'seismic', 'RdYlBu_r', etc.)
    :type cm: str
    :param tit: display the title of the figure when set to True
    :type tit: bool
    :param cbar: display the colorbar when set to True
    :type cbar: bool
    :param vmax: maximum value of the contour levels
    :type vmax: float
    :param vmin: minimum value of the contour levels
    :type vmin: float
    :param normed: when set to True, the colormap is centered around zero.
                   Default is True, except for entropy/temperature plots.
    :type normed: bool
    """

    nphi, ntheta = data.shape

    phi = np.linspace(0., 2.*np.pi, nphi)
    rr, pphi = np.meshgrid(radius, phi)
    xx = rr * np.cos(pphi)
    yy = rr * np.sin(pphi)

    if normRad: # Normalise each radius
        maxS = np.sqrt(np.mean(data**2, axis=0))
        data[:, maxS!=0.] /= maxS[maxS!=0.]

    if tit and label is not None:
        if cbar:
            fig = plt.figure(figsize=(6.5,5.5))
            ax = fig.add_axes([0.01, 0.01, 0.76, 0.9])
        else:
            fig = plt.figure(figsize=(5,5.5))
            ax = fig.add_axes([0.01, 0.01, 0.98, 0.9])
        ax.set_title(label, fontsize=24)
    else:
        if cbar:
            fig = plt.figure(figsize=(6.5,5))
            ax = fig.add_axes([0.01, 0.01, 0.76, 0.98])
        else:
            fig = plt.figure(figsize=(5, 5))
            ax = fig.add_axes([0.01, 0.01, 0.98, 0.98])

    cmap = plt.get_cmap(cm)
    if vmax is not None or vmin is not None:
        normed = False
        cs = np.linspace(vmin, vmax, levels)
        im = ax.contourf(xx, yy, data, cs, cmap=cmap, extend='both')
    else:
        cs = levels
        im = ax.contourf(xx, yy, data, cs, cmap=cmap)
        #im = ax.pcolormesh(xx, yy, data, cmap=cmap, antialiased=True)
    ax.plot(radius[0]*np.cos(phi), radius[0]*np.sin(phi), 'k-', lw=1.5)
    ax.plot(radius[-1]*np.cos(phi), radius[-1]*np.sin(phi), 'k-', lw=1.5)

    ax.axis('off')

    # Add the colorbar at the right place
    pos = ax.get_position()
    l, b, w, h = pos.bounds
    if cbar:
        if tit and label is not None:
            cax = fig.add_axes([0.85, 0.46-0.7*h/2., 0.03, 0.7*h])
        else:
            cax = fig.add_axes([0.85, 0.5-0.7*h/2., 0.03, 0.7*h])
        mir = fig.colorbar(im, cax=cax)

    # Normalise data 
    if normed:
        im.set_clim(-max(abs(data.max()), abs(data.min())),
                     max(abs(data.max()), abs(data.min())))

    return fig, xx, yy

def lineContour(data,radius,ri,ro):

    nphi, ntheta = data.shape

    phi = np.linspace(0., 2.*np.pi, nphi)
    rr, pphi = np.meshgrid(radius, phi)
    xx = rr * np.cos(pphi)
    yy = rr * np.sin(pphi)
    fig2, ax = plt.subplots(1)
    ax.set_aspect('equal')
    fig2 = plt.contour(xx, yy, data, 35, colors='black',extent=None);
    circ1 = Circle((0, 0), ro, facecolor='None', edgecolor='k', lw=1) 
    circ2 = Circle((0, 0), ri, facecolor='None', edgecolor='k', lw=1) 
    ax.add_patch(circ1)
    ax.add_patch(circ2)
    #ax.set_xticks([])
    #ax.set_yticks([])
    ax.axis('off')
    return fig2, xx, yy
