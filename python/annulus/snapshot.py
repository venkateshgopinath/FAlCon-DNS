# -*- coding: utf-8 -*-
import numpy as np
import sys
import matplotlib.pyplot as plt
from .plotlib import equatContour, lineContour
from .cheblib import chebgrid, spec_spat, scanDir, chebforward, chebinverse
from .costf import costf, get_dr
from .npfile import *
import cmath


class Snapshot:

    def __init__(self, filename=None, endian='l', iplot=False, cm='seismic',
                 levels=65, cbar=True, title=True):

        if filename is None:
            filenames = scanDir('snapshot_plot_*')
            if len(filenames) > 0:
                filename = filenames[-1]
            else:
                sys.exit('No snapshot found!')

        infile = npfile(filename, endian=endian)
        self.ra, self.pr, self.radratio = infile.fort_read('Float64')
        self.dt, self.time = infile.fort_read('Float64')
        self.n_m, self.n_r = infile.fort_read('int32')

        self.n_phi_max = 3*self.n_m
        self.ri = self.radratio/(1.-self.radratio)
        self.ro = 1./(1.-self.radratio)
        tanom = np.zeros([self.n_phi_max,self.n_r])
        # to be removed 
        #self.radius = chebgrid(self.n_r-1, self.ri, self.ro)
        # and replaced by
        self.radius = infile.fort_read('Float64')

        # Reading temperature and bring back to physical space
        data = infile.fort_read('Complex64')
        data = data.reshape((self.n_r, self.n_m+1))
        self.tFR = data.T # Transpose for compatibility with Fortran array
        self.temperature = spec_spat(self.tFR, self.n_phi_max)
        #self.tanom = self.temperature - np.mean(self.temperature, axis=0)
         
        for i in range(0,self.n_phi_max):
            for j in range(0,self.n_r):
                tanom[i,j]=self.temperature[i,j]-((np.log(self.radius[j])-np.log(self.ro))/np.log(self.ri/self.ro)) 
        self.tanom = tanom
        # Reading vorticity and bring back to physical space
        data = infile.fort_read('Complex64')
        data = data.reshape((self.n_r, self.n_m+1))
        self.omgFR = data.T
        self.vorticity = spec_spat(self.omgFR, self.n_phi_max)

        # Reading radial velocity and bring back to physical space
        data = infile.fort_read('Complex64')
        data = data.reshape((self.n_r, self.n_m+1))
        self.urFR = data.T
        self.vr = spec_spat(self.urFR, self.n_phi_max)

        # Reading azimuthal velocity and bring back to physical space
        data = infile.fort_read('Complex64')
        data = data.reshape((self.n_r, self.n_m+1))
        self.upFR = data.T
        self.vphi = spec_spat(self.upFR, self.n_phi_max)
        
        # Vorticity recreated from snapshots
        self.Nm_max = self.n_m
        self.Nr_max = self.n_r
        self.vr_c = chebforward(self.urFR,self.Nm_max,self.Nr_max)
        self.vr_b = chebinverse(self.vr_c,self.Nm_max,self.Nr_max,self.Nr_max) 
        self.omgFR_snap = np.zeros((self.Nm_max+1,self.Nr_max))

        self.dupdr = get_dr(self.upFR)
        self.omgFR_snap = -self.dupdr +self.upFR/self.radius
        for m in range(1,self.Nm_max+1):
            self.omgFR_snap[m,:] =  self.omgFR_snap[m,:] - 1j*m*self.urFR[m,:]/self.radius
        
        self.vorticity_snap = spec_spat(self.omgFR_snap, self.n_phi_max)
        ## Reading omgFR from psi and bring back to physical space
        #data = infile.fort_read('Complex64')
        #data = data.reshape((self.n_r, self.n_m+1))
        #self.omgFR_check = data.T
        #self.omg_check = spec_spat(self.upFR, self.n_phi_max)

        infile.close()

        if iplot:
            self.cplot(cm=cm, levels=levels, cbar=cbar, title=title)
            #self.lplot(title=title)

    def cplot(self, cm='jet', levels=65, cbar=True, title=True):


        fig1, xx, yy = equatContour(self.temperature, self.radius, normed=False,
                                   levels=levels, cm=cm, cbar=cbar,tit=title)
        #fig, xx, yy = equatContour(self.vorticity, self.radius, normed=False,
        #                           levels=levels, cm=cm, cbar=cbar,tit=title)
        fig2, xx, yy = equatContour(self.vorticity, self.radius,levels=levels,vmin=0,vmax=0.01, cm=cm, cbar=cbar,tit=title)
        fig2, xx, yy = equatContour(self.vr, self.radius,levels=levels, cm=cm, cbar=cbar,tit=title)
        #fig, xx, yy = equatContour(self.vr, self.radius, 
        #                           levels=levels, cm=cm, cbar=cbar,label='t=1.8e-3',vmin=-2800,vmax=2800, tit=title)
        #fig, xx, yy = equatContour(self.vphi, self.radius, label='Azimuthal velocity',
        #                           levels=levels, cm=cm, cbar=cbar, tit=title)

    def lplot(self, title=True):

        fig2, xx, yy = lineContour(self.tanom,self.radius,self.ri,self.ro)
        plt.title('Temperature')
