# -*- coding: utf-8 -*-
import glob, os
import numpy as np
from scipy.fftpack import dct, idct, fft, ifft

def chebgrid(nr, a, b):
    """
    This function defines a Gauss-Lobatto grid from a to b.

    >>> r_icb = 0.5 ; r_cmb = 1.5; n_r_max=65
    >>> rr = chebgrid(n_r_max, r_icb, r_cmb)

    :param nr: number of radial grid points
    :type nr: int
    :param a: lower limit of the Gauss-Lobatto grid
    :type a: float
    :param b: upper limit of the Gauss-Lobatto grid
    :type b: float
    :returns: the Gauss-Lobatto grid
    :rtype: numpy.ndarray
    """
    rst = (a+b)/(b-a)
    rr = 0.5*(rst+np.cos(np.pi*(1.-np.arange(nr+1.)/nr)))*(b-a)
    return rr

def matder(nr, z1, z2):
    """ 
    This function calculates the derivative in Chebyshev space.

    >>> r_icb = 0.5 ; r_cmb = 1.5; n_r_max=65
    >>> d1 = matder(n_r_max, r_icb, r_cmb)
    >>> # Chebyshev grid and data
    >>> rr = chebgrid(n_r_max, r_icb, r_cmb)
    >>> f = sin(rr)
    >>> # Radial derivative
    >>> df = dot(d1, f)

    :param nr: number of radial grid points
    :type nr: int
    :param z1: lower limit of the Gauss-Lobatto grid
    :type z1: float
    :param z2: upper limit of the Gauss-Lobatto grid
    :type z2: float
    :returns: a matrix of dimension (nr,nr) to calculate the derivatives
    :rtype: numpy.ndarray
    """
    nrp = nr+1
    w1 = np.zeros((nrp, nrp), dtype='Float64')
    zl = z2-z1
    for i in range(nrp):
        for j in range(nrp):
            w1[i, j] = spdel(i, j, nr, zl)
     
    return w1


def intcheb(f, nr, z1, z2):
    """ 
    This function integrates an input function f defined on the Gauss-Lobatto grid.

    >>> print(intcheb(f, 65, 0.5, 1.5))

    :param f: an input array
    :type: numpy.ndarray
    :param nr: number of radial grid points
    :type nr: int
    :param z1: lower limit of the Gauss-Lobatto grid
    :type z1: float
    :param z2: upper limit of the Gauss-Lobatto grid
    :type z2: float
    :returns: the integrated quantity
    :rtype: float
    """
     
    #fn = np.abs(np.real(dct(f,1))*0.5*np.sqrt(2.0/np.real(nr-1))) 
    fn = np.real(chebtransform(nr,f))
    int = 0.
    for i in range(0, nr, 2):
        if i==0 or i==nr-1:
           int = int-0.5*(z2-z1)/(i**2.-1.)*fn[i]
        else:
           int = int-(z2-z1)/(i**2.-1.)*fn[i]
    int = int * np.sqrt(2.0/np.real(nr-1.)) 
    return int


def spdel(kr, jr, nr, zl):
    if kr != nr :
        fac = 1.
        k = kr
        j = jr
    else:
        fac = -1.
        k = 0.
        j = nr-jr
     
    spdel = fac*dnum(k, j, nr)/den(k, j, nr)
    return -spdel*(2./zl)

def dnum(k, j, nr):
    if k == 0:
        if (j == 0 or j == nr):
            dnum = 0.5
            a = nr % 2
            if a == 1:
                dnum = -dnum
            if j == 0:
                dnum = 1./3.*float(nr*nr)+1./6.
            return dnum
     
        dnum = 0.5*(float(nr)+0.5)*((float(nr)+0.5)+(1./np.tan(np.pi*float(j) \
               /float(2.*nr)))**2)+1./8.-0.25/(np.sin(np.pi*float(j)/ \
               float(2*nr))**2) - 0.5*float(nr*nr)
        return dnum
     
    dnum = ff(k+j, nr)+ff(k-j, nr)
    return dnum

def ff(i, nr):
    if i == 0:
        return 0
    ff = float(nr)*0.5/np.tan(np.pi*float(i)/float(2.*nr))
     
    a = i % 2
    if a == 0:
        ff = -ff
    return ff
 
def den(k, j, nr):
    if k == 0:
        den = 0.5*float(nr)
        a = j % 2
        if a == 1:
            den = -den
        if (j == 0 or j == nr):
            den = 1.
        return den
     
    den = float(nr)*np.sin(np.pi*float(k)/float(nr))
    if (j == 0 or j == nr):
        den = 2.*den
    return den


def scanDir(pattern, tfix=None):
    """
    This function sorts the files which match a given input pattern from the oldest
    to the most recent one (in the current working directory)

    >>> dat = scanDir('log.*')
    >>> print(log)

    :param pattern: a classical regexp pattern
    :type pattern: str
    :param tfix: in case you want to add only the files that are more recent than
                 a certain date, use tfix (computer 1970 format!!)
    :type tfix: float
    :returns: a list of files that match the input pattern
    :rtype: list
    """
    dat = [(os.stat(i).st_mtime, i) for i in glob.glob(pattern)]
    dat.sort()
    if tfix is not None:
        out = []
        for i in dat:
            if i[0] > tfix:
                out.append(i[1])
    else:
        out = [i[1] for i in dat]
    return out

def spat_spec(data, nm, np):
    out = fft(data, nm)
    return out/(np)

def spec_spat2(data, n):
    out = ifft(data, n)
    return out.real

def spec_spat(data, n, axis=0):
    out = np.fft.irfft(data, axis=axis, n=n)*n
    return out.real

def chebforward(varFR,Nm_max,Nr_max):
    
    varFC = np.zeros((Nm_max+1,Nr_max),dtype=np.complex128)
    for i in range(0,Nm_max+1):
       varFC[i][:] = dct(varFR[i][:],1)

    return varFC

def chebinverse(varFC,Nm_max,Nr_max_ref,Nr_max_cur):
  
    varFR = np.zeros((Nm_max+1,Nr_max_ref),dtype=np.complex128)
    for i in range(0,Nm_max+1):
       varFR[i][:] = idct(varFC[i][:],1)/(2.0*(Nr_max_cur-1))

    return varFR


def padding(ttFC,Nm_max_ref,Nr_max_ref,Nm_max_cur,Nr_max_cur):


    var_comp = np.zeros((Nm_max_ref+1,Nr_max_ref),dtype=np.complex128)
    if (Nm_max_ref>Nm_max_cur) and (Nr_max_ref>Nr_max_cur) :
        for i in range(0,Nm_max_cur+1):     
           for j in range(0,Nr_max_cur):
                   var_comp[i,j]=ttFC[i,j]

        for i in range(Nm_max_cur+1,Nm_max_ref+1):
           for j in range(Nr_max_cur,Nr_max_ref):    
                   var_comp[i,j]=0.0

    elif (Nm_max_ref==Nm_max_cur) and (Nr_max_ref>Nr_max_cur):
        for i in range(0,Nm_max_cur+1):     
           for j in range(0,Nr_max_cur):
                   var_comp[i,j]=ttFC[i,j]

        for i in range(0,Nm_max_cur+1):     
           for j in range(Nr_max_cur,Nr_max_ref):    
                   var_comp[i,j]=0.0

    elif (Nm_max_ref>Nm_max_cur) and (Nr_max_ref==Nr_max_cur):
        for i in range(0,Nm_max_cur+1):     
           for j in range(0,Nr_max_cur):
                   var_comp[i,j]=ttFC[i,j]

        for i in range(Nm_max_cur+1,Nm_max_ref+1):
           for j in range(0,Nr_max_cur):
                   var_comp[i,j]=0.0
    elif (Nm_max_ref==Nm_max_cur) and (Nr_max_ref==Nr_max_cur):
             var_comp = ttFC   

    return var_comp

def calc_L2norm(varFR_comp,varFR_ref,Nm_max_ref,Nr_max_ref,rmin,rmax):

    var_diff_phi=np.zeros(Nr_max_ref)
    varp=np.zeros(Nr_max_ref)
    rr= chebgrid(Nr_max_ref-1,rmin,rmax)
    for j in range(0,Nr_max_ref-0):
       var_diff_phi[j]=0.0
       for i in range(1,Nm_max_ref+1):
          var_diff_phi[j] = var_diff_phi[j] + abs((varFR_comp[i,j] - varFR_ref[i,j])*(varFR_comp[i,j]-varFR_ref[i,j])) 
       var_diff_phi[j] = 2*var_diff_phi[j] + abs((varFR_comp[0,j] - varFR_ref[0,j])*(varFR_comp[0,j]-varFR_ref[0,j]))
    varp = rr*var_diff_phi
    L2_error = np.sqrt((1.0/(np.pi*(rmax**2.0-rmin**2.0)))*(2.0*np.pi*intcheb(varp, Nr_max_ref, rmin, rmax)))

    return L2_error

def calc_rel_L2norm(varFR_comp,varFR_ref,Nm_max_ref,Nr_max_ref,rmin,rmax):

    var_diff_phi=np.zeros(Nr_max_ref)
    varp=np.zeros(Nr_max_ref)
    denom=np.zeros(Nr_max_ref)
    rr= chebgrid(Nr_max_ref-1,rmin,rmax)
    for j in range(0,Nr_max_ref):
       var_diff_phi[j]=0.0
       denom[j]=0.0
       for i in range(1,Nm_max_ref+1):
          var_diff_phi[j] = var_diff_phi[j] + abs((varFR_comp[i,j] - varFR_ref[i,j])*(varFR_comp[i,j]-varFR_ref[i,j])) 
          denom[j] = denom[j] + abs((varFR_ref[i,j])*(varFR_ref[i,j]))
       var_diff_phi[j] = 2*var_diff_phi[j] + abs((varFR_comp[0,j] - varFR_ref[0,j])*(varFR_comp[0,j]-varFR_ref[0,j]))
       denom[j] = 2*denom[j] + abs((varFR_ref[0,j])*(varFR_ref[0,j]))
    varp = rr*var_diff_phi
    denom = np.real(np.multiply(rr,denom))
    rel_L2_error = np.sqrt((intcheb(varp, Nr_max_ref, rmin, rmax))/(intcheb(denom, Nr_max_ref, rmin, rmax)))
    #rel_L2_error = np.sqrt((2.0*np.pi*intcheb(var_diff_phi, Nr_max_ref, rmin, rmax))/(2.0*np.pi*intcheb(denom, Nr_max_ref, rmin, rmax)))
    #rel_L2_error = np.sqrt((1.0/(np.pi*(rmax**2.0-rmin**2.0)))*(4.0*np.pi*intcheb(varp, Nr_max_ref, rmin, rmax)))

    return rel_L2_error

def calc_combine_error(var1_comp,var1_ref,var2_comp,var2_ref,var3_comp,var3_ref,Nm_max_ref,Nr_max_ref,rmin,rmax):

    var_diff_phi1=np.zeros(Nr_max_ref)
    denom1=np.zeros(Nr_max_ref)
    var_diff_phi2=np.zeros(Nr_max_ref)
    denom2=np.zeros(Nr_max_ref)
    var_diff_phi3=np.zeros(Nr_max_ref)
    denom3=np.zeros(Nr_max_ref)
    rr= chebgrid(Nr_max_ref-1,rmin,rmax)

    for j in range(0,Nr_max_ref):
       var_diff_phi1[j]=0.0
       denom1[j]=0.0
       var_diff_phi2[j]=0.0
       denom2[j]=0.0
       var_diff_phi3[j]=0.0
       denom3[j]=0.0
       for i in range(1,Nm_max_ref+1):

          var_diff_phi1[j] = var_diff_phi1[j] + abs((var1_comp[i,j] - var1_ref[i,j])*(var1_comp[i,j]-var1_ref[i,j])) 
          denom1[j] = denom1[j] + abs((var1_ref[i,j])*(var1_ref[i,j]))

          var_diff_phi2[j] = var_diff_phi2[j] + abs((var2_comp[i,j] - var2_ref[i,j])*(var2_comp[i,j]-var2_ref[i,j])) 
          denom2[j] = denom2[j] + abs((var2_ref[i,j])*(var2_ref[i,j]))

          var_diff_phi3[j] = var_diff_phi3[j] + abs((var3_comp[i,j] - var3_ref[i,j])*(var3_comp[i,j]-var3_ref[i,j])) 
          denom3[j] = denom3[j] + abs((var3_ref[i,j])*(var3_ref[i,j]))

       var_diff_phi1[j] = 2*var_diff_phi1[j] + abs((var1_comp[0,j] - var1_ref[0,j])*(var1_comp[0,j]-var1_ref[0,j]))
       denom1[j] = 2*denom1[j] + abs((var1_ref[0,j])*(var1_ref[0,j]))

       var_diff_phi2[j] = 2*var_diff_phi2[j] + abs((var2_comp[0,j] - var2_ref[0,j])*(var2_comp[0,j]-var2_ref[0,j]))
       denom2[j] = 2*denom2[j] + abs((var2_ref[0,j])*(var2_ref[0,j]))

       var_diff_phi3[j] = 2*var_diff_phi3[j] + abs((var3_comp[0,j] - var3_ref[0,j])*(var3_comp[0,j]-var3_ref[0,j]))
       denom3[j] = 2*denom3[j] + abs((var3_ref[0,j])*(var3_ref[0,j]))

    var_diff_phi1 = np.real(np.multiply(rr,var_diff_phi1))
    denom1 = np.real(np.multiply(rr,denom1))
    n1=2.0*np.pi*intcheb(var_diff_phi1, Nr_max_ref, rmin, rmax)
    d1=2.0*np.pi*intcheb(denom1, Nr_max_ref, rmin, rmax)

    var_diff_phi2 = np.real(np.multiply(rr,var_diff_phi2))
    denom2 = np.real(np.multiply(rr,denom2))
    n2=2.0*np.pi*intcheb(var_diff_phi2, Nr_max_ref, rmin, rmax)
    d2=2.0*np.pi*intcheb(denom2, Nr_max_ref, rmin, rmax)

    var_diff_phi3 = np.real(np.multiply(rr,var_diff_phi3))
    denom3 = np.real(np.multiply(rr,denom3))
    n3=2.0*np.pi*intcheb(var_diff_phi3, Nr_max_ref, rmin, rmax)
    d3=2.0*np.pi*intcheb(denom3, Nr_max_ref, rmin, rmax)

    combine_error = np.sqrt(n1/d1 + n2/d2 + n3/d3)

    return combine_error

def calc_L2norm_two(var1FR_comp,var1FR_ref,var2FR_comp,var2FR_ref,Nm_max_ref,Nr_max_ref,rmin,rmax):

    var_diff_phi=np.zeros(Nr_max_ref)
    rr= chebgrid(Nr_max_ref-1,rmin,rmax)
    for j in range(0,Nr_max_ref):
       var_diff_phi[j]=0.0
       for i in range(0,Nm_max_ref+1):
          var_diff_phi[j] = var_diff_phi[j] + np.abs(var1FR_comp[i,j] - var1FR_ref[i,j])*np.abs(var1FR_comp[i,j]-var1FR_ref[i,j]) +  np.abs(var2FR_comp[i,j] - var2FR_ref[i,j])*np.abs(var2FR_comp[i,j]-var2FR_ref[i,j])
    var_diff_phi = np.real(np.multiply(rr,var_diff_phi))
    L2_error = np.sqrt(2.0*np.pi*intcheb(var_diff_phi, Nr_max_ref, rmin, rmax))

    return L2_error


def calc_rmsvel(var1FR_ref,var2FR_ref,Nm_max_ref,Nr_max_ref,rmin,rmax): # From output variables from code in Fourier-real space

    var1=np.zeros(Nr_max_ref)
    var2=np.zeros(Nr_max_ref)
    rr= chebgrid(Nr_max_ref-1,rmin,rmax)
    for j in range(0,Nr_max_ref):
       var1[j]=0.0
       var2[j]=0.0
       for i in range(0,Nm_max_ref+1):
          var1[j] = var1[j] + np.abs(var1FR_ref[i,j]*var1FR_ref[i,j]) 
          var2[j] = var2[j] + np.abs(var2FR_ref[i,j]*var2FR_ref[i,j]) 
    var1 = np.real(np.multiply(rr,var1))
    var2 = np.real(np.multiply(rr,var2))
    Ekin = (2.0*np.pi*intcheb(var1, Nr_max_ref, rmin, rmax))+(2.0*np.pi*intcheb(var2, Nr_max_ref, rmin, rmax))
    rmsvel = np.sqrt(2.0*Ekin/(np.pi*(rmax**2.0-rmin**2.0)))

    return rmsvel

def calc_rmsvel2(var1_r,var2_r,Nm_max,Nr_max,rmin,rmax): # From variables in physical space and taking Fourier transform

    Np_max=3*Nm_max
    var1_FR=np.zeros((Np_max,Nr_max),dtype=complex) 
    var2_FR=np.zeros((Np_max,Nr_max),dtype=complex) 
    var1FR=np.zeros((Nm_max+1,Nr_max),dtype=complex) 
    var2FR=np.zeros((Nm_max+1,Nr_max),dtype=complex) 

    for i in range(0,Nr_max):
       var1_FR[:,i]=fft(var1_r[:,i])/Np_max 
       var2_FR[:,i]=fft(var2_r[:,i])/Np_max 
     
    for i in range(0,Nm_max+1):
       for j in range(0,Nr_max):
          var1FR[i,j]=var1_FR[i,j]
          var2FR[i,j]=var2_FR[i,j]

    var1=np.zeros(Nr_max)
    var2=np.zeros(Nr_max)
    rr= chebgrid(Nr_max-1,rmin,rmax)
     
    for j in range(0,Nr_max):
       var1[j]=0.0
       var2[j]=0.0
       #for i in range(0,Nm_max+1):
       for i in range(0,Np_max):
          var1[j] = var1[j] + np.abs(var1_FR[i,j]*var1_FR[i,j]) 
          var2[j] = var2[j] + np.abs(var2_FR[i,j]*var2_FR[i,j]) 
          #var1[j] = var1[j] + np.abs(var1FR[i,j]*var1FR[i,j]) 
          #var2[j] = var2[j] + np.abs(var2FR[i,j]*var2FR[i,j]) 
    var1 = np.real(np.multiply(rr,var1))
    var2 = np.real(np.multiply(rr,var2))
    Ekin = (0.5*((2.0*np.pi*intcheb(var1, Nr_max, rmin, rmax))+(2.0*np.pi*intcheb(var2, Nr_max, rmin, rmax))))
    rmsvel = np.sqrt(2.0*Ekin/(np.pi*(rmax**2.0-rmin**2.0)))

    return rmsvel, Ekin

def calc_rmsvel3(var1_r,var2_r,Nm_max,Nr_max,rmin,rmax): # From variables in physical space using np.trapz 

    Np_max=3*Nm_max

    var1=np.zeros(Nr_max)
    var2=np.zeros(Nr_max)
    rr= chebgrid(Nr_max-1,rmin,rmax)
    phi = np.zeros(Np_max)
    for i in range(1,Np_max):
       phi[i]=phi[i-1]+2.0*np.pi/(Np_max)
 

    for j in range(0,Nr_max):
       var1[j] = np.trapz(var1_r[:,j]*var1_r[:,j],phi) 
       var2[j] = np.trapz(var2_r[:,j]*var2_r[:,j],phi) 
    var1 = np.real(np.multiply(rr,var1))
    var2 = np.real(np.multiply(rr,var2))
    Ekin = 0.5*((intcheb(var1, Nr_max, rmin, rmax))+(intcheb(var2, Nr_max, rmin, rmax)))
    rmsvel = np.sqrt(2.0*Ekin/(np.pi*(rmax**2.0-rmin**2.0)))

    return rmsvel, Ekin

def calc_rmsvel_phi(var1FR_ref,var2FR_ref,Nm_max_ref,Nr_max_ref,rmin,rmax):

    var=np.zeros(Nr_max_ref)
    rr= chebgrid(Nr_max_ref-1,rmin,rmax)
    for j in range(0,Nr_max_ref):
       var[j]=0.0
       for i in range(0,Nm_max_ref+1):
          var[j] = var[j] + np.abs(var1FR_ref[i,j])*np.abs(var1FR_ref[i,j]) +  np.abs(var2FR_ref[i,j])*np.abs(var2FR_ref[i,j]) 
       var[j] = np.sqrt(var[j]/(2.0*np.pi*rr[j]))

    return var

def time_avg(var_ref,time,Nr_max_ref,nsnaps):

    var=np.zeros(Nr_max_ref,dtype=complex)
    fac = time[nsnaps-1] - time[0]
    for j in range(0,Nr_max_ref):
       var[j]=0.0
       var[j] = 1./fac * np.trapz(var_ref[:,j],time)   

    return var

def calc_maxnorm(var_comp,var_ref,Nm_max_ref,Nr_max_ref):
       
    #var_c = spec_spat(var_comp, 3*Nm_max_ref)
    #var_r = spec_spat(var_ref, 3*Nm_max_ref)
    #diff = np.zeros(shape=(Nm_max_ref+1,Nr_max_ref))
    #for j in range(1,Nr_max_ref-1):
    #   for i in range(0,Nm_max_ref+1):
       #for i in range(0,1):
    #      diff[i,j] = abs(np.imag(var_comp[i,j]-var_ref[i,j])) 
    #max_error = (diff).max()
    max_error = abs(var_comp - var_ref).max()
    #max_error = abs(var_c - var_r).max()
    #max_error = abs(np.real(var_comp) - np.real(var_ref)).max()
    #max_error = abs(np.imag(var_comp) - np.imag(var_ref)).max()

    return max_error

def calc_maxnorm_omg(var_comp,var_ref,Nm_max_ref,Nr_max_ref):
       
    diff = np.zeros(shape=(Nm_max_ref+1,Nr_max_ref))
    for j in range(1,Nr_max_ref-1):
       for i in range(0,Nm_max_ref+1):
          diff[i,j] = abs(np.imag(var_comp[i,j]-var_ref[i,j])) 
    max_error = diff.max()

def chebtransform(Nr_max,f):

    f00=np.zeros([2*Nr_max-2],dtype=complex)
    ff=np.zeros([2*Nr_max-2],dtype=complex)
    f2=np.zeros([Nr_max],dtype=complex)
    ft=np.zeros([Nr_max],dtype=complex)

    f0 = f[::-1]
    # Pre-processing
    f00[0:Nr_max]=f0[0:Nr_max]
    f00[Nr_max:2*Nr_max-2]=f[1:Nr_max-1]
    
    ff[:] = fft(f00) # Execute dft

    # Post-processing 
    ff=ff/(2*Nr_max-2)
    
    f2[0]=ff[0]
    f2[Nr_max-1]=ff[Nr_max-1]
    f2[1:Nr_max-1]=2*ff[1:Nr_max-1]
    
    fac=np.sqrt(2./(Nr_max-1))
    
    ft=f2/fac
    ft[0]=2*ft[0]
    ft[Nr_max-1]=2*ft[Nr_max-1]
    
    return ft
     
def chebinvtran(Nr_max,fc):
     
    fin=np.zeros([Nr_max],dtype=complex)
    ff=np.zeros([Nr_max],dtype=complex)
    f2=np.zeros([Nr_max],dtype=complex)
    ft=np.zeros([Nr_max],dtype=complex)

    fin = idct(fc,1)
    fac = np.sqrt(2./(Nr_max-1))
    fin = fac * fin * 0.5

    fin = fin[::-1]
    return fin

def chebinvtranD1(Nr_max,ft):

    f2r=np.zeros([2*Nr_max-2],dtype=complex)
    f2c=np.zeros([2*Nr_max-2],dtype=complex)
    df=np.zeros([Nr_max],dtype=complex)
    beta1=np.zeros([Nr_max],dtype=complex)

    f=ft
    fac = np.sqrt(2./(Nr_max-1))

    # Recurrence for the 1st derivative coefficients
    beta1[Nr_max-1] = 0.0
    beta1[Nr_max-2] = 2. * (Nr_max-1) * f[Nr_max-1]
    for i in range(Nr_max-2,0,-1):
       beta1[i-1] = beta1[i+1] + 4. * (i) * f[i] 

    beta1=beta1*fac
    beta1[0]=beta1[0]/2
    beta1[Nr_max-1]=beta1[Nr_max-1]/2

    # Fast inverse Chebyshev transform
    # Pre-processing
    f2c[0]=beta1[0]
    f2c[1:Nr_max-1]=beta1[1:Nr_max-1]/2
    f2c[Nr_max:2*Nr_max-2]=beta1[Nr_max-2:0:-1]/2

    f2r = fft(f2c)

    df[0:Nr_max]=f2r[0:Nr_max]
    df=df[::-1] 

    return df 

def chebinvtranD2(Nr_max,ft):

    f2r=np.zeros([2*Nr_max-2],dtype=complex)
    f2c=np.zeros([2*Nr_max-2],dtype=complex)
    d2f=np.zeros([Nr_max],dtype=complex)
    beta1=np.zeros([Nr_max],dtype=complex)
    beta2=np.zeros([Nr_max],dtype=complex)

    f=ft
    fac = np.sqrt(2./(Nr_max-1))

    # Recurrence for the 1st derivative coefficients
    beta1[Nr_max-1] = 0.0
    beta1[Nr_max-2] = 2. * (Nr_max-1) * f[Nr_max-1]
    for i in range(Nr_max-2,0,-1):
       beta1[i-1] = beta1[i+1] + 4. * (i) * f[i] 

    # Recurrence for the 2nd derivative coefficients
    beta2[Nr_max-1] = 0.0
    beta2[Nr_max-2] = 0.0 
    for i in range(Nr_max-2,0,-1):
       beta2[i-1] = beta2[i+1] + 4. * (i) * beta1[i] 

    beta1=beta1*fac
    beta1[0]=beta1[0]/2
    beta1[Nr_max-1]=beta1[Nr_max-1]/2

    beta2=beta2*fac
    beta2[0]=beta2[0]/2
    beta2[Nr_max-1]=beta2[Nr_max-1]/2

    # Fast inverse Chebyshev transform
    # Pre-processing
    f2c[0]=beta2[0]
    f2c[1:Nr_max-1]=beta2[1:Nr_max-1]/2
    f2c[Nr_max:2*Nr_max-2]=beta2[Nr_max-2:0:-1]/2

    f2r = fft(f2c) # Execute dft 

    d2f[0:Nr_max]=f2r[0:Nr_max]
    d2f=d2f[::-1] 

    return d2f 

   # NOTES
   # Note that when we do a integral in phi direction we get 2*pi as a prefactor due to Parseval's theorem 
