import numpy as np

EPS = np.radians(23+26./60+21.406/3600) # axial tilt J2000
SEPS = np.sin(EPS); CEPS = np.cos(EPS)
MET_MJD_REF = 51910. + 7.428703703703703e-4  # MJD of Fermi MET = 0

def nancay_zero(prof_name):
    # check for a fiducial_point in comment
    for line in file(prof_name):
        if line.strip()[0] != '#': continue
        if 'fiducial_point' in line:
            print ('Found fiducial phase in comments.')
            return float(line.split('=')[-1].split()[0])
    print ('No fiducial phase found; using Fourier transform.')
    return -2

def prof_offset(prof_name,average_level=0):
    """ Compute the zero of phase of a radio profile by determining the 
        position of the fundamental peak."""
    TWOPI = 2*np.pi
    vals = np.asarray([x.strip().split()[1] for x in file(prof_name) if x.strip()[0] != '#']).astype(float)
    # compute zero phase before any averaging
    zero_ph = nancay_zero(prof_name)
    if zero_ph == -2:
        ph = np.linspace(0,TWOPI,len(vals)+1)[:-1] # LEFT bin edges
        a1 = (np.sin(ph)*vals).sum()
        a2 = (np.cos(ph)*vals).sum()
        zero_ph = np.arctan2(a1,a2)/TWOPI
    # average if desired
    for i in range(average_level):
        if len(vals)%2>0: break
        vals = (vals[::2]+vals[1::2])/2
    ph = np.linspace(0,TWOPI,len(vals)+1)[:-1] # LEFT bin edges
    return zero_ph,ph/TWOPI,vals

def cel2ecl(skydir,inverse=False):
    """ If inverse==True, assume SkyDir is in ecliptic and perform
        transformation to celestial."""
    # this implementation could be rather more efficient...
    sign = 1-2*inverse
    mat = np.asarray([[1.,0,0],[0,CEPS,sign*SEPS],[0,-1*(sign)*SEPS,CEPS]])
    ra = np.radians(skydir.ra())
    de = np.radians(90-skydir.dec()) # convert from lat to co-lat
    z_cel = np.cos(de)
    sd = np.sin(de)
    x_cel = sd*np.cos(ra)
    y_cel = sd*np.sin(ra)
    vec_ecl = np.dot(mat,np.asarray([x_cel,y_cel,z_cel]))
    bet = np.arccos(vec_ecl[2])
    lam = np.arctan2(vec_ecl[1],vec_ecl[0])
    if lam < 0: lam += 2*np.pi 
    return np.degrees(lam),90-np.degrees(bet)

def ecl2cel(skydir):
    """ NB -- SkyDir celestial coordinates are treated as ecliptic."""
    return cel2ecl(skydir,inverse=True)
    
def mjd2met(mjd):
    """ Convert MJD (TT) to MET (TT)."""
    return (np.asarray(mjd)-MET_MJD_REF)*86400

def met2mjd(met):
    """ Convert MET (TT) to MJD (TT)."""
    return np.asarray(met,dtype=float)/86400+MET_MJD_REF
   

