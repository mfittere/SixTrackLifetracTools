from read_fortbin import *
import numpy as np
from matplotlib.pyplot import *

def _read_dump(fn):
  ftype=zip('ID turn s x xp y yp z dE/E ktrack'.split(),list((int,)*2+(float,)*8+(int,)))
  return np.loadtxt(fn,dtype=ftype)

class MacroParticle():
  """ converts liftrac input distribution in normalized 
  coordinates into SixTrack or collimation SixTrack 
  coordinates
  emitx,emity,emitz: normalized emittances in SI units
  fndist: distribution in normalized coordinates
  fn: fort.90 file with normalization matrix"""
  def __init__(self, (emitx, emity, emitz),fndist,
             ncore=10000,E0=7000,fn='fort.90'):
    self.fndist=fndist
    self.normdist = np.loadtxt(self.fndist,skiprows=1)# nx,npx,ny,npy,nz,npz,weight,particle_id
    self.current_dir = os.getcwd()
    self.E0 = float(E0)  # energy in GeV
    self.mp = 0.938272046  # mass in GeV
    gamma_rel = self.E0 / self.mp
    self.ncore = int(ncore)
    self.ex = float(emitx) / gamma_rel  # normalized emittance
    self.ey = float(emity) / gamma_rel
    self.ez = float(emitz) # bunchlength*energy spread
    self.header, self.part = read_fortbin(open(os.path.join(self.current_dir, fn),'rb'))
    self.co = self.header['closorb'] # closed orbit in mm,mrad,mm,mrad,mm,1
    self.T = self.Tmatrix() # m,rad,m,rad,m,1
    self.Tinv = np.linalg.inv(self.T) # m,rad,m,rad,m,1
    self.S = self.Smatrix()
    self.Sinv = np.linalg.inv(self.S)
    # print 'REMINDER: emittances in SI * pi && dE*M'
  def Tmatrix(self):
    """normalization matrix M=T R inv(T) in Si units
    T converts from normalized to physical coordinates"""
    R = np.array([list(self.header['rmatrix'][0:6]), \
                  list(self.header['rmatrix'][6:12]), \
                  list(self.header['rmatrix'][12:18]), \
                  list(self.header['rmatrix'][18:24]), \
                  list(self.header['rmatrix'][24:30]), \
                  list(self.header['rmatrix'][30:36])])
    for i in range(6):
        # conversion to SI
        R[i][5] = R[i][5] * 1.0E-3
        R[5][i] = R[5][i] * 1.0E+3
    return R
  def Smatrix(self):
    """matrix of sqrt(eps*)"""
    return np.array([[np.sqrt(self.ex), 0, 0, 0, 0, 0], \
                     [0, np.sqrt(self.ex), 0, 0, 0, 0], \
                     [0, 0, np.sqrt(self.ey), 0, 0, 0], \
                     [0, 0, 0, np.sqrt(self.ey), 0, 0], \
                     [0, 0, 0, 0, np.sqrt(self.ez), 0], \
                     [0, 0, 0, 0, 0, np.sqrt(self.ez)]])
  def toSI(self,verbose=False):
    """return distribution in physical coordinates
       x,x',y,y',z,dp/p in SI units"""
    dist = np.array([ np.dot(self.T,np.dot(self.S,n[:6])) for n in self.normdist ]) # x = self.T*self.S*nx(fndist) SI units x[m] xp[rad] y[m] yp[rad] sig[m] dE/E [1]
    x,px,y,py,sigma,delta=dist.T
    # 0) get the closed orbit
    x0,xp0,y0,yp0,sigma0,delta0=self.co # x[mm] xp[mrad] y[mm] yp[mrad] sig[mm] dE/E [1]
    #    convert to SI units x[m] xp[rad] y[m] yp[rad] sig[m] dE/E [1]
    for var in x0,xp0,y0,yp0,sigma0:
      var=var*1.e-3
    if verbose: print 'closed orbit = ',x0,xp0,y0,yp0,sigma0,delta0
    # 1) add the momentum
    # convert to canonical variables
    xp=px*(1+delta+delta0);yp=py*(1+delta+delta0);
    # 4D closed orbit added automatically by SixTrack
    np.savetxt('%s_SI'%self.fndist,np.array([x,xp,y,yp,sigma,delta]).T,'%4.8e')
    return np.array([x,xp,y,yp,sigma,delta])
  def toCol(self):
    # get amplitudes in SI units x[m] xp[rad] y[m] yp[rad] sig[m] dE/E [1]
    [x,xp,y,yp,sigma,delta] = self.toSI()
    # convert to collimation units x[m] xp[rad] y[m] yp[rad] sig[mm] dE [MeV]
    sigma = sigma*1.e3
    delta = self.E0*1.e3*(1+delta)
    np.savetxt('%s_Col'%self.fndist,np.array([x,xp,y,yp,sigma,delta]).T,'%4.8e')
    return np.array([x,xp,y,yp,sigma,delta])
  def plot_phase_space_fort(self,partid,basedir='.'):
    """plots the phase space of the tracked
    particles saved in fort.* files"""
    head,part=read_allfortdir(basedir=basedir)
    idxs = part['ID']==partid
    for p in 'x y'.split():
      figure('%s%sp'%(p,p))
      plot(part[idxs][p],part[idxs]['%sp'%p],'.')
      xlabel('%s [mm]'%p)
      ylabel('%sp [mrad]'%p)
    figure('zdE/E')
    plot(part[idxs]['z'],part[idxs]['dE/E'],'.')
    xlabel('z [mm]')
    ylabel('dE/E [1]')
  def _six_to_norm(self,x,xp,y,yp,z,delta):
    # subtract closed orbit
    x,xp,y,yp,z,delta = np.array([x,xp,y,yp,z,delta]) - self.co
    # convert to canonical momenta
    delta0 = self.co[5]
    px = xp/(1+delta+delta0);py=yp/(1+delta+delta0)
    # convert to SI units
    x,px,y,py,z = np.array([x,px,y,py,z])*1.e-3
    # normalize
    return tuple(np.dot(self.Sinv,np.dot(self.Tinv,np.array([x,px,y,py,z,delta]))))
  def _dump_norm(self,fn='IP1_DUMP_1'):
    """return coordinates in dump file
    in normalized coordinates using self.T
    and self.S"""
    part = _read_dump(fn)
    idx,turn,s,x,xp,y,yp,z,delta,ktrack = [ part[var] for var in 'ID turn s x xp y yp z dE/E ktrack'.split() ]
    nv = [ (idx[i],turn[i],s[i],)+self._six_to_norm(x[i],xp[i],y[i],yp[i],z[i],delta[i])+(ktrack[i],) for i in range(len(x)) ]
    ftype = zip('ID turn s nx npx ny npy nz npz ktrack'.split(),list((int,)*2+(float,)*8+(int,)))
    return np.array(nv,dtype=ftype)
  def plot_phase_space_dump(self,partid,fn='IP1_DUMP_1',closed_orbit = True):
    """plots phase space of particle saved
    in dumpfile fn generated with DUMP block
    in SixTrack
    Parameter:
    ----------
    closed_orbit: True = closed orbit subtracted"""
    part = _read_dump(fn)
    idxs = part['ID']==partid
    x0,xp0,y0,yp0,z0,delta0 = self.co
    # x-px, y-py
    for p,c in zip('x y'.split(),[(x0,xp0),(y0,yp0)]):
      figure('%s%sp'%(p,p))
      if closed_orbit: plot(part[idxs][p]-c[0],part[idxs]['%sp'%p]-c[1],'.')
      else: plot(part[idxs][p],part[idxs]['%sp'%p],'.')
      xlabel('%s [mm]'%p)
      ylabel('%sp [mrad]'%p)
    # z-pz
    figure('zzp')
    plot(part[idxs]['z'],part[idxs]['dE/E'],'.')
    xlabel('z [mm]')
    ylabel('dE/E [1]')
  def plot_phase_space_dump_norm(self,partid,fn='IP1_DUMP_1'):
    """plots normalized phase space of particle 
    saved in dumpfile fn generated with DUMP block
    in SixTrack. Coordinates are read from fn and
    then normalized with self.T and self.S"""
    part = self._dump_norm(fn)
    idxs = part['ID']==partid
    for p in 'x y z'.split():
      figure('n%s%sp'%(p,p))
      plot(part[idxs]['n%s'%p],part[idxs]['np%s'%p],'.')
      xlabel(r'$n_{%s} [\sigma]$'%p)
      ylabel(r'$n_{p_%s} [\sigma]$'%p)
    # x-y
    figure('axay')
    plot(np.sqrt(part[idxs]['nx']**2+part[idxs]['npx']**2),np.sqrt(part[idxs]['ny']**2+part[idxs]['npy']**2),'.')
    xlabel(r'$\sqrt{n_{x}^2+n_{p_x}^2} \ [\sigma]$')
    ylabel(r'$\sqrt{n_{y}^2+n_{p_y}^2} \ [\sigma]$')

## Usage:
#fndist='distr_xy4-6_z-gauss.dat'
#
## beam parameters
#energy=7000 # TeV
#epsx=2.5 # normalized emittance [um]
#epsy=2.5 # normalized emittance [um]
#sigs=7.5e-2 # bunch length [m]
#dpp=1.1e-4 # energy spread [1]
#epsz=sigs*dpp

#mp = MacroParticle((3.75E-6, 3.75E-6, 8.305E-6),  #emittances in m
#                fndist='distr_xy4-6_z-delta1.dat',
#                ncore=10000,  #number of particles
#                E0=7000,  #Energy in GeV
#                fn='fort.90') #default file for T-matrix extraction
#
#sixcol = mp.sixcol()
