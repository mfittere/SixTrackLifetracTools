import os
import struct
import math
import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt


def write_ascii(lines, fn, status='w'):
    """ write list of lines to ascii file """
    outfile = open(fn, str(status))
    for line in lines:
        line = str(line)
        if line[-1:] != '\n': line = line + '\n'
        outfile.write(line)
    outfile.close()

    
class MacroParticle():
    """ macro-particle: multivariate gaussian-distribution """
    def __init__(self, (emitx, emity, emitz),
                 ncore=10000, 
                 E0=7000, mp=0.938272,
                 fn='fort.90'):
        self.current_dir = os.getcwd()
        self.E0 = float(E0)  # energy in GeV
        self.mp = float(mp)  # mass in GeV
        gamma_rel = self.E0 / self.mp
        self.ncore = int(ncore)
        self.ex = float(emitx) / gamma_rel  # normalized emittance
        self.ey = float(emity) / gamma_rel
        self.ez = float(emitz)
        self.header, self.part = self.read_fortbin(os.path.join(self.current_dir, fn))
        self.T = self.Tmatrix()
        self.SI = self.Smatrix()
        self.SI10 = self.Smatrix10()
        # print 'REMINDER: emittances in SI * pi && dE*M'

    def read_fortbin(self, fn):
        # Function provided by Riccaro De Maria (CERN)
        fmt_head = """\
        head1      1I Fortran header
        title     80s General title of the run
        title2    80s Additional title
        date       8s Date
        time       8s Time
        progname   8s Program name
        partfirst  1I First particle in the file
        partlast   1I Last particle in the file
        parttot    1I Total number of particles
        spacecode  1I Code for dimensionality of phase space (1,2,4 are hor., vert. and longitudinal respectively)
        turnproj   1I Projected number of turns
        qx         1d Horizontal Tune
        qy         1d Vertical Tune
        qs         1d Longitudinal Tune
        closorb    6d Closed Orbit vector
        dispvec    6d Dispersion vector
        rmatrix   36d Six-dimensional transfer map
        mess1     50d 50 additional parameter
        mess2      1I ...\
        """
        #"""
        #seedmax    1d Maximum number of different seeds
        #seednum    1d Actual seed number
        #seedstart  1d Starting value of the seed
        #turnrev    1d Number of turns in the reverse direction (IBM only)
        #lyapcor1   1d Correction-factor for the Lyapunov (sigma=s - v0 t)
        #lyapcor2   1d Correction-factor for the Lyapunov (DeltaP/P0)
        #turnrip    1d Start turn number for ripple prolongation\
        #"""
        fmt_part = """\
        partnum  1I Particle number
        partdist 1d Angular distance in phase space
        x        1d x (mm)
        xp       1d x'(mrad)
        y        1d y (mm)
        yp       1d y'(mrad)
        sig      1d Path-length sigma=s - v0 t
        delta    1d DeltaP/P0
        energy   1d Energy (Mev)\
        """
        def _read(fh, fmt):
            """ read the binary tracking file and returns header -
            see above """
            out = {}
            for line in fmt.splitlines():
                lbl, spec, desc = line.split(None, 2)
                data = fh.read(struct.calcsize(spec))
                obj = struct.unpack(spec, data)
                if len(obj) == 1:
                    obj = obj[0]
                out[lbl] = obj
            return out

        #read the binary tracking file and returns:
        #header - (see above) + particle data
        fh = open(fn, 'rb')
        header = _read(fh, fmt_head)
        partfirst = header['partfirst']
        partlast = header['partlast']
        part = {}
        for i in range(partfirst, partlast + 1):
            part[i] = []
            st = ''
        while fh.read(4) != '':
            turnnum = struct.unpack('I', fh.read(4))
            for i in range(partfirst, partlast + 1):
                pnum1 = struct.unpack('I', fh.read(4))
                orb1 = struct.unpack('8d', fh.read(64))
                part[i].append(list(orb1))
            fh.read(4)
        for i in range(partfirst, partlast + 1):
            part[i] = np.array(part[i])
        fh.close()
        return header, part

    def Tmatrix(self):
        # one turn transport matrix
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
        # sigma matrix
        RT = np.transpose(self.T)
        RI = LA.inv(self.T)
        Sdiag = np.array([[1 / self.ex, 0, 0, 0, 0, 0], \
                          [0, 1 / self.ex, 0, 0, 0, 0], \
                          [0, 0, 1 / self.ey, 0, 0, 0], \
                          [0, 0, 0, 1 / self.ey, 0, 0], \
                          [0, 0, 0, 0, 1 / self.ez, 0], \
                          [0, 0, 0, 0, 0, 1 / self.ez]])
        S = np.dot(LA.inv(RT), np.dot(Sdiag, RI))
        SI = LA.inv(S)
        return SI

    def Smatrix10(self):
        # sigma matrix
        RT = np.transpose(self.T)
        RI = LA.inv(self.T)
        k = 10.0 # see SI matrix below
        Sdiag = np.array([[1 / self.ex / k, 0, 0, 0, 0, 0], \
                          [0, 1 / self.ex / k, 0, 0, 0, 0], \
                          [0, 0, 1 / self.ey / k, 0, 0, 0], \
                          [0, 0, 0, 1 / self.ey / k, 0, 0], \
                          [0, 0, 0, 0, 1 / self.ez / 1., 0], \
                          [0, 0, 0, 0, 0, 1 / self.ez / 1.]])
        S = np.dot(LA.inv(RT), np.dot(Sdiag, RI))
        SI10 = LA.inv(S)
        return SI10

    def GaussianCore(self):
        # multivariate gaussian distribution (coordinates are in SI)
        SI = self.SI
        E0 = self.E0
        mp = self.mp
        cov = (list(SI[0]), list(SI[1]), list(SI[2]), list(SI[3]), list(SI[4]), list(SI[5]))
        cov = list(cov)
        mu = list(np.array(self.header['closorb']) * 1.0E-3)  # in SI
        self.x, self.px, self.y, self.py, self.z, self.delta = \
            np.random.multivariate_normal(mu, cov, self.ncore).T
        # dE -> dP
        p0 = np.sqrt((E0 - mp) * (E0 + mp))
        self.E = (1.0 + self.delta) * E0
        p = np.sqrt((self.E - mp) * (self.E + mp))
        self.dPP = (p - p0) / p0
        self.weight = [1.0]*self.ncore
        #plt.hist(self.px, 100, color='blue', alpha=0.5)
        return self.x, self.px, self.y, self.py, \
               self.z, self.dPP, self.delta, self.E, self.weight

    def Normalized(self, X):
        TI = LA.inv(self.T)
        vect = np.array([X[0], X[1], \
                         X[2], X[3], \
                         X[4], X[5]])
        Nvect = np.dot(TI, vect)
        u1 = Nvect[0] / (self.ex * self.ey)**(0.25)
        u2 = Nvect[1] / (self.ex * self.ey)**(0.25)
        u3 = Nvect[2] / (self.ex * self.ey)**(0.25)
        u4 = Nvect[3] / (self.ex * self.ey)**(0.25)
        u5 = Nvect[4] / (self.ez)**(0.5)
        u6 = Nvect[5] / (self.ez)**(0.5)
        return u1, u2, u3, u4, u5, u6

    def WriteDistribution(self, outfile, *args):
        data = []
        for arg in args:
            for i in range(len(arg[0])):
                st = ''
                for j in range(len(arg)):
                    st = st + str(arg[j][i]) + ' '
                data.append(st)
        write_ascii(data, outfile)


# USAGE:
Bunch = MacroParticle((3.75E-6, 3.75E-6, 8.305E-6),  #emittances in um
                ncore=10000,  #number of particles
                E0=7000,  #Energy in GeV
                mp=0.938272,  #mass of proton in GeV
                fn='fort.90') #default file for T-matrix extraction

Particles_coords = Bunch.GaussianCore() #Physical coordinates in SI: x,px,y,py,z,dP/P, dE/E, E, stat.weight = 1

x = Particles_coords[0]
px = Particles_coords[1]
y = Particles_coords[2]
py = Particles_coords[3]
z = Particles_coords[4]
delta = Particles_coords[6]

Normalized_coords = Bunch.Normalized((x,px,y,py,z,delta)) #in sigmas: ux, upx, uy, upy, uz, udE

Bunch.WriteDistribution('distribution.txt', Particles_coords) #write distribution in ASCII file

#plt.hist(Particles_coords[0], 100)
plt.hist(Normalized_coords[0], 100)
plt.show()



