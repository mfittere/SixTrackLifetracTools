from distribution_generator import *
from numpy import *

fn = 'distr_x4_y5_delta1.dat'
eps=3.5e-6
sigs=0.0755
dpp=1.129E-4
epsz=sigs*dpp

mp = MacroParticle((eps, eps, epsz),  #emittances in m
                fndist=fn,
                ncore=3,  #number of particles
                E0=7000,  #Energy in GeV
                fn='fort.90_da') #default file for T-matrix extraction

bet_ip1=6.0
sig=sqrt(bet_ip1*mp.ex)
sigp=sqrt(1/bet_ip1*mp.ex)

print 'sig(IP1) = %4.2f um'%(sig*1.e6)
print 'sigp(IP1) = %4.2f urad'%(sigp*1.e6)
# -------------- generate input distributions
varsi = mp.toSI()
varcol = mp.toCol()

# check for 'distr_x4_y5_delta1.dat' - note index starts here at 1
# 4.000000 0.000000 5.000000 0.000000 0.000000 1.000000 1 59
# 0.000000 4.000000 0.000000 5.000000 1.000000 0.000000 1 64
fmt='%s '+'%2.2f '*2
print(fmt%(('x/sig =',)+tuple([varsi[0,58]/sig,varsi[0,63]/sig])))
print(fmt%(('xp/sig =',)+tuple([varsi[1,58]/sigp,varsi[1,63]/sigp])))
print(fmt%(('y/sig =',)+tuple([varsi[2,58]/sig,varsi[2,63]/sig])))
print(fmt%(('yp/sig =',)+tuple([varsi[3,58]/sigp,varsi[3,63]/sigp])))
print(fmt%(('sigm/sigs =',)+tuple([varsi[4,58]/sigs,varsi[4,63]/sigs])))
print(fmt%(('delta/dpp =',)+tuple([varsi[5,58]/dpp,varsi[5,63]/dpp])))

