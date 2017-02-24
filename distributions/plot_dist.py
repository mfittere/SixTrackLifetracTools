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

print('sig(IP1) = %4.2f um'%(sig*1.e6))
print('sigp(IP1) = %4.2f urad'%(sigp*1.e6))

# -------------- plot phase space
close('all')
for i in range(1,64):
  mp.plot_phase_space_dump(i,fn='IP1_DUMP_1')
for i in range(1,64):
  mp.plot_phase_space_dump_norm(i,fn='IP1_DUMP_1')
# save figures
for f in 'xxp zzp yyp'.split():
  figure(f)
  savefig('%s.png'%f)
  figure('n%s'%f)
  savefig('%s_norm.png'%f)
figure('axay')
savefig('axay_norm.png')

#plot initial distribution in normalized phase space
init_dist=loadtxt('distr_xy4-6_z-delta1.dat',skiprows=1)
# plot first 64 particles in xy
plot(init_dist[0:64,0],init_dist[0:64,2])
xlabel(r'$n_{x} [\sigma]$'%p)
ylabel(r'$n_{y} [\sigma]$'%p)
 
# plot all particles in xy
plot(init_dist[:,0],init_dist[:,2],'.')
xlabel(r'$n_{x} [\sigma]$'%p)
ylabel(r'$n_{y} [\sigma]$'%p)
