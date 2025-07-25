start.dat contains initial conditions in VMEC coordinates

s, th, ph, vnorm, vparnorm

v0=sqrt(2.d0*E_alpha*ev/(n_d*p_mass))

E_alpha = 3.5e6

double precision :: n_d=4
double precision, parameter  :: pi=3.14159265358979d0
double precision, parameter  :: twopi=6.28318530717958d0
double precision, parameter  :: c=2.9979d10
double precision, parameter  :: e_charge=4.8032d-10
double precision, parameter  :: e_mass=9.1094d-28
double precision, parameter  :: p_mass=1.6726d-24
double precision, parameter  :: ev=1.6022d-12
double precision, parameter  :: sqrt2=dsqrt(2.0d0)
