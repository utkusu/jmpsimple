program try
use global
use randomgen
use emax
USE IFPORT ! for intel fortran only
implicit none


real(dble) coef(Gsize+1)
real(dble) coefoc(Gsizeoc+1)
real(dble) mutype(3)
real(dble) parA(12),parU(7),parW(7),parH(6),beta,sigma(shocksize1,shocksize1),parFV(Gsize+1),parFVv(Gsizeoc+1)
real(dble) parFVmat(Gsize+1,nctype+1)
real(dble) omega1(3),omega2(4),omega3(4),eps(Nmc,shocksize1)
real(dble) ftype(5)
real(dble) solw(Gsize,nperiods-deltamin+2,deltamax-deltamin+2)
real(dble) wcoefficients(Gsize+1,nperiods-deltamin+2,deltamax-deltamin+2,2)
real(dble) solv(Gsizeoc,nperiods)
real(dble) parB(Bsizeexo+1)
real(dble) typevec(1)
real(dble) typeprob(1)
integer k
omega1=(/4.d0,3.0d0,0.0d0/)
omega2=(/2.0d0,1.0d0,20.0d0,20.0d0/)
omega3=(/10.0d0,1.0d0,18.0d0,10.0d0/)

parA=1.0d0*0.001
parU=1.0d0*0.001
parW=0.1d0*0.00010
parH=0.1d0*0.00010
Sigma=0.5d0*0.10
beta=0.9d0
parFV=1.0d0
parFVv=1.0d0
parFVmat=1.0d0
mutype=0.010d0
eps=0.0d0
ftype=(/0.d0,999.999d0, 0.0d0,1.5d0,1.0d0/)
typevec=0.5d0
typeprob=1.0d0
! try the coefochcb
parB=0.8d0
print*, emaxlate(omega1,omega2,omega3,mutype,eps,parA,parU,parW,parH,beta,parFV)
print*, emaxhc(omega1,omega2,omega3,mutype,eps,parA,parU,parW,parH,beta,parFV,1.0d0)
print*, emaxhcx(omega1,omega2,omega3,mutype,eps,parA,parU,parW,parH,beta,parFV,1.0d0)
print*, emaxoclate(omega1,omega2,omega3,mutype,eps,parA,parU,parW,parH,beta,parFVv)
print*, emaxochc(omega1,omega2,omega3,mutype,eps,parA,parU,parW,parH,beta,parFVv,1.0d0)
print*, emaxochcb(omega1,omega2,omega3,mutype,eps,parA,parU,parW,parH,beta,parFVmat,parFVv,parB,typevec,typeprob,1.0d0)

call coeffinal(coef, 3.0d0, 2.0d0, mutype, parA, parU,parW,parH,beta,sigma)


print*, coef
print*, '--------------------------------------------'
call coeflate(coef, 3.0d0, 2.0d0, mutype, parA, parU,parW,parH,beta,sigma,parFV)
print*, coef
print*, '--------------------------------------------'

call coefhc(coef, 3.0d0, 2.0d0, mutype, parA, parU,parW,parH,beta,sigma,parFV,1.0d0)

call coefhcx(coef, 3.0d0, 2.0d0, mutype, parA, parU,parW,parH,beta,sigma,parFV,1.0d0)
print*, coef
print*, '--------------------------------------------'

call coefocfinal(coefoc, 3.0d0, mutype, parA, parU,parW,parH,beta,sigma)

print*, coefoc
print*, '--------------------------------------------'

call coefoclate(coefoc, 3.0d0, mutype, parA, parU,parW,parH,beta,sigma,parFVv)
print*, coefoc
print*, '--------------------------------------------'
call coefochc(coefoc, 3.0d0, mutype, parA, parU,parW,parH,beta,sigma,parFVv,1.0d0)
print*, coefoc
print*, '--------------------------------------------'

call coefochcb(coefoc, 3.0d0, mutype, parA, parU,parW,parH,beta,sigma,parFVmat,parFVv,parB,typevec,typeprob,1.0d0)
print*, coefoc
print*, '--------------------------------------------'
!call wsolver(solw,ftype,parA,parW,parH,parU,beta,Sigma)
!call vsolver(solv,ftype,parA,parW,parH,parU,parB,beta,Sigma, wcoefficients,typevec,typeprob)





end program try
