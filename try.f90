program try
use global
use randomgen
use emax
use opt
USE IFPORT ! for intel fortran only
implicit none


real(dble) coef(Gsize+1)
real(dble) coefoc(Gsizeoc+1)
real(dble) mutype(3)
real(dble) parA(12),parU(7),parW(7),parH(6),beta,sigma(shocksize1,shocksize1),parFV(Gsize+1),parFVv(Gsizeoc+1)
real(dble) parFVmat(Gsize+1,nctype+1)
real(dble) omega1(3),omega2(4),omega3(4),eps(Nmc,shocksize1)
real(dble) ftype(5)
real(dble) solw(Gsize,nperiods-deltamin+2)
real(dble) sollw(Gsize+1,nperiods-deltamin+2)
real(dble) wcoeff(Gsize+1,nperiods-deltamin+2,deltamax-deltamin+1,nttypes)
real(dble) vcoeff(Gsizeoc+1,nperiods,deltamax-deltamin+1,nttypes)
real(dble) solv(Gsizeoc+1,nperiods)
real(dble) parBmat(Bsizeexo+1,nfert)
real(dble) typevec(1)
real(dble) typeprob(1)
integer k,i

! the simulation result collectors
real(dble) SS(6,nperiods,Npaths) 		!< Simulated State space (A1,A2,E,age1,age2,agem)xnperiods,Npaths	
real(dble) outcomes(2,nperiods,Npaths) !< outcomes: wages of the father and mother.
integer choices(1,nperiods,Npaths) 	!< choices : the choice history of h.
integer xchoices(1,nperiods,Npaths) 	!< not in the data, but I will keep an history of the x choices as well.
integer birthhist(Npaths) 			!< the birth timing vec, 0 if one child throughout.
! test score stuff
real(dble) lambdas(Ntestage) 			!< old factor loadings in a two test world 
real(dble) sigmaetas(2,Ntestage) 			!<measurement error variance for each age.
real(dble) testoutcomes(4,Ntestage, Npaths)

! smoothed
real(dble) smpar
real(dble) smchoices(3,nperiods, Npaths)
real(dble) smexperience(nperiods, Npaths)
real(dble) smAs(2,nperiods, Npaths)
real(dble) smtestoutcomes(4,Ntestage, Npaths)
		
			
! stuff for simulation
real(dble) intercepts(3),a1type(2),pa1type(2),llmsvec(SampleSize),rho
integer id


omega1=(/8.d0,5.0d0,0.0d0/)
omega2=(/2.0d0,1.0d0,20.0d0,20.0d0/)
omega3=(/10.0d0,1.0d0,18.0d0,10.0d0/)

parA=0.10d0
parU=1.0d0
parU(4)=100.0d0
parU(1)=-100
parW=0.1d0*0.010
parW(6)=-0.01d0
parH=0.1d0*0.0000020
Sigma=0.5d0*0.10
beta=0.95d0
parFV=0.10d0
parFVv=0.10d0
parFVmat=1.0d0
mutype=1.0d0
eps=0.0d0
ftype=(/1.0d0,1.0d0, 1.0d0,1.0d0,0.0d0/)
typevec=0.5d0
typeprob=1.0d0
wcoeff=0.10d0
vcoeff=0.1d0
! try the coefochcb
parBmat=0.0001d0
intercepts=(/1.0d0,1.0d0,1.0d0/)
a1type=(/0.0d0,1.0d0/)
pa1type=(/0.3d0,0.7d0/)
rho=1.0d0
id=18

smpar=1.0d0
lambdas=1.0d0
sigmaetas=0.1d0

call simhist(SS,outcomes,testoutcomes, choices, xchoices,birthhist,smchoices, smexperience, smAs, smtestoutcomes,omega3,intercepts,parA,parU,parW,parH,beta,Sigma,a1type,pa1type,parBmat,vcoeff,wcoeff,llmsvec,id,rho,lambdas,sigmaetas,smpar)
!print*, choices(1,:,80)
!print*, '------------------------'
!print*, xchoices(1,:,80)


!print*,'------------------'
!print*, testoutcomes(:,:,80)


!print*, '-----------------'
!print*, 'STATE SPACE'
!do i=1,22
	!print*, SS(:,i,8)
	!print*, '----------------'
!end do 

print*,smchoices(:,:,80)
print*,'---------'
print*, smexperience(:,80)

!print*, emaxlate(omega1,omega2,omega3,mutype,eps,parA,parU,parW,parH,beta,parFV)
!print*, emaxhc(omega1,omega2,omega3,mutype,eps,parA,parU,parW,parH,beta,parFV,1.0d0)
!print*, emaxhcx(omega1,omega2,omega3,mutype,eps,parA,parU,parW,parH,beta,parFV,1.0d0)
!print*, emaxoclate(omega1,omega2,omega3,mutype,eps,parA,parU,parW,parH,beta,parFVv)
!print*, emaxochc(omega1,omega2,omega3,mutype,eps,parA,parU,parW,parH,beta,parFVv,1.0d0)
!print*, emaxochcb(omega1,omega2,omega3,mutype,eps,parA,parU,parW,parH,beta,parFVmat,parFVv,parB,typevec,typeprob,1.0d0)

!call coeffinal(coef, 22.0d0, 2.0d0, mutype, parA, parU,parW,parH(5:6),beta,sigma)


!print*, coef
!print*, '--------------------------------------------'
!call coeflate(coef, 3.0d0, 2.0d0, mutype, parA, parU,parW,parH,beta,sigma,parFV)
!print*, coef
!print*, '--------------------------------------------'

!call coefhc(coef, 3.0d0, 2.0d0, mutype, parA, parU,parW,parH,beta,sigma,parFV,1.0d0)

!call coefhcx(coef, 3.0d0, 2.0d0, mutype, parA, parU,parW,parH,beta,sigma,parFV,1.0d0)
!print*, coef
!print*, '--------------------------------------------'

!call coefocfinal(coefoc, 3.0d0, mutype, parA, parU,parW,parH,beta,sigma)

!print*, coefoc
!print*, '--------------------------------------------'

!call coefoclate(coefoc, 3.0d0, mutype, parA, parU,parW,parH,beta,sigma,parFVv)
!print*, coefoc
!print*, '--------------------------------------------'
!call coefochc(coefoc, 3.0d0, mutype, parA, parU,parW,parH,beta,sigma,parFVv,1.0d0)
!print*, coefoc
!print*, '--------------------------------------------'

!call coefochcb(coefoc, 3.0d0, mutype, parA, parU,parW,parH,beta,sigma,parFVmat,parFVv,parB,typevec,typeprob,1.0d0)
!print*, coefoc
!print*, '--------------------------------------------'
!call wsolver(solw,ftype,parA,parW,parH,parU,beta,Sigma)
!call vsolver(solv,ftype,parA,parW,parH,parU,parB,beta,Sigma, wcoefficients,typevec,typeprob)

!jcall wsolver(sollw,3.0d0,ftype,parA,parW,parH(5:6),parU,beta,Sigma,1.0d0)

!print *,sollw(:,4)
!call vsolver(solv,ftype,parA,parW,parH,parU,parBmat,beta,Sigma, wcoeff,typevec,typeprob,1.0d0)


end program try


