!!> MAIN PROGRAM FOR THE SOLUTION AND ESTIMATION OF JMP MODEL
program actserial
use global
use randomgen
use emax
use opt
USE IFPORT ! for intel fortran only
implicit none
!include 'mpif.h'

!real(dble) coef(Gsize1+1)
!real(dble) coeaf(Gsize1+1)
!real(dble) mutype(3)

! mpi stuff
!integer nproc, rank, ier, status(MPI_STATUS_SIZE),i, sender, position, number_sent, number_received, order,tag, rorder
real(dble) buffer(1000)

! parameters of the model
real(dble) parA(10),parU(7),parW(7),parH(4),beta,sigma1(shocksize1,shocksize1),sigma2(shocksize3,shocksize3)
real(dble) ctype(nctype), mtype(nmtype), atype(natype), a1type(na1type)
real(dble) pctype(nctype),pmtype(nmtype),patype(natype),pa1type(na1type)

real(dble) packed(10+7+7+4+1+(shocksize1+shocksize1*(shocksize1-1)/2)+(shocksize3+shocksize3*(shocksize3-1)/2))
integer npack

real(dble) solw(Gsize1+1, nperiods-deltamin+2, deltamax-deltamin+2) !< collects solution coef.
real(dble) solwall(Gsize1+1, nperiods-deltamin+2, deltamax-deltamin+2,nttypes) !< collects solution coef for all types
real(dble) ftype(5)
real(dble) ftypemat(5,nttypes)


!real(dble) omega1(4),omega2(5),omega3(4),eps(Nmc,shocksize1)
!real(dble) ftype(5)
!real(dble) solw(Gsize1,nperiods-deltamin+2,deltamax-deltamin+2)
!real(dble) wcoefficients(Gsize1+1,nperiods-deltamin+2,deltamax-deltamin+2,2)
!real(dble) solv(Gsize1+1,nperiods)
!real(dble) parB(4)
!real(dble) typevec(2)
!real(dble) typeprob(2)

end program actserial
