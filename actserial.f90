!!> MAIN PROGRAM FOR THE SOLUTION AND ESTIMATION OF JMP MODEL
!! ifort global.f90 randomgen.f90 emax.f90 opt.f90 actserial.f90 -llapack -CB -debug -o exec
!program act
!use global
!use randomgen
!use emax
!use opt
!USE IFPORT ! for intel fortran only
!implicit none

!!real(dble) coef(Gsize+1)
!!real(dble) coeaf(Gsize+1)
!!real(dble) mutype(3)


!! parameters of the model
!real(dble) parA(12),parU(7),parW(7),parH(6),beta,sigma1(shocksize1,shocksize1),parB(Bsizeexo+1)
!real(dble) ctype, mtype, atype, a1type(na1type)
!!real(dble) pctype(nctype),pmtype(nmtype),patype(natype),pa1type(na1type)
!real(dble) condprob
!real(dble) packed(parAsize+parWsize+parUsize+parHsize+1+(shocksize1+shocksize1*(shocksize1-1)/2)+Bsizeexo+1+4)
!integer npack
!real(dble) ftypemat(2,na1type*(deltamax-deltamin+1))
!real(dble) ftypematoc(5,nttypesoc)
!real(dble) ocp(nctype+Bsize+2+nctype*nmtype)

!real(dble) solw(Gsize+1, nperiods-deltamin+1) !< collects solution coef.
!real(dble) solwall(Gsize+1, nperiods-deltamin+1, deltamax-deltamin+1,nttypes) !< collects solution coef for all types
!real(dble) solv(Gsize+1,nperiods)
!real(dble) solvall(Gsize+1,nperiods, nttypesoc)
!integer momorder
!integer k,j,l,m
!real(dble) wcoef(Gsize+1, nperiods-deltamin+2, deltamax-deltamin+1,nctype) 		! for the vsolver trial
!real(dble) start,endtime, vtime
!real(dble) rho
!integer nftype
!!real(dble) omega1(4),omega2(5),omega3(4),eps(Nmc,shocksize1)
!!real(dble) ftype(5)
!!real(dble) solw(Gsize,nperiods-deltamin+2,deltamax-deltamin+2)
!!real(dble) wcoefficients(Gsize+1,nperiods-deltamin+2,deltamax-deltamin+2,2)
!!real(dble) solv(Gsize+1,nperiods)
!!real(dble) typevec(2)
!!real(dble) typeprob(2)
!rho=1.0d0
!npack=parAsize+parWsize+parUsize+parHsize+1+(shocksize1+shocksize1*(shocksize1-1)/2)+Bsizeexo+1+4
!nftype=2*na1type*(deltamax-deltamin+1)



	!call initparam(parA,parU,parW, parH, beta, sigma1,parB,ctype, mtype, atype, a1type,condprob)
	!call utku_pack(packed,npack,parA,parU,parW, parH, beta, sigma1,parB,ctype,mtype, atype, condprob)
	!call utku_unpack(packed,npack,(/12,7,7,6,1,15,(Bsizeexo+1),1,1,1,1/),parA,parU,parW, parH, beta, sigma1,parB,ctype,mtype,atype,condprob)
	


!end program act
