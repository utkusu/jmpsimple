!> holds global definitions for the entire project
module global
implicit none
save

integer, parameter:: dble=kind(1d0)		!< double precision
! parameters of the model
integer, parameter:: period=1 				!< parameter for frequency of data, 1 is yearly, 0.5 half a year, 0.25 quarterly
integer, parameter:: abar=7*(1/period)		!< when switch to old child- at abar, the child is old.
integer, parameter:: astar=15*(1/period)	!< childhood ends
integer, parameter:: deltamax=7*(1/period)	!< maximum birth of second child
integer, parameter:: deltamin=2*(1/period)		!< minimum birth timing for second, set to 1 to allow twins
integer, parameter:: finalage=65*(1/period)		!< final age of decision making for mother: retirement
integer, parameter:: fertmaxage=37*(1/period)	!< latest age a mom can give birth.
integer, parameter:: fertminage=20*(1/period)	!< first age a mom can give birth. Will get rid of teenage moms
integer, parameter::nperiods=22 				!< number of solution periods
!-----------------size of of stuff----------------------
integer, parameter:: shocksize1=5		!< how many shocks in emax calculation, with one child
integer, parameter:: shocksize2=5		!< how many shocks in emax calculation, with two children
integer, parameter:: shocksize3=8 		!< how many shocks in emax calculation, with one child, fecund period
integer, parameter:: Gsize=14			!< the number of regressors in the intrapolating function g. NOT INC. INTERCEPT 
integer, parameter:: Gsizeoc=10			!< the number of regressors in the intrapolating function g. for ONE CHILD FAMILIES, intercept sold separately.
integer, parameter:: Ssize=16			!< the number of elements in state space
integer, parameter:: Bsize=2			!< the number of regressors in the birth probability function other than the contraception
integer, parameter:: Bsizeexo=4			!< the number of regressors in the birth probability function other than the contraception

integer, parameter::o1size=3 			!< size of omega1
integer, parameter::o2size=4	!< size of omega1
integer, parameter::o3size=4 			!< size of omega1


real(dble), parameter::emaxscale=1.0d0 !< scale emaxs down  !!! CAN YOU DO THIS??? MAYBE ONLY IF USE THE SAME SCALE for U, too.

! number of types
integer, parameter::nctype=1 			!< number of unobserved types for children
integer, parameter::nmtype=1 			!< number of unobserved types for mom (ability)
integer, parameter::natype=1 			!< number of unobserve alpha types
integer, parameter::na1type=2 			!< number of unobserve alpha1 types
integer, parameter::nttypes=2 			!< number of total types for two child families
integer, parameter::nttypesoc=2 		!< number of total types for two child families


! x and c grids and related utilities
integer,parameter::xgridsize=11 
integer,parameter::cgridsize=1 
integer,parameter::cxgridsize=11

real(dble), parameter::xgrid(xgridsize)=(/0.01d0,0.1d0,0.2d0,0.3d0,0.4d0,0.5d0,0.6d0,0.7d0,0.8d0,0.9d0,0.99d0/)
real(dble), parameter::cgrid(cgridsize)=0.5d0  ! fixing the consumption choice at 0.5 for now.
real(dble), parameter::onescgrid(cgridsize)=1.0d0
real(dble), parameter::onesxgrid(xgridsize)=1.0d0
real(dble), parameter::onescxgrid(cxgridsize)=1.0d0
real(dble), parameter::onescgridwithbirth(cgridsize*2)=1.0d0

! emax interpolation parameters
integer, parameter:: Nmc=100		!<monte carlo integration draw size

! following are the fixed vectors of observed types for the emax approximations
integer, parameter::svecage0m=2
integer, parameter::svecsch0m=2
integer, parameter::svecaqft=2
integer, parameter::svecomegaf=2

integer, parameter::sveca1=2
integer, parameter::sveca2=2
integer, parameter::svecE=5
!integer, parameter::sveclagh=3 ! you don't really need this for the simple one

real(dble), parameter:: vecage0m(svecage0m)=(/20.0d0,24.0d0/)
real(dble), parameter:: vecsch0m(svecsch0m)=(/9.0d0,16.0d0/)
real(dble), parameter:: vecomegaf(svecomegaf)=(/9.0d0,16.0d0/)
real(dble) , parameter:: vecaqft(svecaqft)=(/.1d0,.3d0/)

! temp A values for interpolation: TODO this will need to change every year, depending on the test scores.
real(dble), parameter::veca1(sveca1)=(/1.0d0,2.0d0/)
real(dble), parameter::veca2(sveca2)=(/1.0d0,2.0d0/)

! temp A values for interpolation

!real(dble), parameter:: vecdelta(3)=(/2.0d0,4.0d0,6.0d0/)
!integer, parameter:: nintpfinal=432*375*3 			!< # of interpolating points: above 5 each for E and As and 3 for delta
! # of interpolating points: above 5 each for E and As 
integer, parameter::nintp= svecage0m*svecsch0m*svecomegaf*svecaqft*sveca1*sveca2*svecE
integer, parameter:: nintpoc=svecage0m*svecsch0m*svecomegaf*svecaqft*sveca1*svecE

! auxillary
real(dble), parameter :: nr_small = 1.0d-300
integer, parameter:: gseed(12)=1
real(dble), parameter :: pi=3.141592653589793238462643d0

! lapack parameters
integer, parameter:: blocksize=64


end module global
