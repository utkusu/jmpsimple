!> holds global definitions for the entire project
module global
implicit none
save

integer, parameter:: dble=kind(1d0)		!< double precision
! parameters of the model
integer, parameter:: timeperiod=1 				!< parameter for frequency of data, 1 is yearly, 0.5 half a year, 0.25 quarterly
integer, parameter:: abar=7*(1/timeperiod)		!< when switch to old child- at abar, the child is old.
integer, parameter:: astar=15*(1/timeperiod)	!< childhood ends
integer, parameter:: deltamax=7*(1/timeperiod)	!< maximum birth of second child
integer, parameter:: deltamin=2*(1/timeperiod)		!< minimum birth timing for second, set to 1 to allow twins
integer, parameter:: finalage=65*(1/timeperiod)		!< final age of decision making for mother: retirement
integer, parameter:: fertmaxage=37*(1/timeperiod)	!< latest age a mom can give birth.
integer, parameter:: fertminage=20*(1/timeperiod)	!< first age a mom can give birth. Will get rid of teenage moms
integer, parameter::nperiods=22 				!< number of solution periods
integer, parameter::nfert=6 					!< number of periods where can the mother can pregnant.
!-----------------size of of stuff----------------------
integer, parameter:: shocksize1=5		!< how many shocks in emax calculation, with one child
integer, parameter:: shocksize2=5		!< how many shocks in emax calculation, with two children
integer, parameter:: shocksize3=8 		!< how many shocks in emax calculation, with one child, fecund period
integer, parameter:: Gsize=14			!< the number of regressors in the intrapolating function g. NOT INC. INTERCEPT 
integer, parameter:: Gsizeoc=10			!< the number of regressors in the intrapolating function g. for ONE CHILD FAMILIES, intercept sold separately.
integer, parameter:: Ssize=16			!< the number of elements in state space
integer, parameter:: Bsize=2			!< the number of regressors in the birth probability function other than the contraception
integer, parameter:: Bsizeexo=4			!< the number of regressors in the birth probability function other than the contraception

! estimation stuff 
integer,parameter::Npaths=5 			!<number of simulated paths for each family.
integer,parameter::simvecsize=10 		!<the size of the simulated data vector for each period. obviosly very related to ss.
integer, parameter::SampleSize=10 		!< Size of the estimation sample
integer, parameter:: Ntestage=10 		!< number of ages in which we have test scores. start at 5, end in 14.
integer, parameter:: testminage=5 		!< earliest test scores
integer, parameter:: testmaxage=14 		!< latest test scores
integer, parameter::o1size=3 			!< size of omega1
integer, parameter::o2size=4			!< size of omega2
integer, parameter::o3size=4 			!< size of omega3

! moments vector sizes used in moments()

integer, parameter:: lfpsize=4 			!< number of periods for the participation equations
integer, parameter:: expsize=4 			!< number of mother's ages at which the experieence level distributions are matched
integer, parameter:: tssize=4 			!< number of ages for which the children's test score averages to be calculated.
integer, parameter:: nreglfp=6 			!< number of regressors in labor force participation equations. intercept not included.
integer, parameter:: nregtsdiff=13 		!> number of regressors in tsdiff equations. intercept not included.

integer, parameter:: MomentSize= (2*nreglfp+2) + (3*expsize-1) + (tssize) + (nregtsdiff+1) 		!< number of moments to be matched


! parameter vector sizes

integer, parameter::parAsize=12 	!< size of parA
integer, parameter::parWsize=7 		!< size of parW
integer, parameter::parHsize=6 		!< size of parH
integer, parameter::parUsize=7 		!< size of parU

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

integer, parameter::svecage0m=4
integer, parameter::svecsch0m=4
integer, parameter::svecaqft=4
integer, parameter::svecomegaf=4

integer, parameter::sveca1=4
integer, parameter::sveca2=4
integer, parameter::svecE=5
!integer, parameter::sveclagh=3 ! you don't really need this for the simple one

real(dble), parameter:: vecage0m(svecage0m)=(/20.0d0,24.0d0,28.0d0,30.0d0/)
real(dble), parameter:: vecsch0m(svecsch0m)=(/9.0d0,16.0d0,18.0d0,20.0d0/)
real(dble), parameter:: vecomegaf(svecomegaf)=(/0.090d0,0.10d0, 0.130d0,0.15d0/)
real(dble) , parameter:: vecaqft(svecaqft)=(/.1d0,.3d0,0.6d0,0.8d0/)

! temp A values for interpolation: TODO this will need to change every year, depending on the test scores.
real(dble), parameter::veca1(sveca1)=(/1.0d0,1.20d0,1.3d0,1.4d0/)
real(dble), parameter::veca2(sveca2)=(/1.0d0,1.20d0,1.3d0,1.4d0/)
!real(dble), parameter::veca2(sveca2)=(/1.0d0,2.0d0/)

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

! ------------------------- NEW GLOBALS--------------------------------------


! MPI Globals
! don't put status yet
integer nproc, rank, ier, sender, position, number_sent, number_received, order, tag, rorder, order3, order4

! parameter globals
integer, parameter:: parameterssize=20

! wage parameters
real(dble), parameter::gparW(parWsize)=(/ 0.01 , 0.01 , 0.01 , 0.01 , 0.01 , 0.01 , 0.01 /)*1.0d0
real(dble), parameter::gparH(parHsize)=(/ 0.01 , 0.01 , 0.01 , 0.01 , 0.01 , 0.01  /) *1.0d0
real(dble) , parameter:: gpart2_parA(9)=(/ 0.01 , 0.01 , 0.01 , 0.01 , 0.01 , 0.01 , 0.01, 0.01, 0.01 /)*1.0d0
real(dble) gparB(Bsizeexo+1,nfert)


! parameters of other kind
real(dble), parameter:: grho=1.0d0
real(dble), parameter:: gbeta=0.950d0    ! beta is a free parameter for now


! factor model parameters: read these from a matrix
real(dble) glambdas(Ntestage)
real(dble) gsigmaetas(2,Ntestage)

! smoothing parameter for the smoothing
integer, parameter:: gsmpar=0.5


! intercepts - later on these will be types, too
real(dble) , parameter:: gctype=0.1d0
real(dble) , parameter:: gmtype=0.1d0
real(dble) , parameter:: gatype=0.1d0



! data to read
integer gidmat(SampleSize, MomentSize)
real(dble) gomega3data(osize, SampleSize)
real(dble) llmsmat(nperiods, SampleSize)



! global parameters for moment calculation
integer,parameter:: gexpperiods(expsize)=(/25,30,35,40/) 		! the array that holds the ages of mom for which experience moments calculated
! the array that holds the period numbers for which labor for particapation equations are estimated.
integer,parameter:: glfpperiods(lfpsize)=(/3,5,7,9/) 
! the array that holds the ages for which test score averages are calc. see explanation if confused 
integer,parameter:: gtsperiods(expsize)=(/2,4,6,8/) 				


! ------------ initialize variables---------------------





contains 

	subroutine readdata()
		implicit none
		! read idmat and omega3 -NEED TO MAKE THESE READING FROM THE ACTUAL DATA FILES.
		gidmat=1
		gomega3data=1.0d0
		llmsmat=0.5d0
	end subroutine readdata

	subroutine readsomeparam()
		implicit none
		! read some of the complicated members of the parameters space
		glambdas=1.0d0
		gsigmaetas=1.d0
		
		gparB(:,1)=(/ 0.01, 0.01, 0.01, 0.01, 0.01/)*1.0d0
		gparB(:,2)=(/ 0.01, 0.01, 0.01, 0.01, 0.01/)*1.0d0
		gparB(:,3)=(/ 0.01, 0.01, 0.01, 0.01, 0.01/)*1.0d0
		gparB(:,4)=(/ 0.01, 0.01, 0.01, 0.01, 0.01/)*1.0d0
		gparB(:,5)=(/ 0.01, 0.01, 0.01, 0.01, 0.01/)*1.0d0
		gparB(:,6)=(/ 0.01, 0.01, 0.01, 0.01, 0.01/)*1.0d0
 	end subroutine readsomeparam





end module global
