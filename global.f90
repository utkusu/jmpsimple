!> holds global definitions for the entire project
module global
implicit none
save

integer, parameter:: dble=kind(1d0)		!< double precision
! parameters of the model
integer, parameter:: timeperiod=1 				!< parameter for frequency of data, 1 is yearly, 0.5 half a year, 0.25 quarterly
integer, parameter:: abar=9*(1/timeperiod)		!< when switch to old child- at abar, the child is old.
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
integer,parameter::Npaths=100			!<number of simulated paths for each family.
integer,parameter::simvecsize=10 		!<the size of the simulated data vector for each period. obviosly very related to ss.
integer, parameter::SampleSize=301 		!< Size of the estimation sample
integer, parameter:: Ntestage=10 		!< number of ages in which we have test scores. start at 5, end in 14.
integer, parameter:: Noldtest=5 		!< number of ages used in diffdiff regressions
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
integer, parameter:: nregtsdiff=14 		!> number of regressors in tsdiff equations. intercept not included.

integer, parameter:: MomentSize= (2*nreglfp+2) + (3*expsize-1) + (tssize) + (nregtsdiff+1) 		!< number of moments to be matched
integer, parameter::idmatsize=lfpsize+expsize+expsize-1+tssize+Noldtest 						!< idmat's second dimension

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
integer, parameter::nttypesoc=2 		!< number of total types for one child families


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
integer, parameter:: Nmc=50	!<monte carlo integration draw size
! sizes of interpolating grid parts
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

! temp A values for interpolation. See setvecAs to see how these were set. 
real(dble) veca1(sveca1)
real(dble) veca2(sveca2)


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
integer, parameter:: parsize=20


! wage parameters
! >>>>>>>>>>>>>  put the betaf and betam here <<<<<<<<<<<<<<<<<<<<<<<<<
real(dble), parameter::gparW(parWsize)=(/ 0.0315141d0 , 0.22902d0 , 0.31238d0 , -0.0053421d0 , 0.0568242d0 , -0.026432d0 , -0.4756638d0 /)
real(dble), parameter::gparH(parHsize)=(/ 0.066522d0 , 0.037802d0 , -0.0004797d0 , 0.0549845d0 , -0.0013848d0 , 8.865523d0  /)
! removing age0m, it had a huge negative coefficient in front of age0m x tau2, making the effect of tau2 negative, and dragging the
! A1 and A2 to negative as the mom did work full time.
! THIS WAS THE OLD ONE: real(dble) , parameter:: gpart2_parA(9)=(/ 9.9107d0 , 0.0127697d0 , 3.361842d0 , -0.07951778d0 ,  4.8295807d0 , 0.22579117d0, -1.8741083d0, -0.11462243d0, 0.0000160622d0 /)
real(dble) , parameter:: gpart2_parA(9)=(/ 8.9498d0 , -0.0987d0 , 3.738713d0 , 0.0d0 ,  3.221d0 , 0.11145d0, -1.874183d0, 0.0d0, 0.0000160622d0 /)


real(dble) gparBmat(Bsizeexo+1,nfert)

! parameters of other kind
real(dble), parameter:: grho=1.0d0
real(dble), parameter:: gbeta=0.950d0    ! beta is a free parameter for now


! factor model parameters: read these from a matrix
real(dble) glambdas(Ntestage)
real(dble) gsigmaetas(2,Ntestage)

! smoothing parameter for the smoothing
real(dble) , parameter:: gsmpar=0.0005d0

! estimated shock variances for wages
real(dble), parameter:: gwshock=0.25068398d0
real(dble), parameter:: gwfathershock=0.467246870d0

! intercepts - later on these will be types, too

! the constant estimate from the estimation of logwage regression needs to be entered here
! the code is like this, because it helps to facilitate to make this hetero. later.
!>>>>>>>>>>>>>>>>>> INTERCEPT HERE<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
real(dble) , parameter:: gmtype=-2.340616d0

!real(dble) , parameter:: gatype=0.01d0



! data to read
integer gidmat(SampleSize,idmatsize)
real(dble) gomega3data(o3size, SampleSize)
real(dble) llmsmat(nperiods, SampleSize)



! global parameters for moment calculation
integer,parameter:: gexpperiods(expsize)=(/25,30,35,40/) 		! the array that holds the ages of mom for which experience moments calculated
! the array that holds the period numbers for which labor for particapation equations are estimated.
integer,parameter:: glfpperiods(lfpsize)=(/3,5,7,9/) 
! the array that holds the ages for which test score averages are calc. see explanation if confused 
integer,parameter:: gtsperiods(expsize)=(/2,4,6,8/) 				


real(dble) parameters(parsize), lb(parsize), ub(parsize), targetvec(MomentSize), weightmat(MomentSize,MomentSize)

! just moved this here
real(dble) momentvec(MomentSize) 

integer itercounter, evaliter

real(dble), parameter:: wm=2000.0d0  !< wage multiplier, the hours of work in a period if the mom works full time. 

! grid add-on

integer, parameter:: parmatsize=4 !<size of the matrix of parameters to be evaluated
real(dble) distvec(parmatsize) 
real(dble) parmat(parmatsize,parsize) 


! ------------ initialize variables---------------------

! THIS IS THE AREA WHERE WE READ SOME EXTERNAL FILES. HERE IS A LIST OF FILES THAT ARE REQUIRED
! rest of the parameters are manually set above. 
! 1- idmat.csv (Sample Size x idmatsize): Integer values, used in the moments calculation. 
! 	Indicates which members of the sample enters which calculation.  See explanation in moments routine
! 2- omega3mat.csv (o3size x Sample): holds the omega3 values for the sample. double precision.
! 3- llms.csv: holds llms data for the sample, double. T x Nsample
! 4- factorpar.csv (3 x Ntestage) : holds the factor analysis results. First row hold factor loading for the second measurement equation
! 	2 and 3 holds measurement error variance for the first and second equations, double
! 5- optim.csv (3 x parsize): holds initial values, lower and upper bounds for the optimization in its rows. double
! 6- setii.csv ( MomentSize+1 x MomentSize): first row holds target moments, the rest of it is the optimal weighting matrix for the
!    indirect inference.
! 7- NEW: parBmat (Bsizeexo+1,nfert): columns hold the the parameters of bprob 

contains 

	subroutine readdata()
		implicit none
		integer i,j,l
		! read idmat and omega3 -NEED TO MAKE THESE READING FROM THE ACTUAL DATA FILES.
		
		open(unit=12, file="idmat.csv")

		do i=1, SampleSize
			read(12,*) (gidmat(i,j),j=1,idmatsize)
		end do
		close(12)
		
		open(unit=13, file="omega3mat.csv")
		do i=1, o3size
			read (13,*) (gomega3data(i,j),j=1,SampleSize)
		end do
		close(13)

		open(unit=14, file="llms.csv")
		do i=1, nperiods
			read (14,*) (llmsmat(i,j),j=1,SampleSize)
		end do
		close(14)
	end subroutine readdata

	subroutine readsomeparam()
		implicit none
		! read some of the complicated members of the parameters space
		real(dble) factormat(3,Ntestage)
		! first line lambdas
		! second line variance for first equation, 3rd variance of second.
		integer i,j

		open(unit=12, file="factorpar.csv")
		do i=1, 3
			read(12,*) (factormat(i,j),j=1,Ntestage)
		end do
		close(12)	
		glambdas=factormat(1,:)
		gsigmaetas(1,:)=factormat(2,:)
		gsigmaetas(2,:)=factormat(3,:)
		

		! also assign gparBmat by hand.	
		open(unit=19,file="parBmat.csv")
		do i=1,Bsizeexo+1
			read(19,*) (gparBmat(i,j),j=1,nfert)
		end do
		close(19)
 	end subroutine readsomeparam
	
	! use the the mean test scores and their variance to set up grid points for veca1 and veca2
	! These points are used in the interpolations
	subroutine setvecAs(period, delta, noc)
		implicit none
		real(dble) period, delta
		integer noc
		real(dble) meantestscore(Ntestage)
		real(dble) variancetestscores(Ntestage)
		integer i, j, age1, age2
		! insert data here!!!
		meantestscore=(/ 14.48d0, 17.8d0 , 25.5d0 , 35.38d0, 40.28d0,46.28d0, 50.27d0, 53.67d0, 56.06d0, 58.87d0  /)
		variancetestscores=(/16.0d0, 25.0d0, 64.0d0, 81.0d0, 81.0d0,81.0d0, 100.0d0, 121.0d0, 100.0d0,121.0d0 /)
		
		! to understand this later: if age1
		age1=max(min(nint(period)-testminage+1,Ntestage),1)
		age2=max(min((nint(period)-nint(delta)+1)-testminage+1, Ntestage),1)

		veca1=(/ meantestscore(age1)-2*SQRT(variancetestscores(age1)),meantestscore(age1)-SQRT(variancetestscores(age1)), &
			& meantestscore(age1)+SQRT(variancetestscores(age1)),meantestscore(age1)+2*SQRT(variancetestscores(age1)) /) 
		if (noc==2) then		
			veca2=(/ meantestscore(age2)-2*SQRT(variancetestscores(age2)),meantestscore(age2)-SQRT(variancetestscores(age2)) ,&
				& meantestscore(age2)+SQRT(variancetestscores(age2)),meantestscore(age2)+2*SQRT(variancetestscores(age2)) /) 
		end if
	end subroutine setvecAs

	! set initial parameters and bounds for optimization
	subroutine setoptimstuff()
		implicit none 
		real(dble) parmat(3,parsize)
		! first row initial value, then lb, then u2
		integer i,j
		open(unit=15, file="optim.csv")
		do i=1,3
			read(15,*) (parmat(i,j),j=1,parsize)
		end do
		close(15)	
		parameters=parmat(1,:)
		lb=parmat(2,:)
		ub=parmat(3,:)
	end subroutine setoptimstuff

	! set indirect inference things, first
	subroutine setii()
		implicit none
		! first lline=targetvec the rest is the weight matrix
		real(dble) setiimat(MomentSize+1,MomentSize)
		integer i,j
		open(unit=16, file="setii.csv")
		do i=1,MomentSize+1
			read(16,*) (setiimat(i,j),j=1,MomentSize)
		end do
		close(16)	
		targetvec=setiimat(1,:)
		weightmat=setiimat(2:MomentSize+1,:)

end subroutine setii

! ---- parameter utilities------
! utility to print variable parameters
subroutine printpar(parvec)
	implicit none
	real(dble), intent(in) :: parvec(parsize)
		print*, '------- PARAMETERS----------'
		print*, 'gh1=', parvec(1)
		print*, 'gh2=', parvec(2)
		print*, 'gh3=', parvec(3)
		print*, 'phi=', parvec(4)
		print*, 'alpha=', parvec(5)
		print*, 'alpha2=', parvec(6)
		print*, 'alpha3=', parvec(7)
		print*, 'alpha4=', parvec(8)
		print*, 'alpha5=', parvec(9)
		print*, 'a1type1=', parvec(10)
		print*, 'a1type2=', parvec(11)
		print*, 'pa1type=', parvec(12)
		print*, 'sigma11=', parvec(13)
		print*, 'sigma22=', parvec(14)
		print*, 'sigma33=', parvec(15)
		print*, 'sigma12=', parvec(16)
		print*, 'sigma13=', parvec(17)
		print*, 'sigma23=', parvec(18)
		print*, 'beta=', parvec(19)
		print*, 'muc=', parvec(20)
		print*, '------ End of Parameters -----'
end subroutine printpar

! if called (after optimstuff) sets the parameter vectors
subroutine setparameters(mode)
	implicit none
	integer, intent(in):: mode ! =0 if hand written, anything else if reading from a line by line from givenpar.csv
	real(dble) gh1, gh2, gh3, phi, alpha, alpha2, alpha3, alpha4, alpha5, a1type1 ,a1type2 ,pa1type
	real(dble) 	sigma11 ,sigma22 ,sigma33 ,sigma12 ,sigma23, sigma13 ,beta ,muc
integer i 	
	gh1= 		0.5525d0 
	gh2= 		4.75025d0
	gh3= 		0.33825d0
	phi= 		-9.5971d0 
	alpha= 		-1.24d0
	alpha2=  	0.15d0
	alpha3= 	-9.76d0	
	alpha4= 	-0.76d0
	alpha5= 	-2.57d0
	a1type1= 	0.10374625d0	
	a1type2= 	1.14124625d0
	pa1type= 	0.5124998d0	
	sigma11= 	1.0374625d0
	sigma22= 	1.0374625d0
	sigma33= 	1.0374625d0
	sigma12= 	0.25d0
	sigma13= 	0.25d0
	sigma23= 	0.25d0
	beta= 	 	0.951874d0
	muc= 		7.9375374d0
	if (mode==0) then
		parameters= (/gh1, gh2, gh3, phi, alpha, alpha2, alpha3, alpha4, alpha5, a1type1, a1type2, pa1type, sigma11, sigma22, sigma33, sigma12, sigma23, sigma13, beta, muc/)
	else
		open(unit=17,file="givenpar.csv")
		do i=1, parsize 
			read (17,*)  parameters(i)
		end do 
		close(17)	
	end if 	
end subroutine setparameters	

! slightly change the existing parameter vector, it takes VECTOR values.
subroutine twistparameters(dimvec, valvec)
	implicit none
	integer, intent(in) :: dimvec(:)  ! dimensions of the parameter vector you want to change
	real(dble), intent(in) :: valvec(:)  ! values that needs to be inserted, vector of same size
	integer n, i
	n=size(dimvec)
	do i=1,n
		parameters(dimvec(i))=valvec(i)
	end do
end subroutine twistparameters


! subroutine to write interpolation coef to files
subroutine writeintpar(solwall, solvall)
	implicit none
	real(dble), intent(in):: solvall(Gsizeoc+1,nperiods, nttypes)
	real(dble), intent(in):: solwall(Gsize+1,nperiods-deltamin+2, deltamax-deltamin+1, nttypes)
	!integer, intent(in):: itercounter
	integer i,j,k,l, m

	if (itercounter==1) then
		open(77,file='wcoefALL.txt')
		write(77,*) '---- ITERCOUNTER=     ', itercounter
		do l=1,2
			!write(77,*) "-------------- type order=",l,"------------" 
			do m=1,deltamax-deltamin+1
				!write(77,*) "########   delta type=",m,"------------"
				write(77,70) ((solwall(j,k,m,l),k=1,nperiods-deltamin+2),j=1,Gsize+1)
			end do
		end do
		close(77)

		open(66,file='vcoefALL.txt')
		write(66,*) '---- ITERCOUNTER=     ', itercounter
		do l=1,2
			write(66,60) ((solvall(j,k,l),k=1,nperiods),j=1,Gsizeoc+1)
		end do
		close(66)
	else 
		open(77,file='wcoefALL.txt', position="append")
		write(77,*) '---- ITERCOUNTER=     ', itercounter, '------'
		do l=1,2
			!write(77,*) "-------------- type order=",l,"------------" 
			do m=1,deltamax-deltamin+1
				!write(77,*) "########   delta type=",m,"------------"
				write(77,70) ((solwall(j,k,m,l),k=1,nperiods-deltamin+2),j=1,Gsize+1)
			end do
		end do
		close(77)

		open(66,file='vcoefALL.txt', position="append")
		write(66,*) '---- ITERCOUNTER=     ', itercounter, '------'
		do l=1,2
			write(66,60) ((solvall(j,k,l),k=1,nperiods),j=1,Gsizeoc+1)
		end do
		close(66) 
	end if
	70 format(22F46.9)
	60 format(22F46.9)
end subroutine writeintpar

! writing the moments calculated to a file for each iteration
subroutine writemomentvec(momentvec)
	implicit none
	real(dble), intent(in) :: momentvec(MomentSize)
	integer i,j,k,l, m
	character (len=6):: momentnames(MomentSize)
	momentnames= (/"ft_sch", "ft_afq", "ft_age","ft_bby","ft_exp", "ft_2ch", "ft_con" ,   &
		"pt_sch", "pt_afq", "pt_age","pt_bby","pt_exp", "pt_2ch","pt_con" ,   &
		"mexp25", "vexp25", "mexp30", "vexp30", "mexp35", "vexp35","mexp40", "vexp40", "cv2530","cv3035", "cv3540", &
		"mcom_1", "mcom_2","mcom_3","mcom_4", & 
		"d_tau1", "d_tau2", "dt1sch", "dt1afq", "dt1age", "dt2sch", "dt2afq", "dt2age", &
		"d_agem", "d_schm", "d_afqt", "d_agef", "d_agf2", "d_ages", "consta"/)

	if (itercounter==2) then
		open(88,file='moments.txt')
		write(88,81) momentnames
		write(88,80) itercounter-1, momentvec 
		close(88)
	else 
		open(88,file='moments.txt', position="append")
		write(88,80) itercounter-1, momentvec 
		close(88)

	end if
	80 format(1I5,<MomentSize>G19.9)  
	81 format(<MomentSize+1>G19.9)
	! Note: this <> thing only works in ifort.
end subroutine writemomentvec


!>writing the iteration number, distance, and parameters to a file for each iteration
subroutine writeresults(dist, parvec) 
	implicit none
	real(dble), intent(in) :: dist
	real(dble), intent(in):: parvec(parsize)
	character (len=7):: parnames(parsize+1)
	integer i,j,k,l, m
	parnames= (/"distanc","gh1_sch", "gh2_afq", "gh3_age",  &
		"uti_phi", "u_alpha", "a2_cons", "a3_cexp", "a4_baby", "a5_hbby", &
		"a1_type", "a2_type", "prob_a1", &
		"sigma11", "sigma22", "sigma33", "sigma12", "sigma13", "sigma23", &
		"_beta__", "muchild"/)

	if (itercounter==2) then
		open(99,file='results.txt')
		write(99,91) parnames
		write(99,90) itercounter-1,dist, parvec 
		close(99)
	else 
		open(99,file='results.txt', position="append")
		write(99,90) itercounter-1, dist, parvec 
		close(99)

	end if
	90 format(1I5, <parsize+1>G19.9)  
	91 format(<parsize+1>G19.9)
	! Note: this <> thing only works in ifort.
end subroutine writeresults

end module global
