!> This module holds the functions for random number generation and picking random ss points
module randomgen
use global
implicit none

contains

!#####################---------RANDOM NORMALS------##################################

!> generate an array of n normal random variables using Box Muller method, with mean mu and VARIANCE sigma
!>  opt determines how the seed is assigned.
!>opt=0: reset to default seed, opt=1 : random seed using system clock,
!> opt=2: use system clock but modify the seed by a multiplier. This is going to help create different arrays for different nodes while using mpi
!> new option: opt=3 - provide the seed yourself as an integer array of desired size
function randnv(n,mu,sigma, opt,seedm, userseed)

	implicit none
	real(dble) randnv(n)					!< output, scalar
	integer n, opt							!< number of random variables, option switch
	real(dble) mu, sigma					!< mean and standard deviation
	integer, optional, intent(in):: seedm	!< optional seed multiplier
	integer, optional :: userseed(:)		!<optional user supplied seed, array of size 15.

	integer nseed, clock, i
	integer, allocatable:: seed(:)

	real(dble), dimension(n):: ubig, sn
	real(dble), dimension(n/2):: u1, u2, omega, R, x, y

	!----------setting the seed------------
	! seed array size: it will be 12 in gfortran but elsewhere might be different (seems to be 2 with ifort)
	call random_seed(size = nseed)
	allocate(seed(nseed))
	call system_clock(count=clock)
	select case (opt)
		case(0)
			call random_seed(get=seed)
		case(1)
			seed = clock + 37 * (/ (i - 1, i = 1, nseed) /)
		case(2)
			seed = clock + 37 * seedm* (/ (i - 1, i = 1, nseed) /)
		case(3)
			seed=userseed
	end select
	!------------
	! seed the seed
	call random_seed(put=seed)
	deallocate(seed)
	! generate a n randu array, and chop into two.
	call random_number(ubig)
	u1=ubig(1:(n/2))
	u2=ubig((n/2)+1:n)
	! array box muller
	omega=2*pi*u1
	R=sqrt(-2*log(u2))
	x=R*cos(omega)
	y=R*sin(omega)
	sn=(/x,y/)
	randnv=mu+sqrt(sigma)*sn

end function

!> returns an n x k matrix generated from a normal distribution with mean mu and variance sigma..
!> opt=0: reset to default seed, opt=1 : random seed using system clock, opt=2: use system clock
!> but modify the seed by a multiplier. This is going to help create different arrays for different
!>  nodes while using mpi, opt=3 user supplies the seed array.
function randmnv(n,k,mu,sigma, opt, seedm, userseed)
implicit none
	integer n, k							!< draw size, number of variables
	real(dble), dimension(n,k):: randmnv	!< output (n,k)
	real(dble) mu(k)						!< mean vector of size k
	real(dble), dimension(k,k):: sigma		!< variance covariance matrix
	integer opt 							!< option switch
	integer, optional :: userseed(:)		!< user provided seed
	integer, optional, intent(in):: seedm	!< seed multiplier


	integer i,j, damn
	real(dble), dimension(k,k):: T
	real(dble), dimension(k,n):: Z, laststep
	T=sigma
	call DPOTRF('l', k, T, k, damn )
	! make T proper lower triangule matrix
	do i=1,k
		do j=1,k
			if (j>i) T(i,j)=0.0d0
		end do
	end do
	! generate a matrix N(0,1)'

	do i=1,k
		Z(i,:)=randnv(n,0.0d0,1.0d0, opt,seedm,userseed)
	end do
	! multiply with T
	laststep=matmul(T,Z)
	! add the means88
	do j=1,k
		laststep(j,:)=laststep(j,:) + mu(j)
	end do
	! lastly, transpose
	randmnv=transpose(laststep)
	! CORRECT THIS SO YOU DON'T HAVE TO TRANSPOSE AT THE END
end function

!**********************************************

! ------------------BIRTH PROBABILITY FUNCTION------------------------------


!> Normal CDF via Navarro
function cdf_normal(zin,mu,var)
implicit none
	real(dble) cdf_normal
	real(dble), intent(in) :: zin
	real(dble), intent(in) :: mu,var
	real(dble) :: zabs,p,arg,logpdf,z,std
	real(dble), parameter :: p0=220.2068679123761d0,p1=221.2135961699311d0,p2=112.0792914978709d0, &
		 p3 = 33.91286607838300d0,p4 = 6.373962203531650d0,p5 = .7003830644436881d0, &
		 p6 = .3526249659989109d-01,q0 = 440.4137358247522d0,q1 = 793.8265125199484d0, &
		 q2 = 637.3336333788311d0,q3 = 296.5642487796737d0,q4 = 86.78073220294608d0, &
		 q5=16.06417757920695d0,q6=1.755667163182642d0,q7=.8838834764831844d-1,cutoff = 7.071d0, &
		 logroot2pi = 0.918938533204672780563271317078d0
	z=zin-mu
	zabs=abs(z)  
	if (zabs<nr_small) then
	   cdf_normal=0.5d0 !is this right?
	   return
	end if
	std = sqrt(var)
	if (std<nr_small) then
	   if (zin-mu>=0.0d0) then
		  cdf_normal = 1.0d0
	   else if (zin-mu<0.0d0) then
		  cdf_normal = 0.0d0
	   end if
	end if
	if (z > 37.0d0) then
	   cdf_normal = 1.0d0
	   return
	else if (z < -37.0d0) then
	   cdf_normal = 0.0d0
	   return
	end if
	zabs=zabs/std
	arg = -0.5d0*zabs*zabs
	logpdf = -logroot2pi - log(std) + arg
	if (zabs < cutoff) then
	   p = arg + log(((((((p6*zabs + p5)*zabs + p4)*zabs + p3)*zabs + &
			p2)*zabs + p1)*zabs + p0)) - log((((((((q7*zabs + q6)*zabs + &
			q5)*zabs + q4)*zabs + q3)*zabs + q2)*zabs + q1)*zabs + &
			q0))
	else
	   p = logpdf - log((zabs + 1.0d0/(zabs + 2.0d0/(zabs + 3.0d0/(zabs + 4.0d0/ &
			(zabs + 0.65d0))))))
	end if
	p = exp(p)
	if (z < 0.0d0) then
	   cdf_normal=p
	   return
	else
	   cdf_normal = 1.0d0 - p
	   return
	end if
	return

end function cdf_normal


!> This is a simple birth probability function using standard normal cdf to deliver probabilities. Contraception use determines the intercept
!> of the model
function bprob(choice, variables, parameters)
	implicit none
	integer choice 				!< the choice of contraception, 1: use it, 0: no use
	real(dble) variables(Bsize) 	!< whatever determines the birth probability other than contraception SIZE=Bsize
	real(dble) parameters(Bsize+2) 	!< parameters of the birth SIZE=Bsize + 2 
	real(dble) bprob
	! the first two parameters are the ones for the intercept.
	real(dble) indeks
	indeks=sum(variables*parameters(3:Bsize + 2)) + parameters(1) + parameters(2)*choice
	bprob=cdf_normal(indeks,0.0d0,1.0d0)
end function bprob

!> Exogenous birth probability calculator, for the simplified model.
function bprobexo(variables,parameters)
	implicit none
	real(dble) variables(Bsizeexo) 		! variables of the birth prob.
	real(dble) parameters(Bsizexo + 1) 	! parameters of the model, intercept is the last item
	real(dble) bprobexo
	real(dble) indeks
	indeks = sum(variables*parameters(1:Bsizexo)) + parameters(Bsizexo + 1)
	bprobexo = cdf_normal(indeks,0.0d0,1.0d0)
end function bprobexo
 



! --------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
!                                 OLDER STUFF
!---------------------------------------------------------------------------------
!_________________________________________________________________________________

!#################----PICK AN State Space Point----#########################

!!> aux function to get cumulative distributions, helps: randmem
!function cumx(x)
!implicit none
	!real(dble) x(:)				!< input vector x
	!real(dble) cumx(size(x))	!<cumulative sum of the vector x

	!integer i
	!real(dble) sumto
	!sumto=0
	!do i=1,size(x)
		!sumto=sumto+x(i)
		!cumx(i)=sumto
	!end do
!end function cumx

!!> This is a function to randomly pick an element from array x
!!> where x's has the pdf in w (sum(w)=1). This function
!!> requires a uniform random,u, draw to be sent to it.
!!> it uses inverse cdf method. This is different than randpick functions
!!> because this one allows non uniform discrete distributions
!function randmem(x,w,u)
!implicit none
	!real(dble) randmem		 	!< output
	!real(dble) x(:)				!< input: array to drawn from.
	!real(dble) w(size(x))		!< distribution function
	!real(dble) u				!< uniform (0,1) draw

	!real(dble) cw(size(w)+1)
	!integer i, pick

	!cw=(/0.0d0,cumx(w)/)
	!do i=1,size(x)
		!if ((cw(i)<=u).AND.(cw(i+1)>u))  pick=i
	!end do
	!randmem=x(pick)
!end function randmem

!!> integer version of randmem
!function Irandmem(x,w,u)
!implicit none
	!integer Irandmem		 	!< output
	!integer x(:)				!< input: array to drawn from.
	!real(dble) w(size(x))		!< distribution function
	!real(dble) u				!< uniform (0,1) draw

	!real(dble) cw(size(w)+1)
	!integer i, pick

	!cw=(/0.0d0,cumx(w)/)
	!do i=1,size(x)
		!if ((cw(i)<=u).AND.(cw(i+1)>u))  pick=i
	!end do
	!Irandmem=x(pick)
!end function Irandmem

!!> This is a function to randomly pick a value from the
!!> range (rangemin, rangemax). a U(0,1) random number needs
!!> to be provided.
!function randpick(rangemin, rangemax, u)
!implicit none
	!real(dble) randpick !>output scalar
	!real(dble) rangemin, rangemax !>  the interval
	!real(dble) u			! U(0,1) draw

	!randpick=(rangemin+(rangemax-rangemin+1.0d0)*u)
!end function randpick


!!>  This is same as randpcik, only integer valued
!function Irandpick(rangemin, rangemax, u)
!implicit none
	!integer Irandpick !>output scalar
	!integer rangemin, rangemax !>  the interval
	!real(dble) u			! U(0,1) draw

	!Irandpick=floor(rangemin+(rangemax-rangemin+1)*u)
!end function Irandpick

!!> randomly pick a !FEASIBLE! state space point for families with two kids
!!> For this, I randomly pick a history for h and x up to time t and other state variables
!!> And create a
!subroutine pickSS(pick, t, parvecW, parvecH, sdw,sdh, omega4, allru)
!implicit none
	!real(dble), intent(out):: pick(Ssize)	!< output: S array of state space points
	!integer, intent(in):: t					!< time period as of now, also age of kid1
	!real(dble), intent(in):: parvecH(:), parvecW(:) !< parameter vectors for wages
	!real(dble), intent(in):: omega4(:)			!< type vector
	!real(dble), intent(in):: sdw			!< sd of wife's wage shock
	!real(dble), intent(in):: sdh			!< husband wage shock sd. parameter(sqrt(Sigma(5,5)))
	!real(dble), intent(in)::allru(:)		!< a vector of U(0,1) random numbers. size=10+(t-1)x2

	!! locals
	!! state space elements

	!integer age2, schoolingm, agem
	!real(dble) hprevious, afqt, llms
	!integer ageh, schoolingh
	!integer delta
	!real(dble) Iold1, Iold2, Iyoung1, Iyoung2, experience,y1, y2

	!integer i

	!!integer deltavec(deltamax-deltamin+1)	! potential deltas
	!real(dble) wdelta(deltamax-deltamin+1)	! weights of deltas
	!real(dble) hhist(t-1), xhist(t-1)
	!real(dble) logw0, yhmin, yhmax, yh2min,yh2max,ymin, ymax,y2min, y2max
	!integer ft,pt,cumexp


	!!**************************************d
	!! 2. ---draw delta=age2----
	!! I could potentially allow the picks determined by data, for now equal weights:
	!!deltavec=(/(i,i=deltamin, deltamax)/)
	!!	wdelta= (1.0d0/(deltamax-deltamin+1))

	!delta=Irandpick(deltamin, min(deltamax,t), allru(1))
	!age2=t-delta+1

	!!*****************
	!! 3. Pick schooling. Because we don't allow for women be in school when they have their baby.
	!! schooling will determine the minumum age of the mother.
	!schoolingm=Irandpick(9*(1/period),18*(1/period),allru(2))

	!! 4. pick Mom's age!! min= fertminage or school finishing age, whichever is bigger.
	!! max= fertmaxage+age2 "mom can't have a baby after fertmaxage"
	!agem=Irandpick(max(fertminage,(1/period)*7+schoolingm)+t,fertmaxage+age2, allru(3))
	!! so age0m=agem-t

	!!5. pick hprevious, schoolingh, ageh, afqt, llms

	!schoolingh=Irandpick(9*(1/period),18*(1/period),allru(5))
	!ageh=Irandpick(agem-3,agem+15,allru(6))*(1/period)			! sexist but practical interim solution for dad's age.
	!afqt=allru(7)
	!llms=allru(8)
	!!*********************************
	!!6. Now: More complicated stuff
	!! ---------------first way: generate a random history-------------------
	!! generate purely random histories
	!hhist=(floor(allru(11:(9+t))*3.0d0))/2.0d0			! size: t-1

	!hprevious=hhist(t-1)
	!xhist=floor(allru(10+t:(10+2*t-2))*11)/10.0d0		! size: t-1, but will not use all
	!experience=sum(hhist)

	!! new way: run a loop to create investment levels
	!! initialize investments
	!Iold1=0
	!Iold2=0
	!Iyoung1=0
	!Iyoung2=0

	!do i=1,(t-1)
		!! FIRST KID
		!! --------- second period of childhood---------------
		!if (i>=abar) then
			!Iold1=Iold1			! first period investment stays the same
			!! second period
			!if (i<delta) then	! still one child
				!Iold2=Iold2+1-0.5*hhist(i)
			!elseif (i>=astar) then	! too old to develop
				!Iold2=Iold2
			!else				! with sibling
				!Iold2=Iold2+(1-0.5*hhist(i))*xhist(i)
			!end if
		!else
		!! ******** first period of the childhood***********
			!Iold2=Iold2				! second period investment does not move.
			!if (i<delta) then
				!Iold1=Iold1+1-0.5*hhist(i)
			!else
				!Iold1=Iold1+(1-0.5*hhist(i))*xhist(i)
			!end if
		!end if
		!! SECOND KID
		!! --------- second period of childhood---------------

		!if (i>=delta) then	! ARE YOU EVEN ALIVE???
			!if (i>=delta+abar-1) then
				!Iyoung1=Iyoung1
				!if (i>=astar+delta-1) then	! second kid too old
					!Iyoung2=Iyoung2
				!elseif (i<astar) then	 ! first kid still getting time
					!Iyoung2=Iyoung2+hhist(i)*(1-xhist(i))
				!else 					! first kid is too old, second one is not.
					!Iyoung2=Iyoung2+hhist(i)
				!end if
			!else
			!! --------- first period of childhood---------------
				!Iyoung2=Iyoung2
				!if (i<astar) then	! still have a brother at home
					!Iyoung1=Iyoung1+(1-0.5*hhist(i))*(1-xhist(i))
				!else 				! second gets all the attention
					!Iyoung1=Iyoung1+1-0.5*hhist(i)
				!end if
			!end if
		!else 			! too soon, there is only one kid, so stay at zero.
			!Iyoung1=Iyoung1
			!Iyoung2=Iyoung2
		!end if

	!end do

!!	--------------OLD BUT WRONG WAY----------------------
!!	! first child's inputs
!!	! the relationship between delta and abar is important: can accumulate Iold only upto abar-1. And I2 will only change up to astar.

!!	if (delta<=abar) then
!!		if (t>abar) then
!!			Iold1= sum(1-0.5*hhist(1:delta-1))+sum((1-0.5*hhist(delta:(abar-1)))*xhist(delta:(abar-1)))
!!			Iold2= sum((1-0.5*hhist(abar:min(t-1,astar-1)))*xhist(abar:min(t-1,astar-1)))
!!		else
!!			Iold1=sum(1-0.5*hhist(1:delta-1))+sum((1-0.5*hhist(delta:(t-1)))*xhist(delta:(t-1)))
!!			Iold2=0
!!			if (t==delta) Iold1=sum(1-0.5*hhist(1:delta-1))
!!		end if
!!	else   ! if delta>abar, automatically t>abar, so no need for extra branching for the location of delta wrt t and abar. t=>delta>
!!		Iold1=sum(1-0.5d0*hhist(1:abar-1))
!!		Iold2=sum(1-0.5d0*hhist(abar:delta-1))+sum((1-0.5d0*hhist(delta:min(astar-1,t-1)))*xhist(delta:min(astar-1,t-1)))
!!	end if
!!	! second child's inputs
!!	Iyoung1=sum(1.0d0-0.5d0*hhist(delta:min(t-1,(delta+abar-2)))*(1.0d0-xhist(delta:min(t-1,(delta+abar-2)))))


!!	Iyoung2=sum((1-0.5d0*hhist(delta+abar-1:t-1))*(1-xhist(delta+abar-1:t-1)))
	!! 7 income: with given type and other characteristics, we can figure out the lower and upper boundaries for Y1 and Y2
	!! lower boundaries will come from the case that the history of the person realized and she and her husband receive the worst possible shocks (3xsd)

	!! a.husband wage: receive income from a exogenous process that depends on age and schooling. min and max determined by this and the exogenous shock.
	!! for the minimum: use -3*sd, and +3sd for the maximum
	!yhmin=0
	!yhmax=0
	!yh2min=0
	!yh2max=0
	!do i=ageh-t+1,ageh-1
		!yhmin=yhmin+exp(parvecH(1)+parvecH(2)*schoolingh+parvecH(3)*i+parvecH(4)*i*i-3*sdh)
		!yhmax=yhmax+exp(parvecH(1)+parvecH(2)*schoolingh+parvecH(3)*i+parvecH(4)*i*i+3*sdh)
		!if (i>=ageh-(t-delta)) then
			!yh2min=yhmin+exp(parvecH(1)+parvecH(2)*schoolingh+parvecH(3)*i+parvecH(4)*i*i-3*sdh)
			!yh2max=yhmax+exp(parvecH(1)+parvecH(2)*schoolingh+parvecH(3)*i+parvecH(4)*i*i+3*sdh)
		!end if
	!end do
	!!b. mom's income. more complicated as it depends on the hhist pick. But the idea is the same also with possible llms.
	!logw0=omega4(3)+parvecW(1)*schoolingm+parvecW(2)*afqt+parvecW(3)*agem+parvecW(4)*agem**2
	!ymin=0.0d0
	!ymax=0.0d0
	!y2min=0.0d0
	!y2max=0.0d0
	!cumexp=0
	!pt=0		!previous experience pt
	!ft=0		!previous experience ft
	!do i=1,t-1
		!ymin=ymin+(exp(logw0+parvecW(5)*cumexp+parvecW(6)*ft+parvecW(7)*pt+parvecW(8)*0.0d0-3.0d0*sdw))*hhist(i)
		!ymax=ymax+(exp(logw0+parvecW(5)*cumexp+parvecW(6)*ft+parvecW(7)*pt+parvecW(8)*1.0d0+3.0d0*sdw))*hhist(i)
		!if (i>=delta) then
			!y2min=y2min+(exp(logw0+parvecW(5)*cumexp+parvecW(6)*ft+parvecW(7)*pt+parvecW(8)*0.0d0-3.0d0*sdw))*hhist(i)
			!y2max=y2max+(exp(logw0+parvecW(5)*cumexp+parvecW(6)*ft+parvecW(7)*pt+parvecW(8)*1.0d0+3.0d0*sdw))*hhist(i)
		!end if

		!cumexp=cumexp+hhist(i)
		 !if (hhist(i)==1) then
			 !ft=1
			 !pt=0
		  !elseif  (hhist(i)==0.5) then
			  !pt=1
			  !ft=0
		 !else
			 !ft=0
			 !pt=0
		 !end if
	!end do
	!ymin=ymin+yhmin
	!ymax=ymax+yhmax
	!y2min=y2min+yh2min
	!y2max=y2max+yh2max
	!y1=randpick(ymin,ymax,allru(9))
	!y2=randpick(y2min,y2max,allru(10))

!pick=(/Iold1,Iold2,Iyoung1,Iyoung2,y1,y2,experience,1.0d0*t,1.0d0*age2,1.0d0*agem,1.0d0*hprevious,llms,ageh*1.0d0,1.0d0*schoolingm,afqt,(agem-t+1)*1.0d0,1.0d0*schoolingh,omega4(1),omega4(2),omega4(3)/)
!end subroutine pickSS

end module randomgen

