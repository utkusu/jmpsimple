!> MAIN PROGRAM FOR THE SOLUTION AND ESTIMATION OF JMP MODEL (SIMPLE VERSION)
program act
use global
use randomgen
use emax
use optu
USE IFPORT ! for intel fortran only
implicit none
!external objfunc
include 'mpif.h'
include 'nlopt.f'
real(dble)  dist, fakeg(parsize)
! nlopt stuff
integer*8 opt
integer ires
real(dble) minf
real(dble)  fmat
opt=0; fmat=1.0d0

call MPI_INIT(ier)
call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ier)
call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ier)

if (rank==0) then
	call readdata() ! read some data	
	call readsomeparam() ! read some external parameters
	! read these
	call setoptimstuff() ! read and set initial values and boundaries
	call setii()  ! set target vec and weighting matrix for indirect inference
	! set lower - upper bound
	itercounter=1; evaliter=1
end if
! broadcast all these
call MPI_BCAST(gidmat, SampleSize*idmatsize, MPI_INTEGER, 0, MPI_COMM_WORLD, ier)
call MPI_BCAST(gomega3data, o3size*SampleSize, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ier)
call MPI_BCAST(llmsmat, nperiods*SampleSize, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ier)
call MPI_BCAST(glambdas, Ntestage, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ier)
call MPI_BCAST(gsigmaetas, 2*Ntestage, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ier)
call MPI_BCAST(gparBmat, (Bsizeexo+1)*nfert, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ier)
call MPI_BCAST(parameters, parsize, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ier)
call MPI_BCAST(lb, parsize, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ier)
call MPI_BCAST(ub, parsize, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ier)
call MPI_BCAST(targetvec, MomentSize, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ier)
call MPI_BCAST(weightmat, MomentSize*MomentSize, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ier)




call nlo_create(opt, NLOPT_LD_MMA,parsize)
call nlo_set_lower_bounds(ires, opt, lb)
call nlo_set_upper_bounds(ires, opt, ub)
call nlo_set_min_objective(ires, opt, objfunc, fmat)
call nlo_set_maxeval(ires,opt,40)
!call nlo_set_xtol_rel(ires, opt, 0.0001d0)

call nlo_optimize(ires, opt, parameters, minf)
if (rank==0) then 
	print*, ires
	print*, minf
end if
call nlo_destroy(opt)	



call MPI_FINALIZE(ier)


contains




	subroutine distance(dist,parameters,targetvec,weightmat)
		implicit none
		include 'mpif.h'
		real(dble), intent(in)::parameters(parsize)
		real(dble), intent(in)::weightmat(MomentSize, MomentSize) 
		real(dble), intent(in) :: targetvec(MomentSize)
		real(dble), intent(out) :: dist
		

		! locals
		! first to break up the parameters
		real(dble) part1_parA(3) 	! A production function excluding intercept
		real(dble) parUpart(parUsize-1) ! parameters of utility except alpha1
		real(dble) a1type(na1type)  ! alpha1 typevec
		real(dble) pa1typepart(na1type-1) ! probability of alpha1, last one is 1-sum the rest 
		real(dble) pa1type(na1type) ! probability of alpha1, last one is 1-sum the rest 
		real(dble) sigma(shocksize1,shocksize1) ! variance matrix of shocks. first diagonal then off diagonal. go column wise: (2,1) (3,1) (3,2) for 3
		real(dble) beta

		integer id


		! the stuff needed between steps
		integer nftype
		integer i,j,k,l,m
		real(dble) solw(Gsize+1, nperiods-deltamin+2) ! collects solution coef.
		real(dble) solwall(Gsize+1, nperiods-deltamin+2, deltamax-deltamin+1,nttypes) ! collects solution coef for all types
		real(dble) solv(Gsizeoc+1,nperiods)
		real(dble) solvall(Gsizeoc+1,nperiods, nttypes)
		real(dble) wcoef(Gsize+1, nperiods-deltamin+2, deltamax-deltamin+1,nctype) 		! for the vsolver trial
		real(dble) ftypemat(2,na1type*(deltamax-deltamin+1))
		real(dble) starttime,endtime, soltime
		real(dble) parA(parAsize)
		real(dble) parU(parUsize)
		
		! siumlation stuff
		real(dble) SScollect(6,nperiods,Npaths,SampleSize) 		! Simulated State space (A1,A2,E,age1,age2,agem)xnperiods,Npaths	
		!real(dble) outcomescollect(2,nperiods,Npaths,SampleSize) ! outcomes: wages of the father and mother.
		!integer choicescollect(1,nperiods,Npaths,SampleSize) 	! choices : the choice history of h.
		integer birthhistcollect(Npaths,SampleSize)  		    ! the birth timing vec, 0 if one child throughout.
		real(dble) omega3data(o3size, SampleSize) 	! holds the omega3 data for Sample people
		real(dble) smtestoutcomescollect(4,Ntestage,Npaths,SampleSize)  ! holds test scores for two children 
		real(dble) smchoicescollect(3,nperiods,Npaths,SampleSize)
		real(dble) smexperiencecollect(nperiods,Npaths,SampleSize)
		integer idmat(SampleSize,MomentSize) 	   		!indicates the which sample units are in the jth columnth moment
					
		! one use ones 
		real(dble) SS(6,nperiods,Npaths) 		!< Simulated State space (A1,A2,E,age1,age2,agem)xnperiods,Npaths	
		real(dble) outcomes(2,nperiods,Npaths) !< outcomes: wages of the father and mother.
		integer choices(1,nperiods,Npaths) 	!< choices : the choice history of h.
		integer xchoices(1,nperiods,Npaths) 	!< not in the data, but I will keep an history of the x choices as well.
		integer birthhist(Npaths) 			!< the birth timing vec, 0 if one child throughout.
		real(dble) testoutcomes(4,Ntestage, Npaths)
		real(dble) smchoices(3,nperiods, Npaths)
		real(dble) smexperience(nperiods, Npaths)
		real(dble) smAs(2,nperiods, Npaths)
		real(dble) smtestoutcomes(4,Ntestage, Npaths)
		real(dble) omega3(o3size) 
		
		!calculation
		real(dble) momentvec(MomentSize) 	!

		integer status(MPI_STATUS_SIZE)
		real(dble) difft(1,MomentSize), diff(MomentSize,1), middlestep(1,MomentSize), laststep(1,1)
		
		! delete me =why?
		real(dble) coef(Gsizeoc+1)
	
		! late addition: I forgot to parameterze mu_c
		real(dble) nctype

		! break up parameters for two types
		part1_parA=parameters(1:3)
		parUpart= parameters(4:9)
		a1type=parameters(10:11)
		pa1typepart=parameters(12)
		pa1type=(/pa1typepart,1-sum(pa1typepart)/)
		sigma(1,1)=parameters(13); sigma(2,2)=parameters(14); sigma(3,3)=parameters(15)
		sigma(2,1)=parameters(16); sigma(1,2)=parameters(16)
		sigma(3,1)=parameters(17); sigma(1,3)=parameters(17)
		sigma(3,2)=parameters(18); sigma(2,3)=parameters(18)
		beta=parameters(19)	
		! the following is supposed to come in a distribution in a model, but simple model does not allow heterogeneity.
		! I forgot about this intercept, it needs to estimated inside the model.
		nctype=parameters(20)


		parA=(/part1_parA,gpart2_parA/)
		parU=(/parUpart(1:2), 1.0d0, parUpart(3:6)/) ! put a number in there,a1 will come from a1 distribution but wsolver handles that
		
		! initialize stuff 
		nftype=na1type*(deltamax-deltamin+1)
    	call type_packsimple(ftypemat,a1type)
		! ---------------------------------SOLUTION ALGORITHM---------------------------------------
		
		! ---------------------------------------MASTER----------------------------------
		if (rank==0) then
			starttime=MPI_WTIME()
			print*, 'starting iteration', itercounter
		!------------------- MASTER: W SOLVER -----------------
			! send order numbers to workers
			number_sent=0
			do i=1,min(nproc-1,nftype)
				order=i
				call MPI_SEND(order, 1, MPI_INTEGER,i ,1, MPI_COMM_WORLD, ier)
				number_sent=number_sent+1
			end do 
			! start receiving the wmatrix and place them properly
			number_received=0
			do while (number_received<nftype)
				call MPI_RECV(rorder, 1, MPI_INTEGER, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, status, ier)
				sender=status(MPI_SOURCE)	
				call MPI_RECV(solw,(Gsize+1)*(nperiods-deltamin+2),MPI_DOUBLE_PRECISION,sender,MPI_ANY_TAG,MPI_COMM_WORLD,status,ier)
				order3=mod(rorder,6)
				if (mod(rorder,6)==0) order3=6
				order4=((rorder-1)/6)+1
				solwall(:,:,order3,order4)=solw
				number_received=number_received+1
				if (number_sent<nftype) then
					order=number_sent+1
					call MPI_SEND(order, 1, MPI_INTEGER,sender,2, MPI_COMM_WORLD, ier)
					number_sent=number_sent+1
				else 
					call MPI_SEND(order, 1, MPI_INTEGER,sender ,0, MPI_COMM_WORLD, ier)
					!print*, 'MASTER: finalizing two child solution for worker',sender
					
				end if
			end do
			
			!print*,'--------------------------------------------------------------------------------'
			!print*, '----MASTER: FINISHED THE SOLUTION OF THE WCOEF, NOW SWITCHING TO THE V-----'
			!print*,'--------------------------------------------------------------------------------'
			



			open(77,file='wcoef.txt')
			do l=1,2
				write(77,*) "-------------- type order=",l,"------------" 

				do m=1,deltamax-deltamin+1
					write(77,*) "########   delta type=",m,"------------"
					write(77,70) ((solwall(j,k,m,l),k=1,nperiods-deltamin+2),j=1,Gsize+1)
				end do
			end do
			70 format(22F26.9)
			close(77)

			
			call MPI_BCAST(solwall, (Gsize+1)*(nperiods-deltamin+2)*(deltamax-deltamin+1)*nttypes, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ier)
		
			! ------------MASTER: Solution for V----------------
			number_sent=0
			do i=1,min(nproc-1,na1type)
				order=i
				call MPI_SEND(order, 1, MPI_INTEGER,i ,11, MPI_COMM_WORLD, ier)
				number_sent=number_sent+1
			end do 
			!start receiving the vmatrix and place them properly
			number_received=0
			do while (number_received<nttypesoc)
				call MPI_RECV(rorder, 1, MPI_INTEGER, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, status, ier)
				sender=status(MPI_SOURCE)	
				call MPI_RECV(solv,(Gsizeoc+1)*nperiods,MPI_DOUBLE_PRECISION,sender,MPI_ANY_TAG,MPI_COMM_WORLD,status,ier)
				!print*,'Master recieved one child solution', solv
				solvall(:,:,rorder)=solv
				number_received=number_received+1
				! this version has two types, so it is always less than no processors, there is no second round sending.
				! this just is finalizing workers
				call MPI_SEND(order, 1, MPI_INTEGER,sender ,10, MPI_COMM_WORLD, ier)
			end do
	

		! ----------------------------------workers--------------------------------
		else
				! --------Workers: SOLVING FOR W COEFFICIENTS-------------------
			do
				call MPI_RECV(order, 1, MPI_INTEGER, 0, MPI_ANY_TAG, MPI_COMM_WORLD, status, ier)
				tag=status(MPI_TAG)
				if (tag>0) then
					call wsolver(solw,ftypemat(2,order),(/nctype,nctype,gmtype,parUpart(2),ftypemat(1,order)/),parA,gparW,gparH(5:6),parU, beta,sigma,grho)
					!print*, ftypemat(1,:)
					!print*,  '___worker',rank,'order=',order, 'calculated for', ftypemat(2,order), ftypemat(1,order)
					rorder=order
					call MPI_SEND(rorder, 1, MPI_INTEGER, 0, 1, MPI_COMM_WORLD, ier)
					call MPI_SEND(solw,(Gsize+1)*(nperiods-deltamin+2),MPI_DOUBLE_PRECISION,0,1,MPI_COMM_WORLD,ier)
				else
					EXIT
				end if
			end do

			
			! --------------Workers: SOLVING FOR V COEFFICIENTS-------------
			! receive all the coefficients from the master-from the w stage
			call MPI_BCAST(solwall, (Gsize+1)*(nperiods-deltamin+2)*(deltamax-deltamin+1)*nttypes, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ier)
			if( (rank==1) .OR. (rank==2) ) then
			do
				call MPI_RECV(order, 1, MPI_INTEGER, 0, MPI_ANY_TAG, MPI_COMM_WORLD, status, ier)
				tag=status(MPI_TAG)
				if(tag>10) then
					! get the correct mother type from solwall
					wcoef(:,:,:,1)=solwall(:,:,:,order)
					call vsolver(solv,(/nctype,nctype,gmtype,parUpart(2),a1type(order)/), parA,gparW,gparH(5:6),parU,gparBmat, beta,sigma,wcoef,(/nctype/),(/1.0d0/),grho)
					!print*, solv
					rorder=order   ! returning the order
					call MPI_SEND(rorder, 1, MPI_INTEGER, 0, 1, MPI_COMM_WORLD, ier)
					call MPI_SEND(solv,(Gsizeoc+1)*nperiods,MPI_DOUBLE_PRECISION,0,1,MPI_COMM_WORLD,ier)
					!print*, 'worker',rank, 'sent back order', rorder
				else
					EXIT
				end if
			end do
			end if
		end if
		
		! now that we have the solwall and solvall (the coeffiecients for the interpolation, it is time to simulate data from these.
		! it is only master who should do this.
		if (rank==0) then
			print*, 'MODEL SOLVED, MOVING ONTO SIMULATIONS'
			soltime=MPI_WTIME()
			open(66,file='vcoef.txt')
			do l=1,2
				write(66,*) "-------------- type order=",l,"------------" 
				write(66,60) ((solvall(j,k,l),k=1,nperiods),j=1,Gsizeoc+1)
			end do
			60 format(22F46.9)
			close(66)

			do id=1,SampleSize
				call simhist(SS,outcomes,testoutcomes, choices, xchoices,birthhist,smchoices, smexperience, smAs, smtestoutcomes,gomega3data(:,id),(/gctype,gmtype/),parA,parU,gparW,gparH,beta,sigma,a1type,pa1type,gparBmat,solvall,solwall,llmsmat(:,id),id,grho,glambdas,gsigmaetas, gsmpar)
				SScollect(:,:,:,id)=SS
				birthhistcollect(:,id)=birthhist
				smchoicescollect(:,:,:,id)=smchoices
				smexperiencecollect(:,:,id)=smexperience
				smtestoutcomescollect(:,:,:,id)=smtestoutcomes
				!print*, 'simulations complete for individual', id
			end do
			print*, 'Simulations Complete, Now Calculating Moments'
			call moments(momentvec, SScollect, smtestoutcomescollect,  birthhistcollect, smchoicescollect, smexperiencecollect, gomega3data,glfpperiods, gexpperiods, gtsperiods, gidmat)
			difft(1,:)=momentvec-targetvec
			diff(:,1)=momentvec-targetvec
			middlestep=matmul(difft,weightmat)
			laststep=matmul(middlestep,diff)
			dist=laststep(1,1)
			endtime=MPI_WTIME()
			print*, '-------------- ITERATION SUMMARY ----------------'
			print*, 'This is iteration', itercounter
			print*, 'Number of Evaluated Points', evaliter
			print*, 'Distance is', dist
			print*, 'solution calc took', soltime-starttime
			print*, 'simulation and moments', endtime-soltime
			print*, 'whole thing took', endtime-starttime
			print*, '-------------------------------------------------'
			itercounter=itercounter+1
		end if
	end subroutine distance


	!> make the distance function and its numerical derivatives a function to be evaluated by the optimizer
	subroutine objfunc(val, n, xvec, grad, need_gradient, fmat)
	implicit none
	integer n, need_gradient
	real(dble) val, xvec(n), grad(n), jacobian(n)
	real(dble) dist, fmat 
		call distance(dist, parameters, targetvec, weightmat)
		evaliter=evaliter+1
	val=dist
	if ( need_gradient .NE. 0 ) then
		call diffdistance(jacobian,parameters,targetvec, weightmat, dist)
		grad=jacobian
	end if
	end subroutine objfunc
	! takes the forward numerical derivative of the distance function, using already calculated f0


	subroutine diffdistance(jacobian, parameters,targetvec,weightmat,f0)
		implicit none
		real(dble) jacobian(parsize)
		real(dble) parameters(parsize), targetvec(MomentSize), weightmat(MomentSize,MomentSize)
		real(dble) f0
		real(dble) f1
		real(dble) h(parsize),hbar(parsize),xh(parsize), fakeparameters(parsize)
		integer i
		
		hbar=(0.00000001)*max(abs(parameters),1.0d0)
		xh=parameters+hbar
		h=xh-parameters
		do i=1,parsize
			fakeparameters=parameters; fakeparameters(i)=xh(i)
			call distance(f1,fakeparameters,targetvec,weightmat)
			jacobian(i)=(f1-f0)/h(i)
			if (rank==0) print*, 'jacobian in direction', i
		end do
	end subroutine diffdistance


end program act
