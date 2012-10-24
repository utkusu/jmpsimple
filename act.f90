!> MAIN PROGRAM FOR THE SOLUTION AND ESTIMATION OF JMP MODEL (SIMPLE VERSION)
program act
use global
use randomgen
use emax
use opt
USE IFPORT ! for intel fortran only
implicit none
include 'mpif.h'
real(dble) parameters(parsize), dist, targetvec(MomentSize), weightmat(MomentSize,MomentSize)


call MPI_INIT(ier)
call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ier)
call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ier)

call readdata()	
call readsomeparam()
parameters=0.01d0
targetvec=0.0d0
weightmat=1.0d0

call distance(dist, parameters, targetvec, weightmat)

print*,dist




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
		real(dble) pa1type(na1type-1) ! probability of alpha1, last one is 1-sum the rest 
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
		real(dble) start,endtime, vtime
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
			
		
		! break up parameters for two types
		part1_parA=parameters(1:3)
		parUpart= parameters(5:10)
		a1type=parameters(11:12)
		pa1type=parameters(13)
		sigma(1,1)=parameters(14); sigma(2,2)=parameters(15); sigma(3,3)=parameters(16)
		sigma(2,1)=parameters(17); sigma(1,2)=parameters(17)
		sigma(3,1)=parameters(18); sigma(1,3)=parameters(18)
		sigma(3,2)=parameters(19); sigma(2,3)=parameters(19)
		beta=parameters(20)	

		parA=(/part1_parA,gpart2_parA/)
		parU=(/parUpart(1:2), 1.0d0, parUpart(3:6)/) ! put a number in there,a1 will come from a1 distribution but wsolver handles that
		
		! initialize stuff 
		nftype=na1type*(deltamax-deltamin+1)
    	call type_packsimple(ftypemat,a1type)


		! ---------------------------------SOLUTION ALGORITHM---------------------------------------
		
		! ---------------------------------------MASTER----------------------------------
		if (rank==0) then
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
					print*, 'MASTER: finalizing two child solution for worker',sender
					
				end if
			end do
			
			print*,'--------------------------------------------------------------------------------'
			print*, '----MASTER: FINISHED THE SOLUTION OF THE WCOEF, NOW SWITCHING TO THE V-----'
			print*,'--------------------------------------------------------------------------------'
		
			call MPI_BCAST(solwall, (Gsize+1)*(nperiods-deltamin+2)*(deltamax-deltamin+1)*nttypes, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ier)
			print*, 'master broadcasted last stage results to workers successfully'
		
			! ------------MASTER: Solution for V----------------
			number_sent=0
			do i=1,min(nproc-1,nttypes)
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
					call wsolver(solw,ftypemat(2,order),(/gctype,gctype,gmtype,gatype,ftypemat(1,order)/),parA,gparW,gparH(5:6),parU, beta,sigma,grho)
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
			if(rank<3) then
			do
				call MPI_RECV(order, 1, MPI_INTEGER, 0, MPI_ANY_TAG, MPI_COMM_WORLD, status, ier)
				tag=status(MPI_TAG)
				if(tag>10) then
					! get the correct mother type from solwall
					wcoef(:,:,:,1)=solwall(:,:,:,order)
					call vsolver(solv,(/gctype,gctype,gmtype,gatype,a1type(order)/), parA,gparW,gparH(5:6),parU,gparBmat, beta,sigma,wcoef,(/gctype/),(/1.0d0/),grho)
					rorder=order   ! returning the order
					call MPI_SEND(rorder, 1, MPI_INTEGER, 0, 1, MPI_COMM_WORLD, ier)
					call MPI_SEND(solv,(Gsizeoc+1)*nperiods,MPI_DOUBLE_PRECISION,0,1,MPI_COMM_WORLD,ier)
					print*, 'worker',rank, 'sent back order', rorder
				else
					EXIT
				end if
			end do
			end if
		end if
		
		print*, 'MODEL SOLVED, MOVING ONTO SIMULATIONS'

		! now that we have the solwall and solvall (the coeffiecients for the interpolation, it is time to simulate data from these.
		! it is only master who should do this.
		if (rank==0) then

			do id=1,SampleSize
				call simhist(SS,outcomes,testoutcomes, choices, xchoices,birthhist,smchoices, smexperience, smAs, smtestoutcomes,gomega3data(:,id),(/gctype,gmtype/),parA,parU,gparW,gparH,beta,sigma,a1type,pa1type,gparBmat,solvall,solwall,llmsmat(:,id),id,grho,glambdas,gsigmaetas, gsmpar)
				SScollect(:,:,:,id)=SS
				birthhistcollect(:,id)=birthhist
				smchoicescollect(:,:,:,id)=smchoices
				smexperiencecollect(:,:,id)=smexperience
				smtestoutcomescollect(:,:,:,id)=smtestoutcomes
			end do
			call moments(momentvec, SScollect, smtestoutcomescollect,  birthhistcollect, smchoicescollect, smexperiencecollect, gomega3data,glfpperiods, gexpperiods, gtsperiods, gidmat)
			
			difft(1,:)=momentvec-targetvec
			diff(:,1)=momentvec-targetvec
			middlestep=matmul(difft,weightmat)
			laststep=matmul(middlestep,diff)
			dist=laststep(1,1)
		end if
	end subroutine distance


	!> make the distance function and its numerical derivatives a function to be evaluated by the optimizer
	subroutine objfunc(val, n, xvec, grad, need_gradient, fmat)
	implicit none
	integer n, need_gradient
	real(dble) val, xvec(n), grad(n), fmat
	

	! locals

	! separate out the x  
	val=0.d0
end subroutine objfunc





end program act
