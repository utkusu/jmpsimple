!> MAIN PROGRAM FOR THE SOLUTION AND ESTIMATION OF JMP MODEL
program act
use global
use randomgen
use emax
use opt
USE IFPORT ! for intel fortran only
implicit none
include 'mpif.h'

!real(dble) coef(Gsize+1)
!real(dble) coeaf(Gsize+1)
!real(dble) mutype(3)

! mpi stuff
integer nproc, rank, ier, status(MPI_STATUS_SIZE),i, sender, position, number_sent, number_received, order,tag, rorder

! parameters of the model
real(dble) parA(12),parU(7),parW(7),parH(6),beta,sigma1(shocksize1,shocksize1),parB(Bsizeexo+1)
real(dble) ctype, mtype, atype, a1type(na1type)
real(dble) condprob
real(dble) packed(parAsize+parWsize+parUsize+parHsize+1+(shocksize1+shocksize1*(shocksize1-1)/2)+Bsizeexo+1+4)
integer npack
real(dble) ftypemat(2,na1type*(deltamax-deltamin+1))
!real(dble) ftypematoc(5,nttypesoc)
!real(dble) ocp(nctype+Bsize+2+nctype*nmtype)

real(dble) solw(Gsize+1, nperiods-deltamin+2) !< collects solution coef.
real(dble) solwall(Gsize+1, nperiods-deltamin+2, deltamax-deltamin+1,nttypes) !< collects solution coef for all types
real(dble) solv(Gsizeoc+1,nperiods)
real(dble) solvall(Gsizeoc+1,nperiods, nttypes)
integer k,j,l,m
real(dble) wcoef(Gsize+1, nperiods-deltamin+2, deltamax-deltamin+1,nctype) 		! for the vsolver trial
real(dble) start,endtime, vtime
real(dble) rho
integer nftype
integer order3,order4
!real(dble) omega1(4),omega2(5),omega3(4),eps(Nmc,shocksize1)
!real(dble) ftype(5)
!real(dble) solw(Gsize,nperiods-deltamin+2,deltamax-deltamin+2)
!real(dble) wcoefficients(Gsize+1,nperiods-deltamin+2,deltamax-deltamin+2,2)
!real(dble) solv(Gsize+1,nperiods)
!real(dble) typevec(2)
!real(dble) typeprob(2)
rho=1.0d0
npack=parAsize+parWsize+parUsize+parHsize+1+(shocksize1+shocksize1*(shocksize1-1)/2)+Bsizeexo+1+4
nftype=na1type*(deltamax-deltamin+1)


call MPI_INIT(ier)
call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ier)
call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ier)

if (rank==0) then
	!------------------------------ MASTER ------------------------------------
	start=MPI_Wtime()
	call initparam(parA,parU,parW, parH, beta, sigma1,parB,ctype, mtype, atype, a1type,condprob)
	call utku_pack(packed,npack,parA,parU,parW, parH, beta, sigma1,parB,ctype,mtype, atype, condprob)
	call MPI_BCAST(packed, npack, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ier)
	! pack unobserved types in to matrices
    call type_packsimple(ftypemat,a1type)
    call MPI_BCAST(ftypemat, 2*nftype, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ier)
 		print*, 'MASTER: broadcasted everything, starting stuff'
 	! send everyone one index for which to calculate the wmatrix
 	
 	! -----------------------------Master: Solution of W--------------------------
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
		print*, 'MASTER: received solw for order=',rorder,', total number received=', number_received+1, solw(:,10), '--------------------------------------was for', order3, order4
		number_received=number_received+1
		if (number_sent<nftype) then
			order=number_sent+1
			call MPI_SEND(order, 1, MPI_INTEGER,sender,2, MPI_COMM_WORLD, ier)
			print*, 'MASTER sent order', order,'to',sender
			number_sent=number_sent+1
		else 
			call MPI_SEND(order, 1, MPI_INTEGER,sender ,0, MPI_COMM_WORLD, ier)
			print*, 'MASTER: finalizing two child solution for worker',sender
			
		end if
	end do
print*,'--------------------------------------------------------------------------------'
print*, '----MASTER: FINISHED THE SOLUTION OF THE WCOEF, NOW SWITCHING TO THE V-----'
print*,'--------------------------------------------------------------------------------'
	 !first things first: workers need to the results of the last round
	vtime=MPI_Wtime()
	call MPI_BCAST(solwall, (Gsize+1)*(nperiods-deltamin+2)*(deltamax-deltamin+1)*nttypes, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ier)
	print*, 'master broadcasted last stage results to workers successfully'
	 !sending the first group
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
		solvall(:,:,rorder)=solv
		number_received=number_received+1
		call MPI_SEND(order, 1, MPI_INTEGER,sender ,10, MPI_COMM_WORLD, ier)
		!if (number_sent<nttypesoc) then
			!order=number_sent+1
			!call MPI_SEND(order, 1, MPI_INTEGER,sender,11, MPI_COMM_WORLD, ier)
			!number_sent=number_sent+1
		!else 
			!call MPI_SEND(order, 1, MPI_INTEGER,sender ,10, MPI_COMM_WORLD, ier)
			!print *, 'MASTERV:finalizing sender',sender
		!end if
	end do
	
	 !write vsolve matrix to a (Gsize+1,nperiods)x nttypesoc matrix to print it
	 !out.

	!endtime=MPI_Wtime()	
		!print*, '----------------------- SUMMARY ---------------------'
		!print*, 'w coef calc took', vtime-start
		!print*, 'v coef calc took', endtime-vtime
		!print*, 'whole thing took', endtime-start
		
	open(66,file='vcoef.txt')
	do l=1,2
		write(66,*) "-------------- type order=",l,"------------" 
		write(66,60) ((solvall(j,k,l),k=1,nperiods),j=1,Gsizeoc+1)
	end do
	60 format(22F46.9)
	close(66)
	open(77,file='wcoef.txt')
	do l=1,2
		write(77,*) "-------------- type order=",l,"------------" 
		do m=1,deltamax-deltamin+1
			write(77,*) "########   delta type=",m,"------------"
			write(77,70) ((solwall(j,k,m,l),k=1,nperiods-deltamin+2),j=1,Gsize+1)
		end do
	end do
	70 format(22F46.9)
	close(77)
	
else
	!--------------------------- WORKERS -------------------------------------------
	call MPI_BCAST(packed, npack, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ier)
	call utku_unpack(packed,npack,(/12,7,7,6,1,15,(Bsizeexo+1),1,1,1,1/),parA,parU,parW, parH, beta, sigma1,parB,ctype,mtype,atype,condprob)
    call MPI_BCAST(ftypemat, 2*nftype, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ier)	
	! --------SOLVING FOR W COEFFICIENTS-------------------
	do
		call MPI_RECV(order, 1, MPI_INTEGER, 0, MPI_ANY_TAG, MPI_COMM_WORLD, status, ier)
		tag=status(MPI_TAG)
		if (tag>0) then
			call wsolver(solw,ftypemat(2,order),(/ctype,ctype,mtype,atype,ftypemat(1,order)/),parA,parW,parH(5:6),parU, beta,sigma1,rho)
			!print*, ftypemat(1,:)
			!print*,  '___worker',rank,'order=',order, 'calculated for', ftypemat(2,order), ftypemat(1,order)
			if (order==12) print*, solw(:,8)
			rorder=order
			call MPI_SEND(rorder, 1, MPI_INTEGER, 0, 1, MPI_COMM_WORLD, ier)
			call MPI_SEND(solw,(Gsize+1)*(nperiods-deltamin+2),MPI_DOUBLE_PRECISION,0,1,MPI_COMM_WORLD,ier)
		else
			EXIT
		end if
	end do
	! -------SOLVING FOR V COEFFICIENTS-------------
	! receive all the coefficients from the master-from the w stage
	call MPI_BCAST(solwall, (Gsize+1)*(nperiods-deltamin+2)*(deltamax-deltamin+1)*nttypes, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ier)
	
	if(rank<3) then
	do
		call MPI_RECV(order, 1, MPI_INTEGER, 0, MPI_ANY_TAG, MPI_COMM_WORLD, status, ier)
		tag=status(MPI_TAG)
		if(tag>10) then
			! get the correct mother type from solwall
			wcoef(:,:,:,1)=solwall(:,:,:,order)
			call vsolver(solv,(/ctype,ctype,mtype,atype,a1type(order)/), parA,parW,parH(5:6),parU,parB, beta,sigma1,wcoef,(/ctype/),(/1.0d0/),rho)
			rorder=order   ! returning the order
			call MPI_SEND(rorder, 1, MPI_INTEGER, 0, 1, MPI_COMM_WORLD, ier)
			call MPI_SEND(solv,(Gsizeoc+1)*nperiods,MPI_DOUBLE_PRECISION,0,1,MPI_COMM_WORLD,ier)
			print*, 'worker',rank, 'sent back order', rorder
		else
			EXIT
			print*, 'master told me to shutdown'
		end if
	end do
	end if
end if	
	call MPI_FINALIZE(ier)


end program act
