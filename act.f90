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
real(dble) buffer(1000)

! parameters of the model
real(dble) parA(12),parU(7),parW(7),parH(6),beta,sigma1(shocksize1,shocksize1),parB(Bsizeexo+1)
real(dble) ctype(nctype), mtype(nmtype), atype(natype), a1type(na1type)
!real(dble) pctype(nctype),pmtype(nmtype),patype(natype),pa1type(na1type)
real(dble) condprob
real(dble) packed(parAsize+parWsize+parUsize+parHsize+1+(shocksize1+shocksize1*(shocksize1-1)/2)+Bsizeexo+1+4)
integer npack
real(dble) ftypemat(5,nttypes)
real(dble) ftypematoc(5,nttypesoc)
real(dble) ocp(nctype+Bsize+2+nctype*nmtype)

real(dble) solw(Gsize+1, nperiods-deltamin+2, deltamax-deltamin+1) !< collects solution coef.
real(dble) solwall(Gsize+1, nperiods-deltamin+2, deltamax-deltamin+1,nttypes) !< collects solution coef for all types
real(dble) solv(Gsize+1,nperiods)
real(dble) solvall(Gsize+1,nperiods, nttypesoc)
integer momorder
integer k,j,l,m
real(dble) wcoef(Gsize+1, nperiods-deltamin+2, deltamax-deltamin+1,nctype) 		! for the vsolver trial
real(dble) start,endtime, vtime
real(dble) rho
!real(dble) omega1(4),omega2(5),omega3(4),eps(Nmc,shocksize1)
!real(dble) ftype(5)
!real(dble) solw(Gsize,nperiods-deltamin+2,deltamax-deltamin+2)
!real(dble) wcoefficients(Gsize+1,nperiods-deltamin+2,deltamax-deltamin+2,2)
!real(dble) solv(Gsize+1,nperiods)
!real(dble) typevec(2)
!real(dble) typeprob(2)
rho=1.0d0
npack=parAsize+parWsize+parUsize+parHsize+1+(shocksize1+shocksize1*(shocksize1-1)/2)+Bsizeexo+1+4



call MPI_INIT(ier)
call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ier)
call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ier)

if (rank==0) then
	!------------------------------ MASTER ------------------------------------
	start=MPI_Wtime()
	call initparam(parA,parU,parW, parH, parB, beta, sigma1, sigma2,parB,ctype, mtype, atype, a1type,condprob)
	call utku_pack(packed,npack,parA,parU,parW, parH, beta, sigma1,ctype,mtype, atype)
	call MPI_BCAST(packed, npack, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ier)
	! pack unobserved types in to matrices
    call type_packsimple(ftypemat,a1type)
    call MPI_BCAST(ftypemat, 2*na1type*(deltamax-deltamin+1), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ier)
 		print*, 'MASTER: broadcasted everything, starting stuff'
 	! send everyone one index for which to calculate the wmatrix
 	number_sent=0
	do i=1,nproc-1
		order=i
		call MPI_SEND(order, 1, MPI_INTEGER,i ,1, MPI_COMM_WORLD, ier)
	number_sent=number_sent+1
	end do 
	! start receiving the wmatrix and place them properly
	number_received=0
 	do while (number_received<nttypes)
		call MPI_RECV(rorder, 1, MPI_INTEGER, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, status, ier)
		sender=status(MPI_SOURCE)	
		call MPI_RECV(solw,(Gsize+1)*(nperiods-deltamin+2),MPI_DOUBLE_PRECISION,sender,MPI_ANY_TAG,MPI_COMM_WORLD,status,ier)
		solwall(:,:,:,rorder)=solw
		print*, 'MASTER: received solw for order=',rorder,', total number received=', number_received+1
		number_received=number_received+1
			
		if (number_sent<nttypes) then
			order=number_sent+1
			call MPI_SEND(order, 1, MPI_INTEGER,sender,2, MPI_COMM_WORLD, ier)
			number_sent=number_sent+1
		else 
			call MPI_SEND(order, 1, MPI_INTEGER,sender ,0, MPI_COMM_WORLD, ier)
			print*, 'MASTER: finaling solw for worker',sender
			
		end if
	end do
print*,'--------------------------------------------------------------------------------'
print*, '----MASTER: FINISHED THE SOLUTION OF THE WCOEF, NOW SWITCHING TO THE V-----'
print*,'--------------------------------------------------------------------------------'
	! first things first: workers need to the results of the last round
	vtime=MPI_Wtime()
	call MPI_BCAST(solwall, (Gsize+1)*(nperiods-deltamin+2)*(deltamax-deltamin+2)*nttypes, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ier)
	print*, 'master broadcasted last stage results to workers successfully'
	! sending the first group
	number_sent=0
	do i=1,nproc-1
		order=i
		call MPI_SEND(order, 1, MPI_INTEGER,i ,11, MPI_COMM_WORLD, ier)
	number_sent=number_sent+1
	end do 

	! start receiving the vmatrix and place them properly
	number_received=0
 	do while (number_received<nttypesoc)
		call MPI_RECV(rorder, 1, MPI_INTEGER, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, status, ier)
		sender=status(MPI_SOURCE)	
		call MPI_RECV(solv,(Gsize+1)*nperiods,MPI_DOUBLE_PRECISION,sender,MPI_ANY_TAG,MPI_COMM_WORLD,status,ier)
		solvall(:,:,rorder)=solv
		number_received=number_received+1
			
		if (number_sent<nttypesoc) then
			order=number_sent+1
			call MPI_SEND(order, 1, MPI_INTEGER,sender,11, MPI_COMM_WORLD, ier)
			number_sent=number_sent+1
		else 
			call MPI_SEND(order, 1, MPI_INTEGER,sender ,10, MPI_COMM_WORLD, ier)
		end if
	end do
	
	! write vsolve matrix to a (Gsize+1,nperiods)x nttypesoc matrix to print it
	! out.

	endtime=MPI_Wtime()	
		print*, '----------------------- SUMMARY ---------------------'
		print*, 'w coef calc took', vtime-start
		print*, 'v coef calc took', endtime-vtime
		print*, 'whole thing took', endtime-start
		
	open(66,file='vcoef.txt')
	do l=1,nttypesoc
		write(66,*) "-------------- type order=",l,"------------" 
		write(66,60) ((solvall(j,k,l),k=1,nperiods),j=1,Gsize+1)
	end do
	60 format(22F46.9)
	close(66)
	
	open(77,file='wcoef.txt')
	do l=1,nttypes
		write(77,*) "-------------- type order=",l,"------------" 
		do m=1,deltamax-deltamin+2
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
    call MPI_BCAST(ftypemat, 2*na1type*(deltamax-deltamin+1), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ier)
	
	! --------SOLVING FOR W COEFFICIENTS-------------------
	do
		call MPI_RECV(order, 1, MPI_INTEGER, 0, MPI_ANY_TAG, MPI_COMM_WORLD, status, ier)
		tag=status(MPI_TAG)
		if (tag>0) then
			call wsolver(solw,ftypemat(2,order),ftypemat(1,order),parA,parW,parH,parU, beta,sigma1,rho)
			rorder=order
			call MPI_SEND(rorder, 1, MPI_INTEGER, 0, 1, MPI_COMM_WORLD, ier)
			call MPI_SEND(solw,(Gsize+1)*(nperiods-deltamin+2),MPI_DOUBLE_PRECISION,0,1,MPI_COMM_WORLD,ier)
			print*, 'worker',rank,'sent solw with order',order
		else
			EXIT
		end if
	end do
	! -------SOLVING FOR V COEFFICIENTS-------------
	! receive all the coefficients from the master-from the w stage
	call MPI_BCAST(solwall, (Gsize+1)*(nperiods-deltamin+2)*(deltamax-deltamin+2)*nttypes, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ier)
	do
		call MPI_RECV(order, 1, MPI_INTEGER, 0, MPI_ANY_TAG, MPI_COMM_WORLD, status, ier)
		tag=status(MPI_TAG)
		if (tag>10) then
			!for the order they receive, there are two corresping 3-d matrices(for each child 2 type) and these need to be put back to back.
			! I created the ftypemat and so solwall matrix so that the orders corresponds nicely. first 5*nttypesoc are for the
			! child2=type 1, the second 5*nttypesoc is for 2, etc.
			do i=1,nctype
				wcoef(:,:,:,i)=solwall(:,:,:,order+(i-1)*nttypesoc)
			end do
			call ocunpack(ocp,nocpack,ctype,parB,condprob)
			! from the order of I received, for which kinda mom I am doing this. Because the mom type has the slowest loop speed
			! while forming ftypematoc, if order<=nctype*na1type*na1type (=8), then we have type 1 mom, if it is between 9,16 a type
			! 2 mom, etc. So if I subtract one from order and then divide nctype*natype*na1type and add one, I get the type of mom. 
			! example order=2, (2-1)/8=0 (integer!!) 0+1=1. Order=16, (16-1)/8+1=2
			momorder=((order-1)/(nctype*na1type*natype))+1
			call vsolver(solv,ftypematoc(:,order), parA,parW,parH,parU,parB, beta,sigma2,wcoef,ctype,condprob(:,momorder))
			!if (order==1) then 
				!print*,' OUTSIDE type',ftypematoc(:,order)
				!print*,'OUTSIDE:', solv
			!end if
			! WORKER DETAILED REPORT ON SOLV
				!print *, 'Worker=', rank, 'at order=',order, 'Types=', ftypematoc(:,order), 'mom type', momorder,'using condprob=', condprob(:,momorder)
			rorder=order   ! returning the order
			call MPI_SEND(rorder, 1, MPI_INTEGER, 0, 1, MPI_COMM_WORLD, ier)
			call MPI_SEND(solv,(Gsize+1)*nperiods,MPI_DOUBLE_PRECISION,0,1,MPI_COMM_WORLD,ier)
		else
			EXIT
		end if
	end do

end if	
	call MPI_FINALIZE(ier)


end program act
