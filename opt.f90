!> holds optimization related routines
module opt
	use global 
	use randomgen
	use emax
	USE IFPORT
	implicit none

contains 
	! ---------------------------------------- matrix utilities-------------------------------------------------
	 function get_diagonal(mat) 	! of a square matrix
 		implicit none
		real(dble), dimension(:,:), intent(in) :: mat
		real(dble), dimension(size(mat,1)) :: get_diagonal
		integer :: j
		do j=1,size(mat,1)
 			get_diagonal(j)=mat(j,j)
		end do
	end function get_diagonal

	subroutine put_diagonal(diagv,mat) ! to a SQUARE MATRIX
		real(dble), dimension(:), intent(in) :: diagv
		real(dble), dimension(:,:), intent(inout) :: mat
		integer :: j
		do j=1,size(diagv)
			mat(j,j)=diagv(j)
		end do
	end subroutine put_diagonal

	!> extract the lower off diagonal elements of a symmetric matrix in a particular order. lowe triangle to vector
	!> vec=(a21,a31,a41,a32,a42,a43) for a 4x4 matrix.
	subroutine lt2vec(mat,vec,ndim)
		implicit none
		real(dble), intent(in) ::  mat(:,:) 		!< input
		integer, intent(in),optional :: ndim     !< dim of the matrix, optional
		real(dble), intent(out) :: vec(size(mat,1)*(size(mat,1)-1)/2) !<output
		!locals
		integer i,alt, ust,n,nvec
		n=size(mat,1)
		nvec=size(mat,1)*(size(mat,1)-1)/2
	 	
		if (ndim.ne.n) print*, 'possible error, check your dimensions- this is lttovec'
		! last element of the vector is the mat(n,n-1)
		vec(nvec)=mat(n,n-1)
		alt=1
		ust=n-1
		! then recursively fill in the vec	
		do i=1,n-2
			vec(alt:ust)=mat(i+1:n,i)
			alt=alt+(n-i)		
			ust=ust+(n-i-1)
		end do
		end subroutine lt2vec	
	
	!> place a vector to the lower diagonal elements of a (soon to be symmetric) matrix. order of the elements as lttovec subroutine.
	!>this is the inverse of that routine. These two will be used to go back and forth between vectors variance covariance matrices.
	subroutine vec2lt(vec,mat,ndim)
		implicit none
		real(dble), intent(in) :: vec(:)
		integer, intent(in)::ndim	
		real(dble), intent(out) ::  mat(ndim,ndim)
		integer i, alt, ust
		! check if dimensions match
		if ((ndim*(ndim-1)/2).NE.size(vec)) print*, 'possible error, check your dimensions. This is vectolt'

		mat=0.0d0
		alt=1
		ust=ndim-1
		do i=1,ndim-2
			mat(i+1:ndim,i)=vec(alt:ust)
			alt=alt+(ndim-i)		
			ust=ust+(ndim-i-1)
		end do
		mat(ndim,ndim-1)=vec(size(vec))
	end subroutine vec2lt

	!> take sigma matrix of dim=ndim and convert it to a vector to be fed into the optimizer
	!> this vector would have the size ndim+ndim*(ndim-1)/2 for the diagonal and off diagonal elements.
	!> the ordering is so that the first ndim elements are the variance terms.
	subroutine smat2vec(smat,vec,ndim)
		implicit none
		real(dble),intent(in)::smat(:,:)
		integer, intent(in)::ndim
		real(dble), intent(out) :: vec(ndim+ndim*(ndim-1)/2)
		
		real(dble) offdiag(ndim+ndim*(ndim-1)/2)
		integer i
		! get the diagonal elements
		do i=1,ndim
			vec(i)=smat(i,i)
		end do
		! get the off diag elements
		call lt2vec(smat,offdiag,ndim)

		! and add it to the tail of vec
		vec(ndim+1:size(vec))=offdiag

	end subroutine smat2vec

	!> take a vector, which has the same structure explain in smat2vec, and create a symmetric variance matrix.
	subroutine vec2smat(vec,smat,ndim)
		implicit none
		real(dble), intent(in) :: vec(:)
		integer, intent(in)::ndim
		real(dble), intent(out) :: smat(ndim,ndim)

		real(dble) lt(ndim,ndim),ut(ndim,ndim)
		integer i
		call vec2lt(vec(ndim+1:size(vec)),lt,ndim)
		ut=transpose(lt)
		smat=lt+ut

		do i=1,ndim
			smat(i,i)=vec(i)
		end do

	end subroutine vec2smat
	!> this subroutine packs up the parameters for the mpi send. it could also
	!>be used for the optimizer routine. The order of the parameters in the
	!>array 'packed' is as the input list of the function.
	subroutine utku_pack(packed,npack,parA,parU,parW, parH, beta, sigma1,parB,ctype,mtype, atype, condprob)
		implicit none
		real(dble), intent(in) ::  parA(:),parU(:),parW(:),parH(:),beta,sigma1(shocksize1,shocksize2)
		real(dble), intent(in):: parB(Bsizeexo+1),ctype, mtype, atype,condprob
		integer,intent(in)::npack 						!< size of the pack
		real(dble), intent(out) :: packed(npack) 
		! locals
		real(dble) vecs1(shocksize1+shocksize1*(shocksize1-1)/2)
		! first undo the sigma1 and sigma2
		call smat2vec(sigma1,vecs1,shocksize1)
		packed=(/parA,parU,parW,parH,beta,vecs1,Bsizeexo,ctype,mtype, atype, condprob/)

	end subroutine utku_pack


	!> This undoes what utku_pack does. Takes an array of size npack and creates
	!> parameter values, ready to be evaluated.
	subroutine utku_unpack(packed,npack,ppart,parA,parU,parW, parH, beta, sigma1,parB,ctype,mtype, atype, condprob)
		implicit none
		integer,intent(in)::npack,ppart(:) 			!<npack=size of packed,
		!<ppart=array partiotining packed array containing the sizes of the
		!<subcomponents. for beta, it should have a "1" in the correct spot.
		
		real(dble), intent(in) :: packed(:)
		real(dble), intent(out) :: parA(ppart(1)),parU(ppart(2)),parW(ppart(3)), parH(ppart(4)), beta, sigma1(shocksize1,shocksize1)  
		real(dble), intent(out):: parB(Bsizeexo+1),ctype, mtype, atype,condprob
		
		real(dble) svec1(shocksize1+shocksize1*(shocksize1-1)/2))
		! first basic stuff, don't touch the order of these. notice this would
		! still work if you added or subtracted parameters from this.
		parA=packed(1:ppart(1))
		parU=packed(ppart(1):sum(ppart(1:2)))
		parW=packed(sum(ppart(1:2))+1:sum(ppart(1:3)))
		parH=packed(sum(ppart(1:3))+1:sum(ppart(1:4)))
		beta=packed(sum(ppart(1:4))+1)
			
		! now extract svec1 and svec2 so that you can build sigma1 and sigma2,
		! via vec2smat
		svec1=packed(sum(ppart(1:5))+1:sum(ppart(1:6)))
		call vec2smat(svec1,sigma1, shocksize1)
		parB=packed(sum(ppart(1:6))+1:sum(ppart(1:7)))
		ctype=packed(sum(ppart(1:7)+1))
		mtype=packed(sum(ppart(1:7)+2))
		atype=packed(sum(ppart(1:7)+3))
	    condprob=packed(sum(ppart(1:7)+4))	
	end subroutine utku_unpack
		

	!> creates a family type matrix with possible family types are on the columns of the matrix.
	!> for the case where each type can take 2 values, we have 5x32 matrix
	subroutine type_pack(ftypemat,ctype,mtype,atype,a1type)
		implicit none
		real(dble), intent(in) :: ctype(:),mtype(:),a1type(:),atype(:)
		real(dble), intent(out) :: ftypemat(5,nctype*nctype*nmtype*na1type*natype)
		integer i,j,k,l,m,counter	
		counter=1
		! first counter (and the slowest loop) helps to separate coefficients for second child type.
		!see how wcoef matrix is formed for the vsolver in act.f90                                     						
		! looping thru mom later help for the vsolving, it is easier to tell what kinda mum from order 										
	do j=1,nctype 
		do k=1,nmtype  				
			do i=1,nctype
				do l=1,natype
					do m=1,na1type
						ftypemat(:,counter)=(/ctype(i),ctype(j),mtype(k),atype(l),a1type(m)/)
						counter=counter+1
					end do
				end do
			end do
		end do
	end do
	end subroutine type_pack
	
	!> creates typematrix of size 2 x na1type*(deltamax-deltamin+1) matrix that contains the family types for which coefficients
	!> needs to be calculated. This is for the simple model where the only unobserved heterogeneity is in a1.  
	subroutine type_packsimple(typemat, a1type)
	implicit none
	real(dble), intent(in) :: a1type
	real(dble), intent(out) :: typemat(2,na1type*(deltamax-deltamin+1)) 
	integer i, j, counter
	counter=1
	do j=1,2
		do i=deltamax,2,-1
			typemat(:,counter)=(/a1type(j),i/)  ! first row type values, second row delta's.
			counter=counter+1
		end do
	end do
end subroutine  type_packsimple

	!> creates the family type combinations matrix for one child families. This is going to be 5x16 matrix for the case with
	!unobserved type can take two values. this is the case even though ftype has 4 elements. We are gonna filling the second child
	!type with 999.0d0. this silliness to keep the mutype(=omega4) used in the emaxoc routines working.
	subroutine type_pack_oc(ftypematoc,ctype,mtype,atype,a1type)
		implicit none
		real(dble), intent(in) :: ctype(:),mtype(:),a1type(:),atype(:)
		real(dble), intent(out) :: ftypematoc(5,nctype*nmtype*na1type*natype)
		integer i,j,k,l,m,counter	
		counter=1

		! loop thru types, first mom (helps for vsolving), then first child
		do k=1,nmtype 		! looping thru mom first help for the vsolving, it is easier to tell what kinda mum from order
			do i=1,nctype
				do l=1,natype
					do m=1,na1type
						ftypematoc(:,counter)=(/ctype(i),999.0d0,mtype(k),atype(l),a1type(m)/)
						counter=counter+1
					end do
				end do
			end do
		end do
	end subroutine type_pack_oc

	!> ocpack packs the additional bit of info required for the one child families in vsolver calculation after for each parameter trial (or initialization).
	subroutine ocpack(ocp, nocpack,ctype,parB,condprob)
		implicit none
		real(dble), intent(in)::ctype(nctype),parB(Bsize+2),condprob(nctype,nmtype)
		integer, intent(in)::nocpack
		real(dble), intent(out) :: ocp(nocpack)
	
		! local
		integer i
		real(dble) condprobvec(nctype*nmtype)

		do i=1,nmtype
			condprobvec((i-1)*nctype+1:i*nctype)=condprob(:,i)
		end do
		ocp=(/ctype,parB,condprobvec/)
	end subroutine ocpack

	!> ocunpack unpacks the additional bit of info required for the one child families in vsolver calculationafter for each parameter trial (or initialization).
	subroutine ocunpack(ocp, nocpack,ctype,parB,condprob)
		implicit none
		real(dble), intent(out)::ctype(nctype),parB(Bsize+2),condprob(nctype,nmtype)
		integer, intent(in)::nocpack
		real(dble), intent(in) :: ocp(nocpack)
	
		! local
		integer i
		real(dble) condprobvec(nctype*nmtype)
		condprobvec=ocp(nctype+Bsize+2+1:size(ocp))

		do i=1,nmtype
			condprob(:,i)=condprobvec((i-1)*nctype+1:i*nctype)
		end do
		ctype=ocp(1:nctype)
		parB=ocp(nctype+1:nctype+Bsize+2)
	
	end subroutine ocunpack

	!----------------------OPTIMIZATION RELATED STUFF--------------------------

	!> initialize parameters
	subroutine initparam(parA,parU,parW, parH,parB, beta,sigma1, ctype, mtype, atype, a1type, condprob) 
	implicit none
		real(dble), intent(out) ::  parA(12),parU(7),parW(7),parH(6),beta,sigma1(shocksize1,shocksize2),parB(Bsizeexo+1)
		real(dble), intent(out):: ctype(nctype), mtype(nmtype), atype(natype), a1type(na1type)
		real(dble), intent(out):: condprob(nctype,nmtype) 	! prob of second child of being a certain type, conditional on mother type 
		! locals
		real(dble) var0,varp,varf,varw, varwh, var00, var01, varp0,varp1,varf0,varf1
		real(dble) cov0p, cov0f, covpf, covwwh, cov0w,cov0wh,covpw, covpwh,covfw,covfwh
		real(dble) cov0001,cov00p0,cov00p1,cov00f0,cov00f1,cov01p0,cov01p1,cov01f0,cov01f1,covp0p1,covp0f0,covp0f1,covp1f0,covp1f1,covf0f1			
		! -----------FUNDAMENTAL PARAMETERS---------------
		! gamma0's and gamma's: 
		parA=(/1.0d0,0.1d0,0.1d0,0.1d0,0.1d0, 0.1d0,0.1d0,0.1d0,0.1d0,0.1d0,0.1d0,0.1d0/)
		! gamma, alpha, alpha1,...alpha5	
		parU=(/0.1d0,0.1d0,0.1d0,0.1d0,0.1d0, 0.1d0,0.1d0/)
		! beta0 (4) and beta (3) in wage equation
		parW=(/0.1d0,0.1d0,0.1d0,0.1d0,0.1d0, 0.1d0,0.1d0/)*0.00005
		! beta_ws (4)
		parH=(/0.1d0,0.1d0,0.1d0,0.1d0,0.1d0,0.1d0/)*0.04

		beta=0.95
		condprob=1
		
		! ------------SHOCK STUFF-----------------------
		! ----   VARIANCES
		
		! hc periods
		var0=0.1d0 		! variance h=0
		varp=0.1d0 		! variance h=0.5
		varf=0.1d0 		! variance h=1
		varw=0.1d0 		! variance of wage shock
		varwh=0.1d0 	! variance of father's wage shock
		
		! hcb periods
		!var00=0.1d0 	! variance of h=0, b=0
		!var01=0.1d0 	! variance of h=0, b=1
		!varp0=0.1d0 	! variance of h=0.5, b=0
		!varp1=0.1d0 	! variance of h=0.5, b=1
		!varf0=0.1d0 	! variance of h=1, b=0
		!varf1=0.1d0 	! variance of h=1, b=1
		
		! ---  COVARIANCES
		! hc periods
		cov0p=0.1d0 	! cov(h=0,h=0.5)
		cov0f=0.1d0 	! cov(h=0,h=1)
		covpf=0.1d0 	! cov(h=0.5,h=1)
		covwwh=0.1d0 	! cov(wages of mom and dad)
		! cross w and preference shocks- JUST ZERO FOR NOW
		cov0w=0.0d0 	! cov(wage, h=0)
		cov0wh=0.0d0 	! cov(wageh, h=0)
		covpw=0.0d0 	! cov(wage, h=0.5)
		covpwh=0.0d0 	! cov(wageh, h=0.5)
		covfw=0.0d0 	! cov(wage, h=1)
		covfwh=0.0d0 	! cov(wageh, h=1)
		
		! hcb periods
	
		!cov0001=0.1d0 	! cov(h=0,b=0 and h=0, b=1)
		!cov00p0=0.1d0 	! cov(h=0,b=0 and h=0.5, b=0)	
		!cov00p1=0.1d0 	! cov(h=0,b=0 and h=0.5, b=1)	
		!cov00f0=0.1d0 	! cov(h=0,b=0 and h=1, b=0)	
		!cov00f1=0.1d0 	! cov(h=0,b=0 and h=1, b=1)	
		
		!cov01p0=0.1d0 	! cov(h=0,b=1 and h=0.5, b=0)	
		!cov01p1=0.1d0 	! cov(h=0,b=1 and h=0.5, b=1)	
		!cov01f0=0.1d0 	! cov(h=0,b=1 and h=1, b=0)	
		!cov01f1=0.1d0 	! cov(h=0,b=1 and h=1, b=1)	
	
		!covp0p1=0.1d0 	! cov(h=0.5,b=0 and h=0.5,b=1)
		!covp0f0=0.1d0
		!covp0f1=0.1d0
		
		!covp1f0=0.1d0
		!covp1f1=0.1d0

		!covf0f1=0.1d0

		
		! cross w and preference shocks
		! fuck, let's just say they are zero.	
		! TODO: NEED TO PROPERLY HANDLE PARAMETERS. YOU WILL NEED TO UPDATE THESE AT SOME POINT
		! fill in sigma => might fill this in properly later
		sigma1=0.1d0
		
		! get parB
		parB=0.3d0
		
		ctype=1.0d0
		mtype=1.0d0
		atype=1.0d0
		a1type=(/0.0d0,1.0d0/)
		!ctype=(/0.0d0,1.0d0/)
		!mtype=(/0.d0,1.0d0/)
		!atype=(/1.0d0,1.5d0/)
		!a1type=(/1.0d0,2.0d0/)
	
		! conditional probabilities of child 2 type should be pr(x|y)=pr(ctype=x,mtype=y)/pr(mtype=y)
		! so I need to specify the joint distribution of ctype and mtype. if this does not go anywhere else I can just use the
		! conditional probs as a parameter 
			
	end subroutine initparam	

end module opt

