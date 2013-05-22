!> holds optimization related routines
module optu
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
	packed=(/parA,parU,parW,parH,beta,vecs1,parB,ctype,mtype, atype, condprob/)

end subroutine utku_pack


!> This undoes what utku_pack does. Takes an array of size npack and creates
!> parameter values, ready to be evaluated.
subroutine utku_unpack(packed,npack,ppart,parA,parU,parW, parH, beta, sigma1,parB,ctype,mtype, atype, condprob)
	implicit none
	integer,intent(in)::npack,ppart(:) 			!<npack=size of packed,
	!<ppart=array partiotining packed array containing the sizes of the
	!<subcomponents. for beta, it should have a "1" in the correct spot.
	
	real(dble), intent(in) :: packed(npack)
	real(dble), intent(out) :: parA(ppart(1)),parU(ppart(2)),parW(ppart(3)), parH(ppart(4)), beta, sigma1(shocksize1,shocksize1)  
	real(dble), intent(out):: parB(Bsizeexo+1),ctype, mtype, atype,condprob
	
	real(dble) svec1(shocksize1+shocksize1*(shocksize1-1)/2)
	! first basic stuff, don't touch the order of these. notice this would
	! still work if you added or subtracted parameters from this.
	parA=packed(1:ppart(1))
	parU=packed(ppart(1)+1:sum(ppart(1:2)))
	parW=packed(sum(ppart(1:2))+1:sum(ppart(1:3)))
	parH=packed(sum(ppart(1:3))+1:sum(ppart(1:4)))
	beta=packed(sum(ppart(1:4))+1)
		
	! now extract svec1 and svec2 so that you can build sigma1 and sigma2,
	! via vec2smat
	svec1=packed(sum(ppart(1:5))+1:sum(ppart(1:6)))
	call vec2smat(svec1,sigma1, shocksize1)
	parB=packed(sum(ppart(1:6))+1:sum(ppart(1:7)))
	ctype=packed(sum(ppart(1:7))+1)
	mtype=packed(sum(ppart(1:7))+2)
	atype=packed(sum(ppart(1:7))+3)
	condprob=packed(sum(ppart(1:7))+4)	
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
	real(dble), intent(in) :: a1type(:)
	real(dble), intent(out) :: typemat(2,na1type*(deltamax-deltamin+1)) 
	integer i, j, counter
	counter=1
	do j=1,2
		do i=2, deltamax
			typemat(:,counter)=(/a1type(j),i*1.0d0/)  ! first row type values, second row delta's.
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


! -----------------------------SIMULATION ROUTINES-----------------------------------
! --------------------------AUX FUNCTIIONS FOR SIM----------------------
!> same as wagef but accepts E as a vector, so calculated the wages for a bunch of Es and eps. wagef in emax.f90 is designed to
!> give you wages for one E with multiple shocks so you can calculate the mc integral.
function wagefsim(E,time,llms,logw0,par,eps)
	implicit none
	real(dble) E(:),time,llms,logw0 				 	!< experience,time, llms,logw0, as they appear in extended
	real(dble) eps(:) 							!< emax mc integration random draw, size=Nmc
	real(dble) par(:) 		 					!< parameter vectors, for wage 
	real(dble) wagefsim(size(eps)) 				!< output, wage income for the mother, size=Nmc
	real(dble) logw(size(eps))
	logw=logw0+par(1)*E+par(2)*time+par(3)*llms
	wagefsim=exp(logw+eps)
end function wagefsim

! pick a random element from the vector x, with probabilities given in w (sumw=1), with the uniform random draw u.
function randtype(x,w,u)
implicit none
	integer randtype		 	!< output
	integer x(:)				!< input: array to drawn from. (/1,...,nttypes/)
	real(dble) w(size(x))		!< distribution function ( P(x) function)
	real(dble) u				!< uniform (0,1) draw

	real(dble) cw(size(w)+1)
	integer i, pick

	cw=(/0.0d0,cumx(w)/)
	do i=1,size(x)
		if ((cw(i)<=u).AND.(cw(i+1)>u))  pick=i
	end do
	randtype=x(pick)
end function randtype


!> aux function to get cumulative distributions, helps: randtype
function cumx(x)
implicit none
	real(dble) x(:)				!< input vector x
	real(dble) cumx(size(x))	!<cumulative sum of the vector x

	integer i
	real(dble) sumto
	sumto=0
	do i=1,size(x)
		sumto=sumto+x(i)
		cumx(i)=sumto
	end do
end function cumx

!> simulate histories of choice, experience, test scores, birth from a given
!> parameter trial, approximate solution parameters. Also produces smoothed
!> histories a la Keane and Smith. This is for one person. It needs to be called
!>for each person in the sample with the correct unique sample ID.
subroutine simhist(SS,outcomes,testoutcomes, choices, xchoices,birthhist,smchoices, smexperience, smAs, smtestoutcomes,omega3,intercepts,parA,parU,parW,parH,beta,sigma,a1type,pa1type,parBmat,vcoef,wcoef,llmsvec,id,rho,lambdas,sigmaetas,smpar)
	implicit none
	real(dble), intent(in) :: omega3(:) 		!< observed family type
	real(dble), intent(in) :: intercepts(:) 	!< normaly hetero, but for now constant: ctype,mtype,atype
	real(dble), intent(in) :: parA(:),parU(:),parW(:),parH(:) 	!<model parameters
	real(dble), intent(in) :: sigma(:,:) 						!< variance covariance matrix of shocks
	real(dble), intent(in) :: beta
	real(dble), intent(in) :: a1type(:),pa1type(:) 				!<a1 distribution, vectors of size 2 each.
	real(dble), intent(in) :: parBmat(Bsizeexo+1,nfert) 		!<(4+1,6) birth probability parameters.
	integer, intent(in):: id 							!<unique sample id of the family
	real(dble), intent(in) :: llmsvec(:) 				!< actual local labor market shocks for the sample object.
	real(dble), intent(in) :: vcoef(Gsizeoc+1,nperiods,nttypes),wcoef(Gsize+1,nperiods-deltamin+2,deltamax-deltamin+1,nttypes)
	real(dble), intent(in) :: rho 						!<skill depreciation para, added later.
	real(dble),intent(out):: SS(6,nperiods,Npaths) 		!< Simulated State space (A1,A2,E,age1,age2,agem)xnperiods,Npaths	
	real(dble),intent(out):: outcomes(2,nperiods,Npaths) !< outcomes: wages of the father and mother.
	integer,intent(out):: choices(1,nperiods,Npaths) 	!< choices : the choice history of h.
	integer,intent(out):: xchoices(1,nperiods,Npaths) 	!< not in the data, but I will keep an history of the x choices as well.
	integer,intent(out):: birthhist(Npaths) 			!< the birth timing vec, 0 if one child throughout.
	! test score stuff
	real(dble), intent(in) ::lambdas(Ntestage) 			!< old factor loadings in a two test world 
	real(dble), intent(in) :: sigmaetas(2,Ntestage) 			!<measurement error variance for each age.
	real(dble), intent(out) :: testoutcomes(4,Ntestage, Npaths)

	! smoothed stuff
	real(dble), intent(in) :: smpar
	real(dble), intent(out) :: smchoices(3,nperiods, Npaths)
	real(dble), intent(out) :: smexperience(nperiods, Npaths)
	real(dble), intent(out) :: smAs(2,nperiods, Npaths)
	real(dble) , intent(out):: smtestoutcomes(4,Ntestage, Npaths)
	
		
	! LOCALS
	integer period, mainseed(2)
	integer a1holder(Npaths)
	real(dble) pa1(Npaths)
	real(dble) Ainit,wage(Npaths),wageh(Npaths)
	real(dble) mu(shocksize1)
	real(dble) eps(Npaths,shocksize1)
	
	real(dble) logw0
	real(dble) power, mult

	
	! choice specific utilities WITHIN the period
	real(dble) umat(3,Npaths), umatbig(3,11,Npaths),uc
	! interpolating vectors within a period and for an npath alone, constant parts predicted beforehand.
	real(dble) intcons(4), fvconsw(nttypes), fvconsv(nttypes)
	real(dble) intw(9), intv(5) ,fvw,fvv
	real(dble) A1next, A2next, Enext, Anext(2),inputs(5), income
	real(dble) nextcollect(3,3) ! to collect collected future state space points and carry them to next period.
	real(dble) nextcollectbig(3,3,xgridsize) ! 
	real(dble) next(3,Npaths)  ! to collect collected future state space points and carry them to next period.
	real(dble) pbirth, omegaB(4), birthdraw(Npaths) 
	integer birth, baby
	integer i,j,k,l
	integer coordbig(2)
	real(dble) incrat ! ratio of the income that goes to consumption. I use this in the periods 15,21 to change income for W.

	integer a1order(nttypes)
	integer picker(1)

	! adding Test Scores
	real(dble) sigmafortests(2,2)
	real(dble) etashock(Npaths,4)

	! smoothing shit
	real(dble) smdenom
	real(dble) smxchoices(xgridsize, nperiods, Npaths)
	real(dble) smh, smx
	real(dble) smnext(3,Npaths)
	real(dble) smbig(3,xgridsize)
		
	integer counter
	
	! Iniatilize some stuff
	mu=0.0d0
	choices=-9999
	xchoices=-9999 !if not chosen, be -9999.
	smchoices=-9999
	smxchoices=-9999 !if not chosen, be -9999.
	period=1
	birthhist=0 	! this records the timing of the second birth. if none, stays at zero.	
	
	counter=0
	! first draw Npaths alternative folks to assign a1 values.	
	mainseed=id*(/3,2/)
	call random_seed(put=mainseed)
	call random_number(pa1)
	a1order=(/(i,i=1, nttypes)/)
	
	! a1 holds 
	do i=1, Npaths
		a1holder(i)=randtype(a1order,pa1type,pa1(i))
	end do

	! initialize A
	Ainit=intercepts(1)+parA(1)*omega3(1)+parA(2)*omega3(2)+parA(3)*omega3(3) 
	! initialize wage	
	logw0=intercepts(2)+parW(1)*omega3(1)+parW(2)*omega3(2)+parW(3)*omega3(3)+parW(4)*omega3(3)*omega3(3)
	! birth prob determiner
	omegaB=(/omega3(3),omega3(3)**2,omega3(1),omega3(2)/)

	! these are for the temp term values
	power=(finalage/timeperiod)-omega3(3)
	mult=beta*(1.0d0-beta**power)/(1.0d0-beta)
	
	eps= randmnv(Npaths,5,mu,sigma,3,1,(/1,id*1/))
	! initialize omegas
	SS(1,1,:)=Ainit 	!A1
	SS(2,1,:)=0.0d0 	!A2
	SS(3,1,:)=0.0d0 	!E
	SS(4,1,:)=1.0d0 	!age1
	SS(5,1,:)=0.0d0 	!age2
	SS(6,1,:)=omega3(3) !agem
	smAs=0.0d0
	smexperience=0.0d0
	! initialize test scores
	testoutcomes=-9999.0d0
	smtestoutcomes=-9999.0d0
	wage=wagef(0.d0,1.0d0,llmsvec(1),omega3(1),omega3(2),omega3(3),eps(:,4),intercepts(2),parW(5:7),parW(1:4))
	wageh=wagehfquick(omega3(4),1.0d0,eps(:,5),parH(4:5))	
	outcomes(1,1,:)=wage
	outcomes(2,1,:)=wageh
	! current utilities
	umat(1,:)=parU(2)*SS(1,1,:) 	+     parU(4)*(wageh)**parU(5)+parU(6)
	umat(2,:)=parU(2)*SS(1,1,:)+a1type(a1holder)*0.5d0+ parU(4)*(wageh+0.5d0*wage)**parU(5)+parU(6)+parU(7)*0.5d0
	umat(3,:)=parU(2)*SS(1,1,:)+a1type(a1holder)*1.0d0+ parU(4)*(wageh+wage)**parU(5)+parU(6)+parU(7)
	! --------------------future value calculations----------------------
	! first, birth probabilities next period
	pbirth=bprobexo(omegaB,parBmat(:,1))
	! what I do this period will affect future A and E, the rest of the state space (which are used in the interpolations) are
	! deterministic.
	! first common deterministic parts
	intcons=(/omega3(3)+1.0d0,llmsvec(period),omega3(1),omega3(2)/)
	do k=1,nttypes
		fvconsv(k)=sum(vcoef(4:7,period+1,k)*intcons)+vcoef(Gsizeoc+1,period+1,k)
		fvconsw(k)=sum(wcoef(4:7,period+1,period,k)*intcons)+wcoef(Gsize+1,period+1,period,k)
	end do
	do i=1, Npaths
		do j=1,3
			A1next=pfone((/SS(1,period,i),(1-0.25d0*(j-1)),0.5d0*(wageh(i)+wage(i)*0.5d0*(j-1))/), omega3(1:3), period*1.0d0,parA(4:12),rho) 
			A2next=Ainit
			Enext=SS(3,period,i)+(j-1)*0.5d0
			intw=(/A1next,A2next,Enext,A1next**2,A2next**2,Enext**2, A1next*A2next,A1next*Enext,A2next*Enext/)
			intv=(/A1next,Enext,A1next**2,Enext**2, A1next*Enext/)
			fvw=sum(intw*(/wcoef(1:3,period+1,period,a1holder(i)),wcoef(9:Gsize,period+1,period,a1holder(i))/))
			fvv=sum(intv*(/vcoef(1:2,period+1,a1holder(i)),vcoef(8:Gsizeoc,period+1,a1holder(i))/))
			! now add the expected future value to corresponding umat element
			umat(j,i)=umat(j,i)+beta*(pbirth*(fvw+fvconsw(a1holder(i)))+(1-pbirth)*(fvv+fvconsv(a1holder(i))))
			nextcollect(:,j)=(/A1next,A2next,Enext/)
		end do
		picker=maxloc(umat(:,i))
		choices(1,period,i)=picker(1)
		next(:,i)=nextcollect(:,picker(1))

		! NEW: Smoothed As and Es and choices and later test scores (not this period, though)
		smdenom=sum(exp((umat(:,i)-maxval(umat(:,i)))/smpar))
		smchoices(:,period,i)=exp((umat(:,i)-maxval(umat(:,i)))/smpar)/smdenom
		smh=sum(smchoices(:,period,i)*(/0.d0,0.5d0,1.0d0/))
		smnext(1,i)=pfone((/SS(1,period,i),smh,0.5d0*(wageh(i)+wage(i)*smh)/), omega3(1:3), period*1.0d0,parA(4:12),rho) 
		smnext(2,i)=SS(2,period,i)
		smnext(3,i)=SS(3,period,i)+smh
	end do
	

	!                            -----------------LOOOOOP TO THE FUTURE-----------------
	
	do period=2,nperiods
		! draw the period shocks
		eps= randmnv(Npaths,5,mu,sigma,3,1,(/period,id*period/))
		! update state space vectors as much as you can outside the loop	
		! USE THE NEXT MATRIX TO UPDATE THE STATE SPACE FOR As and E, so that I don't need to calculate them twice.
		! I am updating A2 as well, but if it turns out this is a period where there is no second or birth, it needs to be
		! revert back to zero. Same as for age2
		SS(1,period,:)=next(1,:) 
		SS(2,period,:)=next(2,:)
		SS(3,period,:)=next(3,:)
		SS(4,period,:)=SS(4,period-1,:)+1.0d0 	! update age1
		SS(5,period,:)=SS(5,period-1,:) 		!  DO NOT update age2
		SS(6,period,:)=SS(6,period-1,:)+1.0d0 	! udpate agem
		!SS(3,period,:)=SS(3,period-1,:)+(choices(1,period-1,:)*1.0d0-1.0d0)/2.0d0
		! these are all at vector level, i.e., for all Npaths 
		wage=wagefsim(SS(3,period,:),period*1.0d0,llmsvec(period),logw0,parW(5:7),eps(:,4))	
		wageh=wagehfquick(omega3(4),period*1.0d0,eps(:,5),parH(4:5))	
		outcomes(1,period,:)=wage
		outcomes(2,period,:)=wageh
		smexperience(period,:)=smnext(3,:)
		smAs(:,period,:)=smnext(1:2,:)
		
		if (period>=testminage .AND. period<=testmaxage) then
			sigmafortests=0.0d0
			sigmafortests(1,1)=sigmaetas(1,period-testminage+1)
			sigmafortests(2,2)=sigmaetas(2,period-testminage+1)
			etashock=randmnv(Npaths,4,mu,sigmafortests,3,1,(/100*period+5*period+period,300*id*period+15*period+period+id/))	
		end if

		! get some period specific but not path specific stuff to calculate interpolated fv
		intcons=(/omega3(3)+period*1.0d0,llmsvec(period),omega3(1),omega3(2)/)
		if (period<nperiods) then
			do k=1,nttypes
				fvconsv(k)=sum(vcoef(4:7,period+1,k)*intcons)+vcoef(Gsizeoc+1,period+1,k)
			end do
		end if

		! -------------------------- 1a. PERIOD<7: STILL IN THE FECUND PERIOD ---------------------------------
		if (period<7)  then

			! determine the prob of a child birth. this part is constant all the Npaths
			pbirth=bprobexo(omegaB,parBmat(:,period-1))	
			call random_seed(put=mainseed*((period+1)*10)+period+2) 
			call random_number(birthdraw)
			! FIRST FIGURE OUT WHO HAD BIRTH
			
			do i=1,Npaths
				if (SS(5,period-1,i)<0.001) then   ! they don't have a second child
					if (birthdraw(i)<pbirth) then
						! record the birth time
						birthhist(i)=period 
						! initialize ss for the young one.
						!SS(2,period,i)=Ainit  !Already done this above	
						SS(5,period,i)=1.0d0
					else
						birthhist(i)=0
						SS(2,period,i)=0.0d0 	! reset A2 back to zero.
						SS(5,period,i)=0.0d0 	! reset age2 back to zero
					end if
				end if
			end do

			! now we know who has two child in this period (<8) so we can simulate their decisions. This entails the calculation
			! of choice specific utilities.

			do i=1,Npaths
				! ----------------------------1.a PERIOD<7: ONE CHILD-----------------------------------
			
				if (SS(5,period,i)<0.001) then 		! if still with one child. 
					if (SS(4,period,i)<3) then
						baby=1
					else 
						baby=0
					end if
					uc=parU(2)*SS(1,period,i) 
					umat(1,i)=uc+parU(4)*(wageh(i))**parU(5)+parU(6)*baby
					umat(2,i)=uc+a1type(a1holder(i))*0.5d0+ parU(4)*(wageh(i)+0.5d0*wage(i))**parU(5)+parU(6)*baby+parU(7)*0.5d0*baby
					umat(3,i)=uc+a1type(a1holder(i))*1.0d0+ parU(4)*(wageh(i)+wage(i))**parU(5)+parU(6)*baby+parU(7)*baby
					
					do k=1,nttypes
						fvconsw(k)=sum(wcoef(4:7,period+1,period,k)*intcons)+wcoef(Gsize+1,period+1,period,k)
					end do
						! now calculate future values for the possible 3 decisions, so we can calculate the future value.
					
					do j=1,3
						A1next=pfone((/SS(1,period,i),(1-0.25d0*(j-1)),0.5d0*(wageh(i)+wage(i)*0.5d0*(j-1))/), omega3(1:3), period*1.0d0,parA(4:12),rho) 
						A2next=Ainit
						Enext=SS(3,period,i)+(j-1)*0.5d0
						intw=(/A1next,A2next,Enext,A1next**2,A2next**2,Enext**2, A1next*A2next,A1next*Enext,A2next*Enext/)
						intv=(/A1next,Enext,A1next**2,Enext**2, A1next*Enext/)
						fvw=sum(intw*(/wcoef(1:3,period+1,period,a1holder(i)),wcoef(9:Gsize,period+1,period,a1holder(i))/))
						fvv=sum(intv*(/vcoef(1:2,period+1,a1holder(i)),vcoef(8:Gsizeoc,period+1,a1holder(i))/))
						! now add the expected future value to corresponding umat element
						umat(j,i)=umat(j,i)+beta*(pbirth*(fvw+fvconsw(a1holder(i)))+(1-pbirth)*(fvv+fvconsv(a1holder(i))))
						nextcollect(:,j)=(/A1next,A2next,Enext/)
					end do
					! finally, pick the highest expected utility one.
					picker=maxloc(umat(:,i))
					choices(1,period,i)=picker(1)
					! get the optimized future state space points to use them in future.
					next(:,i)=nextcollect(:,picker(1))
					
					! NEW: Smoothed As and Es and choices and later test scores (not this period, though)
					smdenom=sum(exp((umat(:,i)-maxval(umat(:,i)))/smpar))
					smchoices(:,period,i)=exp((umat(:,i)-maxval(umat(:,i)))/smpar)/smdenom
					smh=sum(smchoices(:,period,i)*(/0.d0,0.5d0,1.0d0/))
					smnext(1,i)=pfone((/smAs(1,period,i),smh,0.5d0*(wageh(i)+wage(i)*smh)/), omega3(1:3), period*1.0d0,parA(4:12),rho) 
					smnext(2,i)=SS(2,period,i)
					smnext(3,i)=smexperience(period,i)+smh
					
					if (id==1) print*, period, '------------'
					if (id==1) print*,smnext(:,1:50)
				!-------------------------- PERIOD<7: 1.b TWO CHILDREN ------------------------------
				else
					do k=1,nttypes
						fvconsw(k)=sum(wcoef(4:7,period+1,birthhist(i),k)*intcons)+wcoef(Gsize+1,period+1,birthhist(i),k)
					end do
					! increase the age2 by one if it is different than zero
					SS(5,period,i)=SS(5,period,i)+1.0d0
					! is there a baby at home
					if (SS(4,period,i)<3 .OR. SS(5,period,i)<3) then
						baby=1
					else 
						baby=0
					end if
					! now we x is a choice, so I will have to use the bigger umatbig
					! period utilities don't depend on x, so first get rid of those
					! because there are two children, I need to use the CES utility function
					uc=uctwo(SS(1,period,i),SS(2,period,i),parU(1)) 
					umatbig(1,:,i)=	uc + 			 parU(4)*(wageh(i))**parU(5)+parU(6)*baby
					umatbig(2,:,i)= uc +a1type(a1holder(i))*0.5d0+ parU(4)*(wageh(i)+0.5d0*wage(i))**parU(5)+parU(6)*baby+parU(7)*0.5d0*baby
					umatbig(3,:,i)= uc +a1type(a1holder(i))*1.0d0+ parU(4)*(wageh(i)+wage(i))**parU(5)+parU(6)*baby+parU(7)*baby

					! now calculate and add future values.
					do j=1,3
						! E does not depend on x, so calculate it beforehand, also income
						Enext=SS(3,period,i)+(j-1)*0.5d0
						income=wageh(i)*0.5d0+wage(i)*(j-1)*0.25 ! it is 0.25 because half the income goes to goods inputs
						do k=1,xgridsize
							inputs=(/SS(1,period,i),SS(2,period,i),(1-0.25*(j-1))*xgrid(k),(1-0.25*(j-1))*(1-xgrid(k)),income/)
							Anext=pftwo(inputs,omega3(1:3),SS(4:5,period,i),parA(4:12), rho)
							A1next=Anext(1)
							A2next=Anext(2)
							intw=(/A1next,A2next,Enext,A1next**2,A2next**2,Enext**2, A1next*A2next,A1next*Enext,A2next*Enext/)
							fvw=sum(intw*(/wcoef(1:3,period+1,birthhist(i),a1holder(i)),wcoef(9:Gsize,period+1,birthhist(i),a1holder(i))/))
							umatbig(j,k,i)=umatbig(j,k,i)+beta*(fvw+fvconsw(a1holder(i)))
							nextcollectbig(:,j,k)=(/A1next,A2next,Enext/)
						end do
					end do
					coordbig=maxloc(umatbig(:,:,i))
					choices(1,period,i)=coordbig(1)
					xchoices(1,period,i)=coordbig(2)
					next(:,i)=nextcollectbig(:,coordbig(1),coordbig(2))
					
					! NEW: Smoothed As and Es and choices and later test scores 
					smdenom=sum(exp((umatbig(:,:,i)-maxval(umatbig(:,:,i)))/smpar))
					
					smbig=exp((umatbig(:,:,i)-maxval(umatbig(:,:,i)))/smpar)/smdenom
					smchoices(:,period,i)=sum(smbig,2) 							! marginal of h
					smxchoices(:,period,i)=sum(smbig,1) 						! marginal of x
					smh=sum(smchoices(:,period,i)*(/0.d0,0.5d0,1.0d0/))
					smx=sum(smxchoices(:,period,i)*xgrid)
					
					income=(wageh(i)+wage(i)*smh)*0.5
					inputs=(/smAs(1,period,i),smAs(2,period,i),(1-0.5*smh)*smx,(1-0.5*smh)*(1-smx),income/)
					smnext(1:2,i)=pftwo(inputs,omega3(1:3),SS(4:5,period,i),parA(4:12), rho)
					smnext(3,i)=smexperience(period,i)+smh
					
					
				end if  ! end of period<8 WV if
					
				if (SS(4,period,i)>=testminage) then
					testoutcomes(1,period-testminage+1,i)=SS(1,period,i)+etashock(i,1)
					testoutcomes(2,period-testminage+1,i)=lambdas(period-testminage+1)*SS(1,period,i)+etashock(i,1)
					smtestoutcomes(1,period-testminage+1,i)=smAs(1,period,i)+etashock(i,1)
					smtestoutcomes(2,period-testminage+1,i)=lambdas(period-testminage+1)*smAs(1,period,i)+etashock(i,1)
				end if

				if (SS(5,period,i)>=testminage) then
					testoutcomes(3,period-testminage+1,i)=SS(2,period,i)+etashock(i,1)
					testoutcomes(4,period-testminage+1,i)=lambdas(period-testminage+1)*SS(2,period,i)+etashock(i,1)
					smtestoutcomes(3,period-testminage+1,i)=smAs(2,period,i)+etashock(i,1)
					smtestoutcomes(4,period-testminage+1,i)=lambdas(period-testminage+1)*smAs(2,period,i)+etashock(i,1)
				end if
			end do

			! --%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%--
			!----------------------------2. 7<period<15 NO CHANCE OF BIRTH, FIRST KID STILL IN CARE---------------------------------
		elseif ((period>=7) .AND. (period<astar)) then
			do i=1,Npaths
				! calculations are different for V and W.
				
				! ---------------------------2.a PERIOD (7,15): ONE CHILD-----------------------------------
				if (SS(5,period,i)<0.001) then 		! if still with one child. 
					if (SS(4,period,i)<3) then
						baby=1
					else 
						baby=0
					end if
					uc=parU(2)*SS(1,period,i) 
					umat(1,i)=uc+parU(4)*(wageh(i))**parU(5)+parU(6)*baby
					umat(2,i)=uc+a1type(a1holder(i))*0.5d0+ parU(4)*(wageh(i)+0.5d0*wage(i))**parU(5)+parU(6)*baby+parU(7)*0.5d0*baby
					umat(3,i)=uc+a1type(a1holder(i))*1.0d0+ parU(4)*(wageh(i)+wage(i))**parU(5)+parU(6)*baby+parU(7)*baby
					
					! now calculate future values for the possible 3 decisions, so we can calculate the future value.
					do j=1,3
						A1next=pfone((/SS(1,period,i),(1-0.25d0*(j-1)),0.5d0*(wageh(i)+wage(i)*0.5d0*(j-1))/), omega3(1:3), 1.0d0*period,parA(4:12),rho) 
						A2next=Ainit
						Enext=SS(3,period,i)+(j-1)*0.5d0
						intv=(/A1next,Enext,A1next**2,Enext**2, A1next*Enext/)
						fvv=sum(intv*(/vcoef(1:2,period+1,a1holder(i)),vcoef(8:Gsizeoc,period+1,a1holder(i))/))
						! now add the expected future value to corresponding umat element
						umat(j,i)=umat(j,i)+beta*(fvv+fvconsv(a1holder(i)))
						nextcollect(:,j)=(/A1next,A2next,Enext/)
					end do
					! finally, pick the highest expected utility one.
					picker=maxloc(umat(:,i))
					choices(1,period,i)=picker(1)
					! get the optimized future state space points to use them in future.
					next(:,i)=nextcollect(:,choices(1,period,i))
					
					! NEW: Smoothed As and Es and choices and later test scores (not this period, though)
					smdenom=sum(exp((umat(:,i)-maxval(umat(:,i)))/smpar))
					smchoices(:,period,i)=exp((umat(:,i)-maxval(umat(:,i)))/smpar)/smdenom
					smh=sum(smchoices(:,period,i)*(/0.d0,0.5d0,1.0d0/))
					smnext(1,i)=pfone((/smAs(1,period,i),smh,0.5d0*(wageh(i)+wage(i)*smh)/), omega3(1:3), period*1.0d0,parA(4:12),rho) 
					smnext(2,i)=SS(2,period,i)
					smnext(3,i)=smexperience(period,i)+smh
				
				!-------------------------- 2.b PERIOD(7,15): TWO CHILDREN ------------------------------
				else
					! increase the age2 by one if it is different than zero
					 SS(5,period,i)=SS(5,period,i)+1.0d0
					
					if (SS(4,period,i)<3 .OR. SS(5,period,i)<3) then
						baby=1
					else 
						baby=0
					end if
					! now we x is a choice, so I will have to use the bigger umatbig
					! period utilities don't depend on x, so first get rid of those
					! because there are two children, I need to use the CES utility function
					uc=uctwo(SS(1,period,i),SS(2,period,i),parU(1)) 
					
					umatbig(1,:,i)=	uc + 			 parU(4)*(wageh(i))**parU(5)+parU(6)*baby
					umatbig(2,:,i)= uc +a1type(a1holder(i))*0.5d0+ parU(4)*(wageh(i)+0.5d0*wage(i))**parU(5)+parU(6)*baby+parU(7)*0.5d0*baby
					umatbig(3,:,i)= uc +a1type(a1holder(i))*1.0d0+ parU(4)*(wageh(i)+wage(i))**parU(5)+parU(6)*baby+parU(7)*baby
					! now calculate and add future values.
				
					do k=1,nttypes
						fvconsw(k)=sum(wcoef(4:7,period+1,birthhist(i),k)*intcons)+wcoef(Gsize+1,period+1,birthhist(i),k)
					end do
					
					do j=1,3
						! E does not depend on x, so calculate it beforehand, also income
						Enext=SS(3,period,i)+(j-1)*0.5d0
						income=wageh(i)*0.5d0+wage(i)*(j-1)*0.25 ! it is 0.25 because half the income goes to goods inputs
						do k=1,xgridsize
							inputs=(/SS(1,period,i),SS(2,period,i),(1-0.25*(j-1))*xgrid(k),(1-0.25*(j-1))*(1-xgrid(k)),income/)
							Anext=pftwo(inputs,omega3(1:3),SS(4:5,period,i),parA(4:12), rho)
							A1next=Anext(1)
							A2next=Anext(2)
							intw=(/A1next,A2next,Enext,A1next**2,A2next**2,Enext**2, A1next*A2next,A1next*Enext,A2next*Enext/)
							fvw=sum(intw*(/wcoef(1:3,period+1,birthhist(i),a1holder(i)),wcoef(9:Gsize,period+1,birthhist(i),a1holder(i))/))
							umatbig(j,k,i)=umatbig(j,k,i)+beta*(fvw+fvconsw(a1holder(i)))
							nextcollectbig(:,j,k)=(/A1next,A2next,Enext/)
						end do
					end do
					coordbig=maxloc(umatbig(:,:,i))
					choices(1,period,i)=coordbig(1)
					xchoices(1,period,i)=coordbig(2)
					next(:,i)=nextcollectbig(:,coordbig(1),coordbig(2))
					
					! NEW: Smoothed As and Es and choices and later test scores 
					smdenom=sum(exp((umatbig(:,:,i)-maxval(umatbig(:,:,i)))/smpar))
					smbig=exp((umatbig(:,:,i)-maxval(umatbig(:,:,i)))/smpar)/smdenom
					smchoices(:,period,i)=sum(smbig,2) 							! marginal of h
					smxchoices(:,period,i)=sum(smbig,1) 						! marginal of x
					smh=sum(smchoices(:,period,i)*(/0.d0,0.5d0,1.0d0/))
					smx=sum(smxchoices(:,period,i)*xgrid)
					
					income=(wageh(i)+wage(i)*smh)*0.5
					inputs=(/smAs(1,period,i),smAs(2,period,i),(1-0.5*smh)*smx,(1-0.5*smh)*(1-smx),income/)
					smnext(1:2,i)=pftwo(inputs,omega3(1:3),SS(4:5,period,i),parA(4:12), rho)
					smnext(3,i)=smexperience(period,i)+smh
				end if  ! end of period<8 WV if
			
			! create test scores	

				if (SS(4,period,i)>=testminage) then
					testoutcomes(1,period-testminage+1,i)=SS(1,period,i)+etashock(i,1)
					testoutcomes(2,period-testminage+1,i)=lambdas(period-testminage+1)*SS(1,period,i)+etashock(i,1)
					smtestoutcomes(1,period-testminage+1,i)=smAs(1,period,i)+etashock(i,1)
					smtestoutcomes(2,period-testminage+1,i)=lambdas(period-testminage+1)*smAs(1,period,i)+etashock(i,1)
				end if

				if (SS(5,period,i)>=testminage) then
					testoutcomes(3,period-testminage+1,i)=SS(2,period,i)+etashock(i,1)
					testoutcomes(4,period-testminage+1,i)=lambdas(period-testminage+1)*SS(2,period,i)+etashock(i,1)
					smtestoutcomes(3,period-testminage+1,i)=smAs(2,period,i)+etashock(i,1)
					smtestoutcomes(4,period-testminage+1,i)=lambdas(period-testminage+1)*smAs(2,period,i)+etashock(i,1)
				end if
			
			end do
			
			!#%#%#%#%#%#%#%#%#%#%#%#%%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%%#%#%#%#%%#%#


			! ---------------------3. Period (15,21) No x choice, because Kid 1 has grown
		elseif ((period>=astar).AND.(period<22)) then 
			do i=1,Npaths
				! calculations are different for V and W.
			
				! ---------------------------3.a PERIOD (16,21): ONE CHILD-----------------------------------
				! 
				if (SS(5,period,i)<0.001) then 		! if still with one child. 
					! there is no chance of having a baby at this point
					baby=0
					uc=parU(2)*SS(1,period,i) 
					! NOTE: important difference: all the income goes into the consumption now, because there is no production of
					! skils
					umat(1,i)=uc+parU(4)*(wageh(i))**parU(5)+parU(6)*baby
					umat(2,i)=uc+a1type(a1holder(i))*0.5d0+ parU(4)*(wageh(i)+0.5d0*wage(i))**parU(5)+parU(6)*baby+parU(7)*0.5d0*baby
					umat(3,i)=uc+a1type(a1holder(i))*1.0d0+ parU(4)*(wageh(i)+wage(i))**parU(5)+parU(6)*baby+parU(7)*baby
					
					! now calculate future values for the possible 3 decisions, so we can calculate the future value.
					do j=1,3
						A1next=SS(1,period,i)
						A2next=SS(2,period,i)
						Enext=SS(3,period,i)+(j-1)*0.5d0
						intv=(/A1next,Enext,A1next**2,Enext**2, A1next*Enext/)
						fvv=sum(intv*(/vcoef(1:2,period+1,a1holder(i)),vcoef(8:Gsizeoc,period+1,a1holder(i))/))
						! now add the expected future value to corresponding umat element
						umat(j,i)=umat(j,i)+beta*(fvv+fvconsv(a1holder(i)))
						nextcollect(:,j)=(/A1next,A2next,Enext/)
					end do
					! finally, pick the highest expected utility one.
					picker=maxloc(umat(:,i))
					choices(1,period,i)=picker(1)
					! get the optimized future state space points to use them in future.
					next(:,i)=nextcollect(:,choices(1,period,i))
					
					! NEW: Smoothed As and Es and choices and later test scores (not this period, though)
					smdenom=sum(exp((umat(:,i)-maxval(umat(:,i)))/smpar))
					smchoices(:,period,i)=exp((umat(:,i)-maxval(umat(:,i)))/smpar)/smdenom
					smh=sum(smchoices(:,period,i)*(/0.d0,0.5d0,1.0d0/))
					smnext(1,i)=smAs(1,period,i) 
					smnext(2,i)=SS(2,period,i)
					smnext(3,i)=smexperience(period,i)+smh
				!-------------------------- 3.b PERIOD(15,21): TWO CHILDREN ------------------------------
				else
					! increase the age2 by one if it is different than zero
					SS(5,period,i)=SS(5,period,i)+1.0d0
					
					! both kids are above baby age.
					baby=0
					! because period is above astar, first kid does not get inputs anymore, no x choice
					uc=uctwo(SS(1,period,i),SS(2,period,i),parU(1)) 
					! NOTE Depending on the age of the second child, we need modify how much of the income goes to the
					! production of skills. This depends on the birth history, so good thing we have birthhist, which recorded
					! the timing of the second birth (if any) as the period number it happened. if period is later than second's
					! kids departure from mother's care (periodage 15 and thereafter)
					
					incrat=0.5d0
					if (period>=astar+birthhist(i)-1) incrat=1.0d0

					umat(1,i)=	uc + 	 parU(4)*(wageh(i)*incrat)**parU(5)+parU(6)*baby
					umat(2,i)= uc +a1type(a1holder(i))*0.5d0+ parU(4)*(incrat*(wageh(i)+0.5d0*wage(i)))**parU(5)+parU(6)*baby+parU(7)*0.5d0*baby
					umat(3,i)= uc +a1type(a1holder(i))*1.0d0+ parU(4)*(incrat*(wageh(i)+wage(i)))**parU(5)+parU(6)*baby+parU(7)*baby
					! now calculate and add future values.

					do k=1,nttypes
						fvconsw(k)=sum(wcoef(4:7,period+1,birthhist(i),k)*intcons)+wcoef(Gsize+1,period+1,birthhist(i),k)
					end do
					
					do j=1,3
						Enext=SS(3,period,i)+(j-1)*0.5d0
						income=wageh(i)*incrat+wage(i)*(j-1)*0.5*incrat ! it is 0.25 because half the income goes to goods inputs
						inputs=(/SS(1,period,i),SS(2,period,i),0.d0,(1-0.25*(j-1))*(1-xgrid(k)),income/)
						Anext=pftwo(inputs,omega3(1:3),SS(4:5,period,i),parA(4:12), rho)
						A1next=SS(1,period,i) 	! don't change A1
						A2next=Anext(2) 		
						if (period>=astar+birthhist(i)-1) A2next=SS(2,period,i)
						intw=(/A1next,A2next,Enext,A1next**2,A2next**2,Enext**2, A1next*A2next,A1next*Enext,A2next*Enext/)
						fvw=sum(intw*(/wcoef(1:3,period+1,birthhist(i),a1holder(i)),wcoef(9:Gsize,period+1,birthhist(i),a1holder(i))/))
						umat(j,i)=umat(j,i)+beta*(fvw+fvconsw(a1holder(i)))
						nextcollect(:,j)=(/A1next,A2next,Enext/)
					end do
					picker=maxloc(umat(:,i))
					choices(1,period,i)=picker(1)
					next(:,i)=nextcollect(:,choices(1,period,i))
					
					! NEW: Smoothed As and Es and choices and later test scores (not this period, though)
					smdenom=sum(exp((umat(:,i)-maxval(umat(:,i)))/smpar))
					smchoices(:,period,i)=exp((umat(:,i)-maxval(umat(:,i)))/smpar)/smdenom
					smh=sum(smchoices(:,period,i)*(/0.d0,0.5d0,1.0d0/))
					smnext(1,i)=smAs(1,period,i) 
					smnext(2,i)=SS(2,period,i)
					smnext(3,i)=smexperience(period,i)+smh
					income=(wageh(i)+wage(i)*smh)*incrat
					inputs=(/smAs(1,period,i),smAs(2,period,i),(1-0.5*smh)*smx,(1-0.5*smh)*(1-smx),income/)
					smnext(1:2,i)=pftwo(inputs,omega3(1:3),SS(4:5,period,i),parA(4:12), rho)
					smnext(1,i)=smAs(1,period,i)
					smnext(3,i)=smexperience(period,i)+smh
				end if  ! end of period<8 WV if
			end do
			!#%#%#%#%#%#%#%#%#%#%#%#%%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%%#%#%#%#%%#%#

			! NO SMOOTHED TEST SCORES OR EVEN STATE SPACE FOR THE FINAL PERIOD.

		
			! ##################################################################3###################3

			!-------------------------4. FINAL PERIOD-----------------------------
		else
			! This is the final period. Future values are not calculated but terminal values are used.
			! I have been using these temporary terminal values for the model, it seems like I will stick to those 
			! for the jmpsimple. they will have three choices as before. 
			do i=1,Npaths
				! calculations are different for V and W.  
			
				! --------------------------- 4.a PERIOD 22: ONE CHILD-----------------------------------
				! 
				if (SS(5,period,i)<0.001) then 		! if still with one child. 
					! there is no chance of having a baby at this point
					baby=0
					uc=parU(2)*SS(1,period,i) 
					! NOTE: important difference: all the income goes into the consumption now, because there is no production of
					! skils
					umat(1,i)=uc+parU(4)*(wageh(i))**parU(5)+parU(6)*baby
					umat(2,i)=uc+a1type(a1holder(i))*0.5d0+ parU(4)*(wageh(i)+0.5d0*wage(i))**parU(5)+parU(6)*baby+parU(7)*0.5d0*baby
					umat(3,i)=uc+a1type(a1holder(i))*1.0d0+ parU(4)*(wageh(i)+wage(i))**parU(5)+parU(6)*baby+parU(7)*baby
						
					! now calculate future values for the possible 3 decisions, so we can calculate the future value.
					do j=1,3
						A1next=SS(1,period,i)
						A2next=SS(2,period,i)
						Enext=SS(3,period,i)+(j-1)*0.5d0
						! now add the expected future value to corresponding umat element
						umat(j,i)=umat(j,i)+beta*umat(j,i)*mult
						nextcollect(:,j)=(/A1next,A2next,Enext/)
					end do
					! finally, pick the highest expected utility one.
					picker=maxloc(umat(:,i))
					choices(1,period,i)=picker(1)
					! get the optimized future state space points to use them in future.
					next(:,i)=nextcollect(:,choices(1,period,i))
				
				!-------------------------- 4.b PERIOD 22: TWO CHILDREN ------------------------------
				else
					! increase the age2 by one if it is different than zero
					SS(5,period,i)=SS(5,period,i)+1.0d0
					
					! both kids are above baby age.
					baby=0
					! because period is above astar, first kid does not get inputs anymore, no x choice
					uc=uctwo(SS(1,period,i),SS(2,period,i),parU(1)) 
					! NOTE Depending on the age of the second child, we need modify how much of the income goes to the
					! production of skills. This depends on the birth history, so good thing we have birthhist, which recorded
					! the timing of the second birth (if any) as the period number it happened. if period is later than second's
					! kids departure from mother's care (periodage 15 and thereafter)
					
					incrat=0.5d0
					if (period>=astar+birthhist(i)-1) incrat=1.0d0

					umat(1,i)=	uc + 	 parU(4)*(wageh(i)*incrat)**parU(5)+parU(6)*baby
					umat(2,i)= uc +a1type(a1holder(i))*0.5d0+ parU(4)*(incrat*(wageh(i)+0.5d0*wage(i)))**parU(5)+parU(6)*baby+parU(7)*0.5d0*baby
					umat(3,i)= uc +a1type(a1holder(i))*1.0d0+ parU(4)*(incrat*(wageh(i)+wage(i)))**parU(5)+parU(6)*baby+parU(7)*baby
					! now calculate and add future values.
					do j=1,3
						Enext=SS(3,period,i)+(j-1)*0.5d0
						income=wageh(i)+wage(i)*(j-1)*0.25 ! it is 0.25 because half the income goes to goods inputs
						A1next=SS(1,period,i) 	! don't change A1
						A2next=SS(2,period,i) 	! don't change A2
						umat(j,i)=umat(j,i)+beta*umat(j,i)*mult
						nextcollect(:,j)=(/A1next,A2next,Enext/)
					end do
					picker=maxloc(umat(:,i))
					choices(1,period,i)=picker(1)
					next(:,i)=nextcollect(:,choices(1,period,i))
				end if  ! end of period<8 WV if
			end do ! end of the Npath loop within final period calc
			!#%#%#%#%#%#%#%#%#%#%#%#%%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%%#%#%#%#%%#%#


		end if ! end of if that divides the the periods

	end do !--------------- END PERIOD LOOP---------------
	
end subroutine simhist



! aux routines for moments
! calculate the mean, variance of a matrix 
subroutine mv(datamat,outputvec)
	implicit none
	real(dble), intent(in) :: datamat(:) 	!< data vec  holding the data
	real(dble), intent(out) :: outputvec(2) 	!< output vec: first means, then variances
	integer n,k,i
	n=size(datamat,1)
	!do i=1,n
		!print*, i, datamat(i)
	!end do
	outputvec(1)=sum(datamat)/n
	outputvec(2)= sum((datamat-outputvec(1))**2)/(n-1)
	!	if (i<k) outputvec(2*k+i)=sum((datamat(:,i)-outputvec(i))*(datamat(:,i+1)-outputvec(i+1)))/(n-1)
end subroutine mv

!> calculate consecutive covariances of variables stored in the the cols of
!>datamat, 2 columns only
subroutine concov(datamat,output)
	implicit none
	real(dble), intent(in) :: datamat(:,:)
	real(dble), intent(out) :: output
	real(dble) :: mean1, mean2
	real(dble) n
	n=size(datamat,1)*1.0d0
	mean1=sum(datamat(:,1))/n
	mean2=sum(datamat(:,2))/n
	output=sum((datamat(:,1)-mean1)*(datamat(:,2)-mean2))/(n-1)
end subroutine concov



!subroutine moments(momentvec, SS,outcomes, choices, testoutcomes, birthhist, smchoices, smexperience, omega3data,lfpperiods, expperiods, tsperiods, idmat)
subroutine moments(momentvec, SS,smtestoutcomes, birthhist, smchoices, smexperience, omega3data,lfpperiods, expperiods, tsperiods, idmat)
	implicit none
	real(dble),intent(in):: SS(6,nperiods,Npaths,SampleSize) 		!< Simulated State space (A1,A2,E,age1,age2,agem)xnperiods,Npaths	
	!real(dble),intent(in):: outcomes(2,nperiods,Npaths,SampleSize) !< outcomes: wages of the father and mother.
	!integer,intent(in):: choices(1,nperiods,Npaths,SampleSize) 	!< choices : the choice history of h.
	integer,intent(in):: birthhist(Npaths,SampleSize)  		    !< the birth timing vec, 0 if one child throughout.
	real(dble), intent(in) :: omega3data(o3size, SampleSize) 	!< holds the omega3 data for Sample people
	real(dble), intent(in) ::smtestoutcomes(4,Ntestage,Npaths,SampleSize)  !< holds test scores for two children 
	real(dble), intent(in) :: smchoices(3,nperiods,Npaths,SampleSize)
	real(dble), intent(in) :: smexperience(nperiods,Npaths,SampleSize)
	integer, intent(in):: expperiods(expsize) 		!< the array that holds the ages of mom for which experience moments calculated
	integer, intent(in):: lfpperiods(lfpsize) 			!< the array that holds the period numbers for which labor for particapation equations are estimated.
	integer, intent(in):: tsperiods(expsize) 		!< the array that holds the ages for which test score averages are calculated. 
	integer, intent(in):: idmat(SampleSize,idmatsize) 	   		!<indicates the which sample units are in the jth columnth moment
														!<calculation
	real(dble), intent(out) :: momentvec(MomentSize) 	!< The set of moments produced by the simulated data
		
	! NOTE  					 ------------ON IDMAT---------------
	! * first size(lfpperiods) of the idmat is for the periods at which participation equations are estimated.
	! * next size(expperiods) is for the mean and variance calculations of experience levels.
	! * following size(expperiods)-1 is for the cov(E_t,E_t+1) calculations. marking the exp pairs. 
	! * Next is the test score averages. We will deals with tsperiods of them. And stick only to one of the tests.
	! idmat has a column for each tsperiod. It takes 4 VALUES. 0=not in the data, 1= first child 2=second child 3=both children
	!	so we can figure out for which kid(s) to add the test score data into the series (to calculate average).
	! * last we have one each period in each testminage, .., testmaxage
	! SO IDMAT HAS A SIZE = SAMPLE SIZE  x ( lfpsize + expsize + expsize-1 + tsperiods + Ntestage )



	! locals 
	integer i,j,k,l, counter, counter2, baby, twochild
	
	! create the regressions matrix for participation equations. the first size(lfpperiods) of the idmat holds the indicators
	! for the involvement of the sample members in a particular period for the participation equations. example: if
	! lfpperiods=(/3,6,9,12/), the first four columns of idmat holds the indicators for the data for these periods. By summing
	! up all the columns, we count the total number sample points that should be entering this regressions, including the time
	! dimension. Each of these will be associated with Npaths simulated observations
	integer regsample
	real(dble) lfpmat(sum(idmat(:,1:lfpsize))*Npaths,nreglfp+1) 	
	real(dble) fulltimevec(sum(idmat(:,1:lfpsize))*Npaths)
	real(dble) parttimevec(sum(idmat(:,1:lfpsize))*Npaths)
	real(dble) parttimecoef(nreglfp+1)
	real(dble) fulltimecoef(nreglfp+1)
	
	! for experience calculations, create a matrix that will hold the relevant experience observations
	real(dble),allocatable:: experience(:), dualexperience(:,:)
	real(dble) expest(expsize*2+expsize-1) 	! hold the mean and variances of experience levels. mean1,var1, mean2, var2,.. then covariances
	real(dble) mvvec(2), covexp
	
	! for average test scores
	real(dble), allocatable:: ts(:)
	real(dble) tsest(tssize)
	integer nspots

	! for test score difference equations
	real(dble), allocatable:: tsdiffmat(:,:)
	real(dble), allocatable:: tsdiffvec(:)
	real(dble) tsdiffest(nregtsdiff+1)
	real(dble) smh, smhnext, A1, A2, A1F, A2F

	! lapack stuff
	real(dble) work(nreglfp+1+(nreglfp+1)*blocksize)
	real(dble) worktsdiff(nregtsdiff+1+(nregtsdiff+1)*blocksize)
	integer info, infotsdiff
		
	! initalize stuff
	regsample=sum(idmat(:,1:lfpsize))*Npaths ! participation equation big sample size
	counter=1
	baby=0
	twochild=0	
	! INDICATOR MATRIX IDMAT SHOULD FOLLOW THIS ORDER
	! STARTING OUT WITH 6A and 6B in jmpsimple document. 

	!---------------------------------PARTICIPATION EQUATIONS------------------------------------------
	! linear regressions with work status, where regressor vector is [schooling AFQT agem baby E 2child 1]
	
	! first fill in 1.0d0 s to the last column. 	
	lfpmat(:,nreglfp+1)=1.d0
	do k=1, lfpsize
		do l=1,SampleSize
			if ( idmat(l,k) == 1 ) then
				do i=1, Npaths
					if ((SS(4,lfpperiods(k),i,l)<3.1) .OR. (SS(5,lfpperiods(k),i,l)<3.1)) baby=1
					if (SS(5,lfpperiods(k),i,l)>0.1) twochild=1
					! fill in the regressor matrix
					lfpmat(counter,1:2)=omega3data(1:2,l)
					lfpmat(counter,3)=omega3data(3,l)+lfpperiods(k)
					lfpmat(counter,4)=baby*1.0d0
					!lfpmat(counter,5)=SS(3,lfpperiods(k),i,l)
					lfpmat(counter,5)=smexperience(lfpperiods(k),i,l)
					lfpmat(counter,6)=twochild*1.0d0
					! fill in the dependent variable	
					!if (choices(1,lfpperiods(k),i,l) > 2.5) fulltimevec(counter)=1.0d0
					!if ( (choices(1,lfpperiods(k),i,l) < 2.5) .AND. (choices(1,lfpperiods(k),i,l) > 1.5) ) parttimevec(counter)=1.0d0
					fulltimevec(counter)=smchoices(3,lfpperiods(k),i,l)
					parttimevec(counter)=smchoices(2,lfpperiods(k),i,l)
					counter=counter+1
				end do
			end if
		end do
	end do

	! now run the lapack to get the regression coefficients

	call DGELS('N', regsample, nreglfp+1,1,lfpmat,regsample, fulltimevec, regsample,work,nreglfp+1+(nreglfp+1)*blocksize,info)
	fulltimecoef=fulltimevec(1:nreglfp+1)
	call DGELS('N', regsample, nreglfp+1,1,lfpmat,regsample, parttimevec, regsample,work,nreglfp+1+(nreglfp+1)*blocksize,info)
	parttimecoef=parttimevec(1:nreglfp+1)
	! NOTE: Do I need to go for the variance of error term? 
	

	!-----------------------------EXPERIENCE-----------------------------
	! EXPERIENCE LEVELS: mean experience levels at age 25, 30, 35, 40. Variance at the same ages. Covariance between age 25,30,
	! 30,35 and 35,40 wages.
	! prepare the vectors first. idmat's following 4 columns mark the people who entered into these calculations. 
	
	! THE EXPPERIODS ARE AGES OF THE MOTHER, NOT PERIODS, SO WE NEED TO PICK THE CORRECT PERIOD with OMEGA3(3)
	do k=1, expsize
		counter=1
		allocate( experience ( sum (idmat(:,lfpsize+k)*Npaths ) ) )
		do l=1, SampleSize
			if (idmat(l,lfpsize+k)==1) then
				do i=1,Npaths
					!experience(counter)=SS(3,expperiods(k)+omega3data(3,l),i,l)
					experience(counter)=smexperience(expperiods(k)-nint(omega3data(3,l)+1),i,l)
					counter=counter+1
				end do
			end if
		end do
		! now we have the correct experience vector, calculate mean and the variance.
		call mv(experience,mvvec) 
		expest((k-1)*2+1:2*k)=mvvec
		deallocate (experience)
	end do

	! now covariance of experiences.
	do k=1, expsize-1
		counter=1
		allocate ( dualexperience( sum(idmat(:,lfpsize+expsize+k))*Npaths,2 ) )
		do l=1, SampleSize
			if (idmat(l,lfpsize+expsize+k)==1) then
				do i=1, Npaths
					!dualexperience(counter,1)=SS(3,expperiods(k)+omega3data(3,l),i,l)
					!dualexperience(counter,2)=SS(3,expperiods(k+1)+omega3data(3,l),i,l)
					dualexperience(counter,1)=smexperience(expperiods(k)-nint(omega3data(3,l)+1),i,l)
					dualexperience(counter,2)=smexperience(expperiods(k+1)-nint(omega3data(3,l)+1),i,l)
					counter=counter+1
				end do
			end if
		end do
		! calculate covariance
		call concov(dualexperience,covexp)
		expest((2*expsize+k))=covexp	
		! deallocate
		deallocate (dualexperience)
	end do
	
	!-------------------------- TEST SCORE LEVELS -----------------------------------
	! test score levels will help to pin down gammahat parameters. match levels at different ages.  
	! NOTE: CONTENT OF THE tsperiods is the test score numbers. For example, if first test score is at 5 years old, and the last one
	! is at 14, "1" means at age 5, "5" means age 9. so if somebody wants to match the mean test scores at 6, 8, 10 and 12, they
	! need to enter tsperiods=(/2,4,6,12/). this is because testoutcome matrices start from has a lower dim than SS, just enough
	! for test scores

	do k=1, tssize
		counter=1
		nspots=0
		! figure out how many places in the test score vector you will need.
		do l=1,SampleSize
			if ( (idmat(l,lfpsize+2*expsize-1+k)>0)) nspots=nspots+1 ! if kids exist add one
			if (idmat(l,lfpsize+2*expsize-1+k)==3) nspots=nspots+1 ! one more if two
		end do
		
		allocate (ts (nspots*Npaths))
		
		do l=1, SampleSize
			if ( (idmat(l,lfpsize+2*expsize-1+k)==1) .OR. (idmat(l,lfpsize+2*expsize-1+k)==3)) then  ! first child
				do i=1, Npaths
					if (smtestoutcomes(1, tsperiods(k), i, l) > -1.0d0 ) ts(counter)=smtestoutcomes(1, tsperiods(k), i, l)
					counter=counter+1
				end do
			end if
			! second child: if mom has a second child with test score associated with the AGE tsperiod(k)
			if  (idmat(l,lfpsize+2*expsize-1+k)==2) then 				
				do i=1, Npaths 
					! check if at age2=tsperiod(k) [which corresponds to period tsperiod(k)+birthhist(i,l)-1 for mom]
					! is a period for the ith simulation has a test score. If it does, it is the one:
					! it is tsperiod+birthist-1 period for mom, the age of the second kid is tsperiod
					if (smtestoutcomes(3, tsperiods(k)+birthhist(i,l)-1, i, l) > -1.0d0 ) ts(counter)=smtestoutcomes(3, tsperiods(k)+birthhist(i,l)-1, i, l)
					counter=counter+1
				end do
			end if
		end do
		! now, calculate the average of the ts vector and put it tsest, then deallocate
		tsest(k)=sum(ts)/nspots
		deallocate(ts)
	end do

	!-------------------------TEST SCORE DIFFERENCE REGRESSIONS----------------------
	! we will go through the number of test periods one by one. Starting from age period 6 (because this is the earliest we can have
	! test scores for a pair of kids. Maximum is testmaxage-2. So idmat should have testmaxage-testminage elements showing us who
	! are in the calculation for each period.
	
	! first count the simulated observations that will go into this
	nspots=0
	
	do k=testminage+1,testmaxage-2
		do l=1, SampleSize
			if ( idmat(l,lfpsize+2*expsize-1+tssize+k)==1) then ! if this mom is in calculation for period k...
				do i=1,Npaths 	! ... go through simulated outcomes for that mom and figure which ones have test scores for both kids at that period.
					! if kid 2 is older than 5, then create a spot
					if (SS(5,k,i,l) > 4.8)  nspots=nspots+1
				end do
			end if
		end do
	end do

	! now we can allocate data matrix and dependent variable vectors
	allocate(tsdiffmat (nspots, nregtsdiff+1) )
	allocate(tsdiffvec (nspots) )
	
	tsdiffmat(:,nregtsdiff+1)=1.0d0
	! another loop to fill in these guys
	counter=1
	do k=testminage+1,testmaxage-2
		do l=1, SampleSize
			if ( idmat(l,lfpsize+2*expsize-1+tssize+k)==1) then ! if this mom is in calculation for period k...
				do i=1,Npaths 	! ... go through simulated outcomes for that mom and pick the right ones
					! 1. fill in the labor supply related parts (smh)
					!if ( smtestoutcomes(3, k-testminage+1,i,l) > -1.0d0 )   then ! if the second kid is old enough to have a testscore (first is kid is already at most testmaxage-2
					if (SS(5,k,i,l)>4.8) then 
						smh=sum(smchoices(:,k,i,l)*(/0.0d0,0.5d0,1.0d0/))
						smhnext=sum(smchoices(:,k+1,i,l)*(/0.0d0,0.5d0,1.0d0/))
						tsdiffmat(counter, 1:2)=(/smh,smhnext/)
						tsdiffmat(counter, 3:5)=smh*omega3data(1:3,l)
						tsdiffmat(counter, 6:8)=smhnext*omega3data(1:3,l)
						tsdiffmat(counter,9:11)=omega3data(1:3,l)
						tsdiffmat(counter,12:13)=SS(4:5,k,i,l)
						! now fillin the tsdiffmat vec
						A1=smtestoutcomes(1,k-testminage+1,i,l)
						A2=smtestoutcomes(3,k-testminage+1,i,l)
						A1f=smtestoutcomes(1,k-testminage+3,i,l)
						A2f=smtestoutcomes(3,k-testminage+3,i,l)
						tsdiffvec(counter)= A1f-A1+A2f-A2
						counter=counter+1
					end if
				end do
			end if
		end do
	end do
	!print*, tsdiffmat(:,3:5)


	! lapack
	call DGELS('N', nspots, nregtsdiff+1,1,tsdiffmat,nspots, tsdiffvec, nspots,worktsdiff,nregtsdiff+1+(nregtsdiff+1)*blocksize,infotsdiff)
	tsdiffest=tsdiffvec(1:nregtsdiff+1)
	! deallocate
	deallocate(tsdiffvec)
	deallocate(tsdiffmat)
	! put all the vecs together
	momentvec=(/fulltimecoef, parttimecoef, expest, tsest,tsdiffest/)

end subroutine moments



!>Subroutine moments(momentvec, SS,outcomes, choices, testoutcomes, birthhist, smchoices, smexperience, omega3data,lfpperiods, expperiods, tsperiods, idmat)
!>This is the alternative, because regression coefficients are estimated for
!>non-stacked samples and then averaged.
subroutine momentsalt(momentvec, SS,smtestoutcomes, birthhist, smchoices, smexperience, omega3data,lfpperiods, expperiods, tsperiods, idmat)
	implicit none
	real(dble),intent(in):: SS(6,nperiods,Npaths,SampleSize) 		!< Simulated State space (A1,A2,E,age1,age2,agem)xnperiods,Npaths	
	!real(dble),intent(in):: outcomes(2,nperiods,Npaths,SampleSize) !< outcomes: wages of the father and mother.
	!integer,intent(in):: choices(1,nperiods,Npaths,SampleSize) 	!< choices : the choice history of h.
	integer,intent(in):: birthhist(Npaths,SampleSize)  		    !< the birth timing vec, 0 if one child throughout.
	real(dble), intent(in) :: omega3data(o3size, SampleSize) 	!< holds the omega3 data for Sample people
	real(dble), intent(in) ::smtestoutcomes(4,Ntestage,Npaths,SampleSize)  !< holds test scores for two children 
	real(dble), intent(in) :: smchoices(3,nperiods,Npaths,SampleSize)
	real(dble), intent(in) :: smexperience(nperiods,Npaths,SampleSize)
	integer, intent(in):: expperiods(expsize) 		!< the array that holds the ages of mom for which experience moments calculated
	integer, intent(in):: lfpperiods(lfpsize) 			!< the array that holds the period numbers for which labor for particapation equations are estimated.
	integer, intent(in):: tsperiods(expsize) 		!< the array that holds the ages for which test score averages are calculated. 
	integer, intent(in):: idmat(SampleSize,idmatsize) 	   		!<indicates the which sample units are in the jth columnth moment
														!<calculation
	real(dble), intent(out) :: momentvec(MomentSize) 	!< The set of moments produced by the simulated data
		
	! NOTE  					 ------------ON IDMAT---------------
	! * first size(lfpperiods) of the idmat is for the periods at which participation equations are estimated.
	! * next 2*size(expperiods) is for the mean and variance calculations of experience levels.
	! * following size(expperiods)-1 is for the cov(E_t,E_t+1) calculations. marking the exp pairs. 
	! * Next is the test score averages. We will deals with tsperiods of them. And stick only to one of the tests.
	! idmat has a column for each tsperiod. It takes 4 VALUES. 0=not in the data, 1= first child 2=second child 3=both children
	!	so we can figure out for which kid(s) to add the test score data into the series (to calculate average).
	! * last we have one each period in each testminage, .., testmaxage
	! SO IDMAT HAS A SIZE = SAMPLE SIZE  x ( lfpsize + expsize + expsize-1 + tsperiods + Ntestage )



	! locals 
	integer i,j,k,l, counter, counter2, baby, twochild
	
	! create the regressions matrix for participation equations. the first size(lfpperiods) of the idmat holds the indicators
	! for the involvement of the sample members in a particular period for the participation equations. example: if
	! lfpperiods=(/3,6,9,12/), the first four columns of idmat holds the indicators for the data for these periods. By summing
	! up all the columns, we count the total number sample points that should be entering this regressions, including the time
	! dimension. Each of these will be associated with Npaths simulated observations
	integer regsample, aregsample
	real(dble) lfpmat(sum(idmat(:,1:lfpsize))*Npaths,nreglfp+1) 	
	! new: Smaller matrix
	real(dble) alfpmat(sum(idmat(:,1:lfpsize)),nreglfp+1) 	
	
	real(dble) fulltimevec(sum(idmat(:,1:lfpsize))*Npaths)
	real(dble) parttimevec(sum(idmat(:,1:lfpsize))*Npaths)
	
	real(dble) afulltimevec(sum(idmat(:,1:lfpsize)))
	real(dble) aparttimevec(sum(idmat(:,1:lfpsize)))
	
	real(dble) parttimecoef(nreglfp+1)
	real(dble) fulltimecoef(nreglfp+1)
	real(dble) parttimecoefvec(Npaths,nreglfp+1)
	real(dble) fulltimecoefvec(Npaths,nreglfp+1)
	
	! for experience calculations, create a matrix that will hold the relevant experience observations
	real(dble),allocatable:: experience(:), dualexperience(:,:)
	real(dble) expest(expsize*2+expsize-1) 	! hold the mean and variances of experience levels. mean1,var1, mean2, var2,.. then covariances
	real(dble) expestmat(Npaths,expsize*2+expsize-1) 	! matrix to hold the calculations for each simulation
	real(dble) mvvec(2), covexp
	
	! for average test scores
	real(dble), allocatable:: ts(:)
	real(dble) tsest(tssize)
	integer nspots

	! for test score difference equations
	real(dble), allocatable:: tsdiffmat(:,:)
	real(dble), allocatable:: tsdiffvec(:)
	real(dble) tsdiffest(nregtsdiff+1), tsdiffestmat(Npaths,nregtsdiff+1)


	real(dble) smh, smhnext, A1, A2, A1F, A2F

	! lapack stuff
	real(dble) work(nreglfp+1+(nreglfp+1)*blocksize)
	real(dble) worktsdiff(nregtsdiff+1+(nregtsdiff+1)*blocksize)
	integer info, infotsdiff


	! initalize stuff
	regsample=sum(idmat(:,1:lfpsize))*Npaths ! participation equation big sample size
	aregsample=sum(idmat(:,1:lfpsize)) 		 ! participation equation small non stacked
	
	counter=1
	baby=0
	twochild=0	
	! INDICATOR MATRIX IDMAT SHOULD FOLLOW THIS ORDER
	! STARTING OUT WITH 6A and 6B in jmpsimple document. 

	!---------------------------------PARTICIPATION EQUATIONS------------------------------------------
	! linear regressions with work status, where regressor vector is [schooling AFQT agem baby E 2child 1]

	do i=1, Npaths
		counter=1
		! first fill in 1.0d0 s to the last column. 	
		alfpmat(:,nreglfp+1)=1.0d0
		do k=1, lfpsize
			do l=1, lfpsize
				if (idmat(l,k)==1) then
					if ((SS(4,lfpperiods(k),i,l)<3.1) .OR. (SS(5,lfpperiods(k),i,l)<3.1)) baby=1
					if (SS(5,lfpperiods(k),i,l)>0.1) twochild=1
					! fill in the regressor matrix
					alfpmat(counter,1:2)=omega3data(1:2,l)
					alfpmat(counter,3)=omega3data(3,l)+lfpperiods(k)
					alfpmat(counter,4)=baby*1.0d0
					!lfpmat(counter,5)=SS(3,lfpperiods(k),i,l)
					alfpmat(counter,5)=smexperience(lfpperiods(k),i,l)
					alfpmat(counter,6)=twochild*1.0d0
					! fill in the dependent variable	
					!if (choices(1,lfpperiods(k),i,l) > 2.5) fulltimevec(counter)=1.0d0
					!if ( (choices(1,lfpperiods(k),i,l) < 2.5) .AND. (choices(1,lfpperiods(k),i,l) > 1.5) ) parttimevec(counter)=1.0d0
					afulltimevec(counter)=smchoices(3,lfpperiods(k),i,l)
					aparttimevec(counter)=smchoices(2,lfpperiods(k),i,l)
					counter=counter+1
				end if
			end do
		end do
		
		call DGELS('N', aregsample, nreglfp+1,1,alfpmat,aregsample, afulltimevec, aregsample,work,nreglfp+1+(nreglfp+1)*blocksize,info)
		fulltimecoefvec(i,:)=afulltimevec(1:nreglfp+1)

		call DGELS('N', aregsample, nreglfp+1,1,alfpmat,aregsample, aparttimevec, aregsample,work,nreglfp+1+(nreglfp+1)*blocksize,info)
		parttimecoefvec(i,:)=aparttimevec(1:nreglfp+1)
	end do
	
	fulltimecoef=sum(fulltimecoefvec,1)/Npaths	
	parttimecoef=sum(parttimecoefvec,1)/Npaths	
	
! NOTE: Do I need to go for the variance of error term? 
	

	!-----------------------------EXPERIENCE-----------------------------
	! EXPERIENCE LEVELS: mean experience levels at age 25, 30, 35, 40. Variance at the same ages. Covariance between age 25,30,
	! 30,35 and 35,40 wages.
	! prepare the vectors first. idmat's following 4 columns mark the people who entered into these calculations. 
	

	! New Way: calculating everything for each path and then averaging them out.
	do i=1, Npaths
		! THE EXPPERIODS ARE AGES OF THE MOTHER, NOT PERIODS, SO WE NEED TO PICK THE CORRECT PERIOD with OMEGA3(3)
		do k=1, expsize
			counter=1
			allocate( experience ( sum (idmat(:,lfpsize+k) ) ) )
			do l=1, SampleSize
				if (idmat(l,lfpsize+k)==1) then
					!experience(counter)=SS(3,expperiods(k)+omega3data(3,l),i,l)
					experience(counter)=smexperience(expperiods(k)-nint(omega3data(3,l)+1),i,l)
					counter=counter+1
				end if
			end do
			! calculate the mean and the variance for the ith path and kth experience level. 
			call mv(experience,mvvec) 
			expestmat(i,(k-1)*2+1:2*k)=mvvec
			deallocate (experience)
		end do
		
		! now the covariances
		do k=1, expsize-1
			counter=1
			allocate ( dualexperience( sum(idmat(:,lfpsize+expsize+k)),2 ) )
			do l=1, SampleSize
				if (idmat(l,lfpsize+expsize+k)==1) then
					!dualexperience(counter,1)=SS(3,expperiods(k)+omega3data(3,l),i,l)
					!dualexperience(counter,2)=SS(3,expperiods(k+1)+omega3data(3,l),i,l)
					dualexperience(counter,1)=smexperience(expperiods(k)-nint(omega3data(3,l)+1),i,l)
					dualexperience(counter,2)=smexperience(expperiods(k+1)-nint(omega3data(3,l)+1),i,l)
					counter=counter+1
				end if
			end do
			! calculate covariance
			call concov(dualexperience,covexp)
			expestmat(i,(2*expsize+k))=covexp	
			! deallocate
			deallocate (dualexperience)
		end do
	end do
	! now averaging it out.
	expest=sum(expestmat,1)/Npaths

	!-------------------------- TEST SCORE LEVELS -----------------------------------
	! test score levels will help to pin down gammahat parameters. match levels at different ages.  
	! NOTE: CONTENT OF THE tsperiods is the test score numbers. For example, if first test score is at 5 years old, and the last one
	! is at 14, "1" means at age 5, "5" means age 9. so if somebody wants to match the mean test scores at 6, 8, 10 and 12, they
	! need to enter tsperiods=(/2,4,6,12/). this is because testoutcome matrices start from has a lower dim than SS, just enough
	! for test scores

	! I WILL NOT CHANGE THIS: (1) It is just the average, so it does not matter how you sum it up. (2) I am not even sure I will be
	! using this as a moment condition.
	do k=1, tssize
		counter=1
		nspots=0
		! figure out how many places in the test score vector you will need.
		do l=1,SampleSize
			if ( (idmat(l,lfpsize+2*expsize-1+k)>0)) nspots=nspots+1 ! if kids exist add one
			if (idmat(l,lfpsize+2*expsize-1+k)==3) nspots=nspots+1 ! one more if two
		end do
		
		allocate (ts (nspots*Npaths))
		
		do l=1, SampleSize
			if ( (idmat(l,lfpsize+2*expsize-1+k)==1) .OR. (idmat(l,lfpsize+2*expsize-1+k)==3)) then  ! first child
				do i=1, Npaths
					if (smtestoutcomes(1, tsperiods(k), i, l) > -1.0d0 ) ts(counter)=smtestoutcomes(1, tsperiods(k), i, l)
					counter=counter+1
				end do
			end if
			! second child: if mom has a second child with test score associated with the AGE tsperiod(k)
			if  (idmat(l,lfpsize+2*expsize-1+k)==2) then 				
				do i=1, Npaths 
					! check if at age2=tsperiod(k) [which corresponds to period tsperiod(k)+birthhist(i,l)-1 for mom]
					! is a period for the ith simulation has a test score. If it does, it is the one:
					! it is tsperiod+birthist-1 period for mom, the age of the second kid is tsperiod
					if (smtestoutcomes(3, tsperiods(k)+birthhist(i,l)-1, i, l) > -1.0d0 ) ts(counter)=smtestoutcomes(3, tsperiods(k)+birthhist(i,l)-1, i, l)
					counter=counter+1
				end do
			end if
		end do
		! now, calculate the average of the ts vector and put it tsest, then deallocate
		tsest(k)=sum(ts)/nspots
		deallocate(ts)
	end do

	!-------------------------TEST SCORE DIFFERENCE REGRESSIONS----------------------
	! we will go through the number of test periods one by one. Starting from age period 6 (because this is the earliest we can have
	! test scores for a pair of kids. Maximum is testmaxage-2. So idmat should have testmaxage-testminage elements showing us who
	! are in the calculation for each period.
	do i=1,Npaths
		nspots=0
		do k=1,Noldtest
			do l=1, SampleSize
				!if (l<10) print*, 'midreport:', l, lfpsize+2*expsize-1+tssize+k, idmat(l,lfpsize+2*expsize-1+tssize+k), SS(5,k+10,i,l)   
				if ( idmat(l,lfpsize+2*expsize-1+tssize+k)==1) then ! if this mom is in calculation for period k...
					! if kid 2 is older than 5, then create a spot
					if (SS(5,k+10,i,l) > 4.8) nspots=nspots+1
				end if
			end do
		end do


		! now we can allocate data matrix and dependent variable vectors
		allocate(tsdiffmat (nspots, nregtsdiff+1) )
		allocate(tsdiffvec (nspots) )
	
	
		tsdiffmat(:,nregtsdiff+1)=1.0d0
		! another loop to fill in these guys
		counter=1
		do k=testminage+1,testmaxage-2
			do l=1, SampleSize
				if ( idmat(l,lfpsize+2*expsize-1+tssize+k)==1) then ! if this mom is in calculation for period k...
					! 1. fill in the labor supply related parts (smh)
					!if ( smtestoutcomes(3, k-testminage+1,i,l) > -1.0d0 )   then ! if the second kid is old enough to have a testscore (first is kid is already at most testmaxage-2
					if (SS(5,k,i,l)>4.8) then 
						smh=sum(smchoices(:,k,i,l)*(/0.0d0,0.5d0,1.0d0/))
						smhnext=sum(smchoices(:,k+1,i,l)*(/0.0d0,0.5d0,1.0d0/))
						tsdiffmat(counter, 1:2)=(/smh,smhnext/)
						tsdiffmat(counter, 3:5)=smh*omega3data(1:3,l)
						tsdiffmat(counter, 6:8)=smhnext*omega3data(1:3,l)
						tsdiffmat(counter,9:11)=omega3data(1:3,l)
						tsdiffmat(counter,12:13)=SS(4:5,k,i,l)
						! now fillin the tsdiffmat vec
						A1=smtestoutcomes(1,k-testminage+1,i,l)
						A2=smtestoutcomes(3,k-testminage+1,i,l)
						A1f=smtestoutcomes(1,k-testminage+3,i,l)
						A2f=smtestoutcomes(3,k-testminage+3,i,l)
						tsdiffvec(counter)= A1f-A1+A2f-A2
						counter=counter+1
					end if
				end if
			end do
		end do

!print*, 'LAPACK PRINT'
!print*, nspots, nregtsdiff+1,  blocksize
		! lapack
		call DGELS('N', nspots, nregtsdiff+1,1,tsdiffmat,nspots, tsdiffvec, nspots,worktsdiff,nregtsdiff+1+(nregtsdiff+1)*blocksize,infotsdiff)
		tsdiffestmat(i,:)=tsdiffvec(1:nregtsdiff+1)
		! deallocate
		deallocate(tsdiffvec)
		deallocate(tsdiffmat)
	end do 
	tsdiffest=sum(tsdiffestmat,1)/Npaths




	!! first count the simulated observations that will go into this
	
	!nspots=0
	!do k=testminage+1,testmaxage-2
		!do l=1, SampleSize
			!if ( idmat(l,lfpsize+2*expsize-1+tssize+k)==1) then ! if this mom is in calculation for period k...
				!do i=1,Npaths 	! ... go through simulated outcomes for that mom and figure which ones have test scores for both kids at that period.
					!! if kid 2 is older than 5, then create a spot
					!if (SS(5,k,i,l) > 4.8)  nspots=nspots+1
				!end do
			!end if
		!end do
	!end do

	!! now we can allocate data matrix and dependent variable vectors
	!allocate(tsdiffmat (nspots, nregtsdiff+1) )
	!allocate(tsdiffvec (nspots) )
	
	!tsdiffmat(:,nregtsdiff+1)=1.0d0
	!! another loop to fill in these guys
	!counter=1
	!do k=testminage+1,testmaxage-2
		!do l=1, SampleSize
			!if ( idmat(l,lfpsize+2*expsize-1+tssize+k)==1) then ! if this mom is in calculation for period k...
				!do i=1,Npaths 	! ... go through simulated outcomes for that mom and pick the right ones
					!! 1. fill in the labor supply related parts (smh)
					!!if ( smtestoutcomes(3, k-testminage+1,i,l) > -1.0d0 )   then ! if the second kid is old enough to have a testscore (first is kid is already at most testmaxage-2
					!if (SS(5,k,i,l)>4.8) then 
						!smh=sum(smchoices(:,k,i,l)*(/0.0d0,0.5d0,1.0d0/))
						!smhnext=sum(smchoices(:,k+1,i,l)*(/0.0d0,0.5d0,1.0d0/))
						!tsdiffmat(counter, 1:2)=(/smh,smhnext/)
						!tsdiffmat(counter, 3:5)=smh*omega3data(1:3,l)
						!tsdiffmat(counter, 6:8)=smhnext*omega3data(1:3,l)
						!tsdiffmat(counter,9:11)=omega3data(1:3,l)
						!tsdiffmat(counter,12:13)=SS(4:5,k,i,l)
						!! now fillin the tsdiffmat vec
						!A1=smtestoutcomes(1,k-testminage+1,i,l)
						!A2=smtestoutcomes(3,k-testminage+1,i,l)
						!A1f=smtestoutcomes(1,k-testminage+3,i,l)
						!A2f=smtestoutcomes(3,k-testminage+3,i,l)
						!tsdiffvec(counter)= A1f-A1+A2f-A2
						!counter=counter+1
					!end if
				!end do
			!end if
		!end do
	!end do
	!!print*, tsdiffmat(:,3:5)


	!! lapack
	!call DGELS('N', nspots, nregtsdiff+1,1,tsdiffmat,nspots, tsdiffvec, nspots,worktsdiff,nregtsdiff+1+(nregtsdiff+1)*blocksize,infotsdiff)
	!tsdiffest=tsdiffvec(1:nregtsdiff+1)
	!! deallocate
	!deallocate(tsdiffvec)
	!deallocate(tsdiffmat)
	! put all the vecs together

	
	momentvec=(/fulltimecoef, parttimecoef, expest, tsest,tsdiffest/)

end subroutine momentsalt







end module optu

