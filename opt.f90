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
		logw=logw0+par(1)*E+par(2)*time+par(3)*llms
		wagef=exp(logw+eps)
	end function wagefsim

	subroutine simhist(simdata,omega3,intercepts,parA,parU,parW,parH,beta,sigma,a1type,pa1type,parBmat,vcoef,wcoef,llmsvec,id,rho)
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
		real(dble), intent(in) :: vcoef(Gsizeoc+1,nperiods),wcoef(Gsize+1,nperiods-deltamin+2,deltamax-deltamin+1,1)
		real(dble), intent(in) :: rho 						!<skill depreciation para, added later.
		real(dble), intent(out) :: simdata(simvecsize,nperiods, Npaths)

		! LOCALS
		integer period, mainseed(2)
		real(dble) a1(Npaths)
		real(dble) Ainit,wage(Npaths),wageh(Npaths)
		real(dble) mu(shocksize1)
		real(dble) eps(Npaths,shocksize1)
		! State space elements to be collected, even if they are not data, I will collect omega1 and omega2, except llms
		real(dble) SS(6,nperiods,Npaths) 		
		! outcomes: wages of the father and mother.
		real(dble) outcomes(2,nperiods,Npaths)
		! choices : the choice history of h.
		integer choices(1,neriods,Npaths) 
		! not in the data, but I will keep an history of the x choices as well.
		integer xchoices(1,nperiods,Npaths)
			
		
		! choice specific utilities WITHIN the period
		real(dble) umat(3,Npaths), umatbig(3,11,Npaths),uc
		! interpolating vectors within a period and for an npath alone, constant parts predicted beforehand.
		real(dble) intcons(5), fvconsw, fvconsv
		real(dble) intw(9), intv(5) ,fvw,fvv
		real(dble) A1next, A2next, Enext, Anext(2),inputs(5), income
		real(dble) nextcollect(3,3) ! to collect collected future state space points and carry them to next period.
		real(dble) nextcollectbig(3,3,xgridsize) ! 
		real(dble) next(3,Npaths)  ! to collect collected future state space points and carry them to next period.
		real(dble) pbirth, omegaB(4), birthdraw, 
		integer birthhist(Npaths)
		integer birth, baby
 		integer i,j,k,l
 		integer coordbig(2)
		real(dble) incrat ! ratio of the income that goes to consumption. I use this in the periods 15,21 to change income for W.


		! Iniatilize some stuff
		mu=0.0d0
		xchoices=0 			! if not chosen, be 0.
		period=1
		birthhist=0 	! this records the timing of the second birth. if none, stays at zero.	
		! first draw Npaths alternative folks to assign a1 values.	
		mainseed=id*(/1,1/)
		call random_seed(put=mainseed)
		call random_number(a1)
		where (a1<pa1type(1))
			a1=a1type(1)
		else 
			a1=a1type(2)
		end where

		! initialize A
		Ainit=intercepts(1)+parA(1)*omega3(1)+parA(2)*omega3(2)+parA(3)*omega3(3) 
		! initialize wage	
		logw0=intercepts(2)+parvecW(1)*schooling+parvecW(2)*afqt+parvecW(3)*age0m+parvecW(4)*age0m*age0m
		! birth prob determiner
		omegaB=(/omega3(3),omega3(3)**2,omega3(1),omega3(2)/)
		
		
		
		! ----------------------------------------FIRST PERIOD---------------------------------------------
		eps= randmnv(Npaths,5,mu,sigma,3,1,(1,id*1))
		! initialize omegas
		SS(1,1,:)=Ainit 	!A1
		SS(2,1,:)=0.0d0 	!A2
		SS(3,1,:)=0.0d0 	!E
		SS(4,1,:)=1.0d0 	!age1
		SS(5,1,:)=0.0d0 	!age2
		SS(6,1,:)=omega3(3) !agem
		
		wage=wagef(0.d0,1.0d0,llmsvec(1),omega3(1),omega3(2),omega3(3),eps(:,4),intercepts(2),parW(5:7),parW(1:4))
		wageh=wagehfquick(omega3(4),1.0d0,eps(:,5),parH(5:6))	
		outcomes(1,1,:)=wage
		outcomes(2,1,:)=wageh
		! current utilities
		umat(1,:)=parU(2)*SS(1,1,:) 	+     parU(4)*(wageh)**par(5)+par6(6)
		umat(2,:)=parU(2)*SS(1,1,:)+a1*0.5d0+ parU(4)*(wageh+0.5d0*wage)**parU(5)+parU(6)+parU(7)*0.5d0
		umat(3,:)=parU(2)*SS(1,1,:)+a1*1.0d0+ parU(4)*(wageh+wage)**parU(5)+parU(6)+parU(7)
		! --------------------future value calculations----------------------
		! first, birth probabilities next period
		pbirth=bprobexo(omegaB,parBmat(:,1))
		! what I do this period will affect future A and E, the rest of the state space (which are used in the interpolations) are
		! deterministic.
		! first common deterministic parts
		intcons=(/omega3(3)+1.0d0,llmsvec(period),omega3(1),omega3(2)/)
		fvconsv=sum(vcoef(4:8,period+1)*intcons)+vcoef(Gsizeoc+1,period+1)
		fvconsw=sum(wcoef(4:8,period+1,period,1)*intcons)+wcoef(Gsize+1,period+1period,1)
	 	do i=1, Npaths
	 		do j=1,3
	 			A1next=pfone((/SS(1,period,i),(1-0.25d0*(j-1)),0.5d0*(wageh(i)+wage(i)*0.5d0*(j-1))/), omega3(1:3), period,parA(4:12),rho) 
	 			A2next=Ainit
	 			Enext=SS(3,period,i)+(j-1)*0.5d0
	 			intw=(/A1next,A2next,Enext,A1next**2,A2next**2,Enext**2, A1next*A2next,A1next*Enext+A2next*Enext/)
	 			intv=(/A1next,Enext,A1next**2,Enext**2, A1next*Enext/)
				fvw=sum(intw*(/wcoef(1:3,period+1,period,1),wcoef(9:Gsize,period+1,period,1)/))
				fvv=sum(intv*(/vcoef(1:2,period+1),wcoef(8:Gsizeoc,period+1)/))
				! now add the expected future value to corresponding umat element
				umat(j,i)=umat(j,i)+beta*(pbirth*(fvw+fvconsw)+(1-pbirth)*(fvv+fvconsv))
				nextcollect(:,j)=(/A1next,A2next,Enext/)
			end do
			choices(1,period,i)=maxloc(umat(:,i))
			next(:,i)=nextcollect(:,choices(1,period,i))
		end do
		
		! -----------------looooop to the future-----------
		
		do period=2,nperiods
			! draw the period shocks
			eps= randmnv(Npaths,5,mu,sigma,3,1,(period,id*period))
			! update state space vectors as much as you can outside the loop	
			! USE THE NEXT MATRIX TO UPDATE THE STATE SPACE FOR As and E, so that I don't need to calculate them twice.
			! I am updating A2 as well, but if it turns out this is a period where there is no second or birth, it needs to be
			! revert back to zero. Same as for age2
			SS(1,period,:)=next(1,:) 
			SS(2,period,:)=next(2,:)
			SS(3,period,:)=next(3,:)
			SS(4,period,:)=SS(4,period-1,:)+1.0d0 	! update age1
			SS(5,period,:)=SS(5,period-1,:)+1.0d0 	! update age2
			SS(6,period,:)=SS(6,period-1,:)+1.0d0 	! udpate agem
			!SS(3,period,:)=SS(3,period-1,:)+(choices(1,period-1,:)*1.0d0-1.0d0)/2.0d0
				
			! these are all at vector level, i.e., for all Npaths 
			wage=wagefsim(SS(3,period,:),period*1.0d0,llmsvec(period),logw0,eps(:,4))	
			wageh=wagehfquick(omega3(4),period*1.0d0,eps(:,5),parH(5:6))	
			outcomes(1,period,:)=wage
			outcomes(2,period,:)=wageh
			
			! determine the prob of a child birth. this part is constant all the Npaths
			pbirth=bprobexo(omegaB,parBmat(:,period-1))	
			
			! get some period specific but not path specific stuff to calculate interpolated fv
			intcons=(/omega3(3)+period*1.0d0,llmsvec(period),omega3(1),omega3(2)/)
			fvconsv=sum(vcoef(4:8,period+1)*intcons)+vcoef(Gsizeoc+1,period+1)
			fvconsw=sum(wcoef(4:8,period+1,period,1)*intcons)+wcoef(Gsize+1,period+1period,1)
			
			! -------------------------- 1a. PERIOD<8: STILL IN THE FECUND PERIOD ---------------------------------
			if (period<8)  then

				! FIRST FIGURE OUT WHO HAD BIRTH
				do i=1,Npaths
					if (SS(5,period-1,i)<0.001) then   ! they don't have a second child
						call random_seed(put=mainseed*period*i) 
						call random_number(birthdraw)
						if (birthdraw<pbirth) then
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

				do i=Npaths
					! ----------------------------1.a PERIOD<8: ONE CHILD-----------------------------------
					if (SS(5,period,i)<0.001) then 		! if still with one child. 
						if (SS(4,period,i)<3) then
							baby=1
						else 
							baby=0
						end if
						uc=parU(2)*SS(1,period,i) 
						umat(1,i)=uc+parU(4)*(wageh(i))**par(5)+par6(6)*baby
						umat(2,i)=uc+a1(i)*0.5d0+ parU(4)*(wageh(i)+0.5d0*wage(i))**parU(5)+parU(6)*baby+parU(7)*0.5d0*baby
						umat(3,i)=uc+a1(i)*1.0d0+ parU(4)*(wageh(i)+wage(i))**parU(5)+parU(6)*baby+parU(7)*baby
						
						! now calculate future values for the possible 3 decisions, so we can calculate the future value.
						do j=1,3
							A1next=pfone((/SS(1,period,i),(1-0.25d0*(j-1)),0.5d0*(wageh(i)+wage(i)*0.5d0*(j-1))/), omega3(1:3), period,parA(4:12),rho) 
							A2next=Ainit
							Enext=SS(3,period,i)+(j-1)*0.5d0
							intw=(/A1next,A2next,Enext,A1next**2,A2next**2,Enext**2, A1next*A2next,A1next*Enext+A2next*Enext/)
							intv=(/A1next,Enext,A1next**2,Enext**2, A1next*Enext/)
							fvw=sum(intw*(/wcoef(1:3,period+1,period,1),wcoef(9:Gsize,period+1,period,1)/))
							fvv=sum(intv*(/vcoef(1:2,period+1),wcoef(8:Gsizeoc,period+1)/))
							! now add the expected future value to corresponding umat element
							umat(j,i)=umat(j,i)+beta*(pbirth*(fvw+fvconsw)+(1-pbirth)*(fvv+fvconsv))
							nextcollect(:,j,i)=(/A1next,A2next,Enext/)
						end do
						! finally, pick the highest expected utility one.
						choices(1,period,i)=maxloc(umat(:,i))
						! get the optimized future state space points to use them in future.
						next(:,i)=nextcollect(:,choices(1,period,i),i)
					
 					!-------------------------- PERIOD<8: 1.b TWO CHILDREN ------------------------------
					else
						if (SS(4,period,i)<3 .OR. SS(5,period,i)<3) then
							baby=1
						else 
							baby=0
						end if
						! now we x is a choice, so I will have to use the bigger umatbig
						! period utilities don't depend on x, so first get rid of those
						! because there are two children, I need to use the CES utility function
						uc=uctwo(SS(1,period,i),SS(2,period,i),parU(1)) 
						umatbig(1,:,i)=	uc + 			 parU(4)*(wageh(i))**par(5)+par6(6)*baby
						umatbig(2,:,i)= uc +a1(i)*0.5d0+ parU(4)*(wageh(i)+0.5d0*wage(i))**parU(5)+parU(6)*baby+parU(7)*0.5d0*baby
						umatbig(3,:,i)= uc +a1(i)*1.0d0+ parU(4)*(wageh(i)+wage(i))**parU(5)+parU(6)*baby+parU(7)*baby
						! now calculate and add future values.
						do j=1,3
							! E does not depend on x, so calculate it beforehand, also income
							Enext=SS(3,period,i)+(j-1)*0.5d0
							income=wageh(i)+wage(i)*(j-1)*0.25 ! it is 0.25 because half the income goes to goods inputs
							do k=1,xgridsize
								inputs=(/SS(1,period,i),SS(2,period,i),(1-0.25*(j-1))*xgrid(k),(1-0.25*(j-1))*(1-xgrid(k)),income/)
								Anext=pftwo(inputs,omega3(1:3),parA(4:12), SS(4:5,period,i),rho)
								A1next=Anext(1)
								A2next=Anext(2)
								intw=(/A1next,A2next,Enext,A1next**2,A2next**2,Enext**2, A1next*A2next,A1next*Enext+A2next*Enext/)
								fvw=sum(intw*(/wcoef(1:3,period+1,period,1),wcoef(9:Gsize,period+1,period,1)/))
								umatbig(j,k,i)=umatbig(j,k,i)+beta*(fvw+fvconsw)
								nextcollectbig(:,j,k)=(/A1next,A2next,Enext/)
							end do
						end do
						coordbig=maxloc(umatbig(:,:,i))
						choices(1,period,i)=coordbig(1)
						xchoices(1,period,i)=coordbig(2)
						next(:,i)=nextcollectbig(:,coordbig(1),coordbig(2))
				
					end if  ! end of period<8 WV if


				end do
								
				! --%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%--
				!----------------------------2. 8<period<15 NO CHANCE OF BIRTH, FIRST KID STILL IN CARE---------------------------------
			elseif ((period>=8) .AND. (period<astar))
				do i=Npaths
					! calculations are different for V and W.
					
					! ---------------------------2.a PERIOD (8,15): ONE CHILD-----------------------------------
					if (SS(5,period,i)<0.001) then 		! if still with one child. 
						if (SS(4,period,i)<3) then
							baby=1
						else 
							baby=0
						end if
						uc=parU(2)*SS(1,period,i) 
						umat(1,i)=uc+parU(4)*(wageh(i))**par(5)+par6(6)*baby
						umat(2,i)=uc+a1(i)*0.5d0+ parU(4)*(wageh(i)+0.5d0*wage(i))**parU(5)+parU(6)*baby+parU(7)*0.5d0*baby
						umat(3,i)=uc+a1(i)*1.0d0+ parU(4)*(wageh(i)+wage(i))**parU(5)+parU(6)*baby+parU(7)*baby
						
						! now calculate future values for the possible 3 decisions, so we can calculate the future value.
						do j=1,3
							A1next=pfone((/SS(1,period,i),(1-0.25d0*(j-1)),0.5d0*(wageh(i)+wage(i)*0.5d0*(j-1))/), omega3(1:3), period,parA(4:12),rho) 
							A2next=Ainit
							Enext=SS(3,period,i)+(j-1)*0.5d0
							intv=(/A1next,Enext,A1next**2,Enext**2, A1next*Enext/)
							fvv=sum(intv*(/vcoef(1:2,period+1),wcoef(8:Gsizeoc,period+1)/))
							! now add the expected future value to corresponding umat element
							umat(j,i)=umat(j,i)+beta*(fvv+fvconsv)
							nextcollect(:,j,i)=(/A1next,A2next,Enext/)
						end do
						! finally, pick the highest expected utility one.
						choices(1,period,i)=maxloc(umat(:,i))
						! get the optimized future state space points to use them in future.
						next(:,i)=nextcollect(:,choices(1,period,i),i)
					
					!-------------------------- 2.b PERIOD(8,15): TWO CHILDREN ------------------------------
					else
						if (SS(4,period,i)<3 .OR. SS(5,period,i)<3) then
							baby=1
						else 
							baby=0
						end if
						! now we x is a choice, so I will have to use the bigger umatbig
						! period utilities don't depend on x, so first get rid of those
						! because there are two children, I need to use the CES utility function
						uc=uctwo(SS(1,period,i),SS(2,period,i),parU(1)) 
						umatbig(1,:,i)=	uc + 			 parU(4)*(wageh(i))**par(5)+par6(6)*baby
						umatbig(2,:,i)= uc +a1(i)*0.5d0+ parU(4)*(wageh(i)+0.5d0*wage(i))**parU(5)+parU(6)*baby+parU(7)*0.5d0*baby
						umatbig(3,:,i)= uc +a1(i)*1.0d0+ parU(4)*(wageh(i)+wage(i))**parU(5)+parU(6)*baby+parU(7)*baby
						! now calculate and add future values.
						do j=1,3
							! E does not depend on x, so calculate it beforehand, also income
							Enext=SS(3,period,i)+(j-1)*0.5d0
							income=wageh(i)+wage(i)*(j-1)*0.25 ! it is 0.25 because half the income goes to goods inputs
							do k=1,xgridsize
								inputs=(/SS(1,period,i),SS(2,period,i),(1-0.25*(j-1))*xgrid(k),(1-0.25*(j-1))*(1-xgrid(k)),income/)
								Anext=pftwo(inputs,omega3(1:3),parA(4:12), SS(4:5,period,i),rho)
								A1next=Anext(1)
								A2next=Anext(2)
								intw=(/A1next,A2next,Enext,A1next**2,A2next**2,Enext**2, A1next*A2next,A1next*Enext+A2next*Enext/)
								fvw=sum(intw*(/wcoef(1:3,period+1,period,1),wcoef(9:Gsize,period+1,period,1)/))
								umatbig(j,k,i)=umatbig(j,k,i)+beta*(fvw+fvconsw)
								nextcollectbig(:,j,k)=(/A1next,A2next,Enext/)
							end do
						end do
						coordbig=maxloc(umatbig(:,:,i))
						choices(1,period,i)=coordbig(1)
						xchoices(1,period,i)=coordbig(2)
						next(:,i)=nextcollectbig(:,coordbig(1),coordbig(2))
					end if  ! end of period<8 WV if
				end do
				!#%#%#%#%#%#%#%#%#%#%#%#%%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%%#%#%#%#%%#%#


				! ---------------------3. Period (16,21) No x choice, because Kid 1 has grown
			elseif (period>=astar)
				do i=Npaths
					! calculations are different for V and W.
				
					! ---------------------------3.a PERIOD (16,21): ONE CHILD-----------------------------------
					! F
					if (SS(5,period,i)<0.001) then 		! if still with one child. 
						! there is no chance of having a baby at this point
						baby=0
						uc=parU(2)*SS(1,period,i) 
						! NOTE: important difference: all the income goes into the consumption now, because there is no production of
						! skils
						umat(1,i)=uc+parU(4)*(wageh(i))**par(5)+par6(6)*baby
						umat(2,i)=uc+a1(i)*0.5d0+ parU(4)*(wageh(i)+0.5d0*wage(i))**parU(5)+parU(6)*baby+parU(7)*0.5d0*baby
						umat(3,i)=uc+a1(i)*1.0d0+ parU(4)*(wageh(i)+wage(i))**parU(5)+parU(6)*baby+parU(7)*baby
						
						! now calculate future values for the possible 3 decisions, so we can calculate the future value.
						do j=1,3
							A1next=SS(1,period,i)
							A2next=SS(2,period,i)
							Enext=SS(3,period,i)+(j-1)*0.5d0
							intv=(/A1next,Enext,A1next**2,Enext**2, A1next*Enext/)
							fvv=sum(intv*(/vcoef(1:2,period+1),wcoef(8:Gsizeoc,period+1)/))
							! now add the expected future value to corresponding umat element
							umat(j,i)=umat(j,i)+beta*(fvv+fvconsv)
							nextcollect(:,j,i)=(/A1next,A2next,Enext/)
						end do
						! finally, pick the highest expected utility one.
						choices(1,period,i)=maxloc(umat(:,i))
						! get the optimized future state space points to use them in future.
						next(:,i)=nextcollect(:,choices(1,period,i),i)
					
					!-------------------------- 3.b PERIOD(15,21): TWO CHILDREN ------------------------------
					else
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

						umat(1,i)=	uc + 	 parU(4)*(wageh(i)*incrat)**par(5)+par6(6)*baby
						umat(2,i)= uc +a1(i)*0.5d0+ parU(4)*(incrat*(wageh(i)+0.5d0*wage(i)))**parU(5)+parU(6)*baby+parU(7)*0.5d0*baby
						umat(3,i)= uc +a1(i)*1.0d0+ parU(4)*(incrat*(wageh(i)+wage(i)))**parU(5)+parU(6)*baby+parU(7)*baby
						! now calculate and add future values.
						do j=1,3
							Enext=SS(3,period,i)+(j-1)*0.5d0
							income=wageh(i)+wage(i)*(j-1)*0.25 ! it is 0.25 because half the income goes to goods inputs
							inputs=(/SS(1,period,i),SS(2,period,i),0.d0,(1-0.25*(j-1))*(1-xgrid(k)),income/)
							Anext=pftwo(inputs,omega3(1:3),parA(4:12), SS(4:5,period,i),rho)
							A1next=SS(1,period,i) 	! don't change A1
							A2next=Anext(2)
							if (period>=astar+birthhist(i)-1) A2next=SS(2,period,i)
							intw=(/A1next,A2next,Enext,A1next**2,A2next**2,Enext**2, A1next*A2next,A1next*Enext+A2next*Enext/)
							fvw=sum(intw*(/wcoef(1:3,period+1,period,1),wcoef(9:Gsize,period+1,period,1)/))
							umatbig(j,k,i)=umatbig(j,k,i)+beta*(fvw+fvconsw)
							nextcollectbig(:,j,k)=(/A1next,A2next,Enext/)
						end do
						coordbig=maxloc(umatbig(:,:,i))
						choices(1,period,i)=coordbig(1)
						xchoices(1,period,i)=coordbig(2)
						next(:,i)=nextcollectbig(:,coordbig(1),coordbig(2))
					end if  ! end of period<8 WV if
				end do
				!#%#%#%#%#%#%#%#%#%#%#%#%%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%%#%#%#%#%%#%#



			
				! ##################################################################3###################3

				!-------------------------4. FINAL PERIOD-----------------------------
			else
				! This is the final period. Future values are not calculated but terminal values are used.
				! I have been using these temporary terminal values for the model, it seems like I will stick to those 
				! for the jmpsimple. they will have three choices as before. 

				


			end if


		end do


		
		
		end subroutine simhist











	!----------------------OPTIMIZATION RELATED STUFF--------------------------

   ! !> initialize parameters
	!subroutine initparam(parA,parU,parW, parH, beta,sigma1, parB, ctype, mtype, atype, a1type, condprob) 
	!implicit none
		!real(dble), intent(out) ::  parA(12),parU(7),parW(7),parH(6),beta,sigma1(shocksize1,shocksize2),parB(Bsizeexo+1)
		!real(dble), intent(out):: ctype, mtype, atype, a1type(na1type)
		!real(dble), intent(out):: condprob ! prob of second child of being a certain type, conditional on mother type 
		!! locals
		!real(dble) var0,varp,varf,varw, varwh, var00, var01, varp0,varp1,varf0,varf1
		!real(dble) cov0p, cov0f, covpf, covwwh, cov0w,cov0wh,covpw, covpwh,covfw,covfwh
		!real(dble) cov0001,cov00p0,cov00p1,cov00f0,cov00f1,cov01p0,cov01p1,cov01f0,cov01f1,covp0p1,covp0f0,covp0f1,covp1f0,covp1f1,covf0f1			
		!! -----------FUNDAMENTAL PARAMETERS---------------
		!! gamma0's and gamma's: 
		!parA=(/1.0d0,0.1d0,0.1d0,0.1d0,0.1d0, 0.1d0,0.1d0,0.1d0,0.1d0,0.1d0,0.1d0,0.1d0/)
		!! gamma, alpha, alpha1,...alpha5	
		!parU=(/0.1d0,0.1d0,0.1d0,0.1d0,0.1d0, 0.1d0,0.1d0/)
		!! beta0 (4) and beta (3) in wage equation
		!parW=(/0.1d0,0.1d0,0.1d0,0.1d0,0.1d0, -0.1d0,0.1d0/)*0.00005
		!! beta_ws (4)
		!parH=(/0.1d0,0.1d0,0.1d0,0.1d0,0.1d0,0.1d0/)*0.04

		!beta=0.95
		!condprob=1
		!parB=0.10d0	
		!! ------------SHOCK STUFF-----------------------
		!! ----   VARIANCES
		
		!! hc periods
		!var0=0.1d0 		! variance h=0
		!varp=0.1d0 		! variance h=0.5
		!varf=0.1d0 		! variance h=1
		!varw=0.1d0 		! variance of wage shock
		!varwh=0.1d0 	! variance of father's wage shock
		
		!! hcb periods
		!!var00=0.1d0 	! variance of h=0, b=0
		!!var01=0.1d0 	! variance of h=0, b=1
		!!varp0=0.1d0 	! variance of h=0.5, b=0
		!!varp1=0.1d0 	! variance of h=0.5, b=1
		!!varf0=0.1d0 	! variance of h=1, b=0
		!!varf1=0.1d0 	! variance of h=1, b=1
		
		!! ---  COVARIANCES
		!! hc periods
		!cov0p=0.1d0 	! cov(h=0,h=0.5)
		!cov0f=0.1d0 	! cov(h=0,h=1)
		!covpf=0.1d0 	! cov(h=0.5,h=1)
		!covwwh=0.1d0 	! cov(wages of mom and dad)
		!! cross w and preference shocks- JUST ZERO FOR NOW
		!cov0w=0.0d0 	! cov(wage, h=0)
		!cov0wh=0.0d0 	! cov(wageh, h=0)
		!covpw=0.0d0 	! cov(wage, h=0.5)
		!covpwh=0.0d0 	! cov(wageh, h=0.5)
		!covfw=0.0d0 	! cov(wage, h=1)
		!covfwh=0.0d0 	! cov(wageh, h=1)
		
		!! hcb periods
	
		!!cov0001=0.1d0 	! cov(h=0,b=0 and h=0, b=1)
		!!cov00p0=0.1d0 	! cov(h=0,b=0 and h=0.5, b=0)	
		!!cov00p1=0.1d0 	! cov(h=0,b=0 and h=0.5, b=1)	
		!!cov00f0=0.1d0 	! cov(h=0,b=0 and h=1, b=0)	
		!!cov00f1=0.1d0 	! cov(h=0,b=0 and h=1, b=1)	
		
		!!cov01p0=0.1d0 	! cov(h=0,b=1 and h=0.5, b=0)	
		!!cov01p1=0.1d0 	! cov(h=0,b=1 and h=0.5, b=1)	
		!!cov01f0=0.1d0 	! cov(h=0,b=1 and h=1, b=0)	
		!!cov01f1=0.1d0 	! cov(h=0,b=1 and h=1, b=1)	
	
		!!covp0p1=0.1d0 	! cov(h=0.5,b=0 and h=0.5,b=1)
		!!covp0f0=0.1d0
		!!covp0f1=0.1d0
		
		!!covp1f0=0.1d0
		!!covp1f1=0.1d0

		!!covf0f1=0.1d0

		
		!! cross w and preference shocks
		!! fuck, let's just say they are zero.	
		!! TODO: NEED TO PROPERLY HANDLE PARAMETERS. YOU WILL NEED TO UPDATE THESE AT SOME POINT
		!! fill in sigma => might fill this in properly later
		!sigma1=0.1d0
		
		
		!ctype=1.0d0
		!mtype=1.0d0
		!atype=1.0d0
		!a1type=(/0.0d0,1.0d0/)
		!!ctype=(/0.0d0,1.0d0/)
		!!mtype=(/0.d0,1.0d0/)
		!!atype=(/1.0d0,1.5d0/)
		!!a1type=(/1.0d0,2.0d0/)
	
		!! conditional probabilities of child 2 type should be pr(x|y)=pr(ctype=x,mtype=y)/pr(mtype=y)
		!! so I need to specify the joint distribution of ctype and mtype. if this does not go anywhere else I can just use the
		!! conditional probs as a parameter 
			
	!end subroutine initparam	




end module opt

