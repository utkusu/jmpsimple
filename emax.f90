module emax
!  module holds the emax calculating routines. 
! TODO: remove the omega4 for the interpolation stuff, so that tr and ftr routines don't get that part of those mss and tss
! matrices, they don't use them anyway. This will require the change of the dimensions in tr routines and removing omega4 from
! filler and mutype from mss matrices.
!TODO: Two child coefficients are wild. One is not, why?

	use global 
	use randomgen
	USE IFPORT
	implicit none

contains 

	!------------------------------------------------AUX FUNCTIONS-----------------------------------------------------------
	!> aux function to create N x n matrix from size N vector x, by stacking x's in columns as in [x x x ... x]
	function sidestack(x,n)
		implicit none
		real(dble) x(:)
		integer n
		real(dble) sidestack(size(x),n)
		integer i
		do i=1,n
			sidestack(:,i)=x
		end do
	end function sidestack


	! ----------------------------- MODEL FUNCTIONS -------------------------------------------------------------------------	
	
	!> cognitive skill production function for one child families
function pfone(inputs, constants,age, parameters,rho)
	implicit none
	real(dble)  pfone 			!<output: A1 for one child family
	real(dble)  inputs(:) 		!<inputs: A_{t-1}, i, g, size=3
	real(dble)  constants(:) 	!<constants: things that stay the same from t=1: schooling, AFQT, age0m 
	real(dble)  parameters(:) 	!<parameters: gammas(9 of them, first early childhood, then late, then goods inputs
	real(dble)  age 			!<age of the child
	real(dble) rho  			!< additional cumulative input coefficient, late addition so separate
	
	if (age<abar) then
		pfone=rho*inputs(1)+(parameters(1)+parameters(2)*constants(1)+parameters(3)*constants(2)+parameters(4)*constants(3))*inputs(2)+parameters(9)*inputs(3)
	else
		pfone=rho*inputs(1)+(parameters(5)+parameters(6)*constants(1)+parameters(7)*constants(1)+parameters(8)*constants(3))*inputs(2)+parameters(9)*inputs(3)
	end if

	if (age>=astar) pfone=inputs(1) 		! if older than a star, just stay put. 
end function pfone
!> cognitive skill production for two child families
function pftwo(inputs, constants,age,parameters,rho)
	implicit none
	real(dble) pftwo(2) 		!<output: A1 and A2 (2)
	real(dble) inputs(:) 		!<inputs: lagged A1, lagged A2, i1, i2, g (5)
	real(dble) constants(:) 	!<things that stay the same:schooling, AFQT, age0m  size=3
	real(dble) parameters(:) 	!<parameters gammas, 9 of them, first early childhood, then late, then goods inputs
	real(dble)  age(:) 			!<ages of the children, age1 age2 (2)
	real(dble) rho  			!< additional cumulative input coefficient, late addition so separate

	! first kid first

	if (age(1)<abar) then
		pftwo(1)=rho*inputs(1)+(parameters(1)+parameters(2)*constants(1)+parameters(3)*constants(2)+parameters(4)*constants(3))*inputs(3)+parameters(9)*inputs(5)
	else
		pftwo(1)= rho*inputs(1)+(parameters(5)+parameters(6)*constants(1)+parameters(7)*constants(1)+parameters(8)*constants(3))*inputs(3)+parameters(9)*inputs(5) 
	end if
	! second kid
	if (age(2)<abar) then
		pftwo(2)=rho*inputs(2)+(parameters(1)+parameters(2)*constants(1)+parameters(3)*constants(2)+parameters(4)*constants(3))*inputs(4)+parameters(9)*inputs(5)
	else
		pftwo(2)= rho*inputs(2)+(parameters(5)+parameters(6)*constants(1)+parameters(7)*constants(1)+parameters(8)*constants(3))*inputs(4)+parameters(9)*inputs(5) 
	end if

	! addition, keep first child's achievement the same if he is older than a*
	if (age(1)>=astar) pftwo(1)=inputs(1)
	if (age(2)>=astar) pftwo(2)=inputs(2) ! not that I will need this
end function pftwo


!> utility function part for children achievement: DON'T USE THIS, just write it
! out there.
function ucone(A1, alpha)
	implicit none
	real(dble) ucone 	!< output
	real(dble) A1 		!< achievement level of the first child
	real(dble) alpha 	!< the parameter
	ucone=(alpha*A1)/emaxscale
end function



!> utility function part for children achievement
function uctwo(A1, A2, phi)
	implicit none
	real(dble) uctwo 	!< output
	real(dble) A1,A2	!< achievement level of the first and second child
	real(dble) phi 		!< the parameter
	uctwo=((0.5*A1**phi) + (0.5*A2**phi) )**(1/phi)
	if ( phi<0.000001 .AND. phi>-0.0000001) uctwo=(A1**(0.5d0))*(A2**(0.5d0)) 	! practically converge to cobbdouglas
	
	uctwo=uctwo/emaxscale
end function uctwo

!> utility function h part
function uh(h,cons,ages,par)
	implicit none
	real(dble) h, cons(:) 			!<labor supply choice, consumption
	real(dble)  ages(:) 		!< vector of children ages, size=1 or 2
	real(dble) 	par(:) 		 	!< vector of parameters, as specified in the uh equation in extended.pdf pg. 5 
	real(dble) uh(size(cons)) 		!< output
	
	integer youngchild 

	youngchild=0
	if (minval(ages)<3)  youngchild=1
	uh=par(1)*h+par(2)*cons**par(3)+par(4)*youngchild+par(5)*h*youngchild
	uh=uh/emaxscale
end function uh

!> utility function h part for ONE CHILD FAMILIES
function uhoc(h,cons,ages,par)
	implicit none
	real(dble) h, cons(:) 			!<labor supply choice, consumption
	real(dble)  ages 		!< vector of children ages, size=1 or 2
	real(dble) 	par(:) 		 	!< vector of parameters, as specified in the uh equation in extended.pdf pg. 5 
	real(dble) uhoc(size(cons)) 		!< output
	
	integer youngchild 

	youngchild=0
	if (ages<3)  youngchild=1
	uhoc=par(1)*h+par(2)*cons**par(3)+par(4)*youngchild+par(5)*h*youngchild
	uhoc=uhoc/emaxscale
end function uhoc


!> wage of the mother
!> produces an array of the size of the provided shock vector, for emax calculations 
! I need to get rid of the lagh and introduc t instead, the parameter name can stay the same, I will just insert age1 as a
! variable instead.
function wagef(E,time,llms,schooling,afqt,age0m,eps,typem,par,par0)
	implicit none
	real(dble) E,time,llms,schooling,afqt,age0m 	!< experience,time, llms,schooling, afqt, age0m, as they appear in extended
	real(dble) eps(:) 							!< emax mc integration random draw, size=Nmc
	real(dble) par(:),par0(:) 					!< parameter vectors, for wage and wage0 equations as in document.
	real(dble) typem 							!< mother's type KEEP THIS FOR NOW + IT IS THE INTERCEPT PARAMETERS
	real(dble) wagef(size(eps)) 				!< output, wage income for the mother, size=Nmc

	real(dble) logw0, logw
	
	logw0=typem+par0(1)*schooling+par0(2)*afqt+par0(3)*age0m+par0(4)*age0m*age0m
	logw=logw0+par(1)*E+par(2)*time+par(3)*llms
	wagef=exp(logw+eps)
end function wagef
!> husband wage income
!> again, produces an array for emax mc

! Change add time to switch to the new spec
function wagehf(schooling, ageh,t, eps,par)
	implicit none
	real(dble) schooling, ageh,t 			!< inputs, ageh=age at time=0, t=time
	real(dble) eps(:) 					!< emax mc random draw size=Nmc
	real(dble) par(:) 					!< parameter vector, size=6
	real(dble) wagehf(size(eps)) 		!< output, wage income for the mother, size=Nmc

	real(dble) logwh

	logwh=par(1)+par(2)*schooling+par(3)*ageh+par(4)*ageh*ageh+t*par(5)+t*t*par(6)
	wagehf=exp(logwh+eps)
end function wagehf

! calculate wages given combined father state space point
function wagehfquick(omegaf,t,eps,parquick)
	implicit none
	real(dble) parquick(:)  	! parameters, of size=2
	real(dble) omegaf,t  		! inputs: omegaf is the combined initial value of the father
	real(dble) eps(:) 					!< emax mc random draw size=Nmc
	real(dble) wagehfquick(size(eps))
	real(dble) logwh
	logwh=omegaf+parquick(1)*t+parquick(2)*t*t
	wagehfquick=exp(logwh+eps)
end function wagehfquick

!-------------------------------------------Future and Terminal Values--------------------------------------------------
subroutine termval()
	implicit none
	
end subroutine termval


! ---------------------------------TEMP terminal values---------------------
! TODO Replace the temp terminal value with the final decision
!> This is the simplest terminal value calculation I can come up with.
!> Assume the last period consumption and achievement will be there
!> for the rest of the periods until 65 and sum up.
!> this assumes twins away: at T*, the older kid is already out of mom's influence, so no x choice
subroutine termvaltemp(TV,uct,uh0,uhp,uhf,beta,agem)
implicit none
	real(dble), intent(out)::TV(Nmc,3)	!< terminal value associated with each h for each draw
	real(dble), intent(in):: uct		!< utility from kids at TP
	real(dble), intent(in):: uhf(:),uh0(:),uhp(:)	!<utility levels from each h for Nmc draws at TP.
	real(dble), intent(in):: beta					!< discount factor
	real(dble), intent(in):: agem 					!< age of the mother at TP.

	real(dble) mult, power
	power=(finalage/timeperiod)-(agem)
	mult=beta*(1.0d0-beta**power)/(1.0d0-beta)
	TV(:,1)=(uh0+uct)*mult
	TV(:,2)=(uhp+uct)*mult
	TV(:,3)=(uhf+uct)*mult
end subroutine termvaltemp

!> fv calculation for provided interpolating coefficients
!> This is for later periods where mother only chooses h
!> to calculate fv, we need to first update the state space for
!> all the choices and monte carlo draw and then use the interpolating
!> coefficients (whatever they are) and do the interpolation.
subroutine fvlate(fv,omega1, omega2,omega3, omega4,fvpar)
	implicit none
	real(dble), intent(out):: fv(3) 			!< output, fv for each h 
	real(dble), intent(in):: omega1(:) 			!< omega1: endogenous dynamic
	real(dble), intent(in):: omega2(:) 			!< omega2: exogenous dynamic
	real(dble), intent(in):: omega3(:) 			!< omega3: static
	real(dble), intent(in):: omega4(:) 			!< omega4: type
	real(dble), intent(in):: fvpar(:) 			!< interpolating parameters
	
	!locals
	real(dble) fss((size(omega1)+size(omega2)+size(omega3)+size(omega4)),3) 	!future state space for each h
	real(dble) tfss(size(fvpar),3) 												!transformed future state space 
	real(dble) omegaall(size(omega1)+size(omega2)+size(omega3)+size(omega4)) 	!combined state space

	omegaall=(/omega1,omega2,omega3,omega4/)
	
	
	! create alternative specific future state spaces
	fss(:,1)=omegaall
	fss(:,2)=omegaall
	fss(:,3)=omegaall
	! update ages
	fss(4:6,:)=fss(4:6,:)+1
	! TODO: NOT UPDATING LLMS MEANS myobic expectations, keep this for now.
	! update experience level
	fss(3,2)=fss(3,2)+0.5d0
	fss(3,3)=fss(3,2)+1.0d0
	
	! NO NEED TO UPDATE achievement because both kids are older than a*
	
	call ftr_late(tfss,fss,3)  
	fv=matmul(fvpar,tfss)

end subroutine fvlate

!> fv calculation for provided interpolating coefficients.
!> This is for periods where mother only chooses h and c only.
!> First child is already grown, only second child receives inputs.
!> To calculate fv, we need to first update the state space for
!> all the choices (and monte carlo draw if necessary) and then use the interpolating
!> coefficients (whatever they are) and do the interpolation.

subroutine fvhc(fv,omega1, omega2,omega3, omega4,wage,wageh,fvpar,parA,rho)
	implicit none
	real(dble), intent(out):: fv(Nmc,cgridsize*3)	!< output, fv for each h,c and mc ! 
	real(dble), intent(in):: omega1(:) 				!< omega1: endogenous dynamic
	real(dble), intent(in):: omega2(:) 				!< omega2: exogenous dynamic
	real(dble), intent(in):: omega3(:) 				!< omega3: static
	real(dble), intent(in):: omega4(:) 				!< omega4: type
	real(dble), intent(in):: fvpar(:) 				!< interpolating parameters	
	real(dble), intent(in):: parA(:) 				!< parameters of the production function
	real(dble), intent(in):: wage(:),wageh(:) 		!< wages of parents, Nmcx1 each
	real(dble), intent(in) :: rho 					!< skill depriciation parameter, previously fixed at 1
	!locals
	
	real(dble) fss(Ssize,3*cgridsize,Nmc) 	!future state space for each h,c,MCdraw
	real(dble) filler((size(omega1)+size(omega2)+size(omega3)+size(omega4)),3*cgridsize) 	! filler matrix
	real(dble) tfss(Gsize+1,3*cgridsize,Nmc) 										!transformed future state space for each h,c,MCdraw
	real(dble) omegaall(size(omega1)+size(omega2)+size(omega3)+size(omega4)) 	!combined state space
	integer i,j,k
	real(dble) :: inputs(5),constants(3), ages(2), parameters(9),As(2)
	real(dble) logw0 					! initial logwage of the mother
	
	omegaall=(/omega1,omega2,omega3,omega4/)
		
	! fill up the filler matrix with the current space
	do i=1,3*cgridsize
		filler(:,i)=omegaall
	end do
 	! update the filler matrix with exogenous changes to it.
	filler(4:6,:)=filler(4:6,:)+1   ! kids and mom age

	! update experience level

	filler(3,cgridsize+1:2*cgridsize)=(omega1(3)+0.5d0)*onescgrid
	filler(3,2*cgridsize+1:3:cgridsize)=(omega1(3)+1.0d0)*onescgrid
!	filler(3,6:10)=(omega1(3)+0.5d0)*onescgrid
!	filler(3,11:15)=(omega1(3)+1.0d0)*onescgrid
		
	inputs=(/omega1(1),omega1(2),0.d0,0.d0,0.d0/)
	constants=omega3(1:3)
	ages=(/omega2(1),omega2(2)/)   	! with no child types, fill in the omega4 with the common intercept.
	parameters=parA(4:12)
	
	do k=1,Nmc
		!use the filler matrix to populate fss 3d matrix
		fss(:,:,k)=filler
		do j=1,cgridsize
			do i=1,3
				! inputs to second kid: first i=h, then g=(1-c)Y 
				inputs(4)=1.0d0-0.5d0*((i-1.0d0)/2.0d0)
				inputs(5)=(wage(k)*(i-1.0d0)/2+wageh(k))*(1-cgrid(j))
				As=pftwo(inputs, constants, ages, parameters,rho)
				fss(2,(i-1)*cgridsize+j,k)=As(2)
			end do
		end do
		! here we can project the state space matrix (for each h and c) 
		! for a given MC draw to the space of interpolating coefficients.
		! But for now I will just use the state space vectors as they are to calculate 
		! FVs. Replace fss with tfss (transformed future state space when you decide on the 
		! interpolating matrix

	end do
		call ftr_hc(tfss,fss,3*cgridsize,Nmc) 
		do k=1,Nmc
			fv(k,:)=matmul(fvpar,tfss(:,:,k)) 
		end do
end subroutine fvhc

!> fv calculation for provided interpolating coefficients.
!> This is for periods where mother chooses h,c,x.
!> To calculate fv, we need to first update the state space for
!> all the choices (and monte carlo draw if necessary) and then use the interpolating
!> coefficients (whatever they are) and do the interpolation.


! WARNING  THIS IS WRONG, writing an alt version for this
!subroutine fvhcx(fv,omega1, omega2,omega3, omega4,wage,wageh,fvpar,parA)
	!implicit none
	!real(dble), intent(out)::fv(Nmc,cxgridsize*3)	!< output, fv for each h,c and mc
	!real(dble), intent(in):: omega1(:) 				!< omega1: endogenous dynamic
	!real(dble), intent(in):: omega2(:) 				!< omega2: exogenous dynamic
	!real(dble), intent(in):: omega3(:) 				!< omega3: static
	!real(dble), intent(in):: omega4(:) 				!< omega4: type
	!real(dble), intent(in):: fvpar(:) 				!< interpolating parameters	
	!real(dble), intent(in):: parA(:) 				!< parameters of the production function
	!real(dble), intent(in):: wage(:),wageh(:) 		!< wages of parents, Nmcx1 each

	!!locals
	
	!real(dble) fss((size(omega1)+size(omega2)+size(omega3)+size(omega4)),3*cxgridsize,Nmc) 	!future state space for each h,c,MCdraw
	!real(dble) filler((size(omega1)+size(omega2)+size(omega3)+size(omega4)),3*cxgridsize) 	! filler matrix
	!real(dble) tfss(Gsize+1,3*cxgridsize,Nmc) 										!transformed future state space for each h,c,MCdraw
	!real(dble) omegaall(size(omega1)+size(omega2)+size(omega3)+size(omega4)) 	!combined state space
	!integer i,j,k,l
	!real(dble) :: inputs(5),constants(3), ages(2), parameters(7),As(2)
	!real(dble) logw0 					! initial logwage of the mother
	
	!omegaall=(/omega1,omega2,omega3,omega4/)
		
	!! fill up the filler matrix with the current space
	!do i=1,3*cxgridsize
		!filler(:,i)=omegaall
	!end do
     !! update the filler matrix with exogenous changes to it.
	!filler(4:6,:)=filler(4:6,:)+1   ! kids and mom age
	!filler(9,:)=filler(9,:)+1 	! father's age.
	!! update last experience
	!filler(7,:)=(/0.0d0*onescxgrid,0.50d0*onescxgrid,1.0d0*onescxgrid/)

	!! update experience level
	!filler(3,6:10)=(omega1(3)+0.5d0)*onescgrid
	!filler(3,11:15)=(omega1(3)+1.0d0)*onescgrid
		
	!inputs=(/omega1(1),omega1(2),0.d0,0.d0,0.d0/)
	!constants=(/omega4(1),omega4(2),logw0/)
	!ages=(/omega2(1),omega2(2)/)
	!parameters=parA(4:10)
	
	!logw0=omega4(3)+parA(1)*omega3(1)+parA(2)*omega3(2)+parA(3)*omega3(3)
	!do k=1,Nmc
		!!use the filler matrix to populate fss 3d matrix
		!fss(:,:,k)=filler
		!do j=1,cgridsize
			!do l=1,xgridsize
				!do i=1,3
					!! inputs to first child, i1=(1-0.5*h)*x
					!inputs(3)=(1-0.5d0*((i-1.0d0)/2.0d0))*xgrid(l)			
					!! inputs to second kid
					!inputs(4)=(1-0.5d0*((i-1.0d0)/2.0d0))*(1-xgrid(l))
					!! goods inputs
					!inputs(5)=(wage(k)*(i-1.0d0)/2+wageh(k))*(1-cgrid(j))
					!As=pftwo(inputs, constants, ages, parameters)
					!fss(1,(i-1)*cgridsize+(j-1)*xgridsize+l,k)=As(1)
					!fss(2,(i-1)*cgridsize+(j-1)*xgridsize+l,k)=As(2)
				!end do
			!end do
		!end do
		!! here we can project the state space matrix (for each h and c) 
		!! for a given MC draw to the space of interpolating coefficients.
		!! But for now I will just use the state space vectors as they are to calculate 
		!! FVs. Replace fss with tfss (transformed future state space when you decide on the 
		!! interpolating matrix
	!end do
		!! ftr_hc works for this one, too.	
		!call ftr_hc(tfss,fss,3*cxgridsize,Nmc)
		!do k=1,Nmc
			!fv(k,:)=matmul(fvpar,tfss(:,:,k))
		!end do
!end subroutine fvhcx

!> this subroutines tries to improve upon the othe fvhcx routine, by calculating the components of the fv one at a time and summing
!>them up instead of dealing with big matrices.
subroutine fvhcxalt(fv,omega1, omega2, omega3, omega4,wage,wageh,fvpar, parA,rho)
	implicit none
	real(dble), intent(out)::fv(Nmc,cxgridsize*3)	!< output, fv for each h,c and mc
	real(dble), intent(in):: omega1(:) 				!< omega1: endogenous dynamic
	real(dble), intent(in):: omega2(:) 				!< omega2: exogenous dynamic
	real(dble), intent(in):: omega3(:) 				!< omega3: static
	real(dble), intent(in):: omega4(:) 				!< omega4: type
	real(dble), intent(in):: fvpar(:) 				!< interpolating parameters	
	real(dble), intent(in):: parA(:) 				!< parameters of the production function
	real(dble), intent(in):: wage(:),wageh(:) 		!< wages of parents, Nmcx1 each
	real(dble), intent(in):: rho 					!< depreciation of skills

	real(dble) omegaind(5)
	integer i,j,k,l
	real(dble) inputs(5), constants(3), age(2), parameters(9), A(2)


	! omegaind: the part of the future value that does not depend on the shocks or the choice. This is just one number, that enters
	! the entire fv matrix additively
	! omegaind=(agem, llms, schm,afqtm,omegaf)
	omegaind=(/omega2(3:4),omega3(1:2),omega3(4)/)
	! increase the ages by one
	! kid's ages can't be in the interpolation because emax is calculated at a specific time and delta, hence these are constant
	! across observations of the regression.
	omegaind(1)=omegaind(1)+1 			! update mom's ag
	! NOTE: not updating llms => myobic expectations, might change this later	
	! now calculate the XB related omegaind, by adding it to intercept of the fvpar
	fv=fvpar(Gsize+1)+sum(omegaind*fvpar(4:8))

	! then we will add the parts which only depend on the experience level, that is, to the choice of h.
	fv(:,1:cxgridsize)=fv(:,1:cxgridsize)+fvpar(3)*omega1(3)+fvpar(11)*omega1(3)**2
	fv(:,cxgridsize+1:2*cxgridsize)=fv(:,cxgridsize+1:2*cxgridsize)+fvpar(3)*(omega1(3)+0.50d0)+fvpar(11)*(omega1(3)+0.5d0)**2
	fv(:,2*cxgridsize+1:3*cxgridsize)=fv(:,2*cxgridsize+1:3*cxgridsize)+fvpar(3)*(omega1(3)+1.0d0)+fvpar(11)*(omega1(3)+1.0d0)**2

	! Now, we have to deal with the parts related to both wage shock and choices. These are going to affect wages and both
	! achievement levels. This requires a big loop over the entire fv
	! matrix.
	
	! before the loop: create the inputs to the child pf.
	inputs=(/omega1(1),omega1(2),0.0d0,0.0d0,0.0d0/)
	constants=omega3(1:3)
	age=(/omega2(1),omega2(2)/)
	parameters=parA(4:12)

	do i=1,Nmc
		do j=1,3
			do k=1,cgridsize
				do l=1,xgridsize
					inputs(3)=(1.0d0-0.5d0*((j-1.0d0)*0.5d0))*xgrid(l)
					inputs(4)=(1.0d0-0.5d0*((j-1.0d0)*0.5d0))*(1-xgrid(l))
					inputs(5)=(wage(i)*(j-1.0d0)*0.5d0+wageh(i))*(1-cgrid(k))
					A=pftwo(inputs, constants,age,parameters,rho)
					! now we can add the parts of the fv that are related to these.
					fv(i,(j-1)*cxgridsize+(k-1)*xgridsize+l)=fv(i,(j-1)*cxgridsize+(k-1)*xgridsize+l) &
						& + fvpar(1)*A(1) + fvpar(2)*A(2) + fvpar(9)*A(1)**2 + fvpar(10)*A(2)**2 + fvpar(12)*A(1)*A(2) &
						& + fvpar(13)*A(1)*(omega1(3)+(j-1)*0.5d0)+ fvpar(14)*A(1)*(omega1(3)+(j-1)*0.5d0)
				end do
			end do
		end do
	end do

	! the end	

end subroutine fvhcxalt

!------------------------------FUTURE VALUES OF ONE CHILD FAMILIES-----------------------


! TERMVAL: for now use the temp one



!> fv calculation for provided interpolating coefficients
!> This is for later periods where mother only chooses h
!> to calculate fv, we need to first update the state space for
!> all the choices and monte carlo draw and then use the interpolating
!> coefficients (whatever they are) and do the interpolation.
!> this is identical to fvlate, I still created this function for the sake of consistency.
subroutine fvoclate(fv,omega1, omega2,omega3, omega4,fvpar)
	implicit none
	real(dble), intent(out):: fv(3) 			!< output, fv for each h 
	real(dble), intent(in):: omega1(:) 			!< omega1: endogenous dynamic
	real(dble), intent(in):: omega2(:) 			!< omega2: exogenous dynamic
	real(dble), intent(in):: omega3(:) 			!< omega3: static
	real(dble), intent(in):: omega4(:) 			!< omega4: type
	real(dble), intent(in):: fvpar(:) 			!< interpolating parameters
	
	!locals
	real(dble) fss((size(omega1)+size(omega2)+size(omega3)+size(omega4)),3) 	!future state space for each h
	real(dble) tfss(Gsizeoc,3) 												!transformed future state space 
	real(dble) omegaall(size(omega1)+size(omega2)+size(omega3)+size(omega4)) 	!combined state space

	omegaall=(/omega1,omega2,omega3,omega4/)
	
	! create alternative specific future state spaces
	fss(:,1)=omegaall
	fss(:,2)=omegaall
	fss(:,3)=omegaall
	! update ages DIFFERENCE: update only the age of the first child, second stays at 0
	fss(4,:)=fss(4,:)+1   ! kid's age
	fss(6,:)=fss(6,:)+1 	! mom's age
	! update experience level
	fss(3,2)=fss(3,2)+0.5d0
	fss(3,3)=fss(3,2)+1.0d0
	
	! NO NEED TO UPDATE achievement because both kids are older than a*
	
	! at this point, I need to take ths fss matrix, and transform for predictions
	! for now, just use the state space as it is.
	! two child transformer will work for now.
	call ftroc_late(tfss,fss,3)	
	fv=matmul(fvpar,tfss)
end subroutine fvoclate


!> future value calculation for one child families while there is no option 
!> of having a second child. Mom chooses h and c. This is almost the same as
!> the two child counterpart, fvhc
subroutine fvochc(fv,omega1, omega2,omega3, omega4,wage,wageh,fvpar,parA,rho)
	implicit none
	real(dble), intent(out):: fv(Nmc,cgridsize*3)	!< output, fv for each h,c and mc
	real(dble), intent(in):: omega1(:) 				!< omega1: endogenous dynamic
	real(dble), intent(in):: omega2(:) 				!< omega2: exogenous dynamic
	real(dble), intent(in):: omega3(:) 				!< omega3: static
	real(dble), intent(in):: omega4(:) 				!< omega4: type
	real(dble), intent(in):: fvpar(:) 				!< interpolating parameters	
	real(dble), intent(in):: parA(:) 				!< parameters of the production function
	real(dble), intent(in):: wage(:),wageh(:) 		!< wages of parents, Nmcx1 each
	real(dble), intent(in) :: rho 					!<skill depreciation=1 for now

	!locals
	real(dble) fss(Ssize,3*cgridsize,Nmc) 	!future state space for each h,c,MCdraw
	real(dble) filler((size(omega1)+size(omega2)+size(omega3)+size(omega4)),3*cgridsize) 	! filler matrix
	real(dble) omegaall(size(omega1)+size(omega2)+size(omega3)+size(omega4)) 	!combined state space
	real(dble) tfss(Gsizeoc+1,3*cgridsize,Nmc) 										!transformed future state space for each h,c,MCdraw
	integer i,j,k
	real(dble) :: inputs(3),constants(3), age, parameters(9),As
	
	omegaall=(/omega1,omega2,omega3,omega4/)
		
	! fill up the filler matrix with the current space
	do i=1,3*cgridsize
		filler(:,i)=omegaall
	end do
 	! update the filler matrix with exogenous changes to it.
 	! difference for one child families, second kid stays at 0.
	filler(4,:)=filler(4,:)+1   ! kid's age
	filler(6,:)=filler(6,:)+1 	! mom's age

	! update experience level
	filler(3,cgridsize+1:2*cgridsize)=(omega1(3)+0.5d0)*onescgrid
	filler(3,2*cgridsize+1:3*cgridsize)=(omega1(3)+1.0d0)*onescgrid
	!filler(3,6:10)=(omega1(3)+0.5d0)*onescgrid
	!filler(3,11:15)=(omega1(3)+1.0d0)*onescgrid
		
	inputs=(/omega1(1),0.d0,0.d0/)
	age=omega2(1)
	constants=omega3(1:3)
	parameters=parA(4:12)
	 
	do k=1,Nmc
		!use the filler matrix to populate fss 3d matrix
		fss(:,:,k)=filler
		do j=1,cgridsize
			do i=1,3
				! inputs to the first kid: first i=h, then g=(1-c)Y 
				inputs(2)=1.0d0-0.5d0*((i-1.0d0)/2.0d0)
				inputs(3)=(wage(k)*(i-1.0d0)/2+wageh(k))*(1-cgrid(j))
				As=pfone(inputs, constants,age,parameters,rho)
				fss(2,(i-1)*cgridsize+j,k)=As
			end do
		end do
		! here we can project the state space matrix (for each h and c) 
		! for a given MC draw to the space of interpolating coefficients.
		! But for now I will just use the state space vectors as they are to calculate 
		! FVs. Replace fss with tfss (transformed future state space when you decide on the 
		! interpolating matrix
		
	!	fv(k,:)=matmul(fvpar,fss(:,:,k))
	end do
		call ftroc_hc(tfss,fss,3*cgridsize,Nmc)
		do k=1,Nmc
			fv(k,:)=matmul(fvpar,tfss(:,:,k))
		end do

end subroutine fvochc

!> calculating future values for families with one kids, during the time where they can have a second child.
!> this version of this subroutine uses the similar ideas as the ones above and employs pretty big matrices unnecessarily to keep
!> things flexible. Once we are settled down with a regression specification for approximation, more efficient subroutines should be
!>used, as in fvochcbalt.


subroutine fvochcb(fv,omega1, omega2,omega3, omega4,wage,wageh,fvpar,parA,parB,typevec,typeprob)
	implicit none
	real(dble), intent(out):: fv(Nmc,cgridsize*6)	!< output, fv for each h,c and mc
	real(dble), intent(in):: omega1(:) 				!< omega1: endogenous dynamic
	real(dble), intent(in):: omega2(:) 				!< omega2: exogenous dynamic
	real(dble), intent(in):: omega3(:) 				!< omega3: static
	real(dble), intent(in):: omega4(:) 				!< omega4: type
	real(dble), intent(in):: fvpar(:,:) 			!< interpolating parameters, now in a matrix form
													! first one for the one child regime and the remainder size(typevec)+1 is for
													! two.
	real(dble), intent(in):: parA(:) 				!< parameters of the production function
	real(dble), intent(in):: wage(:),wageh(:) 		!< wages of parents, Nmcx1 each
	real(dble), intent(in) :: typevec(:)			!< the array of values type of the second child can take
	real(dble), intent(in) :: typeprob(:) 			!< conditional(on mom's type) probability of each type for second child
	real(dble), intent(in) :: parB(Bsize+2) 		!< parameters of the birth probability function
fv=0
end subroutine fvochcb
!  -------------------------ALTERNATIVE VERSION----------------------------
!>future value for the mother while she is contemplating about having a second child
!> or not. state space updating is a bit more complex and fvpar depends on the age of
!> the first born.
!> This version is more efficient(hopefully) but less flexible. I will choose the final, depending on the final
!> choice of the interpolating function


subroutine fvochcbalt(fv,omega1, omega2,omega3, omega4,wage,wageh,fvpar,fvparv,parA,parB,typevec,typeprob,rho)
	implicit none
	real(dble), intent(out):: fv(Nmc,cgridsize*3)	!< output, fv for each h,c and b IF NO BIRTH DECISION *3
	real(dble), intent(in):: omega1(:) 				!< omega1: endogenous dynamic
	real(dble), intent(in):: omega2(:) 				!< omega2: exogenous dynamic
	real(dble), intent(in):: omega3(:) 				!< omega3: static
	real(dble), intent(in):: omega4(:) 				!< omega4: type
	real(dble), intent(in):: fvpar(Gsize+1,nctype) 	!< interpolating parameters, now in a matrix form, size (Gsize,nctype), 
	real(dble), intent(in) :: fvparv(Gsizeoc+1) 	!< interpolating parameters for one child regime.
	real(dble), intent(in):: parA(:) 				!< parameters of the production function
	real(dble), intent(in):: wage(:),wageh(:) 		!< wages of parents, Nmcx1 each
	real(dble), intent(in) :: typevec(:)			!< the array of values type of the second child can take, =nctype
	real(dble), intent(in) :: typeprob(:) 			!< conditional(on mom's type) probability of each type for second child
	real(dble), intent(in) :: parB(Bsizeexo+1) 		!< parameters of the birth probability function
	real(dble), intent(in) :: rho 					!< skill depreciation parameter
	! locals
	
	integer i,j,k,l,m
	real(dble) :: inputs(3),constants(3), age, parameters(9)
	real(dble) logw0 					! initial logwage of the mother

	real(dble) omegaind(5)  ! size of shock-decision independent part of the state space	 	
	! new locals
	real(dble) fvbig(Nmc,cgridsize*3,nctype+1) ! first element of the third dimension is for staying with one child.
	real(dble) As, Ainit
	real(dble) p0,p1, omegaB(Bsizeexo)

	! The idea is to calculate the fv in parts. fvpar coincides with the layout of the transformed state space, which also covers
	! the same omegas in the same way in omega1, omega2, omega3, omega4
	! first calculate the future state space points that does not depend on the choices

	! update: I will kick out age1, age2, age0m from the interpolating function. interpolating function described in jmpsimple.pdf
	
	! 1: insert the omegaind part=(agem, llms, schm,afqtm,omegaf)
	omegaind=(/omega2(3:4),omega3(1:2),omega3(4)/)
	
	! kid's ages can't be in the interpolation because emax is calculated at a specific time and delta, hence these are constant
	! across observations of the regression.
	omegaind(1)=omegaind(1)+1 			! update mom's age
	! NOTE: not updating llms => myobic expectations, might change this later	

	! 2: extract the relevant part from the fvpar vector and multiply with the transformed omegaind to fill in the first part
	! of the FV.  Add these to the intercept to iniatialize fvbig.
	fvbig(:,:,1)=sum(omegaind*(fvparv(4:8)))+fvparv(Gsizeoc+1)
	do i=1,nctype
		fvbig(:,:,i+1)=sum(omegaind*(fvpar(4:8,i)))+fvpar(Gsize+1,i)
	end do 

	! now subtract the 1*fvpar(6,1) from the one child family. NO need to do this anymore, age 2 is not in the interpolation.
	!fvbig(:,:,1)=fvbig(:,:,1)-fvpar(6,1) 


	! then add h_t-1 component+and the experience derived component NO MORE h_t-1 component for jmpsimple, but we have experience
	! squared. (cgridsize*2 if with b choice) (see jmpext september 2012 commits)

	! first staying in one child regime
	fvbig(:,1:cgridsize,1)=fvbig(:,1:cgridsize,1)+fvparv(2)*omega1(3)+fvparv(9)*omega1(3)**2
	fvbig(:,cgridsize+1:2*cgridsize,1)=fvbig(:,cgridsize+1:2*cgridsize,1)+fvparv(2)*(omega1(3)+0.5d0)+fvparv(9)*(omega1(3)+0.5d0)**2
	fvbig(:,2*cgridsize+1:3*cgridsize,1)=fvbig(:,2*cgridsize+1:3*cgridsize,1)+fvparv(2)*(omega1(3)+1.0d0)+fvparv(9)*(omega1(3)+1.0d0)**2

	! then two children, with the possible different types
	do i=1,nctype
		fvbig(:,1:cgridsize,i+1)=fvbig(:,1:cgridsize,i+1)+fvpar(3,i)*omega1(3)+fvpar(11,i)*omega1(3)**2
		fvbig(:,cgridsize+1:2*cgridsize,i+1)=fvbig(:,cgridsize+1:2*cgridsize,i+1)+fvpar(3,i)*(omega1(3)+0.5d0)+fvpar(11,i)*(omega1(3)+0.5d0)**2
		fvbig(:,2*cgridsize+1:3*cgridsize,i+1)=fvbig(:,2*cgridsize+1:3*cgridsize,i+1)+fvpar(3,i)*(omega1(3)+1.0d0)+fvpar(11,i)*(omega1(3)+1.0d0)**2
	end do	
	! lastly add the component of FV that is derived from A's

	inputs=(/omega1(1),0.d0,0.d0/)
	constants=omega3(1:3)
	age=omega2(1)
	parameters=parA(4:12)
	
	! first: if the family stayed one child, A2 will stay zero.
	do j=1,Nmc
		do k=1,3
			do l=1,cgridsize
				!do m=1,1
					inputs(2)=1.0d0-0.5d0*((k-1.0d0)/2.0d0)
					inputs(3)=(wage(j)*(k-1.0d0)/2+wageh(j))*(1-cgrid(l))
					As=pfone(inputs, constants,age,parameters,rho)
					! first staying one child
					fvbig(j,(k-1)*cgridsize+l,1)=fvbig(j,(k-1)*cgridsize+l,1)+As*fvpar(1,i)+As*(omega1(3)+(k-1.0d0)*0.5d0)*fvpar(10,i)
					! then becoming two child family
					do i=1,nctype ! don't have types for jmpsimple, but keep the compatible for the future.
						! how is my kid gonna be? depends on her type. with no UH, typevec is just the intercept of pf.
						Ainit=typevec(i)+parA(1)*omega3(1)+parA(2)*omega3(2)+parA(3)*omega3(3) 
						! then add the parts with As (and interactions with E)
						fvbig(j,(k-1)*cgridsize+l,i+1)=fvbig(j,(k-1)*cgridsize+l,i+1) + As*fvpar(1,i) + Ainit*fvpar(2,i) + As*Ainit*fvpar(12,i) &
							& + As*(omega1(3)+(k-1.0d0)*0.5d0)*fvpar(13,i) + Ainit*(omega1(3)+(k-1.0d0)*0.5d0)*fvpar(14,i) 
					end do
				!end do
			end do
		end do
	end do
	! calculate the birth probability. With p1, we gonna have a second child and move to fvbig(:,:,2) , and with prob (1-p1) we are
	! going to get the fvbig(:,:,1), so we need to weigh those in and then sum them up(along the third dimension) to calculate the
	! expected future value.
	omegaB=(/omega3(3),omega3(3)*omega3(3),omega3(1),omega3(2)/)
	p1=bprobexo(omegaB,parB)

	! NOTE: I am still retaining type code here but without unobserved heterogeneity I will just pass the intercept as the typevec
	! and typeprob=1 as intercept

	!first: one child
	fvbig(:,:,1)=(1-p1)*fvbig(:,:,1)

	do i=2,nctype+1
		fvbig(:,:,i)=fvbig(:,:,i)*typeprob(i-1)*p1
	end do
	! AND FINALLY
	fv=sum(fvbig,3)
end subroutine fvochcbalt
!--------------------------------------------EMAX CALCULATION----------------------------------------------------------

!--------------------------------------------TWO KIDS FIRST---------------------------------



!> final period emax calculation
!> choice set is minimal (just h) and future emax's are not needed.

function emaxfinal(omega1,omega2,omega3,omega4,eps,parA,parU,parW,parH,beta)
	implicit none
	real(dble), intent(in):: omega1(:) 						!< omega1: endogenous dynamic
	real(dble), intent(in):: omega2(:) 						!< omega2: exogenous dynamic
	real(dble), intent(in):: omega3(:) 						!< omega3: static
	real(dble), intent(in):: omega4(:) 						!< omega4: types
	real(dble), intent(in):: eps(Nmc,shocksize1) 			!< random draw for MC integration
	real(dble) parA(:), parU(:),parW(:),parH(:) 			!< parameter vectors except for Sigma and beta
	real(dble) beta 										!< discount factor
	real(dble) emaxfinal 									!< output

	real(dble) uc
	real(dble) wage(Nmc), wageh(Nmc)
	real(dble), dimension(Nmc)::uh0,uhp,uhf
	real(dble), dimension(Nmc,3):: TV, umat

	! calculate utilities
	! uc is shock independent
	uc=uctwo(omega1(1),omega1(2),parU(1))
	! uh is shock dependent and need to calculate wage income for each mc draw, which will equal to consumption because there is no child anymore.
	wage=wagef(omega1(3),omega2(1), omega2(4), omega3(1), omega3(2), omega3(3),eps(:,4),omega4(3),parW(5:7),parW(1:4))  
	! same with husband income
	wageh=wagehfquick(omega3(4), omega2(1), eps(:,5),parH)  
	! uh's for each h choice 
	uh0=uh(0.0d0,(wage+wageh),omega2(1:2),parU(3:7))
	uhp=uh(0.5d0,(wage+wageh),omega2(1:2),parU(3:7))
	uhf=uh(1.0d0,(wage+wageh),omega2(1:2),parU(3:7))
	! terminal values calcuted, for now using a very simple thing
	call termvaltemp(TV,uc,uh0,uhp,uhf,beta,omega2(3)) ! TODO TEMP TERM VALUES
	! distribute choice specific current returns to a Nmcx3 matrix
	umat(:,1)=uc+uh0+eps(:,1)
	umat(:,2)=uc+uhp+eps(:,2)
	umat(:,3)=uc+uhf+eps(:,3)

	! add up the TV
	umat=umat+beta*TV

	! now, take maxval of each row and then take the mean
	emaxfinal=(sum(maxval(umat,2)))/Nmc
end function emaxfinal

!> emax calculations for families where both kids are out of the nest and 
!> mother is choosing only the labor supply
!> the difference with emaxfinal is: FV vs TV
function emaxlate(omega1,omega2,omega3,omega4,eps,parA,parU,parW,parH,beta,parFV)
	implicit none
	real(dble), intent(in):: omega1(:) 						!< omega1: endogenous dynamic
	real(dble), intent(in):: omega2(:) 						!< omega2: exogenous dynamic
	real(dble), intent(in):: omega3(:) 						!< omega3: static
	real(dble), intent(in):: omega4(:) 						!< omega4: types
	real(dble), intent(in):: eps(Nmc,shocksize1) 			!< random draw for MC integration
	real(dble) parA(:), parU(:),parW(:),parH(:), parFV(:)	!< parameter vectors except for Sigma and beta
	real(dble) beta 										!< discount factor
	real(dble) emaxlate 									!< output

	real(dble) uc
	real(dble) wage(Nmc), wageh(Nmc)
	real(dble), dimension(Nmc)::uh0,uhp,uhf
	real(dble), dimension(Nmc,3):: FV, umat

	! calculate utilities
	! uc, current period utility from children, is shock independent
	uc=uctwo(omega1(1),omega1(2),parU(1))
	! uh is shock dependent and need to calculate wage income for each mc draw, which will equal to consumption because there is no child anymore.
	wage=wagef(omega1(3),omega2(1), omega2(4), omega3(1), omega3(2), omega3(3),eps(:,4),omega4(3),parW(5:7),parW(1:4))  
	! same with husband income
	wageh=wagehfquick(omega3(4), omega2(1), eps(:,5),parH) 
	! uh's for each h choice 
	! income=consumption
	uh0=uh(0.0d0,(wage+wageh),omega2(1:2),parU(3:7))
	uhp=uh(0.5d0,(wage+wageh),omega2(1:2),parU(3:7))
	uhf=uh(1.0d0,(wage+wageh),omega2(1:2),parU(3:7))
	! future values needed.
	call fvlate(FV,omega1, omega2,omega3,omega4,parFV)
	! distribute choice specific current returns to a Nmcx3 matrix
	umat(:,1)=uc+uh0+eps(:,1)
	umat(:,2)=uc+uhp+eps(:,2)
	umat(:,3)=uc+uhf+eps(:,3)

	! add up the TV
	umat=umat+beta*FV
	
	! now, take maxval of each row and then take the mean
	emaxlate=(sum(maxval(umat,2)))/Nmc

end function emaxlate



!> function to calculate emax for the periods when firstchild is not affected by
!>the mother but second child is. Mother choose her labor supply and consumption
!>ratio, hence the name hc (h and c)  
function emaxhc(omega1,omega2,omega3,omega4,eps,parA,parU,parW,parH,beta,parFV,rho)
	implicit none
	real(dble), intent(in):: omega1(:) 						!< omega1: endogenous dynamic
	real(dble), intent(in):: omega2(:) 						!< omega2: exogenous dynamic
	real(dble), intent(in):: omega3(:) 						!< omega3: static
	real(dble), intent(in):: omega4(:) 						!< omega4: types
	real(dble), intent(in):: eps(Nmc,shocksize1) 			!< random draw for MC integration
	real(dble) parA(:), parU(:),parW(:),parH(:), parFV(:)	!< parameter vectors except for Sigma and beta
	real(dble) beta 										!< discount factor
	real(dble) emaxhc 									!< output
	real(dble) rho
	!locals

	real(dble) uc
	real(dble) wage(Nmc), wageh(Nmc)
	real(dble), dimension(Nmc,cgridsize)::uh0,uhp,uhf
	integer i
	real(dble) fv(Nmc,cgridsize*3)    
	real(dble) logw0
	real(dble) umat(Nmc,cgridsize*3)
	
	
	! calculate utilities
	! uc, current period utility from children, is shock independent
	uc=uctwo(omega1(1),omega1(2),parU(1))
	! uh is shock dependent and need to calculate wage income for each mc draw, which will equal to consumption because there is no child anymore.
	wage=wagef(omega1(3),omega2(1), omega2(4), omega3(1), omega3(2), omega3(3),eps(:,4),omega4(3),parW(5:7),parW(1:4))  
	! same with husband income
	wageh=wagehfquick(omega3(4), omega2(1), eps(:,5),parH) 

	! calculate utility levels associated with each of the h,c bundles.
	do i=1,cgridsize
		uh0(:,i)=uh(0.0d0,(wageh)*cgrid(i),omega2(1:2),parU(3:7))
		uhp(:,i)=uh(0.50d0,(wage*0.50d0+wageh)*cgrid(i),omega2(1:2),parU(3:7))
		uhf(:,i)=uh(1.0d0,(wage+wageh)*cgrid(i),omega2(1:2),parU(3:7))
	end do 
	
	! get the future values associated with each h,c
	call fvhc(fv,omega1, omega2,omega3,omega4,wage,wageh,parFV,parA,rho)
	! now collect everything together with the shocks, too!!!
	umat(:,1:cgridsize)=uh0+uc+sidestack(eps(:,1),cgridsize)
	umat(:,cgridsize+1:2*cgridsize)=uhp+uc+sidestack(eps(:,2),cgridsize)
	umat(:,2*cgridsize+1:3*cgridsize)=uhf+uc+sidestack(eps(:,3),cgridsize)
	! add the fv
	umat=umat+beta*fv
	! now, take maxval of each row and then take the mean
	emaxhc=(sum(maxval(umat,2)))/Nmc
end function emaxhc

!> function to calculate emax for the periods when both kids are taking both
!> time and goods inputs. Mother choose her labor supply, consumption
!> ratio c and the ratio of time that goes to child 1. hence the name hcx (h and c and x)  
function emaxhcx(omega1,omega2,omega3,omega4,eps,parA,parU,parW,parH,beta,parFV,rho)
	implicit none
	real(dble), intent(in):: omega1(:) 						!< omega1: endogenous dynamic
	real(dble), intent(in):: omega2(:) 						!< omega2: exogenous dynamic
	real(dble), intent(in):: omega3(:) 						!< omega3: static
	real(dble), intent(in):: omega4(:) 						!< omega4: types
	real(dble), intent(in):: eps(Nmc,shocksize1) 			!< random draw for MC integration
	real(dble) parA(:), parU(:),parW(:),parH(:), parFV(:)	!< parameter vectors except for Sigma and beta
	real(dble) beta 										!< discount factor
	real(dble) rho 											!< skill depreciation parameter, for now=1
	real(dble) emaxhcx 									!< output
	
	!locals

	real(dble) uc
	real(dble) wage(Nmc), wageh(Nmc)
	real(dble), dimension(Nmc,cxgridsize)::uh0,uhp,uhf
	integer i,j
	real(dble) fv(Nmc,cxgridsize*3)
	real(dble) logw0
	real(dble) umat(Nmc,cxgridsize*3)

		
	! calculate utilities
	! uc, current period utility from children, is shock independent
	uc=uctwo(omega1(1),omega1(2),parU(1))
	! uh is shock dependent and need to calculate wage income for each mc draw, which will equal to consumption because there is no child anymore.
	wage=wagef(omega1(3),omega2(1), omega2(4), omega3(1), omega3(2), omega3(3),eps(:,4),omega4(3),parW(5:7),parW(1:4))  
	! same with husband income
	wageh=wagehfquick(omega3(4), omega2(1), eps(:,5),parH)  

	! calculate utility levels associated with each of the h,c,x bundles. And fill Nmc x cgridsize*xgridsize matrices
	do i=1,cgridsize  
		do j=1,xgridsize
			uh0(:,(i-1)*xgridsize+j)=uh(0.0d0,(wageh)*cgrid(i),omega2(1:2),parU(3:7))
			uhp(:,(i-1)*xgridsize+j)=uh(0.50d0,(wage*0.50d0+wageh)*cgrid(i),omega2(1:2),parU(3:7))
			uhf(:,(i-1)*xgridsize+j)=uh(1.0d0,(wage+wageh)*cgrid(i),omega2(1:2),parU(3:7))
		end do
	end do 
	
	! get the future values associated with each h,c
	! first calculate logw0
	call fvhcxalt(fv,omega1, omega2,omega3, omega4,wage,wageh,parFV,parA,rho)
	
	! now collect everything together
	! So build the ingridients for a Nmcx cxgridsize (*3 matrix), including the shocks!!! 	
	umat(:,1:cxgridsize)=uh0+uc+sidestack(eps(:,1),cxgridsize)
	umat(:,cxgridsize+1:2*cxgridsize)=uhp+uc+sidestack(eps(:,2),cxgridsize)
	umat(:,2*cxgridsize+1:3*cxgridsize)=uhf+uc+sidestack(eps(:,3),cxgridsize)
	
	! add the fv
	umat=umat+beta*fv
	! now, take maxval of each row and then take the mean
	emaxhcx=sum(maxval(umat,2))/Nmc
end function emaxhcx


! -----------------------EMAXs for one child families---------------------------


!> function to calculate the final period emax for the families with one child.
!> It is almost exactly the same as the emaxfinal but for only one child.
!> terminal values to updated later
function emaxocfinal(omega1,omega2,omega3,omega4,eps,parA,parU,parW,parH,beta)
	implicit none
	real(dble), intent(in):: omega1(:) 						!< omega1: endogenous dynamic
	real(dble), intent(in):: omega2(:) 						!< omega2: exogenous dynamic
	real(dble), intent(in):: omega3(:) 						!< omega3: static
	real(dble), intent(in):: omega4(:) 						!< omega4: types
	real(dble), intent(in):: eps(Nmc,shocksize1) 			!< random draw for MC integration
	real(dble) parA(:), parU(:),parW(:),parH(:) 			!< parameter vectors except for Sigma and beta
	real(dble) beta 										!< discount factor
	real(dble) emaxocfinal 									!< output

	real(dble) uc
	real(dble) wage(Nmc), wageh(Nmc)
	real(dble), dimension(Nmc)::uh0,uhp,uhf
	real(dble), dimension(Nmc,3):: TV, umat

	! calculate utilities
	! uc is shock independent
	uc=parU(2)*omega1(1)/emaxscale
	
	! uh is shock dependent and need to calculate wage income for each mc draw, which will equal to consumption because there is no child anymore.
	wage=wagef(omega1(3),omega2(1), omega2(4), omega3(1), omega3(2), omega3(3),eps(:,4),omega4(3),parW(5:7),parW(1:4))  
	! same with husband income
	wageh=wagehfquick(omega3(4), omega2(1), eps(:,5),parH) 
		! uh's for each h choice 
	uh0=uh(0.0d0,(wageh),omega2(1:2),parU(3:7))
	uhp=uh(0.5d0,(wage*0.50d0+wageh),omega2(1:2),parU(3:7))
	uhf=uh(1.0d0,(wage+wageh),omega2(1:2),parU(3:7))
	! terminal values calcuted, for now using a very simple thing
	call termvaltemp(TV,uc,uh0,uhp,uhf,beta,omega2(3))
	! distribute choice specific current returns to a Nmcx3 matrix
	
	umat(:,1)=uc+uh0+eps(:,1)
	umat(:,2)=uc+uhp+eps(:,2)
	umat(:,3)=uc+uhf+eps(:,3)
	! add up the TV
	umat=umat+beta*TV

	! now, take maxval of each row and then take the mean
	emaxocfinal=(sum(maxval(umat,2)))/Nmc

end function emaxocfinal

!> emax calculation for periods where mother WITH A SINGLE child only choose the labor supply
!> almost the same as 
function emaxoclate(omega1,omega2,omega3,omega4,eps,parA,parU,parW,parH,beta,parFV)
	implicit none
	real(dble), intent(in):: omega1(:) 						!< omega1: endogenous dynamic
	real(dble), intent(in):: omega2(:) 						!< omega2: exogenous dynamic
	real(dble), intent(in):: omega3(:) 						!< omega3: static
	real(dble), intent(in):: omega4(:) 						!< omega4: types
	real(dble), intent(in):: eps(Nmc,shocksize1) 			!< random draw for MC integration
	real(dble) parA(:), parU(:),parW(:),parH(:), parFV(:)	!< parameter vectors except for Sigma and beta
	real(dble) beta 										!< discount factor
	real(dble) emaxoclate 									!< output

	real(dble) uc
	real(dble) wage(Nmc), wageh(Nmc)
	real(dble), dimension(Nmc)::uh0,uhp,uhf
	real(dble), dimension(Nmc,3):: FV, umat

	! calculate utilities
	! uc, current period utility from children, is shock independent
	uc=parU(2)*omega1(1)/emaxscale
	 
	! uh is shock dependent and need to calculate wage income for each mc draw, which will equal to consumption because there is no child anymore.
	wage=wagef(omega1(3),omega2(1), omega2(4), omega3(1), omega3(2), omega3(3),eps(:,4),omega4(3),parW(5:7),parW(1:4))  
	! same with husband income
	wageh=wagehfquick(omega3(4), omega2(1), eps(:,5),parH) 

	! uh's for each h choice 
	uh0=uhoc(0.0d0,(wageh),omega2(1),parU(3:7))
	uhp=uhoc(0.5d0,(wage*0.5d0+wageh),omega2(1),parU(3:7))
	uhf=uhoc(1.0d0,(wage+wageh),omega2(1),parU(3:7))
	! future values calcuted and added.
	call fvoclate(FV,omega1, omega2,omega3,omega4,parFV)
	! distribute choice specific current returns to a Nmcx3 matrix
	umat(:,1)=uc+uh0+eps(:,1)
	umat(:,2)=uc+uhp+eps(:,2)
	umat(:,3)=uc+uhf+eps(:,3)

	! add up the FV
	umat=umat+beta*FV
	
	! now, take maxval of each row and then take the mean
	emaxoclate=(sum(maxval(umat,2)))/Nmc

end function emaxoclate

!> emax for single child mothers, taking care of their child and working so that
!> the choice set is h,c. It will be very similar to the two child version above.
!> the main difference is A2 does not enter anywhere and stays the same at 0.
function emaxochc(omega1,omega2,omega3,omega4,eps,parA,parU,parW,parH,beta,parFV,rho)
	implicit none
	real(dble), intent(in):: omega1(:) 						!< omega1: endogenous dynamic
	real(dble), intent(in):: omega2(:) 						!< omega2: exogenous dynamic
	real(dble), intent(in):: omega3(:) 						!< omega3: static
	real(dble), intent(in):: omega4(:) 						!< omega4: types
	real(dble), intent(in):: eps(Nmc,shocksize1) 			!< random draw for MC integration
	real(dble) parA(:), parU(:),parW(:),parH(:), parFV(:)	!< parameter vectors except for Sigma and beta
	real(dble) beta 										!< discount factor
	real(dble) emaxochc 									!< output
	real(dble) rho 									!<skill depreciation assumed =1 so far.
	!locals

	real(dble) uc
	real(dble) wage(Nmc), wageh(Nmc)
	real(dble), dimension(Nmc,cgridsize)::uh0,uhp,uhf
	integer i
	real(dble) fv(Nmc,cgridsize*3)
	real(dble) logw0
	real(dble) umat(Nmc,cgridsize*3)
	! calculate utilities
	! uc, current period utility from children, is shock independent
	uc=parU(2)*omega1(1)/emaxscale

	 
	! uh is shock dependent and need to calculate wage income for each mc draw, which will equal to consumption because there is no child anymore.
	wage=wagef(omega1(3),omega2(1), omega2(4), omega3(1), omega3(2), omega3(3),eps(:,4),omega4(3),parW(5:7),parW(1:4))  
	! same with husband income
	wageh=wagehfquick(omega3(4), omega2(1), eps(:,5),parH) 

	! calculate utility levels associated with each of the h,c bundles.
	do i=1,cgridsize
		uh0(:,i)=uhoc(0.0d0,(wageh)*cgrid(i),omega2(1),parU(3:7))
		uhp(:,i)=uhoc(0.50d0,(wage*0.5d0+wageh)*cgrid(i),omega2(1),parU(3:7))
		uhf(:,i)=uhoc(1.0d0,(wage+wageh)*cgrid(i),omega2(1),parU(3:7))
	end do 
	
	! get the future values associated with each h,c
	! first calculate logw0
	call fvochc(fv,omega1, omega2,omega3, omega4,wage,wageh,parFV,parA,rho)
	
	! now collect everything together
	umat(:,1:cgridsize)=uh0+uc+sidestack(eps(:,1),cgridsize)
	umat(:,cgridsize+1:2*cgridsize)=uhp+uc+sidestack(eps(:,2),cgridsize)
	umat(:,2*cgridsize+1:3*cgridsize)=uhf+uc+sidestack(eps(:,3),cgridsize)
	! add the fv
	umat=umat+beta*fv
	! now, take maxval of each row and then take the mean
	emaxochc=(sum(maxval(umat,2)))/Nmc


end function emaxochc


!> the final emax of the single mother's. This is when the second kid is still 
!> a possibility. So the mother chooses not only h and c but also b, her contraceptive
!> choice. 

function emaxochcb(omega1,omega2,omega3,omega4,eps,parA,parU,parW,parH,beta,parFV,parFVv,parB,typevec,typeprob,rho)
	implicit none
	
	real(dble) emaxochcb 									!< output
	real(dble), intent(in):: omega1(:) 						!< omega1: endogenous dynamic
	real(dble), intent(in):: omega2(:) 						!< omega2: exogenous dynamic
	real(dble), intent(in):: omega3(:) 						!< omega3: static
	real(dble), intent(in):: omega4(:) 						!< omega4: types
	real(dble), intent(in):: eps(Nmc,shocksize1) 			!< random draw for MC integration for h and b (NOTE: NO B FOR THE jmpsimple)
															!< shocks for contraception decision.
	real(dble),intent(in)::parA(:), parU(:),parW(:),parH(:) 	!< parameter vectors except for Sigma and beta
	real(dble),intent(in)::parFV(Gsize+1,nctype) 						!< interpolating function parameter, second dimension for child types
	real(dble), intent(in) :: parFVv(Gsizeoc+1) 					!< interpolating function parameters for one child regime.
	real(dble),intent(in)::beta 				!< discount factor
	real(dble),intent(in)::typevec(:)			!< the array of values type of the second child can take
	real(dble),intent(in)::typeprob(:) 			!< conditional(on mom's type) probability of each type for second child
	real(dble),intent(in)::parB(Bsizeexo+1) 		!< parameters of the birth probability function 
	real(dble), intent(in) :: rho
	!locals

	real(dble) uc
	real(dble) wage(Nmc), wageh(Nmc)
	real(dble), dimension(Nmc,cgridsize)::uh0,uhp,uhf  	! removed fertility decision
	integer i,j
	real(dble) fv(Nmc,cgridsize*3)
	real(dble) logw0
	real(dble) umat(Nmc,cgridsize*3)
	
	! calculate utilities
	! uc, current period utility from children, is shock independent
	uc=parU(2)*omega1(1)/emaxscale
	! uh is shock dependent and need to calculate wage income for each mc draw, which will equal to consumption because there is no child anymore.
	wage=wagef(omega1(3),omega2(1), omega2(4), omega3(1), omega3(2), omega3(3),eps(:,4),omega4(3),parW(5:7),parW(1:4))  
	! same with husband income
	wageh=wagehfquick(omega3(4), omega2(1), eps(:,5),parH) ! 
 	
 	! CAREFUL WITH THE SHOCKS	
	! calculate the utility levels associated with each of h,c and b bundles. For each h, create a Nmc x cgridsize x 2 matrix
	! Removing b choice will remoce the last bit of x2 there. See sept 2012 jmpext commits for the version with b in it.
	!do i=1,cgridsize
		!do j=1,2
			!uh0(:,(i-1)*2+j)=uhoc(0.0d0,(wageh)*cgrid(i),omega2(1),parU(3:7))
			!uhp(:,(i-1)*2+j)=uhoc(0.50d0,(wage*0.5d0+wageh)*cgrid(i),omega2(1),parU(3:7))
			!uhf(:,(i-1)*2+j)=uhoc(1.0d0,(wage+wageh)*cgrid(i),omega2(1),parU(3:7))
		!end do
		!! here let's add the shocks
		!uh0(:,(i-1)*2+1:(i-1)*2+2)=uh0(:,(i-1)*2+1:(i-1)*2+2)+eps(:,1:2)
		!uhp(:,(i-1)*2+1:(i-1)*2+2)=uhp(:,(i-1)*2+1:(i-1)*2+2)+eps(:,3:4)
		!uhf(:,(i-1)*2+1:(i-1)*2+2)=uhf(:,(i-1)*2+1:(i-1)*2+2)+eps(:,5:6)
	!end do
	do i=1,cgridsize
		uh0(:,i)=uh(0.0d0,(wageh)*cgrid(i),omega2(1:2),parU(3:7))
		uhp(:,i)=uh(0.50d0,(wage*0.50d0+wageh)*cgrid(i),omega2(1:2),parU(3:7))
		uhf(:,i)=uh(1.0d0,(wage+wageh)*cgrid(i),omega2(1:2),parU(3:7))
	end do 

	! now collect everything together
	umat(:,1:cgridsize)=uh0+uc+sidestack(eps(:,1),cgridsize)
	umat(:,cgridsize+1:2*cgridsize)=uhp+uc+sidestack(eps(:,2),cgridsize)
	umat(:,2*cgridsize+1:3*cgridsize)=uhf+uc+sidestack(eps(:,3),cgridsize)

	! get the future values
	call fvochcbalt(fv,omega1, omega2,omega3, omega4,wage,wageh,parFV,parFVv,parA,parB,typevec,typeprob,rho)

	
	! add the fv
	umat=umat+beta*fv
	! now, take maxval of each row and then take the mean
	emaxochcb=(sum(maxval(umat,2)))/Nmc

end function emaxochcb



!                                          END OF EMAX  FUNCTIONS!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! --------------------------------------------------------------------------------------------------------------

! --------------------------------------------SOLVERS-----------------------------------------------------------
! EXPLANATION: following functions solve for the emax interpolating coefficients for different time periods.

!> This routine calculates interpolating coefficients for the final period of two child families.
!> This problem is identical regardless of the delta type of the family, so I will use one interpolating 
!> regression for all. THIS IS ALTERNATIVE TO SOLVING THE COEF FOR EACH DELTA SEPARATELY
!subroutine coef_finalalt(coef, mutype,parA,parU,parW,parH,beta,sigma)
	!implicit none
	!real(dble), intent(out) :: coef(Gsize+1)
	!real(dble), intent(in) :: mutype(:)
	!real(dble), intent(in):: parA(:), parU(:),parW(:),parH(:) 			!< parameter vectors except for Sigma and beta
	!real(dble), intent(in):: beta 										!< discount factor
	!real(dble), intent(in):: sigma(shocksize1,shocksize1) 									!< variance covariance matrix

	!! locals
	!real(dble) vemax(nintpfinal)            ! to collect emaxs, vector emax
	!real(dble) mss(nintpfinal,Ssize) 		! to collect SS elements, matrix state space
	!real(dble) tmss(nintpfinal, Gsize+1) 	! transformed ss
	!real(dble) llmsvec(nintpfinal) 		 	! llms vector to drawn from a uniform (0,0.1) :unemploymen rate	
	!integer i,j,k,l,m,n,p,r,q,s
	!real(dble) veca1(5), veca2(5), vecE(5), veclagh(3)
	!real(dble) omega1(o1size), omega2(o2size), omega3(4)
	!integer counter
	!real(dble) mu(shocksize1)
	!real(dble) eps(Nmc,shocksize1)

	!! lapack stuff
     !real(dble) work(Gsize+Gsize*blocksize)
     !integer info
	!! setting up the A1, A2, and E vecs and veclagh
	!veca1=(/1.0d0,2.0d0,3.0d0,4.0d0,5.0d0/)
	!veca2=(/1.0d0,2.0d0,3.0d0,4.0d0,5.0d0/)
	!vecE=(/1.0d0,2.0d0,3.0d0,4.0d0,5.0d0/)*4.d0
	!veclagh=(/0.0d0,0.50d0,1.0d0/)
	!! get a set of llms, which are fixed
	!call random_seed(put=(/22,1/))
	!call random_number(llmsvec)

	!mu=0.0d0
	!counter=1

	!do n=1,3 										! delta
	!do p=1,5 										! A1
		!do q=1,5 									! A2 
			!do r=1,5 								! E
				!do s=1,3 							! lagh	
					!do i=1,4  						! age0m
						!do j=1,3 					! age0f
							!do k=1,3 				! schoolingm
								!do l=1,4 			! afqt
									!do m=1,3 		! schooligf
										!omega1=(/veca1(p),veca2(q),vecE(r),veclagh(s)/)	
										!omega2=(/22.0d0,22.0d0-vecdelta(n)+1,22+vecage0m(i)-1,llmsvec(counter),vecage0f(j)+22-1/)
										!omega3=(/vecsch0m(k), vecaqft(l),vecage0m(i),vecsch0f(m)/)
										!eps=randmnv(Nmc,shocksize1,mu,sigma, 3,1,(/22,counter/))
										!!eps=randmnv(Nmc,shocksize1,mu,sigma, 3,1,gseed)
										!vemax(counter)=emaxfinal(omega1,omega2,omega3,mutype,eps,parA,parU,parW,parH,beta)
										!mss(counter,:)=(/omega1,omega2,omega3,mutype/)	
										!print*, 'counter at', counter, 'calculated', vemax(counter)
										!counter=counter+1
									!end do
								!end do
							!end do
						!end do
					!end do
				!end do
			!end do
		!end do
	!end do
	!end do
	!call tr_late(tmss, mss,nintpfinal)
	!call DGELS('N',nintpfinal,Gsize,1,tmss, nintpfinal, vemax, nintpfinal,work, Gsize+Gsize*blocksize,info)
	!coef=vemax(1:Gsize+1)
!end subroutine coef_finalalt

subroutine coeffinal(coef,period,delta, mutype,parA,parU,parW,parH,beta,sigma)
	implicit none
	real(dble), intent(out) :: coef(Gsize+1)
	real(dble), intent(in) :: mutype(:) 								!< mu type of the family =omega4
	real(dble), intent(in):: parA(:), parU(:),parW(:),parH(:) 			!< parameter vectors except for Sigma and beta
	real(dble), intent(in):: beta 										!< discount factor
	real(dble), intent(in):: sigma(shocksize1,shocksize1) 				!< variance covariance matrix
	real(dble), intent(in)::period 										!< what period we are in, or age of the first born
	real(dble), intent(in)::delta 										!< delta type of the family
	! locals
	real(dble) vemax(nintp)            ! to collect emaxs, vector emax
	real(dble) mss(nintp,Ssize) 		! to collect SS elements, matrix state space
	real(dble) llmsvec(nintp) 		 	! llms vector to drawn from a uniform (0,0.1) :unemploymen rate	
	integer i,j,k,l,m,n,p,r,q,s
	real(dble) vecE(5)
	real(dble) omega1(o1size) ,omega2(o2size), omega3(o3size)
	integer counter
	real(dble) mu(shocksize1)
	real(dble) eps(Nmc,shocksize1)
	real(dble) tmss(nintp,Gsize+1)

	! lapack stuff
 	real(dble) work(Gsize+Gsize*blocksize)
 	integer info
	! setting up the A1, A2, and E vecs and veclagh
	vecE=(/1.0d0,2.0d0,3.0d0,4.0d0,5.0d0/)*(period-1)/5.d0
	
	call setvecAs(period, delta,2)
	! get a set of llms, which are fixed
	call random_seed(put=(/nint(1000*period+delta),1/))
	call random_number(llmsvec)
	llmsvec=llmsvec*0.1d0 			!< adjust llms to have 0-10% unemployment rate 

	mu=0.0d0
	counter=1

	do p=1,sveca1 											! A1
		do q=1,sveca2										! A2 
			do r=1,svecE									! E
					do i=1,svecage0m  						! age0m
							do k=1,svecsch0m 				! schoolingm
								do l=1,svecaqft 			! afqt
									do m=1,svecomegaf		! father constants
										omega1=(/veca1(p),veca2(q),vecE(r)/)	
										omega2=(/period,(period-delta+1),period+vecage0m(i)-1,llmsvec(counter)/)
										omega3=(/vecsch0m(k), vecaqft(l),vecage0m(i),vecomegaf(m)/)
										eps=randmnv(Nmc,shocksize1,mu,sigma, 3,1,(/nint(period*10+delta),counter/))
										vemax(counter)=emaxfinal(omega1,omega2,omega3,mutype,eps,parA,parU,parW,parH,beta)
										mss(counter,:)=(/omega1,omega2,omega3,mutype/)	
										!print *, 'counter at', counter, 'calculated', vemax(counter)
										counter=counter+1
									end do
								end do
							end do
						!end do
					end do
				!end do
			end do
		end do
	end do

	call tr_late(tmss,mss,nintp) 			! transform the state space 
	
	!open(12, file = 'a.txt')
	!do k=1,nintp
	!write(12, 900)  (tmss(k,j) , j=1,Gsize+1)
	!end do
	!900 format(16f30.10)
	!close(12)
	!open(13, file = 'b.txt')
	!write(13, "(f30.10)")  (vemax(k) , k=1,nintp)
	!close(13)


	
	call DGELS('N',nintp,Gsize+1,1,tmss, nintp, vemax, nintp,work, Gsize+(Gsize)*blocksize,info)
	coef=vemax(1:Gsize+1)
	if (info .NE. 0)	print*, 'final period DGELS exploded in period, delta', period, delta
	end subroutine coeffinal


subroutine coeflate(coef,period,delta, mutype,parA,parU,parW,parH,beta,sigma,parFV)
	implicit none
	real(dble), intent(out) :: coef(Gsize+1)
	real(dble), intent(in) :: mutype(:)
	real(dble), intent(in):: parA(:), parU(:),parW(:),parH(:), parFV(:)	!< parameter vectors except for Sigma and beta
	real(dble), intent(in):: beta 										!< discount factor
	real(dble), intent(in):: sigma(shocksize1,shocksize1) 				!< variance covariance matrix
	real(dble), intent(in)::period 										!< what period we are in, or age of the first born
	real(dble), intent(in)::delta 										!< delta type of the family
	! locals
	real(dble) vemax(nintp)            ! to collect emaxs, vector emax
	real(dble) mss(nintp,Ssize) 		! to collect SS elements, matrix state space
	real(dble) llmsvec(nintp) 		 	! llms vector to drawn from a uniform (0,0.1) :unemploymen rate	
	integer i,j,k,l,m,n,p,r,q,s
	real(dble) vecE(5)                                                                     
	real(dble) omega1(o1size) ,omega2(o2size), omega3(o3size)
	integer counter
	real(dble) mu(shocksize1)
	real(dble) eps(Nmc,shocksize1)
	real(dble) tmss(nintp,Gsize+1)

	! lapack stuff
 	real(dble) work(Gsize+Gsize*blocksize)
 	integer info
	! setting up the A1, A2, and E vecs and veclagh
	vecE=(/1.0d0,2.0d0,3.0d0,4.0d0,5.0d0/)*(period-1)/5.d0
	call setvecAs(period, delta,2)
	
	! get a set of llms, which are fixed
	call random_seed(put=(/nint(1000*period+delta),1/))
	call random_number(llmsvec)
	llmsvec=llmsvec*0.1d0  		! unemployment rate max 10%
	mu=0.0d0
	counter=1
	do p=1,sveca1 											! A1
		do q=1,sveca2										! A2 
			do r=1,svecE									! E
					do i=1,svecage0m  						! age0m
						!do j=1,svecage0f 					! age0f
							do k=1,svecsch0m 				! schoolingm
								do l=1,svecaqft 			! afqt
									do m=1,svecomegaf		! schooligf
										omega1=(/veca1(p),veca2(q),vecE(r)/)	
										omega2=(/period,(period-delta+1),period+vecage0m(i)-1,llmsvec(counter)/)
										omega3=(/vecsch0m(k), vecaqft(l),vecage0m(i),vecomegaf(m)/)
										eps=randmnv(Nmc,shocksize1,mu,sigma, 3,1,(/nint(period*10+delta),counter/))
										!eps=randmnv(Nmc,shocksize1,mu,sigma, 3,1,gseed)
										vemax(counter)=emaxlate(omega1,omega2,omega3,mutype,eps,parA,parU,parW,parH,beta,parFV)
										mss(counter,:)=(/omega1,omega2,omega3,mutype/)	
										!print *, 'counter at', counter, 'calculated', vemax(counter)
										counter=counter+1
									end do
								end do
							end do
						!end do
					end do
				!end do
			end do
		end do
	end do

	call tr_late(tmss,mss,nintp) 			! transform the state space 
	call DGELS('N',nintp,Gsize+1,1,tmss, nintp, vemax, nintp,work, Gsize+(Gsize)*blocksize,info)
	if (info .NE. 0)	print*, 'late DGELS exploded in period, delta', period, delta
	coef=vemax(1:Gsize+1)
	
	!open(12, file = 'a.txt')
	!do k=1,nintp
	!write(12, 900)  (tmss(k,j) , j=1,Gsize-1)
	!end do
	!900 format(15f9.5)
	!close(12)
	!open(13, file = 'b.txt')
	!write(13, "(f10.5)")  (vemax(k) , k=1,nintp)
	!close(13)
end subroutine coeflate


subroutine coefhc(coef,period,delta, mutype,parA,parU,parW,parH,beta,sigma,parFV,rho)
	implicit none
	real(dble), intent(out) :: coef(Gsize+1)
	real(dble), intent(in) :: mutype(:)
	real(dble), intent(in):: parA(:), parU(:),parW(:),parH(:), parFV(:)	!< parameter vectors except for Sigma and beta
	real(dble), intent(in):: beta 										!< discount factor
	real(dble), intent(in):: sigma(shocksize1,shocksize1) 				!< variance covariance matrix
	real(dble), intent(in)::period 										!< what period we are in, or age of the first born
	real(dble), intent(in)::delta 										!< delta type of the family
	real(dble), intent(in) :: rho
	! locals
	real(dble) vemax(nintp)            ! to collect emaxs, vector emax
	real(dble) mss(nintp,Ssize) 		! to collect SS elements, matrix state space
	real(dble) llmsvec(nintp) 		 	! llms vector to drawn from a uniform (0,0.1) :unemploymen rate	
	integer i,j,k,l,m,n,p,r,q,s
	real(dble) vecE(5)
	real(dble) omega1(o1size) ,omega2(o2size), omega3(o3size)
	integer counter
	real(dble) mu(shocksize1)
	real(dble) eps(Nmc,shocksize1)
	real(dble) tmss(nintp,Gsize+1)

	! lapack stuff
 	real(dble) work(Gsize+Gsize*blocksize)
 	integer info
	! setting up the A1, A2, and E vecs and veclagh
	vecE=(/1.0d0,2.0d0,3.0d0,4.0d0,5.0d0/)*(period-1)/5.d0
	call setvecAs(period, delta,2)

	! get a set of llms, which are fixed
	call random_seed(put=(/nint(1000*period+delta),1/))
	call random_number(llmsvec)
	llmsvec=llmsvec*0.1d0    ! 10% unemployment rate max.
	mu=0.0d0
	counter=1
	do p=1,sveca1 											! A1
		do q=1,sveca2										! A2 
			do r=1,svecE									! E
				!do s=1,sveclagh 							! lagh	
					do i=1,svecage0m  						! age0m
						!do j=1,svecage0f 					! age0f
							do k=1,svecsch0m 				! schoolingm
								do l=1,svecaqft 			! afqt
									do m=1,svecomegaf		! schooligf => now omegaf
										omega1=(/veca1(p),veca2(q),vecE(r)/)	
										omega2=(/period,(period-delta+1),period+vecage0m(i)-1,llmsvec(counter)/)
										omega3=(/vecsch0m(k), vecaqft(l),vecage0m(i),vecomegaf(m)/)
										eps=randmnv(Nmc,shocksize1,mu,sigma, 3,1,(/nint(period*10+delta),counter/))
										!eps=randmnv(Nmc,shocksize1,mu,sigma, 3,1,gseed)
										vemax(counter)=emaxhc(omega1,omega2,omega3,mutype,eps,parA,parU,parW,parH,beta,parFV,rho)
										mss(counter,:)=(/omega1,omega2,omega3,mutype/)	
									!	print *, 'counter at', counter, 'calculated', vemax(counter)
										counter=counter+1
									end do
								end do
							end do
						!end do
					end do
				!end do
			end do
		end do
	end do

	call tr_late(tmss,mss,nintp) 			! transform the state space 
	call DGELS('N',nintp,Gsize+1,1,tmss, nintp, vemax, nintp,work, Gsize+(Gsize)*blocksize,info)
	coef=vemax(1:Gsize+1)
	if (info .NE. 0)	print*, 'hc DGELS exploded in period,delta', period, delta
end subroutine coefhc

subroutine coefhcx(coef,period,delta, mutype,parA,parU,parW,parH,beta,sigma,parFV,rho)
	implicit none
	real(dble), intent(out) :: coef(Gsize+1)
	real(dble), intent(in) :: mutype(:)
	real(dble), intent(in):: parA(:), parU(:),parW(:),parH(:), parFV(:)	!< parameter vectors except for Sigma and beta
	real(dble), intent(in):: beta 										!< discount factor
	real(dble), intent(in):: sigma(shocksize1,shocksize1) 				!< variance covariance matrix
	real(dble), intent(in)::period 										!< what period we are in, or age of the first born
	real(dble), intent(in)::delta 										!< delta type of the family
	real(dble), intent(in) :: rho
	! locals
	real(dble) vemax(nintp)            ! to collect emaxs, vector emax
	real(dble) mss(nintp,Ssize) 		! to collect SS elements, matrix state space
	real(dble) llmsvec(nintp) 		 	! llms vector to drawn from a uniform (0,0.1) :unemploymen rate	
	integer i,j,k,l,m,n,p,r,q,s
	real(dble) vecE(svecE)
	real(dble) omega1(o1size) ,omega2(o2size), omega3(o3size)
	integer counter
	real(dble) mu(shocksize1)
	real(dble) eps(Nmc,shocksize1)
	real(dble) tmss(nintp,Gsize+1)

	! lapack stuff
 	real(dble) work(Gsize+Gsize*blocksize)
 	integer info
	! setting up the A1, A2, and E vecs and veclagh
	vecE=(/1.0d0,2.0d0,3.0d0,4.0d0,5.0d0/)*(period-1)/5.d0
	call setvecAs(period, delta,2)
	! get a set of llms, which are fixed
	call random_seed(put=(/nint(1000*period+delta),1/))
	call random_number(llmsvec)
	llmsvec=llmsvec*0.1d0    ! 10% unemployment rate max.
	mu=0.0d0
	counter=1
	do p=1,sveca1 											! A1
		do q=1,sveca2										! A2 
			do r=1,svecE									! E
				!do s=1,sveclagh 							! lagh	
					do i=1,svecage0m  						! age0m
						!do j=1,svecage0f 					! age0f
							do k=1,svecsch0m 				! schoolingm
								do l=1,svecaqft 			! afqt
									do m=1,svecomegaf		! schooligf
										omega1=(/veca1(p),veca2(q),vecE(r)/)	
										omega2=(/period,(period-delta+1),period+vecage0m(i)-1,llmsvec(counter)/)
										omega3=(/vecsch0m(k), vecaqft(l),vecage0m(i),vecomegaf(m)/)
										eps=randmnv(Nmc,shocksize1,mu,sigma, 3,1,(/nint(period*10+delta),counter/))
										!eps=randmnv(Nmc,shocksize1,mu,sigma, 3,1,gseed)
										vemax(counter)=emaxhcx(omega1,omega2,omega3,mutype,eps,parA,parU,parW,parH,beta,parFV,rho)
										mss(counter,:)=(/omega1,omega2,omega3,mutype/)	
										!print *, 'counter at', counter, 'calculated', vemax(counter)
										counter=counter+1
									end do
								end do
							!end do
						end do
					end do
				!end do
			end do
		end do
	end do
	 

	call tr_late(tmss,mss,nintp) 			! transform the state space 
	!if ((nint(period)==7) .AND. (nint(delta)==4)) then
		!open(12, file = 'a.txt')
		!do k=1,nintp
			!write(12, 900)  (tmss(k,j) , j=1,Gsize+1)
		!end do
		!900 format(16f30.10)
		!close(12)
		!open(13, file = 'b.txt')
		!write(13, "(f30.10)")  (vemax(k) , k=1,nintp)
		!close(13)
	!end if 
call DGELS('N',nintp,Gsize+1,1,tmss, nintp, vemax, nintp,work, Gsize+(Gsize)*blocksize,info)
	coef=vemax(1:Gsize+1)

	if (info .NE. 0)	print*, 'hcx DGELS exploded in period,delta', period,delta
end subroutine coefhcx
!---------------------------------------------------------------------------------------------------
!---------------------                  ONE CHILD VERSIONS            ------------------------------
!---------------------------------------------------------------------------------------------------


!> WITH ONE CHILD, so no need for A2 or age2, get rid of those parts of the loops


!> one child version of coefoc
subroutine coefocfinal(coef,period, mutype,parA,parU,parW,parH,beta,sigma)
	implicit none
	real(dble), intent(out) :: coef(Gsizeoc+1)
	real(dble), intent(in) :: mutype(:) 								!< mu1,0.0d0,mum
	real(dble), intent(in):: parA(:), parU(:),parW(:),parH(:) 			!< parameter vectors except for Sigma and beta
	real(dble), intent(in):: beta 										!< discount factor
	real(dble), intent(in):: sigma(shocksize1,shocksize1) 				!< variance covariance matrix
	real(dble), intent(in)::period 										!< what period we are in, or age of the first born
	! locals
	real(dble) vemax(nintpoc)            ! to collect emaxs, vector emax
	real(dble) mss(nintpoc,Ssize) 		! to collect SS elements, matrix state space
	real(dble) llmsvec(nintpoc) 		 	! llms vector to drawn from a uniform (0,0.1) :unemploymen rate	
	integer i,j,k,l,m,n,p,r,q,s
	real(dble) vecE(5)
	real(dble) omega1(o1size) ,omega2(o2size), omega3(o3size)
	integer counter
	real(dble) mu(shocksize1)
	real(dble) eps(Nmc,shocksize1)
	real(dble) tmss(nintpoc,Gsizeoc+1)

	! lapack stuff
 	real(dble) work(Gsizeoc+Gsizeoc*blocksize)
 	integer info
	! setting up the A1, A2, and E vecs and veclagh
	vecE=(/1.0d0,2.0d0,3.0d0,4.0d0,5.0d0/)*(period-1)/5.d0

	call setvecAs(period, 0.0d0,1)
! get a set of llms, which are fixed
	call random_seed(put=(/nint(1000*period+1),1/))
	call random_number(llmsvec)
	llmsvec=llmsvec*0.1d0    ! 10% unemployment rate max.
	
	mu=0.0d0
	counter=1
	do p=1,sveca1 											! A1
!		do q=1,sveca2										! A2 
			do r=1,svecE									! E
				!do s=1,sveclagh 							! lagh	
					do i=1,svecage0m  						! age0m
						!do j=1,svecage0f 					! age0f
							do k=1,svecsch0m 				! schoolingm
								do l=1,svecaqft 			! afqt
									do m=1,svecomegaf		! schooligf
										omega1=(/veca1(p),0.0d0,vecE(r)/)	
										omega2=(/period,0.0d0,period+vecage0m(i)-1,llmsvec(counter)/)
										omega3=(/vecsch0m(k), vecaqft(l),vecage0m(i),vecomegaf(m)/)
										eps=randmnv(Nmc,shocksize1,mu,sigma, 3,1,(/nint(period*10+1),counter/))
										!eps=randmnv(Nmc,shocksize1,mu,sigma, 3,1,gseed)
										!if  ((rank==1) .AND. (counter==1)) then
											!print*, '--------------CORE------------------'
											!print*, omega1
											!print*, omega2
											!print*, omega3
											!print*, mu
											!print*, eps(1:10,1)
											!print*, '--------------CORE------------------'
										!end if
										vemax(counter)=emaxocfinal(omega1,omega2,omega3,mutype,eps,parA,parU,parW,parH,beta)
										mss(counter,:)=(/omega1,omega2,omega3,mutype/)	
										!print *, 'counter at', counter, 'calculated', vemax(counter)
										counter=counter+1
									end do
								end do
							end do
						!end do
					end do
				!end do
			end do
		!end do
	end do
	call troc_late(tmss,mss,nintpoc) 			! transform the state space 
	!if (rank==1) then
		!open(22, file = 'a.txt')
		!do k=1,nintpoc
			!write(22, 900)  (tmss(k,j) , j=1,Gsizeoc+1)
		!end do
		!close(22)
		!open(23, file = 'b.txt')
		!do k=1,nintpoc
			!write(23, 900)  vemax(k)
		!end do
		!900 format(16f20.10)
		!close(23)
	!end if
close(rank)
	
	call DGELS('N',nintpoc,Gsizeoc+1,1, tmss, nintpoc, vemax, nintpoc,work, Gsizeoc+(Gsizeoc)*blocksize,info)
	coef=vemax(1:Gsizeoc+1)
	if (info .NE. 0)	print*, 'OC: final DGELS exploded in period', period
end subroutine coefocfinal

subroutine coefoclate(coef,period, mutype,parA,parU,parW,parH,beta,sigma,parFV)
	implicit none
	real(dble), intent(out) :: coef(Gsizeoc+1)
	real(dble), intent(in) :: mutype(:)
	real(dble), intent(in):: parA(:), parU(:),parW(:),parH(:), parFV(:)	!< parameter vectors except for Sigma and beta
	real(dble), intent(in):: beta 										!< discount factor
	real(dble), intent(in):: sigma(shocksize1,shocksize1) 				!< variance covariance matrix
	real(dble), intent(in)::period 										!< what period we are in, or age of the first born
	! locals
	real(dble) vemax(nintpoc)            ! to collect emaxs, vector emax
	real(dble) mss(nintpoc,Ssize) 		! to collect SS elements, matrix state space
	real(dble) llmsvec(nintpoc) 		 	! llms vector to drawn from a uniform (0,0.1) :unemploymen rate	
	integer i,j,k,l,m,n,p,r,q,s
	real(dble) vecE(5)
	real(dble) omega1(o1size) ,omega2(o2size), omega3(o3size)
	integer counter
	real(dble) mu(shocksize1)
	real(dble) eps(Nmc,shocksize1)
	real(dble) tmss(nintpoc,Gsizeoc+1)

	! lapack stuff
 	real(dble) work(Gsizeoc+Gsizeoc*blocksize)
 	integer info
	! setting up the A1, A2, and E vecs and veclagh
	vecE=(/1.0d0,2.0d0,3.0d0,4.0d0,5.0d0/)*(period-1)/5.d0

	call setvecAs(period, 0.0d0,1)
	! get a set of llms, which are fixed
	call random_seed(put=(/nint(1000*period+1),1/))
	call random_number(llmsvec)
	llmsvec=llmsvec*0.1d0    ! 10% unemployment rate max.
	mu=0.0d0
	counter=1
	do p=1,sveca1 											! A1
		!do q=1,sveca2										! A2 
			do r=1,svecE									! E
				!do s=1,sveclagh 							! lagh	
					do i=1,svecage0m  						! age0m
						!do j=1,svecage0f 					! age0f
							do k=1,svecsch0m 				! schoolingm
								do l=1,svecaqft 			! afqt
									do m=1,svecomegaf		! schooligf
										omega1=(/veca1(p),0.0d0,vecE(r)/)	
										omega2=(/period,0.0d0,period+vecage0m(i)-1,llmsvec(counter)/)
										omega3=(/vecsch0m(k), vecaqft(l),vecage0m(i),vecomegaf(m)/)
										eps=randmnv(Nmc,shocksize1,mu,sigma, 3,1,(/nint(period*10+1),counter/))
										!eps=randmnv(Nmc,shocksize1,mu,sigma, 3,1,gseed)
										vemax(counter)=emaxoclate(omega1,omega2,omega3,mutype,eps,parA,parU,parW,parH,beta,parFV)
										mss(counter,:)=(/omega1,omega2,omega3,mutype/)	
										!print *, 'counter at', counter, 'calculated', vemax(counter)
										counter=counter+1
									end do
								end do
							end do
						!end do
					end do
				!end do
			end do
		!end do
	end do

	call troc_late(tmss,mss,nintpoc) 			! transform the state space FIXME
	call DGELS('N',nintpoc,Gsizeoc+1,1,tmss, nintpoc, vemax, nintpoc,work, Gsizeoc+(Gsizeoc)*blocksize,info)
	! add a zero for the second child A2
	if (info .NE. 0)	print*, 'OC: late DGELS exploded in period', period
	coef=vemax(1:Gsizeoc+1)
end subroutine coefoclate


subroutine coefochc(coef,period, mutype,parA,parU,parW,parH,beta,sigma,parFV,rho)
	implicit none
	real(dble), intent(out) :: coef(Gsizeoc+1)
	real(dble), intent(in) :: mutype(:)
	real(dble), intent(in):: parA(:), parU(:),parW(:),parH(:), parFV(:)	!< parameter vectors except for Sigma and beta
	real(dble), intent(in):: beta 										!< discount factor
	real(dble), intent(in):: sigma(shocksize1,shocksize1) 				!< variance covariance matrix
	real(dble), intent(in)::period 										!< what period we are in, or age of the first born
	real(dble), intent(in) :: rho
	! locals
	real(dble) vemax(nintpoc)            ! to collect emaxs, vector emax
	real(dble) mss(nintpoc,Ssize) 		! to collect SS elements, matrix state space
	real(dble) llmsvec(nintpoc) 		 	! llms vector to drawn from a uniform (0,0.1) :unemploymen rate	
	integer i,j,k,l,m,n,p,r,q,s
	real(dble) vecE(5)
	real(dble) omega1(o1size) ,omega2(o2size), omega3(o3size)
	integer counter
	real(dble) mu(shocksize1)
	real(dble) eps(Nmc,shocksize1)
	real(dble) tmss(nintpoc,Gsizeoc+1)

	! lapack stuff
 	real(dble) work(Gsizeoc+Gsizeoc*blocksize)
 	integer info
	! setting up the A1, A2, and E vecs and veclagh
	vecE=(/1.0d0,2.0d0,3.0d0,4.0d0,5.0d0/)*(period-1)/5.d0

	call setvecAs(period, 0.0d0,1)
	! get a set of llms, which are fixed
	call random_seed(put=(/nint(1000*period+1),1/))
	call random_number(llmsvec)
	llmsvec=llmsvec*0.1d0    ! 10% unemployment rate max.
	
	mu=0.0d0
	counter=1
	do p=1,sveca1 											! A1
		!do q=1,sveca2										! A2 
			do r=1,svecE									! E
				!do s=1,sveclagh 							! lagh	
					do i=1,svecage0m  						! age0m
						!do j=1,svecage0f 					! age0f
							do k=1,svecsch0m 				! schoolingm
								do l=1,svecaqft 			! afqt
									do m=1,svecomegaf		! schooligf
										omega1=(/veca1(p),0.0d0,vecE(r)/)	
										omega2=(/period,0.0d0,period+vecage0m(i)-1,llmsvec(counter)/)
										omega3=(/vecsch0m(k), vecaqft(l),vecage0m(i),vecomegaf(m)/)
										eps=randmnv(Nmc,shocksize1,mu,sigma, 3,1,(/nint(period*10+1),counter/))
										!eps=randmnv(Nmc,shocksize1,mu,sigma, 3,1,gseed)
										vemax(counter)=emaxochc(omega1,omega2,omega3,mutype,eps,parA,parU,parW,parH,beta,parFV,rho)
										mss(counter,:)=(/omega1,omega2,omega3,mutype/)	
										!print *, 'counter at', counter, 'calculated', vemax(counter)
										counter=counter+1
									end do
								end do
							end do
						!end do
					end do
				!end do
			end do
		!end do
	end do

	call troc_late(tmss,mss,nintpoc) 			! transform the state space 
	call DGELS('N',nintpoc,Gsizeoc+1,1,tmss, nintpoc, vemax, nintpoc,work, Gsizeoc+(Gsizeoc)*blocksize,info)
	if (info .NE. 0)	print*, 'OC: hc DGELS exploded in period', period
	coef=vemax(1:Gsizeoc+1)
end subroutine coefochc




!> subroutine to pick some ss points and calculate emaxs and get regression coefficients for one child families fertile period
subroutine coefochcb(coef,period, mutype,parA,parU,parW,parH,beta,sigma,parFV,parFVv,parBmat,typevec,typeprob,rho)
	implicit none
	real(dble), intent(out) :: coef(Gsizeoc+1)
	real(dble), intent(in) :: mutype(:)
	real(dble), intent(in):: parA(:), parU(:),parW(:),parH(:), parFV(Gsize+1,nctype),parFVv(Gsizeoc+1)	!< parameter vectors except for Sigma and beta																	 
	real(dble), intent(in):: beta 										!< discount factor 
	real(dble), intent(in) :: sigma(shocksize1,shocksize1)
	real(dble), intent(in)::period 										!< what period we are in, or age of the first born
	real(dble), intent(in) :: typevec(nctype),typeprob(nctype),parBmat(Bsizeexo+1,nfert)
	real(dble), intent(in) :: rho
	! locals
	real(dble) vemax(nintpoc)            ! to collect emaxs, vector emax
	real(dble) mss(nintpoc,Ssize) 		! to collect SS elements, matrix state space
	real(dble) llmsvec(nintpoc) 		 	! llms vector to drawn from a uniform (0,0.1) :unemploymen rate	
	integer i,j,k,l,m,n,p,r,q,s
	real(dble) vecE(5)
	real(dble) omega1(o1size) ,omega2(o2size), omega3(o3size)
	integer counter
	real(dble) mu(shocksize1)
	real(dble) eps(Nmc,shocksize1)
	real(dble) tmss(nintpoc,Gsizeoc+1)

	! lapack stuff
 	real(dble) work(Gsizeoc+Gsizeoc*blocksize)
 	integer info
	! setting up the A1, A2, and E vecs and veclagh
	vecE=(/1.0d0,2.0d0,3.0d0,4.0d0,5.0d0/)*(period-1)/5.d0
	! TODO for earlier periods, I have too many E points, and at period=1, it actually explodes, need to fix this.
	if (nint(period)==1) vecE=(/1.0d0,2.0d0,3.0d0,4.0d0,5.0d0/)*(2-1)/5.d0
	call setvecAs(period, 0.0d0,1)
	! get a set of llms, which are fixed
	call random_seed(put=(/nint(1000*period+1),1/))
	call random_number(llmsvec)
	llmsvec=llmsvec*0.1d0    ! 10% unemployment rate max.
	
	mu=0.0d0
	counter=1
	do p=1,sveca1 											! A1
		!do q=1,sveca2										! A2 
			do r=1,svecE									! E
				!do s=1,sveclagh 							! lagh	
					do i=1,svecage0m  						! age0m
						!do j=1,svecage0f 					! age0f
							do k=1,svecsch0m 				! schoolingm
								do l=1,svecaqft 			! afqt
									do m=1,svecomegaf		! schooligf
										omega1=(/veca1(p),0.0d0,vecE(r)/)	
										omega2=(/period,0.0d0,period+vecage0m(i)-1,llmsvec(counter)/)
										omega3=(/vecsch0m(k), vecaqft(l),vecage0m(i),vecomegaf(m)/)
										eps=randmnv(Nmc,shocksize1,mu,sigma, 3,1,(/nint(period*10+1),counter/))
										!eps=randmnv(Nmc,shocksize1,mu,sigma, 3,1,gseed)
										vemax(counter)=emaxochcb(omega1,omega2,omega3,mutype,eps,parA,parU,parW,parH,beta,parFV,parFVv,parBmat(:,nint(period)),typevec,typeprob,rho)
										mss(counter,:)=(/omega1,omega2,omega3,mutype/)	
										!print *, 'counter at', counter, 'calculated', vemax(counter)
										counter=counter+1
									end do
								end do
							end do
						!end do
					end do
				!end do
			end do
		!end do
	end do
	call troc_late(tmss,mss,nintpoc) 			! transform the state space troc_late
	!if (nint(period)==1) then
		!print*, rank,	size(tmss,1), size(tmss,2)
		!do k=1,nintpoc
			!write(rank, 900)  (tmss(k,j) , j=1,Gsizeoc+1)
		!end do
	!end if
	!900 format(16f20.10)
	
	call DGELS('N',nintpoc,Gsizeoc+1,1,tmss, nintpoc, vemax, nintp,work, Gsizeoc+(Gsizeoc)*blocksize,info)
	if (info .NE. 0) then
		print*, 'OC: hcb DGELS exploded in period', period
	end if
	coef=vemax(1:Gsizeoc+1)
end subroutine coefochcb

! ####################################END OF COEF_ ROUTINES##############################################
!>Main subroutine that solves the model for an unobserved type of a two child family for a given delta. It produces a set of interpolating
!>coefficients for each time period for each type of delta. Difference between the old one and this is that this is for only one
!>delta.
subroutine wsolver(solw,delta,ftype,parA,parW,parH,parU,beta,Sigma,rho)
	implicit none
	real(dble), intent(out):: solw(Gsize+1, nperiods-deltamin+2) !< collects solution coef.
	real(dble), intent(in) ::  ftype(5) 				!< family unobserved type. mu1,mu2,mu_m,alpha,alpha1
	real(dble), intent(in) :: parA(:),parW(:),parH(:),parU(:),beta,Sigma(:,:), rho
	real(dble), intent(in) :: delta 
	! locals
	integer period
	real(dble) coef(Gsize+1)
	real(dble) coefnew(Gsize+1)
	real(dble) paractualU(size(parU))
	! -------------------------ACTION-----------------------------------
	! 
	! fill up the sol with -10^7
	solw=-10000000.0d0
	period=22
	! insert the heterogenous parameters to their places
	paractualU=parU	
	paractualU(2:3)=ftype(4:5)
	
	! the big loop to get all the coefficients
	call coeffinal(coef,22.0d0,delta*1.0d0,ftype(1:3),parA,paractualU,parW,parH,beta,sigma)
	
	! store these guys
	 solw(:,22)=coef
	do period=21, delta, -1
		! if the second kid is also out of the zone of mother's influence, use coeflate
		if (period>=astar+delta-1) then 			
			call coeflate(coefnew,period*1.0d0,delta*1.0d0, ftype(1:3),parA,paractualU,parW,parH,beta,sigma,coef)
		!print*, '------------------ delta= ', delta, ' --------------------'
		!print*, 'Calculated coef for late period'
		!print*, 'period is ', period
		!print*, 'coef is ', coefnew
!		 second kid is influenced by the mother but not the first one.
		else if ((period<astar+delta-1).AND.(astar<=period)) then
			call coefhc(coefnew, period*1.0d0, delta*1.0d0, ftype(1:3),parA, paractualU, parW,parH,beta,sigma, coef,rho)
		!print*, '------------------ delta= ', delta, ' --------------------'
		!print*, 'Calculated coef for hc period'
		!print*, 'period is ', period
		!print*, 'coef is ', coefnew
		!! period less than 15 but >=delta, so mom also chooses x
		else
			call coefhcx(coefnew, period*1.0d0, delta*1.0d0, ftype(1:3),parA, paractualU, parW,parH,beta,sigma, coef,rho)
		!print*, '------------------ delta= ', delta, ' --------------------'
		!print*, 'Calculated coef for hcx period'
		!print*, 'period is ', period
		!print*, 'coef is ', coefnew
		end if
		
		! now store these coefficients properly
		solw(:,period)=coefnew
		! and update the coef vector for use in the previous period
		coef=coefnew
	end do
end subroutine wsolver

!> Main routine that will get the interpolating coefficients for one child families. Main difference of it from wsolver is the fact
!> that it needs the output from the wsolver for each type of second child.
subroutine vsolver(solv,ftype,parA,parW,parH,parU,parBmat,beta,Sigma, wcoef,typevec,typeprob,rho)
	implicit none
	real(dble), intent(out):: solv(Gsizeoc+1, nperiods) !< collects solution coef.
	real(dble), intent(in) ::  ftype(5) 				!< family unobserved type. mu1,mu2,mu_m,alpha,alpha1
	real(dble), intent(in) :: parA(:),parW(:),parH(:),parU(:),beta,Sigma(:,:),parBmat(Bsizeexo+1,nfert),typevec(:), typeprob(:),rho
	real(dble), intent(in) :: wcoef(Gsize+1, nperiods-deltamin+2, deltamax-deltamin+1,nctype)

	! THE DIFFERENCE: ftype has zero for the second child's A. It is as (mu1,0.0d0, mum, alpha, alpha1)
	! locals
	
	integer period, delta
	real(dble) coef(Gsizeoc+1)
	real(dble) coefnew(Gsizeoc+1)
	real(dble) paractualU(size(parU))
	real(dble) parFV(Gsize+1, nctype+1)
	! fill up the sol with -10^7
	solv=-10000000.0d0
	period=22
	! insert the heterogenous parameters to their places
	paractualU=parU	
	paractualU(2:3)=ftype(4:5)
	! first get the final coefficient.
	call coefocfinal(coef,22.0d0,ftype(1:3),parA,paractualU,parW,parH,beta,sigma)
	! store these guys
	solv(:,22)=coef
	do period=21,1, -1
		! if the kid is  out of the zone of mother's influence, use coefoclate
		if (period>=astar) then 			
			call coefoclate(coefnew,period*1.0d0, ftype(1:3),parA,paractualU,parW,parH,beta,sigma,coef)
			!print*, 'Calculated coef for late period'
			!print*, 'period is ', period
			!print*, 'late period=',period,'coef is ', coefnew
		! first kid is affected by the mom, but mom is out of the fecund period. 
		else if ((period<astar).AND.(deltamax<=period)) then
			call coefochc(coefnew, period*1.0d0, ftype(1:3),parA, paractualU, parW,parH,beta,sigma,coef,rho)
			!print*, 'Calculated coef for hc  period'
			!print*, 'period is ', period
			!print*, 'coef is ', coefnew
		else
		! period less than astar, and mother still decides on contraception so a second baby is still in cards.
		! this one is complicated because it has to properly take not only the coefficients estimated in the earlier iterations in
		! the loop but also coefficients estimated by the wsolver, because mothers can switch to a two child regime. And this
		! requires the knowledge w coefficients for all possible types of second child.
	
		! one child regime fv parameters
		! and parameters of fv for two child regime, with all possible child types for the second one.
			call coefochcb(coefnew,period*1.0d0, ftype(1:3),parA,paractualU,parW,parH,beta,sigma,wcoef(:,period+1,period,:),coef,parBmat,typevec,typeprob,rho)
			
			!print*, 'Calculated coef for hcB  period'
			!print*, 'period is ', period
			!print*, 'coef is ', coefnew
		end if
		
		! now store these coefficients properly
		solv(:,period)=coefnew
		! and update the coef vector for use in the previous period
		coef=coefnew
	end do
end subroutine vsolver 


! ##############      END OF SOLVERS     #########################
!********************************************************************

!-------------------INTERPOLATION TRANSFORMING ROUTINES---------------------------
! tr_x : these subroutines are for the coef routines, they prepare the ss and make them the transformed ss , in preparation for the
! regressions.
! ftr_x: these subroutines are for fv routines. They transform the current state so that fv can be predicted via regression
! interpolations.


!> transform the state space for the regression !!!!NOT USED AT THE MOMENT!!!!
subroutine tr_final(tmss,mss,n)
	implicit none
	integer, intent(in):: n
	real(dble), intent(in) :: mss(n,Ssize)
	real(dble), intent(out) :: tmss(n,Gsize+1)

	! for now: use  
	tmss(:,1:4)=mss(:,1:4)
	tmss(:,5:Gsize-1)=mss(:,7:11)
	tmss(:,Gsize)=mss(:,13)	
	tmss(:,Gsize+1)=mss(:,5)-mss(:,6)+1.0d0	
end subroutine tr_final

!> called tr_late but it will work for periods, too. It takes
!> matrix of ss variables (at which the emax's are calculated) and converts 
!> it to the tmss matrix for regression. 
subroutine tr_late(tmss, mss,n)
	implicit none
	integer, intent(in):: n 				!< the number of rows in the matrix to be transformed
	real(dble), intent(in) :: mss(n,Ssize) 	!< matrix holding a subset of ss in its rows
	real(dble), intent(out) :: tmss(n,Gsize+1) !<regression matrix 
		
	! for now: use a subset of SS
	tmss(:,1:3)=mss(:,1:3) 				!< A1, A2, E
	tmss(:,4:7)=mss(:,6:9) 		!< age of the children are constant for a period and delta, so skip those 
 	tmss(:,8)=mss(:,11)
 	tmss(:,9:11)=mss(:,1:3)**2
 	tmss(:,12)=mss(:,1)*mss(:,2)
 	tmss(:,13)=mss(:,1)*mss(:,3)
 	tmss(:,14)=mss(:,2)*mss(:,3)
	tmss(:,Gsize+1)=1.0d0 				!< NEW: add the intercept
end subroutine tr_late

!> ONE CHILD VERSION of tr_late but it will work for periods, too. It takes
!> matrix of ss variables (at which the emax's are calculated) and converts 
!> it to the tmss matrix for regression. One child means no A2 in the state space so interpolating function is much smaller
subroutine troc_late(tmss, mss,n)
	implicit none
	integer, intent(in):: n 				!< the number of rows in the matrix to be transformed
	real(dble), intent(in) :: mss(n,Ssize) 	!< matrix holding a subset of ss in its rows
	real(dble), intent(out) :: tmss(n,Gsizeoc+1) !<regression matrix 

	tmss(:,1)=mss(:,1) 				! A1
	tmss(:,2)=mss(:,3) 				! E
	tmss(:,3:6)=mss(:,6:9) 		!< age of the children are constant for a period and delta, so skip those 
 	tmss(:,7)=mss(:,11) 		! also skip age0m, because agem is already in.
 	tmss(:,8)=mss(:,1)**2
 	tmss(:,9)=mss(:,3)**2
 	tmss(:,10)=mss(:,1)*mss(:,3)
	tmss(:,Gsizeoc+1)=1.0d0 				!< NEW: add the intercept
end subroutine troc_late


!> for fvlate to prepare the fss matrix for interpolation
!> does the same thing as tr_late but for fv calculations where SS vectors column vectors.
subroutine ftr_late(tfss,fss,n)
	implicit none
	integer,intent(in):: n
	real(dble), intent(out) :: tfss(Gsize+1,n)
	real(dble), intent(in) ::  fss(SSize,n)
	tfss(1:3,:)=fss(1:3,:)
	tfss(4:7,:)=fss(6:9,:) 		! age of the children are constant for a period and delta, so skip those 
 	tfss(8,:)=fss(11,:) 		!age0m is perfectly correlated with mother's current age for a given t and delta, so skip it.
 	tfss(9:11,:)=fss(1:3,:)**2
 	tfss(12,:)=fss(1,:)*fss(2,:)
 	tfss(13,:)=fss(1,:)*fss(3,:)
 	tfss(14,:)=fss(2,:)*fss(3,:)
	! NEW: add a ones vector to include the intercept
	tfss(Gsize+1,:)=1.0d0	
end subroutine ftr_late

!> sanme ftr_late, for one child families. interpolating vector does not have A2 in it.
subroutine ftroc_late(tfss,fss,n)
	implicit none
	integer,intent(in):: n
	real(dble), intent(out) :: tfss(Gsizeoc+1,n)
	real(dble), intent(in) ::  fss(SSize,n)
	tfss(1,:)=fss(1,:)
	tfss(2,:)=fss(3,:)
	tfss(3:6,:)=fss(6:9,:) 		!< age of the children are constant for a period and delta, so skip those 
 	tfss(7,:)=fss(11,:)
 	tfss(8,:)=fss(1,:)**2
 	tfss(9,:)=fss(3,:)**2
 	tfss(10,:)=fss(1,:)*fss(3,:)
	! NEW: add a ones vector to include the intercept
	tfss(Gsizeoc+1,:)=1.0d0	
end subroutine ftroc_late

!> for the fvhc to convert the big matrix of fss values to interpolating matrix
!> the main difference between ftr_late and this one is the fact that fss matrix is
!> much more complicated and bigger in the hc periods, because of the production of A.
!> ftr_hc is essentially does the same thing but in one more dimension
subroutine ftr_hc(tfss,fss,n,m)
	integer,intent(in):: n,m
	real(dble), intent(in) :: fss(SSize,n,m)
	real(dble), intent(out) :: tfss(GSize+1,n,m)
	tfss(1:3,:,:)=fss(1:3,:,:)
	tfss(4:7,:,:)=fss(6:9,:,:) 		!< age of the children are constant for a period and delta, so skip those 
 	tfss(8,:,:)=fss(11,:,:)
 	tfss(9:11,:,:)=fss(1:3,:,:)**2
 	tfss(12,:,:)=fss(1,:,:)*fss(2,:,:)
 	tfss(13,:,:)=fss(1,:,:)*fss(3,:,:)
 	tfss(14,:,:)=fss(2,:,:)*fss(3,:,:)
	! new: add ones for the intercept
	tfss(Gsize+1,:,:)=1.0d0
end subroutine ftr_hc

!> one child version of the ftr_hc: one less dimension due to missing A2
subroutine ftroc_hc(tfss,fss,n,m)
	integer,intent(in):: n,m
	real(dble), intent(in) :: fss(SSize,n,m)
	real(dble), intent(out) :: tfss(GSizeoc+1,n,m)
	tfss(1,:,:)=fss(1,:,:)
	tfss(2,:,:)=fss(3,:,:)
	tfss(3:6,:,:)=fss(6:9,:,:) 		!< age of the children are constant for a period and delta, so skip those 
 	tfss(7,:,:)=fss(11,:,:)
 	tfss(8,:,:)=fss(1,:,:)**2
 	tfss(9,:,:)=fss(3,:,:)**2
 	tfss(10,:,:)=fss(1,:,:)*fss(3,:,:)
	! new: add ones for the intercept
	tfss(Gsizeoc+1,:,:)=1.0d0
end subroutine ftroc_hc


end module

