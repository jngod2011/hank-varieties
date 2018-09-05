SUBROUTINE IterateTransitionStickyPi

USE Parameters
USE Globals
USE Procedures

IMPLICIT NONE

INTEGER 	:: it,ii,infix,insmooth
REAL(8) 	:: lldK,ldiffK,ldiffY,ldiffT,lstepB,lstepK,lstepT,lqa,lqb,lqc,lminmargcost,ltaxincr,ltaxincr_implied,ltaxincr1,lpvgovbc,lpvtaxableinc,lworldbond,lpiupdate
REAL(8), DIMENSION(Ttransition) :: lcapital,linflation,lfirmdiscount,lgovbond,lexcessdemand,lbondmarket
REAL(8), DIMENSION(Ttransition) :: lcapital1,lcapital_xold,lcapital_fold,lbond1,lbond_xold,lbond_fold

iteratingtransition = .true.

lminmargcost = 0.01

IF(fsptransition==.false.) THEN
	IF(Display>=1 .and. stickytransition==.true.) write(*,*)' Solving for sticky price transition without ZLB'
	IF(Display>=1 .and. zlbtransition==.true.) write(*,*)' Solving for sticky price transition with ZLB'
ELSE IF (fsptransition==.true.) THEN
	IF(Display>=1 .and. stickytransition==.true.) write(*,*)' Solving for sticky price transition without ZLB (with FSP)'
	IF(Display>=1 .and. zlbtransition==.true.) write(*,*)' Solving for sticky price transition with ZLB (with FSP)'
END IF


!guess capital demand and inflation
IF(fsptransition==.false. .and. UseFlexTransitionAsGuess==0) THEN

	!initial capital pins down labor and KN ratio in first period of transition
	equmTRANS(1)%capital = equmINITSS%capital
	!construct sequence of guesses of KN ratio: assume log linear in KN ratio (constant if a temporary transition)
	lldK = (log(equmFINALSS%capital)  - log(equmTRANS(1)%capital)) / real(Ttransition)
	DO it = 2,Ttransition
		equmTRANS(it)%capital = equmTRANS(1)%capital  *exp(lldK*it)
	END DO

	!inflation guess zero
	equmTRANS(:)%pi = 0.0

ELSE IF(fsptransition==.false. .and. UseFlexTransitionAsGuess==1 ) THEN
	equmTRANS%capital = irfpointer%equmFLEX%capital
	equmTRANS%pi = irfpointer%equmFLEX%pi

ELSE IF	(fsptransition==.true.) THEN
	IF(stickytransition==.true.) THEN
		equmTRANS%capital = irfpointer%equmSTICKY%capital
		equmTRANS%pi = irfpointer%equmSTICKY%pi
	ELSE IF(zlbtransition==.true.) THEN
		equmTRANS%capital = irfpointer%equmZLB%capital
		equmTRANS%pi = irfpointer%equmZLB%pi
	END IF
END IF

!government expenditure path (in baseline these are implied by tax revenues)
IF	(fsptransition==.true.) THEN
	IF(stickytransition==.true.) equmTRANS%govexp = irfpointer%equmSTICKY%govexp		
	IF(zlbtransition==.true.) equmTRANS%govexp = irfpointer%equmZLB%govexp		
	DO it = 1,Ttransition
		IF(cumdeltatrans(it)>=fspointer%govstart .and. cumdeltatrans(it)<fspointer%govend) THEN
			equmTRANS(it)%govexp = equmTRANS(it)%govexp + fspointer%govamount
		END  IF
	END DO
ELSE IF	(fsptransition==.false.) THEN
	equmTRANS(:)%govexp = equmINITSS%govexp
END IF

!labor tax rate
IF(fsptransition==.false.) THEN 
	equmTRANS(:)%labtax = equmINITSS%labtax
	equmTRANS(:)%govbond = equmINITSS%govbond
ELSE IF(fsptransition==.true.) THEN
	
	!tax increase guess
	ltaxincr = 0.003
	equmTRANS(:)%labtax = equmINITSS%labtax
	DO it = 1,Ttransition
		IF(cumdeltatrans(it)>=fspointer%labtaxstart .and. cumdeltatrans(it)<fspointer%labtaxend) THEN
			equmTRANS(it)%labtax = equmINITSS%labtax +ltaxincr
		END  IF
	END DO
		
END IF

!real and nominal interest rates
equmTRANS(:)%rnom = max(equmINITSS%rnom + phitaylor*equmTRANS(:)%pi + equmTRANS(:)%mpshock, 0.0)	!taylor rule
equmTRANS(:)%rb = equmTRANS(:)%rnom - equmTRANS(:)%pi  !fisher equn
equmTRANS(:)%rborr = equmTRANS(:)%rb + equmTRANS(:)%borrwedge

!solve phillips curve backwards for marginal costs
IF (FirmDiscountRate==1) lfirmdiscount = equmTRANS(:)%rho
IF (FirmDiscountRate==2) lfirmdiscount = equmINITSS%rb

!marginal costs
!final period of transition
it = Ttransition
lqa = equmTRANS(it)%elast/theta

lqb = (equmFINALSS%pi-equmTRANS(it)%pi)/deltatransvec(it) - (equmTRANS(it)%pi/deltatransvec(it)) *(1.0-alpha)*frisch/(1.0+alpha*frisch) &
 		-lfirmdiscount(it)*equmFINALSS%pi - (equmTRANS(it)%elast-1.0)/theta &
		+ ((1.0+ frisch)/(1.0+alpha*frisch))*(equmFINALSS%tfp-equmTRANS(it)%tfp)/(equmTRANS(it)%tfp*deltatransvec(it)) &
		+ (alpha*(1.0+ frisch)/(1.0+alpha*frisch))*(equmFINALSS%capital-equmTRANS(it)%capital)/(equmTRANS(it)%capital*deltatransvec(it)) &
		- ((1.0-alpha)*(1.0+ frisch)/(1.0+alpha*frisch))*(equmFINALSS%labtax-equmTRANS(it)%labtax)/(equmTRANS(it)%labtax*deltatransvec(it))

lqc = ((1.0-alpha)*frisch/(1.0+alpha*frisch))*equmFINALSS%mc*equmTRANS(it)%pi/deltatransvec(it)

IF (lqb**2 - 4.0*lqa*lqc >=0.0)	THEN
	equmTRANS(it)%mc = (-lqb + sqrt(lqb**2 - 4.0*lqa*lqc)) / (2.0*lqa)
	equmTRANS(it)%mc = max(lminmargcost,equmTRANS(it)%mc)
ELSE
	equmTRANS(it)%mc = lminmargcost
END IF

!solve backwards
DO it = Ttransition-1,1,-1
	lqa = equmTRANS(it)%elast/theta
	lqb = (equmTRANS(it+1)%pi-equmTRANS(it)%pi)/deltatransvec(it) - (equmTRANS(it)%pi/deltatransvec(it)) *(1.0-alpha)*frisch/(1.0+alpha*frisch) &
	 		-lfirmdiscount(it)*equmTRANS(it+1)%pi - (equmTRANS(it)%elast-1.0)/theta &
			+ ((1.0+ frisch)/(1.0+alpha*frisch))*(equmTRANS(it+1)%tfp-equmTRANS(it)%tfp)/(equmTRANS(it)%tfp*deltatransvec(it)) &
			+ (alpha*(1.0+ frisch)/(1.0+alpha*frisch))*(equmTRANS(it+1)%capital-equmTRANS(it)%capital)/(equmTRANS(it)%capital*deltatransvec(it)) &
			- ((1.0-alpha)*(1.0+ frisch)/(1.0+alpha*frisch))*(equmTRANS(it+1)%labtax-equmTRANS(it)%labtax)/(equmTRANS(it)%labtax*deltatransvec(it))
			

	lqc = ((1.0-alpha)*frisch/(1.0+alpha*frisch))*equmTRANS(it+1)%mc*equmTRANS(it)%pi/deltatransvec(it)

	IF (lqb**2 - 4.0*lqa*lqc >=0.0)	THEN
		equmTRANS(it)%mc = (-lqb + sqrt(lqb**2 - 4.0*lqa*lqc)) / (2.0*lqa)
		equmTRANS(it)%mc = max(lminmargcost,equmTRANS(it)%mc)
	ELSE
		equmTRANS(it)%mc = lminmargcost
	END IF

END DO


equmTRANS(:)%gap = equmTRANS(:)%elast*equmTRANS(:)%mc / (equmTRANS(:)%elast-1.0) - 1.0
equmTRANS(:)%labor = ((1.0-equmTRANS(:)%labtax)*equmTRANS(:)%mc*(1.0-alpha)* equmTRANS(:)%tfp *(equmTRANS(:)%capital**alpha)/chi)**(frisch/(1.0 + frisch*alpha))
equmTRANS(:)%KNratio = equmTRANS(:)%capital/equmTRANS(:)%labor
equmTRANS(:)%KYratio = (equmTRANS(:)%KNratio**(1.0-alpha)) / equmTRANS(:)%tfp
equmTRANS(:)%output = equmTRANS(:)%tfp * (equmTRANS(:)%capital**alpha) * (equmTRANS(:)%labor**(1.0-alpha)) + rhousing*housefrac*(1.0-equmTRANS(:)%fundlev)*equmTRANS(:)%capital/(1.0-housefrac)	
equmTRANS(:)%rcapital = equmTRANS(:)%mc * alpha / equmTRANS(:)%KYratio
equmTRANS(:)%wage = equmTRANS(:)%mc*(1.0-alpha)* equmTRANS(:)%tfp * (equmTRANS(:)%KNratio**alpha)
equmTRANS(:)%netwage = (1.0-equmTRANS(:)%labtax)*equmTRANS(:)%wage
equmTRANS(:)%taxrev = equmTRANS(:)%labtax*equmTRANS(:)%wage*equmTRANS(:)%labor - lumptransfer
IF(fsptransition==.false.) equmTRANS(:)%govexp = equmTRANS(:)%taxrev + equmTRANS(:)%rb*equmINITSS%govbond
equmTRANS(:)%priceadjust = (theta/2.0)*(equmTRANS(:)%pi**2)*equmTRANS(:)%output


!real bond return and borrowing rate
equmTRANS(:)%fundbond = -equmTRANS(:)%capital*equmTRANS(:)%fundlev
equmTRANS(1)%worldbond = equmINITSS%worldbond
DO it = 1,Ttransition-1
	CALL WorldBondFunction2(equmTRANS(it)%rb,lworldbond,equmINITSS%worldbond,equmINITSS%rb,bondelast)
	equmTRANS(it+1)%worldbond = equmTRANS(it)%worldbond + bondadjust*deltatransvec(it)*(lworldbond-equmTRANS(it)%worldbond)
END DO

equmTRANS(:)%bond = -equmTRANS(:)%worldbond - equmTRANS(:)%govbond - equmTRANS(:)%fundbond


!solve backward for investment
it = Ttransition
equmTRANS(it)%investment = (equmFINALSS%capital-equmTRANS(it)%capital)/deltatransvec(it) + deprec*equmTRANS(it)%capital
DO it = Ttransition-1,1,-1
	equmTRANS(it)%investment = (equmTRANS(it+1)%capital-equmTRANS(it)%capital)/deltatransvec(it) + deprec*equmTRANS(it)%capital
END DO

equmTRANS(:)%ra = (equmTRANS(:)%rcapital - deprec + (1.0-equmTRANS(:)%mc)/equmTRANS(:)%KYratio - equmTRANS(:)%priceadjust/equmTRANS(:)%capital - equmTRANS(:)%fundlev*equmTRANS(:)%rb) / (1.0-equmTRANS(:)%fundlev)

ii = 1
lstepK = 0.1
lpiupdate = 0.01
	 
ldiffK = 1.0
ldiffY = 1.0
DO WHILE (ii<=maxitertranssticky .and. max(ldiffK,ldiffY)>toltransition )
	!solve for transtion
	CALL Transition
	
! 	!compute implied govt bond demand
! 	IF(fsptransition==.true.) THEN 
! 		lgovbond(1) = equmINITSS%govbond
! 		DO it = 1,Ttransition-1
! 			lgovbond(it+1) = (1.0 + deltatransvec(it)*equmTRANS(it)%rb)*lgovbond(it) &
! 							+ deltatransvec(it)*equmTRANS(it)%labtax*equmTRANS(it)%wage*equmTRANS(it)%labor &
! 							- deltatransvec(it)*equmTRANS(it)%govexp &
! 							- deltatransvec(it)*statsTRANS(it)%Efsp
! 		END DO
! 		equmTRANS(:)%govbond = lgovbond
! 	ELSE IF(fsptransition==.false.) THEN
! 		lgovbond  = equmTRANS(:)%govbond
! 	END IF
! 	
! 	!compute labor tax increase that would balance the budget
! 	IF(fsptransition==.true.) THEN 
! 		lpvgovbc = 0.0
! 		lpvtaxableinc = 0.0
! 		DO it = Ttransition,1,-1
! 			lpvgovbc = exp(-deltatransvec(it)*equmTRANS(it)%rb) *(lpvgovbc + deltatransvec(it)*(equmTRANS(it)%labtax*equmTRANS(it)%wage*equmTRANS(it)%labor - equmTRANS(it)%govexp - statsTRANS(it)%Efsp))
! 			IF(cumdeltatrans(it)>=fspointer%labtaxstart .and. cumdeltatrans(it)>=fspointer%labtaxend) THEN
! 				lpvtaxableinc = exp(-deltatransvec(it)*equmTRANS(it)%rb) *(lpvtaxableinc + deltatransvec(it)*equmTRANS(it)%wage*equmTRANS(it)%labor)
! 			ELSE
! 				lpvtaxableinc = exp(-deltatransvec(it)*equmTRANS(it)%rb) *lpvtaxableinc
! 			END  IF
! 			
! 		END DO
! 		ltaxincr_implied  = -lpvgovbc/lpvtaxableinc + ltaxincr
! 	END IF
	
	!computed implied equilibrium quantities	
	lcapital = (1.0-housefrac)*statsTRANS(:)%Ea / (1.0 - equmTRANS(:)%fundlev)
	ldiffK = maxval(abs(lcapital/equmTRANS(:)%capital - 1.0))
	
	lbondmarket = statsTRANS(:)%Eb + equmTRANS(:)%worldbond + equmTRANS(:)%fundbond + equmTRANS(:)%govbond
	DO it = 1,Ttransition-1
		lexcessdemand(it) = equmTRANS(it)%rb*lbondmarket(it) - (lbondmarket(it+1)-lbondmarket(it))/deltatransvec(it)
	END DO
	it = Ttransition
	lexcessdemand(it) = equmTRANS(it)%rb*lbondmarket(it) + lbondmarket(it)/deltatransvec(it)
	
	ldiffY = maxval(abs(lexcessdemand/equmTRANS(:)%output - 1.0))

	IF(fsptransition==.true.) ldiffT = lpvgovbc
	IF (Display>=1) write(*,"(A,I,A)") '  Transition iter ',ii, ':'
	IF (Display>=1 .and. fsptransition==.false.) write(*,"(A,E10.3,A,E10.3,A,E10.3)") '   K err',ldiffK, ',  Y err',ldiffY

	!update capital, labor and bond supply
	IF (ii<=maxitertranssticky .and. max(ldiffK,ldiffY)>toltransition ) THEN		
		!updae capital
		CALL PartialUpdate(Ttransition-1,lstepK,equmTRANS(2:Ttransition)%capital,lcapital(2:Ttransition),lcapital1(2:Ttransition))
		lcapital1(1)  = equmTRANS(1)%capital
		lcapital_xold = equmTRANS(:)%capital
		lcapital_fold = lcapital
		equmTRANS(:)%capital = lcapital1

		!update inflation
		equmTRANS(:)%pi = equmTRANS(:)%pi + lpiupdate*lexcessdemand
		
		
		!impose it settles down
		IF(ImposeGuessSettlesDown==1) THEN
			infix = 1
			equmTRANS(Ttransition-infix+1:Ttransition)%pi = equmFINALSS%pi
			equmTRANS(Ttransition-infix+1:Ttransition)%capital = equmFINALSS%capital
			insmooth = Ttransition-10
 			CALL MASmooth(insmooth-infix+1,equmTRANS(Ttransition-insmooth:Ttransition-infix)%pi,deltatransvec(Ttransition-insmooth:Ttransition-infix),3,equmTRANS(Ttransition-insmooth:Ttransition-infix)%pi)
 			CALL MASmooth(insmooth-infix+1,equmTRANS(Ttransition-insmooth:Ttransition-infix)%capital,deltatransvec(Ttransition-insmooth:Ttransition-infix),3,equmTRANS(Ttransition-insmooth:Ttransition-infix)%capital)
			
		END IF
		
! 		!update labor tax rate
! 		IF(fsptransition==.true.) THEN
! 			ltaxincr1 = ltaxincr
! 			ltaxincr = lstepT*ltaxincr_implied + (1.0-lstepT)*ltaxincr1
! 			write(*,*) 'ltaxincr_implied',ltaxincr_implied,'ltaxincr',ltaxincr
! 			equmTRANS(:)%labtax = equmINITSS%labtax
! 			DO it = 1,Ttransition
! 				IF(cumdeltatrans(it)>=fspointer%labtaxstart .and. cumdeltatrans(it)<fspointer%labtaxend) THEN
! 					equmTRANS(it)%labtax = equmINITSS%labtax +ltaxincr
! 					equmTRANS(it)%labtax = min(max(equmTRANS(it)%labtax,0.0_8),1.0/(1.0+frisch))
! 				END  IF
! 			END DO
! 		END IF
	ElSE
		!run distribution stats with full
		iteratingtransition = .false.
		CALL Transition
		equmTRANS(:)%pi = equmTRANS(:)%pi
		equmTRANS(:)%capital = lcapital
		
	END IF
	

	!real and nominal interest rates
	equmTRANS(:)%rnom = max(equmINITSS%rnom + phitaylor*equmTRANS(:)%pi + equmTRANS(:)%mpshock, 0.0)	!taylor rule
	equmTRANS(:)%rb = equmTRANS(:)%rnom - equmTRANS(:)%pi  !fisher equn
	equmTRANS(:)%rborr = equmTRANS(:)%rb + equmTRANS(:)%borrwedge

	!solve phillips curve backwards for marginal costs
	IF (FirmDiscountRate==1) lfirmdiscount = equmTRANS(:)%rho
	IF (FirmDiscountRate==2) lfirmdiscount = equmINITSS%rb

	!marginal costs
	!final period of transition
	it = Ttransition
	lqa = equmTRANS(it)%elast/theta

	lqb = (equmFINALSS%pi-equmTRANS(it)%pi)/deltatransvec(it) - (equmTRANS(it)%pi/deltatransvec(it)) *(1.0-alpha)*frisch/(1.0+alpha*frisch) &
	 		-lfirmdiscount(it)*equmFINALSS%pi - (equmTRANS(it)%elast-1.0)/theta &
			+ ((1.0+ frisch)/(1.0+alpha*frisch))*(equmFINALSS%tfp-equmTRANS(it)%tfp)/(equmTRANS(it)%tfp*deltatransvec(it)) &
			+ (alpha*(1.0+ frisch)/(1.0+alpha*frisch))*(equmFINALSS%capital-equmTRANS(it)%capital)/(equmTRANS(it)%capital*deltatransvec(it)) &
			- ((1.0-alpha)*(1.0+ frisch)/(1.0+alpha*frisch))*(equmFINALSS%labtax-equmTRANS(it)%labtax)/(equmTRANS(it)%labtax*deltatransvec(it))

	lqc = ((1.0-alpha)*frisch/(1.0+alpha*frisch))*equmFINALSS%mc*equmTRANS(it)%pi/deltatransvec(it)

	IF (lqb**2 - 4.0*lqa*lqc >=0.0)	THEN
		equmTRANS(it)%mc = (-lqb + sqrt(lqb**2 - 4.0*lqa*lqc)) / (2.0*lqa)
		equmTRANS(it)%mc = max(lminmargcost,equmTRANS(it)%mc)
	ELSE
		equmTRANS(it)%mc = lminmargcost
	END IF

	!solve backwards
	DO it = Ttransition-1,1,-1
		lqa = equmTRANS(it)%elast/theta
		lqb = (equmTRANS(it+1)%pi-equmTRANS(it)%pi)/deltatransvec(it) - (equmTRANS(it)%pi/deltatransvec(it)) *(1.0-alpha)*frisch/(1.0+alpha*frisch) &
		 		-lfirmdiscount(it)*equmTRANS(it+1)%pi - (equmTRANS(it)%elast-1.0)/theta &
				+ ((1.0+ frisch)/(1.0+alpha*frisch))*(equmTRANS(it+1)%tfp-equmTRANS(it)%tfp)/(equmTRANS(it)%tfp*deltatransvec(it)) &
				+ (alpha*(1.0+ frisch)/(1.0+alpha*frisch))*(equmTRANS(it+1)%capital-equmTRANS(it)%capital)/(equmTRANS(it)%capital*deltatransvec(it)) &
				- ((1.0-alpha)*(1.0+ frisch)/(1.0+alpha*frisch))*(equmTRANS(it+1)%labtax-equmTRANS(it)%labtax)/(equmTRANS(it)%labtax*deltatransvec(it))


		lqc = ((1.0-alpha)*frisch/(1.0+alpha*frisch))*equmTRANS(it+1)%mc*equmTRANS(it)%pi/deltatransvec(it)

		IF (lqb**2 - 4.0*lqa*lqc >=0.0)	THEN
			equmTRANS(it)%mc = (-lqb + sqrt(lqb**2 - 4.0*lqa*lqc)) / (2.0*lqa)
			equmTRANS(it)%mc = max(lminmargcost,equmTRANS(it)%mc)
		ELSE
			equmTRANS(it)%mc = lminmargcost
		END IF

	END DO


	equmTRANS(:)%gap = equmTRANS(:)%elast*equmTRANS(:)%mc / (equmTRANS(:)%elast-1.0) - 1.0
	equmTRANS(:)%labor = ((1.0-equmTRANS(:)%labtax)*equmTRANS(:)%mc*(1.0-alpha)* equmTRANS(:)%tfp *(equmTRANS(:)%capital**alpha)/chi)**(frisch/(1.0 + frisch*alpha))
	equmTRANS(:)%KNratio = equmTRANS(:)%capital/equmTRANS(:)%labor
	equmTRANS(:)%KYratio = (equmTRANS(:)%KNratio**(1.0-alpha)) / equmTRANS(:)%tfp
	equmTRANS(:)%output = equmTRANS(:)%tfp * (equmTRANS(:)%capital**alpha) * (equmTRANS(:)%labor**(1.0-alpha)) + rhousing*housefrac*(1.0-equmTRANS(:)%fundlev)*equmTRANS(:)%capital/(1.0-housefrac)	
	equmTRANS(:)%rcapital = equmTRANS(:)%mc * alpha / equmTRANS(:)%KYratio
	equmTRANS(:)%wage = equmTRANS(:)%mc*(1.0-alpha)* equmTRANS(:)%tfp * (equmTRANS(:)%KNratio**alpha)
	equmTRANS(:)%netwage = (1.0-equmTRANS(:)%labtax)*equmTRANS(:)%wage
	equmTRANS(:)%taxrev = equmTRANS(:)%labtax*equmTRANS(:)%wage*equmTRANS(:)%labor - lumptransfer
	IF(fsptransition==.false.) equmTRANS(:)%govexp = equmTRANS(:)%taxrev + equmTRANS(:)%rb*equmINITSS%govbond
	equmTRANS(:)%priceadjust = (theta/2.0)*(equmTRANS(:)%pi**2)*equmTRANS(:)%output


	!real bond return and borrowing rate
	equmTRANS(:)%fundbond = -equmTRANS(:)%capital*equmTRANS(:)%fundlev
	equmTRANS(1)%worldbond = equmINITSS%worldbond
	DO it = 1,Ttransition-1
		CALL WorldBondFunction2(equmTRANS(it)%rb,lworldbond,equmINITSS%worldbond,equmINITSS%rb,bondelast)
		equmTRANS(it+1)%worldbond = equmTRANS(it)%worldbond + bondadjust*deltatransvec(it)*(lworldbond-equmTRANS(it)%worldbond)
	END DO

	equmTRANS(:)%bond = -equmTRANS(:)%worldbond - equmTRANS(:)%govbond - equmTRANS(:)%fundbond

	!solve backward for investment
	it = Ttransition
	equmTRANS(it)%investment = (equmFINALSS%capital-equmTRANS(it)%capital)/deltatransvec(it) + deprec*equmTRANS(it)%capital
	DO it = Ttransition-1,1,-1
		equmTRANS(it)%investment = (equmTRANS(it+1)%capital-equmTRANS(it)%capital)/deltatransvec(it) + deprec*equmTRANS(it)%capital
	END DO

	equmTRANS(:)%ra = (equmTRANS(:)%rcapital - deprec + (1.0-equmTRANS(:)%mc)/equmTRANS(:)%KYratio - equmTRANS(:)%priceadjust/equmTRANS(:)%capital - equmTRANS(:)%fundlev*equmTRANS(:)%rb) / (1.0-equmTRANS(:)%fundlev)

	
	ii = ii+1	
END DO

IF(fsptransition==.false.) THEN 
	IF(stickytransition==.true.) THEN
		irfpointer%equmSTICKY = equmTRANS
		irfpointer%statsSTICKY = statsTRANS
		irfpointer%solnSTICKY = solnTRANS
	ELSE IF	(zlbtransition==.true.) THEN
		irfpointer%equmZLB = equmTRANS
		irfpointer%statsZLB = statsTRANS
		irfpointer%solnZLB = solnTRANS
	END IF
ELSE IF(fsptransition==.true.) THEN 
	IF(stickytransition==.true.) THEN
		irfpointer_fs%equmSTICKY = equmTRANS
		irfpointer_fs%statsSTICKY = statsTRANS
		irfpointer_fs%solnSTICKY = solnTRANS
	ELSE IF	(zlbtransition==.true.) THEN
		irfpointer_fs%equmZLB = equmTRANS
		irfpointer_fs%statsZLB = statsTRANS
		irfpointer_fs%solnZLB = solnTRANS
	END IF
END IF

END SUBROUTINE IterateTransitionStickyPi