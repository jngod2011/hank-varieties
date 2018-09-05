SUBROUTINE IterateTransitionStickyB

USE Parameters
USE Globals
USE Procedures

IMPLICIT NONE

INTEGER 	:: it,ii,infix,insmooth
REAL(8) 	:: lldB,lldK,ldiffB,ldiffK,ldiffT,lstepT,lqa,lqb,lqc,lminmargcost,ltaxincr,ltaxincr_implied,ltaxincr1,lpvgovbc,lpvtaxableinc,lpvprofit,lpvfirmrb
REAL(8), DIMENSION(Ttransition) :: lcapital,lbond,lfirmdiscount,lgovbond,lrb,lrb1
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


!guess capital and bond demands
IF(fsptransition==.false. .and. UseFlexTransitionAsGuess==0) THEN

	!initial capital pins down labor and KN ratio in first period of transition
	equmTRANS(1)%capital = equmINITSS%capital
	!construct sequence of guesses of KN ratio: assume log linear in KN ratio (constant if a temporary transition)
	lldK = (log(equmFINALSS%capital)  - log(equmTRANS(1)%capital)) / real(Ttransition)
	DO it = 2,Ttransition
		equmTRANS(it)%capital = equmTRANS(1)%capital  *exp(lldK*it)
	END DO

	!initial bond are given by steady state
	equmTRANS(1)%bond = equmINITSS%bond
	!construct sequence of guesses of B: assume log linear in (B-borrlim) (constant if a temporary transition)
	lldB = (log(equmFINALSS%bond)  - log(equmTRANS(1)%bond)) / real(Ttransition)
	DO it = 2,Ttransition
		equmTRANS(it)%bond = equmTRANS(1)%bond  *exp(lldB*it)
	END DO

ELSE IF(fsptransition==.false. .and. UseFlexTransitionAsGuess==1 ) THEN
	equmTRANS%capital = irfpointer%equmFLEX%capital
	equmTRANS%bond = irfpointer%equmFLEX%bond

ELSE IF	(fsptransition==.true.) THEN
	IF(stickytransition==.true.) THEN
		equmTRANS%capital = irfpointer%equmSTICKY%capital
		equmTRANS%bond = irfpointer%equmSTICKY%bond
	ELSE IF(zlbtransition==.true.) THEN
		equmTRANS%capital = irfpointer%equmZLB%capital
		equmTRANS%bond = irfpointer%equmZLB%bond
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

!real bond return and borrowing rate
equmTRANS(:)%fundbond = -equmTRANS(:)%capital*equmTRANS(:)%fundlev
equmTRANS(:)%worldbond = -equmTRANS(:)%bond - equmTRANS(:)%govbond - equmTRANS(:)%fundbond

it = Ttransition
CALL WorldBondInverse2( (equmFINALSS%worldbond-equmTRANS(it)%worldbond)/(bondadjust*deltatransvec(it)) + equmTRANS(it)%worldbond &
					,equmTRANS(it)%rb,equmINITSS%worldbond,equmINITSS%rb,bondelast)
DO it = Ttransition-1,1,-1
	CALL WorldBondInverse2( (equmTRANS(it+1)%worldbond-equmTRANS(it)%worldbond)/(bondadjust*deltatransvec(it)) + equmTRANS(it)%worldbond &
						,equmTRANS(it)%rb,equmINITSS%worldbond,equmINITSS%rb,bondelast)
END DO

equmTRANS(:)%rborr = equmTRANS(:)%rb + equmTRANS(:)%borrwedge

!inflation and nominal interest rates
IF(zlbtransition==.false.) THEN
	IF(BackwardTermInTaylorRule==0) THEN
		equmTRANS(:)%pi = (equmTRANS(:)%rb - equmINITSS%rnom - equmTRANS(:)%mpshock) / (phitaylor-1.0) !taylor rule
		equmTRANS(:)%rnom = equmTRANS(:)%rb + equmTRANS(:)%pi !fisher equn
	ELSE IF(BackwardTermInTaylorRule==1) THEN
		equmTRANS(1)%pi = (equmTRANS(1)%rb - equmINITSS%rnom - taylorpers*deltatransvec(1)*equmTRANS(1)%mpshock) / (phitaylor*taylorpers*deltatransvec(1)-1.0)
		equmTRANS(1)%rnom = equmTRANS(1)%rb + equmTRANS(1)%pi !fisher equn
		DO it = 2,Ttransition
			equmTRANS(it)%pi = (equmTRANS(it)%rb - taylorpers*deltatransvec(it)*equmINITSS%rnom &
								- (1.0-taylorpers*deltatransvec(it))*equmTRANS(it-1)%rnom  &
								- taylorpers*deltatransvec(it)*equmTRANS(it)%mpshock) / (phitaylor*taylorpers*deltatransvec(it)-1.0) !taylor rule
			equmTRANS(it)%rnom = equmTRANS(it)%rb + equmTRANS(it)%pi !fisher equn
		END DO
	END IF
	
ELSE IF(zlbtransition==.true.) THEN
	DO it = 1,Ttransition
		IF(equmTRANS(it)%rb > (equmINITSS%rnom + equmTRANS(it)%mpshock) / phitaylor) THEN !ZLB does not bind
			equmTRANS(it)%pi = (equmTRANS(it)%rb - equmINITSS%rnom - equmTRANS(it)%mpshock) / (phitaylor-1.0) !taylor rule
			equmTRANS(it)%rnom = equmTRANS(it)%rb + equmTRANS(it)%pi !fisher equn		
		ELSE IF(equmTRANS(it)%rb <= (equmINITSS%rnom + equmTRANS(it)%mpshock)/phitaylor) THEN !ZLB binds
			equmTRANS(it)%pi = -equmTRANS(it)%rb
			equmTRANS(it)%rnom = 0.0
		END IF
	END DO
END IF

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
equmTRANS(:)%priceadjust = (theta/2.0)*(equmTRANS(:)%pi**2)*equmTRANS(:)%capital/equmTRANS(:)%KYratio
equmTRANS(:)%profit = (1.0-equmTRANS(:)%mc)*equmTRANS(:)%capital/equmTRANS(:)%KYratio - equmTRANS(:)%priceadjust


!solve backward for investment
it = Ttransition
equmTRANS(it)%investment = (equmFINALSS%capital-equmTRANS(it)%capital)/deltatransvec(it) + deprec*equmTRANS(it)%capital
DO it = Ttransition-1,1,-1
	equmTRANS(it)%investment = (equmTRANS(it+1)%capital-equmTRANS(it)%capital)/deltatransvec(it) + deprec*equmTRANS(it)%capital
END DO

!dividends and illiquid return
IF(DividendSmoothing ==0) THEN
	IF(FixProfitsOutOfSteadyState==0) equmTRANS(:)%dividend 	= equmTRANS(:)%profit
	IF(FixProfitsOutOfSteadyState==1) equmTRANS(:)%dividend 	= equmINITSS%profit
	
	equmTRANS(:)%divrate 	= equmTRANS(:)%dividend/equmTRANS(:)%capital 
	equmTRANS(:)%ra = (equmTRANS(:)%rcapital - deprec + equmTRANS(:)%divrate - equmTRANS(:)%fundlev*equmTRANS(:)%rb) / (1.0-equmTRANS(:)%fundlev)
	equmTRANS(:)%intfirmbond = 0.0

ELSE IF(DividendSmoothing ==1) THEN
	lpvprofit = equmFINALSS%profit/equmFINALSS%rb
	lpvfirmrb = 1.0/equmFINALSS%rb
	DO it = Ttransition,1
		lpvprofit = exp(-equmTRANS(it)%rb*deltatransvec(it))*(lpvprofit + equmTRANS(it)%profit*deltatransvec(it))
		lpvfirmrb = exp((-lfirmdiscount(it)/firmgamma + (1.0/firmgamma-1.0)*equmTRANS(it)%rb)*deltatransvec(it)) * lpvfirmrb
	END DO
	it = 1
! 	equmTRANS(it)%dividend = (lpvprofit/lpvfirmrb) / (1.0 + (lfirmdiscount(it) - equmTRANS(it)%rb)*deltatransvec(it))
	equmTRANS(it)%dividend = (lpvprofit/lpvfirmrb) * (1.0 - (lfirmdiscount(it) - equmTRANS(it)%rb)*deltatransvec(it)/firmgamma)
	DO it = 2,Ttransition
! 		equmTRANS(it)%dividend = equmTRANS(it-1)%dividend / (1.0 + (lfirmdiscount(it) - equmTRANS(it)%rb)*deltatransvec(it))
		equmTRANS(it)%dividend = equmTRANS(it-1)%dividend * (1.0 - (lfirmdiscount(it) - equmTRANS(it)%rb)*deltatransvec(it)/firmgamma)
	END DO
	equmTRANS(:)%divrate 	= equmTRANS(:)%dividend/equmTRANS(:)%capital 
	equmTRANS(:)%ra = (equmTRANS(:)%rcapital - deprec + equmTRANS(:)%divrate - equmTRANS(:)%fundlev*equmTRANS(:)%rb) / (1.0-equmTRANS(:)%fundlev)
	equmTRANS(1)%intfirmbond = 0.0
	DO it = 1, Ttransition-1
		equmTRANS(it+1)%intfirmbond = (equmTRANS(it)%intfirmbond + deltatransvec(it)*(equmTRANS(it)%profit-equmTRANS(it)%dividend)) / (1.0 - equmTRANS(it)%rb*deltatransvec(it))
	END DO

ELSE IF(DividendSmoothing ==2) THEN

	equmTRANS(1)%dividend =  (1.0-divsmooth*cumdeltatrans(Ttransition)/(1.0 - exp(-divsmooth*cumdeltatrans(Ttransition)))) * equmINITSS%profit &
	 						+ (divsmooth/(1.0 - exp(-divsmooth*cumdeltatrans(Ttransition)))) * SUM(equmTRANS(:)%profit*deltatransvec)
	DO it = 1,Ttransition-1
		equmTRANS(it+1)%dividend = equmINITSS%dividend + (equmTRANS(it)%dividend - equmINITSS%dividend) / (1.0 + divsmooth*deltatransvec(it))
	END DO
	equmTRANS(:)%divrate 	= equmTRANS(:)%dividend/equmTRANS(:)%capital 
	equmTRANS(:)%ra = (equmTRANS(:)%rcapital - deprec + equmTRANS(:)%divrate - equmTRANS(:)%fundlev*equmTRANS(:)%rb) / (1.0-equmTRANS(:)%fundlev)
	equmTRANS(1)%intfirmbond = 0.0
	DO it = 1, Ttransition-1
		equmTRANS(it+1)%intfirmbond = (equmTRANS(it)%intfirmbond + deltatransvec(it)*(equmTRANS(it)%profit-equmTRANS(it)%dividend))
	END DO
	
END IF
! write(*,*) ' lpvprofit is ',lpvprofit
! DO it = 1,Ttransition
! 	write(*,*) it, ' mc is ',equmTRANS(it)%mc
! 	write(*,*) it, ' price adjust 1 is ',(theta/2.0)*(equmTRANS(it)%pi**2)*equmTRANS(it)%output
! 	write(*,*) it, ' price adjust 2 is ',(theta/2.0)*(equmTRANS(it)%pi**2)*equmTRANS(it)%capital/equmTRANS(it)%KYratio
! 	write(*,*) it, ' profit is ',equmTRANS(it)%profit
! 	write(*,*) it, ' dividend is ',equmTRANS(it)%dividend
! 	write(*,*) it, ' divrate is ',equmTRANS(it)%divrate
! 	write(*,*) it, ' ra is ',equmTRANS(it)%ra
! 	write(*,*) it, ' intfirmbond is ',equmTRANS(it)%intfirmbond
! END DO


ii = 1	 
ldiffK = 1.0
ldiffB = 1.0
DO WHILE (ii<=maxitertranssticky .and. max(ldiffK,ldiffB)>toltransition )
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
	lbond = statsTRANS(:)%Eb
	lcapital = (1.0-housefrac)*statsTRANS(:)%Ea / (1.0 - equmTRANS(:)%fundlev)
	ldiffK= maxval(abs(lcapital/equmTRANS(:)%capital - 1.0))
	ldiffB= maxval(abs(lbond/equmTRANS(:)%bond - 1.0))
	IF(fsptransition==.true.) ldiffT = lpvgovbc
	IF (Display>=1) write(*,"(A,I,A)") '  Transition iter ',ii, ':'
	IF (Display>=1 .and. fsptransition==.false.) write(*,"(A,E10.3,A,E10.3,A,E10.3)") '   K err',ldiffK, ',  B err',ldiffB

	!update capital, labor and bond supply
	IF (ii<maxitertranssticky .and. max(ldiffK,ldiffB)>toltransition ) THEN
		CALL PartialUpdate(Ttransition-1,stepstickytransK,equmTRANS(2:Ttransition)%capital,lcapital(2:Ttransition),lcapital1(2:Ttransition))
		lcapital1(1)  = equmTRANS(1)%capital
		lcapital_xold = equmTRANS(:)%capital
		lcapital_fold = lcapital
		equmTRANS(:)%capital = lcapital1
		
		IF(UpdateUsingR==0) THEN
			CALL PartialUpdate(Ttransition-1,stepstickytransB,equmTRANS(2:Ttransition)%bond,lbond(2:Ttransition),lbond1(2:Ttransition))
			lbond1(1)  = equmTRANS(1)%bond
			lbond_xold = equmTRANS(:)%bond
			lbond_fold = lbond
			equmTRANS(:)%bond = lbond1
		ELSEIF(UpdateUsingR==1) THEN
			equmTRANS(:)%fundbond = -lcapital*equmTRANS(:)%fundlev
			equmTRANS(:)%worldbond = -lbond - equmTRANS(:)%govbond - equmTRANS(:)%fundbond - equmTRANS(:)%intfirmbond

			it = Ttransition
			CALL WorldBondInverse2( (equmFINALSS%worldbond-equmTRANS(it)%worldbond)/(bondadjust*deltatransvec(it)) + equmTRANS(it)%worldbond &
								,lrb(it),equmINITSS%worldbond,equmINITSS%rb,bondelast)
			DO it = Ttransition-1,1,-1
				CALL WorldBondInverse2( (equmTRANS(it+1)%worldbond-equmTRANS(it)%worldbond)/(bondadjust*deltatransvec(it)) + equmTRANS(it)%worldbond &
									,lrb(it),equmINITSS%worldbond,equmINITSS%rb,bondelast)
			END DO

			CALL PartialUpdate(Ttransition,stepstickytransB,equmTRANS%rb,lrb,lrb1)
			equmTRANS(:)%rb = lrb1

		END IF

		!impose it settles down
		IF(ImposeGuessSettlesDown==1) THEN
			infix = 1
			insmooth = Ttransition-10

			equmTRANS(Ttransition-infix+1:Ttransition)%capital = equmFINALSS%capital
 			CALL MASmooth(insmooth-infix+1,equmTRANS(Ttransition-insmooth:Ttransition-infix)%capital,deltatransvec(Ttransition-insmooth:Ttransition-infix),3,equmTRANS(Ttransition-insmooth:Ttransition-infix)%capital)

			IF(UpdateUsingR==0) THEN
				equmTRANS(Ttransition-infix+1:Ttransition)%bond = equmFINALSS%bond
	 			CALL MASmooth(insmooth-infix+1,equmTRANS(Ttransition-insmooth:Ttransition-infix)%bond,deltatransvec(Ttransition-insmooth:Ttransition-infix),3,equmTRANS(Ttransition-insmooth:Ttransition-infix)%bond)
			ELSEIF(UpdateUsingR==1) THEN
				equmTRANS(Ttransition-infix+1:Ttransition)%rb = equmFINALSS%rb
	 			CALL MASmooth(insmooth-infix+1,equmTRANS(Ttransition-insmooth:Ttransition-infix)%rb,deltatransvec(Ttransition-insmooth:Ttransition-infix),3,equmTRANS(Ttransition-insmooth:Ttransition-infix)%rb)
			END IF
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
		equmTRANS(:)%capital = lcapital
		equmTRANS(:)%bond = lbond
		equmTRANS(:)%rb = lrb
		
	END IF
	
	!real bond return and borrowing rate
	IF(UpdateUsingR==0) THEN
		equmTRANS(:)%fundbond = -equmTRANS(:)%capital*equmTRANS(:)%fundlev
		equmTRANS(:)%worldbond = -equmTRANS(:)%bond - equmTRANS(:)%govbond - equmTRANS(:)%fundbond - equmTRANS(:)%intfirmbond

		it = Ttransition
		CALL WorldBondInverse2( (equmFINALSS%worldbond-equmTRANS(it)%worldbond)/(bondadjust*deltatransvec(it)) + equmTRANS(it)%worldbond &
							,equmTRANS(it)%rb,equmINITSS%worldbond,equmINITSS%rb,bondelast)
		DO it = Ttransition-1,1,-1
			CALL WorldBondInverse2( (equmTRANS(it+1)%worldbond-equmTRANS(it)%worldbond)/(bondadjust*deltatransvec(it)) + equmTRANS(it)%worldbond &
								,equmTRANS(it)%rb,equmINITSS%worldbond,equmINITSS%rb,bondelast)
		END DO
	ELSEIF(UpdateUsingR==1) THEN

		equmTRANS(1)%worldbond = equmINITSS%worldbond
		DO it = 1,Ttransition-1
			CALL WorldBondFunction2( equmTRANS(it)%rb,equmTRANS(it+1)%worldbond,equmINITSS%worldbond,equmINITSS%rb,bondelast)
			equmTRANS(it+1)%worldbond = equmTRANS(it)%worldbond + bondadjust*deltatransvec(it)*(equmTRANS(it+1)%worldbond-equmTRANS(it)%worldbond)
		END DO

		equmTRANS(:)%fundbond = -equmTRANS(:)%capital*equmTRANS(:)%fundlev
		equmTRANS(:)%bond = -equmTRANS(:)%worldbond - equmTRANS(:)%govbond - equmTRANS(:)%fundbond - equmTRANS(:)%intfirmbond

	END IF
	equmTRANS(:)%rborr = equmTRANS(:)%rb + equmTRANS(:)%borrwedge

	!inflation and nominal interest rates
	IF(zlbtransition==.false.) THEN
		IF(BackwardTermInTaylorRule==0) THEN
			equmTRANS(:)%pi = (equmTRANS(:)%rb - equmINITSS%rnom - equmTRANS(:)%mpshock) / (phitaylor-1.0) !taylor rule
			equmTRANS(:)%rnom = equmTRANS(:)%rb + equmTRANS(:)%pi !fisher equn			
		ELSE IF(BackwardTermInTaylorRule==1) THEN
			equmTRANS(1)%pi = (equmTRANS(1)%rb - equmINITSS%rnom - taylorpers*deltatransvec(1)*equmTRANS(1)%mpshock) / (phitaylor*taylorpers*deltatransvec(1)-1.0)
			equmTRANS(1)%rnom = equmTRANS(1)%rb + equmTRANS(1)%pi !fisher equn
			DO it = 2,Ttransition
				equmTRANS(it)%pi = (equmTRANS(it)%rb - taylorpers*deltatransvec(it)*equmINITSS%rnom &
									- (1.0-taylorpers*deltatransvec(it))*equmTRANS(it-1)%rnom  &
									- taylorpers*deltatransvec(it)*equmTRANS(it)%mpshock) / (phitaylor*taylorpers*deltatransvec(it)-1.0) !taylor rule
				equmTRANS(it)%rnom = equmTRANS(it)%rb + equmTRANS(it)%pi !fisher equn
			END DO
		END IF
	ELSE IF(zlbtransition==.true.) THEN
		DO it = 1,Ttransition
			IF(equmTRANS(it)%rb > (equmINITSS%rnom + equmTRANS(it)%mpshock) / phitaylor) THEN !ZLB does not bind
				equmTRANS(it)%pi = (equmTRANS(it)%rb - equmINITSS%rnom - equmTRANS(it)%mpshock) / (phitaylor-1.0) !taylor rule
				equmTRANS(it)%rnom = equmTRANS(it)%rb + equmTRANS(it)%pi !fisher equn		
			ELSE IF(equmTRANS(it)%rb <= (equmINITSS%rnom + equmTRANS(it)%mpshock)/phitaylor) THEN !ZLB binds
				equmTRANS(it)%pi = -equmTRANS(it)%rb
				equmTRANS(it)%rnom = 0.0
			END IF
		END DO
	END IF

	!solve phillips curve backwards for marginal costs
	IF (FirmDiscountRate==1) lfirmdiscount = equmTRANS(:)%rho
	IF (FirmDiscountRate==2) lfirmdiscount = equmINITSS%rb
	
	!update marginal costs
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
	equmTRANS(:)%priceadjust = (theta/2.0)*(equmTRANS(:)%pi**2)*equmTRANS(:)%capital/equmTRANS(:)%KYratio
	equmTRANS(:)%profit = (1.0-equmTRANS(:)%mc)*equmTRANS(:)%capital/equmTRANS(:)%KYratio - equmTRANS(:)%priceadjust

	!solve backward for investment
	it = Ttransition
	equmTRANS(it)%investment = (equmFINALSS%capital-equmTRANS(it)%capital)/deltatransvec(it) + deprec*equmTRANS(it)%capital
	DO it = Ttransition-1,1,-1
		equmTRANS(it)%investment = (equmTRANS(it+1)%capital-equmTRANS(it)%capital)/deltatransvec(it) + deprec*equmTRANS(it)%capital
	END DO

	!dividends and illiquid return
	IF(DividendSmoothing ==0) THEN
		IF(FixProfitsOutOfSteadyState==0) equmTRANS(:)%dividend 	= equmTRANS(:)%profit
		IF(FixProfitsOutOfSteadyState==1) equmTRANS(:)%dividend 	= equmINITSS%profit
	
		equmTRANS(:)%divrate 	= equmTRANS(:)%dividend/equmTRANS(:)%capital 
		equmTRANS(:)%ra = (equmTRANS(:)%rcapital - deprec + equmTRANS(:)%divrate - equmTRANS(:)%fundlev*equmTRANS(:)%rb) / (1.0-equmTRANS(:)%fundlev)
		equmTRANS(:)%intfirmbond = 0.0

	ELSE IF(DividendSmoothing ==1) THEN
		lpvprofit = equmFINALSS%profit/equmFINALSS%rb
		lpvfirmrb = 1.0/equmFINALSS%rb
		DO it = Ttransition,1
			lpvprofit = exp(-equmTRANS(it)%rb*deltatransvec(it))*(lpvprofit + equmTRANS(it)%profit*deltatransvec(it))
			lpvfirmrb = exp((-lfirmdiscount(it)/firmgamma + (1.0/firmgamma-1.0)*equmTRANS(it)%rb)*deltatransvec(it)) * lpvfirmrb
		END DO
		it = 1
	! 	equmTRANS(it)%dividend = (lpvprofit/lpvfirmrb) / (1.0 + (lfirmdiscount(it) - equmTRANS(it)%rb)*deltatransvec(it))
		equmTRANS(it)%dividend = (lpvprofit/lpvfirmrb) * (1.0 - (lfirmdiscount(it) - equmTRANS(it)%rb)*deltatransvec(it)/firmgamma)
		DO it = 2,Ttransition
	! 		equmTRANS(it)%dividend = equmTRANS(it-1)%dividend / (1.0 + (lfirmdiscount(it) - equmTRANS(it)%rb)*deltatransvec(it))
			equmTRANS(it)%dividend = equmTRANS(it-1)%dividend * (1.0 - (lfirmdiscount(it) - equmTRANS(it)%rb)*deltatransvec(it)/firmgamma)
		END DO
		equmTRANS(:)%divrate 	= equmTRANS(:)%dividend/equmTRANS(:)%capital 
		equmTRANS(:)%ra = (equmTRANS(:)%rcapital - deprec + equmTRANS(:)%divrate - equmTRANS(:)%fundlev*equmTRANS(:)%rb) / (1.0-equmTRANS(:)%fundlev)
		equmTRANS(1)%intfirmbond = 0.0
		DO it = 1, Ttransition-1
			equmTRANS(it+1)%intfirmbond = (equmTRANS(it)%intfirmbond + deltatransvec(it)*(equmTRANS(it)%profit-equmTRANS(it)%dividend)) / (1.0 - equmTRANS(it)%rb*deltatransvec(it))
		END DO

	ELSE IF(DividendSmoothing ==2) THEN

		equmTRANS(1)%dividend =  (1.0-divsmooth*cumdeltatrans(Ttransition)/(1.0 - exp(-divsmooth*cumdeltatrans(Ttransition)))) * equmINITSS%profit &
		 						+ (divsmooth/(1.0 - exp(-divsmooth*cumdeltatrans(Ttransition)))) * SUM(equmTRANS(:)%profit*deltatransvec)
		DO it = 1,Ttransition-1
			equmTRANS(it+1)%dividend = equmINITSS%dividend + (equmTRANS(it)%dividend - equmINITSS%dividend) / (1.0 + divsmooth*deltatransvec(it))
		END DO
		equmTRANS(:)%divrate 	= equmTRANS(:)%dividend/equmTRANS(:)%capital 
		equmTRANS(:)%ra = (equmTRANS(:)%rcapital - deprec + equmTRANS(:)%divrate - equmTRANS(:)%fundlev*equmTRANS(:)%rb) / (1.0-equmTRANS(:)%fundlev)
		equmTRANS(1)%intfirmbond = 0.0
		DO it = 1, Ttransition-1
			equmTRANS(it+1)%intfirmbond = (equmTRANS(it)%intfirmbond + deltatransvec(it)*(equmTRANS(it)%profit-equmTRANS(it)%dividend))
		END DO

	END IF
! 	write(*,*) ' lpvprofit is ',lpvprofit
! 	DO it = 1,Ttransition
! 		write(*,*) it, ' mc is ',equmTRANS(it)%mc
! 		write(*,*) it, ' price adjust 1 is ',(theta/2.0)*(equmTRANS(it)%pi**2)*equmTRANS(it)%output
! 		write(*,*) it, ' price adjust 2 is ',(theta/2.0)*(equmTRANS(it)%pi**2)*equmTRANS(it)%capital/equmTRANS(it)%KYratio
! 		write(*,*) it, ' profit is ',equmTRANS(it)%profit
! 		write(*,*) it, ' dividend is ',equmTRANS(it)%dividend
! 		write(*,*) it, ' divrate is ',equmTRANS(it)%divrate
! 		write(*,*) it, ' ra is ',equmTRANS(it)%ra
! 		write(*,*) it, ' intfirmbond is ',equmTRANS(it)%intfirmbond
! 	END DO

	
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

END SUBROUTINE IterateTransitionStickyB