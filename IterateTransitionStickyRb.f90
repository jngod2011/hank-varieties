SUBROUTINE IterateTransitionStickyRb

USE Parameters
USE Globals
USE Procedures

IMPLICIT NONE

INTEGER 	:: it,ii,infix,insmooth,it1
REAL(8) 	:: lldK,ldiffB,ldiffK,lqa,lqb,lqc,lminmargcost,lpvgovbc,lpvtaxableinc,lpvprofit,lpvfirmrb,lpvlumpincr,linitlumpincr,lpvfirma,lpvfirmb
REAL(8), DIMENSION(Ttransition) :: lcapital,lcapital1,lbond,lbond1,lfirmdiscount,lrb,lrb1,lfundbond,lworldbond,lrfirm,lrgov

iteratingtransition = .true.

lminmargcost = 0.01

IF(fsptransition==.false.) THEN
	IF(Display>=1 .and. stickytransition==.true.) write(*,*)' Solving for sticky price transition without ZLB'
	IF(Display>=1 .and. zlbtransition==.true.) write(*,*)' Solving for sticky price transition with ZLB'
ELSE IF (fsptransition==.true.) THEN
	IF(Display>=1 .and. stickytransition==.true.) write(*,*)' Solving for sticky price transition without ZLB (with FSP)'
	IF(Display>=1 .and. zlbtransition==.true.) write(*,*)' Solving for sticky price transition with ZLB (with FSP)'
END IF


!guess capital demand and liquid return
IF(fsptransition==.false. .and. UseFlexTransitionAsGuess==0) THEN

	!initial capital pins down labor and KN ratio in first period of transition
	equmTRANS(1)%capital = equmINITSS%capital
	!construct sequence of guesses of KN ratio: assume log linear in KN ratio (constant if a temporary transition)
	lldK = (log(equmFINALSS%capital)  - log(equmTRANS(1)%capital)) / real(Ttransition)
	DO it = 2,Ttransition
		equmTRANS(it)%capital = equmTRANS(1)%capital  *exp(lldK*it)
	END DO

	!construct sequence of guesses of Rb
	IF (GuessZeroInflation==1) THEN
		equmTRANS(:)%pi = equmINITSS%pi
		IF(forwardguide==.false. .or. ForwardGuideFixNomPreShock==0) THEN
		
			IF(zlbtransition==.false.) equmTRANS(:)%rnom = equmINITSS%rnom +phitaylor*equmTRANS(:)%pi + equmTRANS(:)%mpshock
			IF(zlbtransition==.true.) equmTRANS(:)%rnom = max(equmINITSS%rnom +phitaylor*equmTRANS(:)%pi + equmTRANS(:)%mpshock, 0.0)

		ELSE IF(forwardguide==.true. .and. ForwardGuideFixNomPreShock==1) THEN
			it1 = MINLOC(cumdeltatrans, 1, MASK = cumdeltatrans>=ForwardGuideShockQtrs)
			equmTRANS(1:it1-1)%rnom = equmINITSS%rnom
			IF(zlbtransition==.false.) equmTRANS(it1:Ttransition)%rnom = equmINITSS%rnom +phitaylor*equmTRANS(it1:Ttransition)%pi + equmTRANS(it1:Ttransition)%mpshock
			IF(zlbtransition==.true.) equmTRANS(it1:Ttransition)%rnom = max(equmINITSS%rnom +phitaylor*equmTRANS(it1:Ttransition)%pi + equmTRANS(:)%mpshock, 0.0)			
		END IF	
		equmTRANS(:)%rb = equmTRANS(:)%rnom - equmTRANS(:)%pi
	ELSEIF (GuessZeroInflation==0) THEN
		equmTRANS(:)%rb = equmINITSS%rb		
	END IF
	
ELSE IF(fsptransition==.false. .and. UseFlexTransitionAsGuess==1 ) THEN
	equmTRANS%capital = irfpointer%equmFLEX%capital
	equmTRANS%rb = irfpointer%equmFLEX%rb

ELSE IF	(fsptransition==.true.) THEN
	IF(stickytransition==.true.) THEN
		equmTRANS%capital = irfpointer%equmSTICKY%capital
		equmTRANS%rb = irfpointer%equmSTICKY%rb
	ELSE IF(zlbtransition==.true.) THEN
		equmTRANS%capital = irfpointer%equmZLB%capital
		equmTRANS%rb = irfpointer%equmZLB%rb
	END IF
END IF
equmTRANS(:)%rborr = equmTRANS(:)%rb + equmTRANS(:)%borrwedge

!world bond
equmTRANS(1)%worldbond = equmINITSS%worldbond
DO it = 1,Ttransition-1
	CALL WorldBondFunction2( equmTRANS(it)%rb,equmTRANS(it+1)%worldbond,equmINITSS%worldbond,equmINITSS%rb,bondelast)
	equmTRANS(it+1)%worldbond = equmTRANS(it)%worldbond + bondadjust*deltatransvec(it)*(equmTRANS(it+1)%worldbond-equmTRANS(it)%worldbond)
END DO

!fund bond
equmTRANS(:)%fundbond = -equmTRANS(:)%capital*equmTRANS(:)%fundlev

!inflation and nominal interest rates
IF(GuessZeroInflation==0) THEN
	IF(zlbtransition==.false.) THEN
		IF(BackwardTermInTaylorRule==0) THEN
			IF(forwardguide==.false. .or. ForwardGuideFixNomPreShock==0) THEN
				equmTRANS(:)%pi = (equmTRANS(:)%rb - equmINITSS%rnom - equmTRANS(:)%mpshock) / (phitaylor-1.0) !taylor rule
				equmTRANS(:)%rnom = equmTRANS(:)%rb + equmTRANS(:)%pi !fisher equn
			ELSE IF(forwardguide==.true. .and. ForwardGuideFixNomPreShock==1) THEN
				it1 = MINLOC(cumdeltatrans, 1, MASK = cumdeltatrans>=ForwardGuideShockQtrs)
				
				equmTRANS(1:it1-1)%rnom = equmINITSS%rnom
				equmTRANS(1:it1-1)%pi = equmINITSS%rnom-equmTRANS(1:it1-1)%rb
				
				equmTRANS(it1:Ttransition)%pi = (equmTRANS(it1:Ttransition)%rb - equmINITSS%rnom - equmTRANS(it1:Ttransition)%mpshock) / (phitaylor-1.0) !taylor rule
				equmTRANS(it1:Ttransition)%rnom = equmTRANS(it1:Ttransition)%rb + equmTRANS(it1:Ttransition)%pi !fisher equn
					
			END IF		


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
END IF

!solve phillips curve backwards for marginal costs
IF (FirmDiscountRate==1) lfirmdiscount = equmTRANS(:)%rho
IF (FirmDiscountRate==2) lfirmdiscount = equmINITSS%rb
IF (FirmDiscountRate==3) lfirmdiscount = equmINITSS%ra
IF (FirmDiscountRate==4) lfirmdiscount = equmTRANS(:)%rb
IF (FirmDiscountRate==5) lfirmdiscount = equmINITSS%ra

!marginal costs
!final period of transition
it = Ttransition
lqa = equmTRANS(it)%elast/theta

lqb = (equmFINALSS%pi-equmTRANS(it)%pi)/deltatransvec(it) - (equmFINALSS%pi/deltatransvec(it)) *(1.0-alpha)*frisch/(wageflex+alpha*frisch) &
 		-lfirmdiscount(it)*equmFINALSS%pi - (equmTRANS(it)%elast-1.0)/theta &
		+ equmFINALSS%pi * ((wageflex + frisch)/(wageflex +alpha*frisch))*(equmFINALSS%tfp-equmTRANS(it)%tfp)/(equmTRANS(it)%tfp*deltatransvec(it)) &
		+ equmFINALSS%pi * (alpha*(wageflex + frisch)/(wageflex+alpha*frisch))*(equmFINALSS%capital-equmTRANS(it)%capital)/(equmTRANS(it)%capital*deltatransvec(it)) &
		- ((1.0-alpha)*frisch*wageflex/(wageflex+alpha*frisch))*(equmFINALSS%labtax-equmTRANS(it)%labtax)/((1.0-equmTRANS(it)%labtax)*deltatransvec(it))

lqc = ((1.0-alpha)*frisch/(wageflex +alpha*frisch))*equmFINALSS%mc*equmFINALSS%pi/deltatransvec(it)

!
! IF(StickyWages==0) THEN
! 	lqb = (equmFINALSS%pi-equmTRANS(it)%pi)/deltatransvec(it) - (equmTRANS(it)%pi/deltatransvec(it)) *(1.0-alpha)*frisch/(1.0+alpha*frisch) &
! 	 		-lfirmdiscount(it)*equmFINALSS%pi - (equmTRANS(it)%elast-1.0)/theta &
! 			+ ((1.0+ frisch)/(1.0+alpha*frisch))*(equmFINALSS%tfp-equmTRANS(it)%tfp)/(equmTRANS(it)%tfp*deltatransvec(it)) &
! 			+ (alpha*(1.0+ frisch)/(1.0+alpha*frisch))*(equmFINALSS%capital-equmTRANS(it)%capital)/(equmTRANS(it)%capital*deltatransvec(it)) &
! 			- ((1.0-alpha)*(1.0+ frisch)/(1.0+alpha*frisch))*(equmFINALSS%labtax-equmTRANS(it)%labtax)/(equmTRANS(it)%labtax*deltatransvec(it))
!
! 	lqc = ((1.0-alpha)*frisch/(1.0+alpha*frisch))*equmFINALSS%mc*equmTRANS(it)%pi/deltatransvec(it)
! ELSE IF (StickyWages==1) THEN
! 	lqb = (equmFINALSS%pi-equmTRANS(it)%pi)/deltatransvec(it) - (equmTRANS(it)%pi/deltatransvec(it)) *(1.0-alpha)/alpha &
! 	 		-lfirmdiscount(it)*equmFINALSS%pi - (equmTRANS(it)%elast-1.0)/theta &
! 			+ (1.0/alpha)*(equmFINALSS%tfp-equmTRANS(it)%tfp)/(equmTRANS(it)%tfp*deltatransvec(it)) &
! 			+ (equmFINALSS%capital-equmTRANS(it)%capital)/(equmTRANS(it)%capital*deltatransvec(it))
!
! 	lqc = ((1.0-alpha)/alpha)*equmFINALSS%mc*equmTRANS(it)%pi/deltatransvec(it)
!
!
! END IF

IF (lqb**2 - 4.0*lqa*lqc >=0.0)	THEN
	equmTRANS(it)%mc = (-lqb + sqrt(lqb**2 - 4.0*lqa*lqc)) / (2.0*lqa)
	equmTRANS(it)%mc = max(lminmargcost,equmTRANS(it)%mc)
ELSE
	equmTRANS(it)%mc = lminmargcost
END IF

!solve backwards
DO it = Ttransition-1,1,-1
	lqa = equmTRANS(it)%elast/theta
	
	lqb = (equmTRANS(it+1)%pi-equmTRANS(it)%pi)/deltatransvec(it) - (equmTRANS(it+1)%pi/deltatransvec(it)) *(1.0-alpha)*frisch/(wageflex+alpha*frisch) &
	 		-lfirmdiscount(it)*equmTRANS(it+1)%pi - (equmTRANS(it)%elast-1.0)/theta &
			+ equmTRANS(it+1)%pi * ((wageflex + frisch)/(wageflex +alpha*frisch))*(equmTRANS(it+1)%tfp-equmTRANS(it)%tfp)/(equmTRANS(it)%tfp*deltatransvec(it)) &
			+ equmTRANS(it+1)%pi * (alpha*(wageflex + frisch)/(wageflex+alpha*frisch))*(equmTRANS(it+1)%capital-equmTRANS(it)%capital)/(equmTRANS(it)%capital*deltatransvec(it)) &
			- ((1.0-alpha)*frisch*wageflex/(wageflex+alpha*frisch))*(equmTRANS(it+1)%labtax-equmTRANS(it)%labtax)/((1.0-equmTRANS(it)%labtax)*deltatransvec(it))

	lqc = ((1.0-alpha)*frisch/(wageflex +alpha*frisch))*equmTRANS(it+1)%mc*equmTRANS(it+1)%pi/deltatransvec(it)
	
	IF (lqb**2 - 4.0*lqa*lqc >=0.0)	THEN
		equmTRANS(it)%mc = (-lqb + sqrt(lqb**2 - 4.0*lqa*lqc)) / (2.0*lqa)
		equmTRANS(it)%mc = max(lminmargcost,equmTRANS(it)%mc)
	ELSE
		equmTRANS(it)%mc = lminmargcost
	END IF

END DO

equmTRANS(:)%gap = equmTRANS(:)%elast*equmTRANS(:)%mc / (equmTRANS(:)%elast-1.0) - 1.0
equmTRANS(:)%tfpadj = (equmTRANS(:)%tfp**((1.0+utilelast)/utilelastalpha)) * ((equmTRANS(:)%mc*alpha/equmINITSS%rcapital)**(alpha*utilelast/utilelastalpha))
equmTRANS(:)%labor = ( ((1.0-equmTRANS(:)%labtax)/chi)** (frisch*wageflex*utilelastalpha/(wageflex*utilelastalpha + frisch*alpha)) ) * &
						((equmTRANS(:)%mc*(1.0-alpha)*equmTRANS(:)%tfpadj *(equmTRANS(:)%capital**(alpha/utilelastalpha)))/(equmINITSS%wage**(1.0-wageflex)) ) ** (frisch*utilelastalpha/(wageflex*utilelastalpha + frisch*alpha))
equmTRANS(:)%KNratio = equmTRANS(:)%capital/equmTRANS(:)%labor
equmTRANS(:)%wage = equmTRANS(:)%mc*(1.0-alpha)* equmTRANS(:)%tfpadj * (equmTRANS(:)%KNratio**(alpha/utilelastalpha))
equmTRANS(:)%netwage = (1.0-equmTRANS(:)%labtax)*equmTRANS(:)%wage
equmTRANS(:)%caputil = ((equmTRANS(:)%mc*alpha*equmTRANS(:)%tfp/equmINITSS%rcapital) * equmTRANS(:)%KNratio**(alpha-1.0)) ** (utilelast/utilelastalpha)

equmTRANS(:)%KYratio = (equmTRANS(:)%KNratio**(1.0-alpha)) / (equmTRANS(:)%tfp* (equmTRANS(:)%caputil**alpha))
equmTRANS(:)%rcapital = ((equmINITSS%rcapital**utilelast) * equmTRANS(:)%mc * alpha / equmTRANS(:)%KYratio ) ** (1.0/(1.0+utilelast))
equmTRANS(:)%priceadjust = (theta/2.0)*(equmTRANS(:)%pi**2)*equmTRANS(:)%capital/equmTRANS(:)%KYratio
equmTRANS(:)%profit = (1.0-equmTRANS(:)%mc-operatecost)*equmTRANS(:)%capital/equmTRANS(:)%KYratio - equmTRANS(:)%priceadjust

equmTRANS(:)%deprec = equmINITSS%deprec + (utilelast*equmINITSS%rcapital/(1.0+ utilelast)) * ((equmTRANS(:)%rcapital/equmINITSS%rcapital)**(1.0+utilelast) -1.0)

!solve backward for investment
it = Ttransition
equmTRANS(it)%investment = (equmFINALSS%capital-equmTRANS(it)%capital)/deltatransvec(it) + equmTRANS(it)%deprec*equmTRANS(it)%capital
DO it = Ttransition-1,1,-1
	equmTRANS(it)%investment = (equmTRANS(it+1)%capital-equmTRANS(it)%capital)/deltatransvec(it) + equmTRANS(it)%deprec*equmTRANS(it)%capital
END DO

!dividends and illiquid return
IF(DividendSmoothing ==0) THEN
	IF(FixProfitsOutOfSteadyState==0) equmTRANS(:)%dividend 	= equmTRANS(:)%profit*(1.0-corptax)
	IF(FixProfitsOutOfSteadyState==1) equmTRANS(:)%dividend 	= equmINITSS%profit*(1.0-corptax)

	IF(DividendFundLumpSum==0) equmTRANS(:)%divrate = equmTRANS(:)%dividend/equmTRANS(:)%capital
	IF(DividendFundLumpSum==1) equmTRANS(:)%divrate = 0.0

	equmTRANS(:)%ra = (equmTRANS(:)%rcapital*equmTRANS(:)%caputil - equmTRANS(:)%deprec + equmTRANS(:)%divrate - equmTRANS(:)%fundlev*equmTRANS(:)%rb) / (1.0-equmTRANS(:)%fundlev)
	equmTRANS(:)%intfirmbond = 0.0

ELSE IF(DividendSmoothing ==1) THEN
	lpvprofit = (1.0-corptax)*equmFINALSS%profit/equmFINALSS%rb
	lpvfirmrb = 1.0/equmFINALSS%rb
	DO it = Ttransition,1
		lpvprofit = exp(-equmTRANS(it)%rb*deltatransvec(it))*(lpvprofit + (1.0-corptax)*equmTRANS(it)%profit*deltatransvec(it))
		lpvfirmrb = exp((-lfirmdiscount(it)/firmgamma + (1.0/firmgamma-1.0)*equmTRANS(it)%rb)*deltatransvec(it)) * lpvfirmrb
	END DO
	it = 1
! 	equmTRANS(it)%dividend = (lpvprofit/lpvfirmrb) / (1.0 + (lfirmdiscount(it) - equmTRANS(it)%rb)*deltatransvec(it))
	equmTRANS(it)%dividend = (lpvprofit/lpvfirmrb) * (1.0 - (lfirmdiscount(it) - equmTRANS(it)%rb)*deltatransvec(it)/firmgamma)
	DO it = 2,Ttransition
! 		equmTRANS(it)%dividend = equmTRANS(it-1)%dividend / (1.0 + (lfirmdiscount(it) - equmTRANS(it)%rb)*deltatransvec(it))
		equmTRANS(it)%dividend = equmTRANS(it-1)%dividend * (1.0 - (lfirmdiscount(it) - equmTRANS(it)%rb)*deltatransvec(it)/firmgamma)
	END DO
	IF(DividendFundLumpSum==0) equmTRANS(:)%divrate = equmTRANS(:)%dividend/equmTRANS(:)%capital
	IF(DividendFundLumpSum==1) equmTRANS(:)%divrate = 0.0
	 
	equmTRANS(:)%ra = (equmTRANS(:)%rcapital*equmTRANS(:)%caputil - equmTRANS(:)%deprec + equmTRANS(:)%divrate - equmTRANS(:)%fundlev*equmTRANS(:)%rb) / (1.0-equmTRANS(:)%fundlev)
	equmTRANS(1)%intfirmbond = 0.0
	DO it = 1, Ttransition-1
		equmTRANS(it+1)%intfirmbond = (equmTRANS(it)%intfirmbond + deltatransvec(it)*((1.0-corptax)*equmTRANS(it)%profit-equmTRANS(it)%dividend)) / (1.0 - equmTRANS(it)%rb*deltatransvec(it))
	END DO

ELSE IF(DividendSmoothing ==2) THEN
	lpvfirma = 0.0
	lpvfirmb = 0.0
	DO it=Ttransition,1,-1
		lpvfirmb = lpvfirmb + deltatransvec(it)*((1.0-corptax)*equmTRANS(it)%profit-equmINITSS%dividend)
		lpvfirma = (lpvfirma + deltatransvec(it)) / (1.0 + deltatransvec(it)*divsmooth)
	END DO
	equmTRANS(1)%dividend = (lpvfirmb+lpvfirma*equmINITSS%dividend)/lpvfirma
	DO it = 1,Ttransition-1
		equmTRANS(it+1)%dividend = (equmTRANS(it)%dividend + divsmooth*deltatransvec(it)*equmINITSS%dividend) / (1.0+divsmooth*deltatransvec(it))
	END DO
	IF(DividendFundLumpSum==0) equmTRANS(:)%divrate = equmTRANS(:)%dividend/equmTRANS(:)%capital
	IF(DividendFundLumpSum==1) equmTRANS(:)%divrate = 0.0
	 
	equmTRANS(:)%ra = (equmTRANS(:)%rcapital*equmTRANS(:)%caputil - equmTRANS(:)%deprec + equmTRANS(:)%divrate - equmTRANS(:)%fundlev*equmTRANS(:)%rb) / (1.0-equmTRANS(:)%fundlev)
	equmTRANS(Ttransition)%intfirmbond = 0.0
	DO it = Ttransition-1,2,-1
		equmTRANS(it)%intfirmbond = equmTRANS(it+1)%intfirmbond - deltatransvec(it)*((1.0-corptax)*equmTRANS(it)%profit-equmTRANS(it)%dividend)
	END DO
	equmTRANS(1)%intfirmbond = 0.0
	
ELSE IF(DividendSmoothing ==3) THEN

	lpvfirma = 0.0
	lpvfirmb = 0.0
	DO it=Ttransition,1,-1
		lpvfirmb = (lpvfirmb + deltatransvec(it)*((1.0-corptax)*equmTRANS(it)%profit-equmINITSS%dividend)) / (1.0 + deltatransvec(it)*equmTRANS(it)%rb)
		lpvfirma = (lpvfirma + deltatransvec(it)) / (1.0 + deltatransvec(it)*(equmTRANS(it)%rb+divsmooth))
	END DO
	equmTRANS(1)%dividend = (lpvfirmb+lpvfirma*equmINITSS%dividend)/lpvfirma
	DO it = 1,Ttransition-1
		equmTRANS(it+1)%dividend = (equmTRANS(it)%dividend + divsmooth*deltatransvec(it)*equmINITSS%dividend) / (1.0+divsmooth*deltatransvec(it))
	END DO
	IF(DividendFundLumpSum==0) equmTRANS(:)%divrate = equmTRANS(:)%dividend/equmTRANS(:)%capital
	IF(DividendFundLumpSum==1) equmTRANS(:)%divrate = 0.0
	 
	equmTRANS(:)%ra = (equmTRANS(:)%rcapital*equmTRANS(:)%caputil - equmTRANS(:)%deprec + equmTRANS(:)%divrate - equmTRANS(:)%fundlev*equmTRANS(:)%rb) / (1.0-equmTRANS(:)%fundlev)
	equmTRANS(Ttransition)%intfirmbond = 0.0
	DO it = Ttransition-1,2,-1
		equmTRANS(it)%intfirmbond = (equmTRANS(it+1)%intfirmbond - deltatransvec(it)*((1.0-corptax)*equmTRANS(it)%profit-equmTRANS(it)%dividend)) / (1.0+deltatransvec(it)*equmTRANS(it)%rb)
	END DO
	equmTRANS(1)%intfirmbond = 0.0
ELSE IF(DividendSmoothing ==4) THEN

! 		lrfirm = 0.0
! 	lrfirm = equmINITSS%rb
	lrfirm = equmTRANS(:)%rb

	lpvfirma = 0.0
	lpvfirmb = 0.0
	DO it=Ttransition,1,-1
		lpvfirmb = (lpvfirmb + deltatransvec(it)*((1.0-corptax)*equmTRANS(it)%profit-equmINITSS%dividend)) / (1.0 + deltatransvec(it)*lrfirm(it))
		lpvfirma = (lpvfirma + deltatransvec(it)) / (1.0 + deltatransvec(it)*(lrfirm(it)+divsmooth))
	END DO
	equmTRANS(1)%dividend = (lpvfirmb+lpvfirma*equmINITSS%dividend)/lpvfirma
	DO it = 1,Ttransition-1
		equmTRANS(it+1)%dividend = (equmTRANS(it)%dividend + divsmooth*deltatransvec(it)*equmINITSS%dividend) / (1.0+divsmooth*deltatransvec(it))
	END DO

	equmTRANS(Ttransition)%intfirmbond = 0.0
	DO it = Ttransition-1,2,-1
		equmTRANS(it)%intfirmbond = (equmTRANS(it+1)%intfirmbond - deltatransvec(it)*((1.0-corptax)*equmTRANS(it)%profit-equmTRANS(it)%dividend)) / (1.0+deltatransvec(it)*lrfirm(it))
	END DO
	equmTRANS(1)%intfirmbond = 0.0

	equmTRANS(:)%dividend = equmTRANS(:)%dividend + (equmTRANS(:)%rb-lrfirm(:))*equmTRANS(:)%intfirmbond
	IF(DividendFundLumpSum==0) equmTRANS(:)%divrate = equmTRANS(:)%dividend/equmTRANS(:)%capital 
	IF(DividendFundLumpSum==1) equmTRANS(:)%divrate = 0.0
	
	equmTRANS(:)%ra = (equmTRANS(:)%rcapital*equmTRANS(:)%caputil - equmTRANS(:)%deprec + equmTRANS(:)%divrate - equmTRANS(:)%fundlev*equmTRANS(:)%rb) / (1.0-equmTRANS(:)%fundlev)

END IF
	
	
!value of equity component of investmemt fund
IF(DividendFundLumpSum==0) THEN
	equmTRANS(:)%equity = 0.0
	illequitydrop = 1.0	
ELSE IF(DividendFundLumpSum==1) THEN
	it = Ttransition
	equmTRANS(it)%equity = (equmFINALSS%equity + equmTRANS(it)%dividend*deltatransvec(it)) / (1.0+deltatransvec(it)*equmTRANS(it)%ra)		
	DO it = Ttransition-1,1,-1
		equmTRANS(it)%equity = (equmTRANS(it+1)%equity + equmTRANS(it)%dividend*deltatransvec(it)) / (1.0+deltatransvec(it)*equmTRANS(it)%ra)		
	END DO
	illequitydrop = ((1.0-equmTRANS(1)%fundlev)*equmTRANS(1)%capital + equmTRANS(1)%equity) / ((1.0-equmINITSS%fundlev)*equmINITSS%capital + equmINITSS%equity)
END IF

equmTRANS(:)%output = equmTRANS(:)%tfpadj * (equmTRANS(:)%capital**(alpha/utilelastalpha)) * (equmTRANS(:)%labor**((1.0-alpha)*(1.0+utilelast)/utilelastalpha)) + rhousing*housefrac*((1.0-equmTRANS(:)%fundlev)*equmTRANS(:)%capital + equmTRANS(:)%equity)/(1.0-housefrac)	

	
!government budget constraint,expenditures and tax rates
IF(fsptransition==.false.) THEN 
	IF (AdjGovBudgetConstraint==1) THEN !adjust spending
		equmTRANS(:)%govbond = equmINITSS%govbond
		equmTRANS(:)%labtax = equmINITSS%labtax
		equmTRANS(:)%lumptransfer = equmINITSS%lumptransfer		
		equmTRANS(:)%taxrev = equmTRANS(:)%labtax*equmTRANS(:)%wage*equmTRANS(:)%labor - equmTRANS(:)%lumptransfer + corptax*equmTRANS(:)%profit
		IF(TaxRetainedEarnings==1) THEN
			equmTRANS(:)%taxrev = equmTRANS(:)%taxrev + (1.0-corptax)*equmTRANS(:)%profit-equmTRANS(:)%dividend 
			equmTRANS(:)%intfirmbond = 0.0
		END IF
		equmTRANS(:)%govexp = equmTRANS(:)%taxrev + equmTRANS(:)%rb*equmINITSS%govbond

	ELSE IF(AdjGovBudgetConstraint==2) THEN  !adjust lump sum taxes
		equmTRANS(:)%govbond = equmINITSS%govbond
		equmTRANS(:)%govexp = equmINITSS%govexp
		equmTRANS(:)%labtax = equmINITSS%labtax
		equmTRANS(:)%lumptransfer = equmTRANS(:)%labtax*equmTRANS(:)%wage*equmTRANS(:)%labor + corptax*equmTRANS(:)%profit + equmTRANS(:)%rb*equmINITSS%govbond - equmTRANS(:)%govexp
		equmTRANS(:)%taxrev = equmTRANS(:)%labtax*equmTRANS(:)%wage*equmTRANS(:)%labor - equmTRANS(:)%lumptransfer + corptax*equmTRANS(:)%profit
		IF(TaxRetainedEarnings==1) THEN
			equmTRANS(:)%lumptransfer = equmTRANS(:)%lumptransfer + (1.0-corptax)*equmTRANS(:)%profit-equmTRANS(:)%dividend 
			equmTRANS(:)%intfirmbond = 0.0
		END IF
				
	ELSE IF(AdjGovBudgetConstraint==3) THEN !adjust debt
		IF(GovExpConstantFracOutput==0) equmTRANS(:)%govexp = equmINITSS%govexp
		IF(GovExpConstantFracOutput==1) equmTRANS(:)%govexp = equmTRANS(:)%output*equmINITSS%govexp/equmINITSS%output

		equmTRANS(:)%lumptransfer = equmINITSS%lumptransfer
		equmTRANS(:)%taxrev = equmTRANS(:)%labtax*equmTRANS(:)%wage*equmTRANS(:)%labor - equmTRANS(:)%lumptransfer + corptax*equmTRANS(:)%profit
		IF(TaxRetainedEarnings==1) THEN
			equmTRANS(:)%taxrev = equmTRANS(:)%taxrev + (1.0-corptax)*equmTRANS(:)%profit-equmTRANS(:)%dividend 
			equmTRANS(:)%intfirmbond = 0.0
		END IF
	
		!compute required increase in lumptransfer
 		lrgov = equmTRANS(:)%rb
! 		lrgov = equmINITSS%rb
	
		lpvgovbc = equmFINALSS%govbond
		lpvlumpincr = 0.0
		DO it = Ttransition,1,-1
			lpvgovbc = (lpvgovbc + deltatransvec(it)*(equmTRANS(it)%govexp - equmTRANS(it)%taxrev))/(1.0+deltatransvec(it)*lrgov(it))
 			IF(cumdeltatrans(it)>=taxincrstart) lpvlumpincr = (lpvlumpincr + deltatransvec(it))/(1.0+deltatransvec(it)*(lrgov(it)+taxincrdecay))
			IF(cumdeltatrans(it)<taxincrstart) lpvlumpincr = lpvlumpincr/(1.0+deltatransvec(it)*lrgov(it))
		END DO	
	
		linitlumpincr = (equmINITSS%govbond-lpvgovbc) / lpvlumpincr
		DO it = 1,Ttransition
			IF(cumdeltatrans(it)>=taxincrstart) equmTRANS(it)%lumptransfer = equmTRANS(it)%lumptransfer + linitlumpincr*exp(-taxincrdecay*(cumdeltatrans(it)-taxincrstart))
		END DO
	
		equmTRANS(:)%taxrev = equmTRANS(:)%labtax*equmTRANS(:)%wage*equmTRANS(:)%labor - equmTRANS(:)%lumptransfer + corptax*equmTRANS(:)%profit
		IF(TaxRetainedEarnings==1) equmTRANS(:)%taxrev = equmTRANS(:)%taxrev + (1.0-corptax)*equmTRANS(:)%profit-equmTRANS(:)%dividend 
			
		equmTRANS(Ttransition)%govbond = equmFINALSS%govbond
		DO it = Ttransition-1,2,-1
			equmTRANS(it)%govbond = (equmTRANS(it+1)%govbond - deltatransvec(it)*(equmTRANS(it)%taxrev-equmTRANS(it)%govexp)) / (1.0+deltatransvec(it)*lrgov(it))
		END DO
		equmTRANS(1)%govbond = equmINITSS%govbond
	
		equmTRANS(:)%lumptransfer = equmTRANS(:)%lumptransfer + (equmTRANS(:)%rb-lrgov(:))*equmTRANS(:)%govbond
		equmTRANS(:)%taxrev = equmTRANS(:)%labtax*equmTRANS(:)%wage*equmTRANS(:)%labor - equmTRANS(:)%lumptransfer + corptax*equmTRANS(:)%profit
		IF(TaxRetainedEarnings==1) equmTRANS(:)%taxrev = equmTRANS(:)%taxrev + (1.0-corptax)*equmTRANS(:)%profit-equmTRANS(:)%dividend 
	
	
	END IF
ELSE IF(fsptransition==.true.) THEN	
! 	!tax increase guess
! 	ltaxincr = 0.003
! 	equmTRANS(:)%labtax = equmINITSS%labtax
! 	DO it = 1,Ttransition
! 		IF(cumdeltatrans(it)>=fspointer%labtaxstart .and. cumdeltatrans(it)<fspointer%labtaxend) THEN
! 			equmTRANS(it)%labtax = equmINITSS%labtax +ltaxincr
! 		END  IF
! 	END DO
		
END IF

IF(RetainedEarningsInBondMarket==1) equmTRANS(:)%bond = -equmTRANS(:)%worldbond - equmTRANS(:)%govbond - equmTRANS(:)%fundbond - equmTRANS(:)%intfirmbond
IF(RetainedEarningsInBondMarket==0) equmTRANS(:)%bond = -equmTRANS(:)%worldbond - equmTRANS(:)%govbond - equmTRANS(:)%fundbond

ii = 1	 
ldiffK = 1.0
ldiffB = 1.0
DO WHILE (ii<=maxitertranssticky .and. max(ldiffK,ldiffB)>toltransition )
	!solve for transtion
	CALL Transition
	
	!computed implied equilibrium quantities
	lbond = statsTRANS(:)%Eb
	lcapital = ((1.0-housefrac)*statsTRANS(:)%Ea - equmTRANS(:)%equity)/ (1.0 - equmTRANS(:)%fundlev)
	IF(RetainedEarningsAsCapital==1) lcapital = lcapital + equmTRANS(:)%intfirmbond / (1.0 - equmTRANS(:)%fundlev)
	IF(ConvergenceRelToOutput==0) THEN
		ldiffK= maxval(abs(lcapital/equmTRANS(:)%capital - 1.0))
		ldiffB= maxval(abs(lbond/equmTRANS(:)%bond - 1.0))
	ELSEIF(ConvergenceRelToOutput==1) THEN
		ldiffK= maxval(abs(lcapital-equmTRANS(:)%capital)/equmINITSS%output)
		ldiffB= maxval(abs(lbond-equmTRANS(:)%bond)/equmINITSS%output)
	END IF
	IF (Display>=1) write(*,"(A,I,A)") '  Transition iter ',ii, ':'
	IF (Display>=1 .and. fsptransition==.false.) write(*,"(A,E10.3,A,E10.3,A,E10.3)") '   K err',ldiffK, ',  B err',ldiffB
	
	!update capital and interest rate
	IF (ii<maxitertranssticky .and. max(ldiffK,ldiffB)>toltransition ) THEN
		CALL PartialUpdate(Ttransition-1,stepstickytransK,equmTRANS(2:Ttransition)%capital,lcapital(2:Ttransition),lcapital1(2:Ttransition))
		lcapital1(1)  = equmTRANS(1)%capital
		equmTRANS(:)%capital = lcapital1
		
		IF(UpdateRbFromMarketClearing==0) THEN
			lfundbond = -lcapital*equmTRANS(:)%fundlev
			IF(RetainedEarningsInBondMarket==1) lworldbond = -lbond - equmTRANS(:)%govbond - lfundbond - equmTRANS(:)%intfirmbond
			IF(RetainedEarningsInBondMarket==0) lworldbond = -lbond - equmTRANS(:)%govbond - lfundbond

			it = Ttransition
			CALL WorldBondInverse2( (equmFINALSS%worldbond-lworldbond(it))/(bondadjust*deltatransvec(it)) + lworldbond(it) &
								,lrb(it),equmINITSS%worldbond,equmINITSS%rb,bondelast)
			DO it = Ttransition-1,1,-1
				CALL WorldBondInverse2( (lworldbond(it+1)-lworldbond(it))/(bondadjust*deltatransvec(it)) + lworldbond(it) &
									,lrb(it),equmINITSS%worldbond,equmINITSS%rb,bondelast)
			END DO
	
			CALL PartialUpdate(Ttransition,stepstickytransB,equmTRANS(:)%rb,lrb,lrb1)
			equmTRANS(:)%rb = lrb1
		ELSE IF(UpdateRbFromMarketClearing==1) THEN
			lrb1(1:Ttransition-1) = equmTRANS(1:Ttransition-1)%rb - stepstickytransRb*(lbond(2:Ttransition)-equmTRANS(2:Ttransition)%bond)/equmINITSS%output
			lrb1(Ttransition) = equmINITSS%rb
			equmTRANS(:)%rb = lrb1
		
		END IF

		!impose it settles down
		IF(ImposeGuessSettlesDown==1) THEN
			infix = 1
			insmooth = Ttransition-10

			equmTRANS(Ttransition-infix+1:Ttransition)%capital = equmFINALSS%capital
 			CALL MASmooth(insmooth-infix+1,equmTRANS(Ttransition-insmooth:Ttransition-infix)%capital,deltatransvec(Ttransition-insmooth:Ttransition-infix),3,equmTRANS(Ttransition-insmooth:Ttransition-infix)%capital)

			equmTRANS(Ttransition-infix+1:Ttransition)%rb = equmFINALSS%rb
 			CALL MASmooth(insmooth-infix+1,equmTRANS(Ttransition-insmooth:Ttransition-infix)%rb,deltatransvec(Ttransition-insmooth:Ttransition-infix),3,equmTRANS(Ttransition-insmooth:Ttransition-infix)%rb)
		END IF

	ElSE
		!run distribution stats with full
		iteratingtransition = .false.
		CALL Transition
		equmTRANS(:)%capital = lcapital
		equmTRANS(:)%bond = lbond
		equmTRANS(:)%rb = lrb
		
	END IF
	
	equmTRANS(:)%rborr = equmTRANS(:)%rb + equmTRANS(:)%borrwedge
	
	!world bond
	equmTRANS(1)%worldbond = equmINITSS%worldbond
	DO it = 1,Ttransition-1
		CALL WorldBondFunction2( equmTRANS(it)%rb,equmTRANS(it+1)%worldbond,equmINITSS%worldbond,equmINITSS%rb,bondelast)
		equmTRANS(it+1)%worldbond = equmTRANS(it)%worldbond + bondadjust*deltatransvec(it)*(equmTRANS(it+1)%worldbond-equmTRANS(it)%worldbond)
	END DO

	!fund bond
	equmTRANS(:)%fundbond = -equmTRANS(:)%capital*equmTRANS(:)%fundlev

	!inflation and nominal interest rates
	IF(zlbtransition==.false.) THEN
		IF(BackwardTermInTaylorRule==0) THEN
			IF(forwardguide==.false. .or. ForwardGuideFixNomPreShock==0) THEN
				equmTRANS(:)%pi = (equmTRANS(:)%rb - equmINITSS%rnom - equmTRANS(:)%mpshock) / (phitaylor-1.0) !taylor rule
				equmTRANS(:)%rnom = equmTRANS(:)%rb + equmTRANS(:)%pi !fisher equn
			ELSE IF(forwardguide==.true. .and. ForwardGuideFixNomPreShock==1) THEN
				it1 = MINLOC(cumdeltatrans, 1, MASK = cumdeltatrans>=ForwardGuideShockQtrs)
				
				equmTRANS(1:it1-1)%rnom = equmINITSS%rnom
				equmTRANS(1:it1-1)%pi = equmINITSS%rnom-equmTRANS(1:it1-1)%rb
				
				equmTRANS(it1:Ttransition)%pi = (equmTRANS(it1:Ttransition)%rb - equmINITSS%rnom - equmTRANS(it1:Ttransition)%mpshock) / (phitaylor-1.0) !taylor rule
				equmTRANS(it1:Ttransition)%rnom = equmTRANS(it1:Ttransition)%rb + equmTRANS(it1:Ttransition)%pi !fisher equn
					
			END IF		
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
	IF (FirmDiscountRate==3) lfirmdiscount = equmINITSS%ra
	IF (FirmDiscountRate==4) lfirmdiscount = equmTRANS(:)%rb
	IF (FirmDiscountRate==5) lfirmdiscount = equmTRANS(:)%ra


	!marginal costs
	!final period of transition
	it = Ttransition
	lqa = equmTRANS(it)%elast/theta

	lqa = equmTRANS(it)%elast/theta

	lqb = (equmFINALSS%pi-equmTRANS(it)%pi)/deltatransvec(it) - (equmFINALSS%pi/deltatransvec(it)) *(1.0-alpha)*frisch/(wageflex+alpha*frisch) &
	 		-lfirmdiscount(it)*equmFINALSS%pi - (equmTRANS(it)%elast-1.0)/theta &
			+ equmFINALSS%pi * ((wageflex + frisch)/(wageflex +alpha*frisch))*(equmFINALSS%tfp-equmTRANS(it)%tfp)/(equmTRANS(it)%tfp*deltatransvec(it)) &
			+ equmFINALSS%pi * (alpha*(wageflex + frisch)/(wageflex+alpha*frisch))*(equmFINALSS%capital-equmTRANS(it)%capital)/(equmTRANS(it)%capital*deltatransvec(it)) &
			- ((1.0-alpha)*frisch*wageflex/(wageflex+alpha*frisch))*(equmFINALSS%labtax-equmTRANS(it)%labtax)/((1.0-equmTRANS(it)%labtax)*deltatransvec(it))

	lqc = ((1.0-alpha)*frisch/(wageflex +alpha*frisch))*equmFINALSS%mc*equmFINALSS%pi/deltatransvec(it)

	IF (lqb**2 - 4.0*lqa*lqc >=0.0)	THEN
		equmTRANS(it)%mc = (-lqb + sqrt(lqb**2 - 4.0*lqa*lqc)) / (2.0*lqa)
		equmTRANS(it)%mc = max(lminmargcost,equmTRANS(it)%mc)
	ELSE
		equmTRANS(it)%mc = lminmargcost
	END IF

	!solve backwards
	DO it = Ttransition-1,1,-1
		lqa = equmTRANS(it)%elast/theta
		
		lqb = (equmTRANS(it+1)%pi-equmTRANS(it)%pi)/deltatransvec(it) - (equmTRANS(it+1)%pi/deltatransvec(it)) *(1.0-alpha)*frisch/(wageflex+alpha*frisch) &
		 		-lfirmdiscount(it)*equmTRANS(it+1)%pi - (equmTRANS(it)%elast-1.0)/theta &
				+ equmTRANS(it+1)%pi * ((wageflex + frisch)/(wageflex +alpha*frisch))*(equmTRANS(it+1)%tfp-equmTRANS(it)%tfp)/(equmTRANS(it)%tfp*deltatransvec(it)) &
				+ equmTRANS(it+1)%pi * (alpha*(wageflex + frisch)/(wageflex+alpha*frisch))*(equmTRANS(it+1)%capital-equmTRANS(it)%capital)/(equmTRANS(it)%capital*deltatransvec(it)) &
				- ((1.0-alpha)*frisch*wageflex/(wageflex+alpha*frisch))*(equmTRANS(it+1)%labtax-equmTRANS(it)%labtax)/((1.0-equmTRANS(it)%labtax)*deltatransvec(it))

		lqc = ((1.0-alpha)*frisch/(wageflex +alpha*frisch))*equmTRANS(it+1)%mc*equmTRANS(it+1)%pi/deltatransvec(it)
		
		IF (lqb**2 - 4.0*lqa*lqc >=0.0)	THEN
			equmTRANS(it)%mc = (-lqb + sqrt(lqb**2 - 4.0*lqa*lqc)) / (2.0*lqa)
			equmTRANS(it)%mc = max(lminmargcost,equmTRANS(it)%mc)
		ELSE
			equmTRANS(it)%mc = lminmargcost
		END IF

	END DO

	equmTRANS(:)%gap = equmTRANS(:)%elast*equmTRANS(:)%mc / (equmTRANS(:)%elast-1.0) - 1.0
	equmTRANS(:)%tfpadj = (equmTRANS(:)%tfp**((1.0+utilelast)/utilelastalpha)) * ((equmTRANS(:)%mc*alpha/equmINITSS%rcapital)**(alpha*utilelast/utilelastalpha))
	equmTRANS(:)%labor = ( ((1.0-equmTRANS(:)%labtax)/chi)** (frisch*wageflex*utilelastalpha/(wageflex*utilelastalpha + frisch*alpha)) ) * &
							((equmTRANS(:)%mc*(1.0-alpha)*equmTRANS(:)%tfpadj *(equmTRANS(:)%capital**(alpha/utilelastalpha)))/(equmINITSS%wage**(1.0-wageflex)) ) ** (frisch*utilelastalpha/(wageflex*utilelastalpha + frisch*alpha))
	equmTRANS(:)%KNratio = equmTRANS(:)%capital/equmTRANS(:)%labor
	equmTRANS(:)%wage = equmTRANS(:)%mc*(1.0-alpha)* equmTRANS(:)%tfpadj * (equmTRANS(:)%KNratio**(alpha/utilelastalpha))
	equmTRANS(:)%netwage = (1.0-equmTRANS(:)%labtax)*equmTRANS(:)%wage
	equmTRANS(:)%caputil = ((equmTRANS(:)%mc*alpha*equmTRANS(:)%tfp/equmINITSS%rcapital) * equmTRANS(:)%KNratio**(alpha-1.0)) ** (utilelast/utilelastalpha)

	equmTRANS(:)%KYratio = (equmTRANS(:)%KNratio**(1.0-alpha)) / (equmTRANS(:)%tfp* (equmTRANS(:)%caputil**alpha))
	equmTRANS(:)%rcapital = ((equmINITSS%rcapital**utilelast) * equmTRANS(:)%mc * alpha / equmTRANS(:)%KYratio ) ** (1.0/(1.0+utilelast))
	equmTRANS(:)%priceadjust = (theta/2.0)*(equmTRANS(:)%pi**2)*equmTRANS(:)%capital/equmTRANS(:)%KYratio
	equmTRANS(:)%profit = (1.0-equmTRANS(:)%mc-operatecost)*equmTRANS(:)%capital/equmTRANS(:)%KYratio - equmTRANS(:)%priceadjust

	equmTRANS(:)%deprec = equmINITSS%deprec + (utilelast*equmINITSS%rcapital/(1.0+ utilelast)) * ((equmTRANS(:)%rcapital/equmINITSS%rcapital)**(1.0+utilelast) -1.0)

	!solve backward for investment
	it = Ttransition
	equmTRANS(it)%investment = (equmFINALSS%capital-equmTRANS(it)%capital)/deltatransvec(it) + equmTRANS(it)%deprec*equmTRANS(it)%capital
	DO it = Ttransition-1,1,-1
		equmTRANS(it)%investment = (equmTRANS(it+1)%capital-equmTRANS(it)%capital)/deltatransvec(it) + equmTRANS(it)%deprec*equmTRANS(it)%capital
	END DO

	!dividends and illiquid return
	IF(DividendSmoothing ==0) THEN
		IF(FixProfitsOutOfSteadyState==0) equmTRANS(:)%dividend 	= equmTRANS(:)%profit*(1.0-corptax)
		IF(FixProfitsOutOfSteadyState==1) equmTRANS(:)%dividend 	= equmINITSS%profit*(1.0-corptax)

		IF(DividendFundLumpSum==0) equmTRANS(:)%divrate = equmTRANS(:)%dividend/equmTRANS(:)%capital
		IF(DividendFundLumpSum==1) equmTRANS(:)%divrate = 0.0

		equmTRANS(:)%ra = (equmTRANS(:)%rcapital*equmTRANS(:)%caputil - equmTRANS(:)%deprec + equmTRANS(:)%divrate - equmTRANS(:)%fundlev*equmTRANS(:)%rb) / (1.0-equmTRANS(:)%fundlev)
		equmTRANS(:)%intfirmbond = 0.0

	ELSE IF(DividendSmoothing ==1) THEN
		lpvprofit = (1.0-corptax)*equmFINALSS%profit/equmFINALSS%rb
		lpvfirmrb = 1.0/equmFINALSS%rb
		DO it = Ttransition,1
			lpvprofit = exp(-equmTRANS(it)%rb*deltatransvec(it))*(lpvprofit + (1.0-corptax)*equmTRANS(it)%profit*deltatransvec(it))
			lpvfirmrb = exp((-lfirmdiscount(it)/firmgamma + (1.0/firmgamma-1.0)*equmTRANS(it)%rb)*deltatransvec(it)) * lpvfirmrb
		END DO
		it = 1
	! 	equmTRANS(it)%dividend = (lpvprofit/lpvfirmrb) / (1.0 + (lfirmdiscount(it) - equmTRANS(it)%rb)*deltatransvec(it))
		equmTRANS(it)%dividend = (lpvprofit/lpvfirmrb) * (1.0 - (lfirmdiscount(it) - equmTRANS(it)%rb)*deltatransvec(it)/firmgamma)
		DO it = 2,Ttransition
	! 		equmTRANS(it)%dividend = equmTRANS(it-1)%dividend / (1.0 + (lfirmdiscount(it) - equmTRANS(it)%rb)*deltatransvec(it))
			equmTRANS(it)%dividend = equmTRANS(it-1)%dividend * (1.0 - (lfirmdiscount(it) - equmTRANS(it)%rb)*deltatransvec(it)/firmgamma)
		END DO
		IF(DividendFundLumpSum==0) equmTRANS(:)%divrate = equmTRANS(:)%dividend/equmTRANS(:)%capital
		IF(DividendFundLumpSum==1) equmTRANS(:)%divrate = 0.0
	 
		equmTRANS(:)%ra = (equmTRANS(:)%rcapital*equmTRANS(:)%caputil - equmTRANS(:)%deprec + equmTRANS(:)%divrate - equmTRANS(:)%fundlev*equmTRANS(:)%rb) / (1.0-equmTRANS(:)%fundlev)
		equmTRANS(1)%intfirmbond = 0.0
		DO it = 1, Ttransition-1
			equmTRANS(it+1)%intfirmbond = (equmTRANS(it)%intfirmbond + deltatransvec(it)*((1.0-corptax)*equmTRANS(it)%profit-equmTRANS(it)%dividend)) / (1.0 - equmTRANS(it)%rb*deltatransvec(it))
		END DO

	ELSE IF(DividendSmoothing ==2) THEN
		lpvfirma = 0.0
		lpvfirmb = 0.0
		DO it=Ttransition,1,-1
			lpvfirmb = lpvfirmb + deltatransvec(it)*((1.0-corptax)*equmTRANS(it)%profit-equmINITSS%dividend)
			lpvfirma = (lpvfirma + deltatransvec(it)) / (1.0 + deltatransvec(it)*divsmooth)
		END DO
		equmTRANS(1)%dividend = (lpvfirmb+lpvfirma*equmINITSS%dividend)/lpvfirma
		DO it = 1,Ttransition-1
			equmTRANS(it+1)%dividend = (equmTRANS(it)%dividend + divsmooth*deltatransvec(it)*equmINITSS%dividend) / (1.0+divsmooth*deltatransvec(it))
		END DO
		IF(DividendFundLumpSum==0) equmTRANS(:)%divrate = equmTRANS(:)%dividend/equmTRANS(:)%capital
		IF(DividendFundLumpSum==1) equmTRANS(:)%divrate = 0.0
	 
		equmTRANS(:)%ra = (equmTRANS(:)%rcapital*equmTRANS(:)%caputil - equmTRANS(:)%deprec + equmTRANS(:)%divrate - equmTRANS(:)%fundlev*equmTRANS(:)%rb) / (1.0-equmTRANS(:)%fundlev)
		equmTRANS(Ttransition)%intfirmbond = 0.0
		DO it = Ttransition-1,2,-1
			equmTRANS(it)%intfirmbond = equmTRANS(it+1)%intfirmbond - deltatransvec(it)*((1.0-corptax)*equmTRANS(it)%profit-equmTRANS(it)%dividend)
		END DO
		equmTRANS(1)%intfirmbond = 0.0
	
	ELSE IF(DividendSmoothing ==3) THEN

		lpvfirma = 0.0
		lpvfirmb = 0.0
		DO it=Ttransition,1,-1
			lpvfirmb = (lpvfirmb + deltatransvec(it)*((1.0-corptax)*equmTRANS(it)%profit-equmINITSS%dividend)) / (1.0 + deltatransvec(it)*equmTRANS(it)%rb)
			lpvfirma = (lpvfirma + deltatransvec(it)) / (1.0 + deltatransvec(it)*(equmTRANS(it)%rb+divsmooth))
		END DO
		equmTRANS(1)%dividend = (lpvfirmb+lpvfirma*equmINITSS%dividend)/lpvfirma
		DO it = 1,Ttransition-1
			equmTRANS(it+1)%dividend = (equmTRANS(it)%dividend + divsmooth*deltatransvec(it)*equmINITSS%dividend) / (1.0+divsmooth*deltatransvec(it))
		END DO
		IF(DividendFundLumpSum==0) equmTRANS(:)%divrate = equmTRANS(:)%dividend/equmTRANS(:)%capital
		IF(DividendFundLumpSum==1) equmTRANS(:)%divrate = 0.0
	 
		equmTRANS(:)%ra = (equmTRANS(:)%rcapital*equmTRANS(:)%caputil - equmTRANS(:)%deprec + equmTRANS(:)%divrate - equmTRANS(:)%fundlev*equmTRANS(:)%rb) / (1.0-equmTRANS(:)%fundlev)
		equmTRANS(Ttransition)%intfirmbond = 0.0
		DO it = Ttransition-1,2,-1
			equmTRANS(it)%intfirmbond = (equmTRANS(it+1)%intfirmbond - deltatransvec(it)*((1.0-corptax)*equmTRANS(it)%profit-equmTRANS(it)%dividend)) / (1.0+deltatransvec(it)*equmTRANS(it)%rb)
		END DO
		equmTRANS(1)%intfirmbond = 0.0
	ELSE IF(DividendSmoothing ==4) THEN

	! 		lrfirm = 0.0
	! 	lrfirm = equmINITSS%rb
		lrfirm = equmTRANS(:)%rb

		lpvfirma = 0.0
		lpvfirmb = 0.0
		DO it=Ttransition,1,-1
			lpvfirmb = (lpvfirmb + deltatransvec(it)*((1.0-corptax)*equmTRANS(it)%profit-equmINITSS%dividend)) / (1.0 + deltatransvec(it)*lrfirm(it))
			lpvfirma = (lpvfirma + deltatransvec(it)) / (1.0 + deltatransvec(it)*(lrfirm(it)+divsmooth))
		END DO
		equmTRANS(1)%dividend = (lpvfirmb+lpvfirma*equmINITSS%dividend)/lpvfirma
		DO it = 1,Ttransition-1
			equmTRANS(it+1)%dividend = (equmTRANS(it)%dividend + divsmooth*deltatransvec(it)*equmINITSS%dividend) / (1.0+divsmooth*deltatransvec(it))
		END DO

		equmTRANS(Ttransition)%intfirmbond = 0.0
		DO it = Ttransition-1,2,-1
			equmTRANS(it)%intfirmbond = (equmTRANS(it+1)%intfirmbond - deltatransvec(it)*((1.0-corptax)*equmTRANS(it)%profit-equmTRANS(it)%dividend)) / (1.0+deltatransvec(it)*lrfirm(it))
		END DO
		equmTRANS(1)%intfirmbond = 0.0

		equmTRANS(:)%dividend = equmTRANS(:)%dividend + (equmTRANS(:)%rb-lrfirm(:))*equmTRANS(:)%intfirmbond
		IF(DividendFundLumpSum==0) equmTRANS(:)%divrate = equmTRANS(:)%dividend/equmTRANS(:)%capital 
		IF(DividendFundLumpSum==1) equmTRANS(:)%divrate = 0.0
	
		equmTRANS(:)%ra = (equmTRANS(:)%rcapital*equmTRANS(:)%caputil - equmTRANS(:)%deprec + equmTRANS(:)%divrate - equmTRANS(:)%fundlev*equmTRANS(:)%rb) / (1.0-equmTRANS(:)%fundlev)

	END IF
	
	
	!value of equity component of investmemt fund
	IF(DividendFundLumpSum==0) THEN
		equmTRANS(:)%equity = 0.0
		illequitydrop = 1.0	
	ELSE IF(DividendFundLumpSum==1) THEN
		it = Ttransition
		equmTRANS(it)%equity = (equmFINALSS%equity + equmTRANS(it)%dividend*deltatransvec(it)) / (1.0+deltatransvec(it)*equmTRANS(it)%ra)		
		DO it = Ttransition-1,1,-1
			equmTRANS(it)%equity = (equmTRANS(it+1)%equity + equmTRANS(it)%dividend*deltatransvec(it)) / (1.0+deltatransvec(it)*equmTRANS(it)%ra)		
		END DO
		illequitydrop = ((1.0-equmTRANS(1)%fundlev)*equmTRANS(1)%capital + equmTRANS(1)%equity) / ((1.0-equmINITSS%fundlev)*equmINITSS%capital + equmINITSS%equity)
	END IF

	equmTRANS(:)%output = equmTRANS(:)%tfpadj * (equmTRANS(:)%capital**(alpha/utilelastalpha)) * (equmTRANS(:)%labor**((1.0-alpha)*(1.0+utilelast)/utilelastalpha)) + rhousing*housefrac*((1.0-equmTRANS(:)%fundlev)*equmTRANS(:)%capital + equmTRANS(:)%equity)/(1.0-housefrac)	

	
	!government budget constraint,expenditures and tax rates
	IF(fsptransition==.false.) THEN 
		IF (AdjGovBudgetConstraint==1) THEN !adjust spending
			equmTRANS(:)%govbond = equmINITSS%govbond
			equmTRANS(:)%labtax = equmINITSS%labtax
			equmTRANS(:)%lumptransfer = equmINITSS%lumptransfer		
			equmTRANS(:)%taxrev = equmTRANS(:)%labtax*equmTRANS(:)%wage*equmTRANS(:)%labor - equmTRANS(:)%lumptransfer + corptax*equmTRANS(:)%profit
			IF(TaxRetainedEarnings==1) THEN
				equmTRANS(:)%taxrev = equmTRANS(:)%taxrev + (1.0-corptax)*equmTRANS(:)%profit-equmTRANS(:)%dividend 
				equmTRANS(:)%intfirmbond = 0.0
			END IF
			equmTRANS(:)%govexp = equmTRANS(:)%taxrev + equmTRANS(:)%rb*equmINITSS%govbond

		ELSE IF(AdjGovBudgetConstraint==2) THEN  !adjust lump sum taxes
			equmTRANS(:)%govbond = equmINITSS%govbond
			equmTRANS(:)%govexp = equmINITSS%govexp
			equmTRANS(:)%labtax = equmINITSS%labtax
			equmTRANS(:)%lumptransfer = equmTRANS(:)%labtax*equmTRANS(:)%wage*equmTRANS(:)%labor + corptax*equmTRANS(:)%profit + equmTRANS(:)%rb*equmINITSS%govbond - equmTRANS(:)%govexp
			equmTRANS(:)%taxrev = equmTRANS(:)%labtax*equmTRANS(:)%wage*equmTRANS(:)%labor - equmTRANS(:)%lumptransfer + corptax*equmTRANS(:)%profit
			IF(TaxRetainedEarnings==1) THEN
				equmTRANS(:)%lumptransfer = equmTRANS(:)%lumptransfer + (1.0-corptax)*equmTRANS(:)%profit-equmTRANS(:)%dividend 
				equmTRANS(:)%intfirmbond = 0.0
			END IF
				
		ELSE IF(AdjGovBudgetConstraint==3) THEN !adjust debt
			IF(GovExpConstantFracOutput==0) equmTRANS(:)%govexp = equmINITSS%govexp
			IF(GovExpConstantFracOutput==1) equmTRANS(:)%govexp = equmTRANS(:)%output*equmINITSS%govexp/equmINITSS%output

			equmTRANS(:)%lumptransfer = equmINITSS%lumptransfer
			equmTRANS(:)%taxrev = equmTRANS(:)%labtax*equmTRANS(:)%wage*equmTRANS(:)%labor - equmTRANS(:)%lumptransfer + corptax*equmTRANS(:)%profit
			IF(TaxRetainedEarnings==1) THEN
				equmTRANS(:)%taxrev = equmTRANS(:)%taxrev + (1.0-corptax)*equmTRANS(:)%profit-equmTRANS(:)%dividend 
				equmTRANS(:)%intfirmbond = 0.0
			END IF
	
			!compute required increase in lumptransfer
	 		lrgov = equmTRANS(:)%rb
! 			lrgov = equmINITSS%rb
	
			lpvgovbc = equmFINALSS%govbond
			lpvlumpincr = 0.0
			DO it = Ttransition,1,-1
				lpvgovbc = (lpvgovbc + deltatransvec(it)*(equmTRANS(it)%govexp - equmTRANS(it)%taxrev))/(1.0+deltatransvec(it)*lrgov(it))
	 			IF(cumdeltatrans(it)>=taxincrstart) lpvlumpincr = (lpvlumpincr + deltatransvec(it))/(1.0+deltatransvec(it)*(lrgov(it)+taxincrdecay))
				IF(cumdeltatrans(it)<taxincrstart) lpvlumpincr = lpvlumpincr/(1.0+deltatransvec(it)*lrgov(it))
			END DO	
	
			linitlumpincr = (equmINITSS%govbond-lpvgovbc) / lpvlumpincr
			DO it = 1,Ttransition
				IF(cumdeltatrans(it)>=taxincrstart) equmTRANS(it)%lumptransfer = equmTRANS(it)%lumptransfer + linitlumpincr*exp(-taxincrdecay*(cumdeltatrans(it)-taxincrstart))
			END DO

			equmTRANS(:)%taxrev = equmTRANS(:)%labtax*equmTRANS(:)%wage*equmTRANS(:)%labor - equmTRANS(:)%lumptransfer + corptax*equmTRANS(:)%profit
			IF(TaxRetainedEarnings==1) equmTRANS(:)%taxrev = equmTRANS(:)%taxrev + (1.0-corptax)*equmTRANS(:)%profit-equmTRANS(:)%dividend 
			
			equmTRANS(Ttransition)%govbond = equmFINALSS%govbond
			DO it = Ttransition-1,2,-1
				equmTRANS(it)%govbond = (equmTRANS(it+1)%govbond - deltatransvec(it)*(equmTRANS(it)%taxrev-equmTRANS(it)%govexp)) / (1.0+deltatransvec(it)*lrgov(it))
			END DO
			equmTRANS(1)%govbond = equmINITSS%govbond
	
			equmTRANS(:)%lumptransfer = equmTRANS(:)%lumptransfer + (equmTRANS(:)%rb-lrgov(:))*equmTRANS(:)%govbond
			equmTRANS(:)%taxrev = equmTRANS(:)%labtax*equmTRANS(:)%wage*equmTRANS(:)%labor - equmTRANS(:)%lumptransfer + corptax*equmTRANS(:)%profit
			IF(TaxRetainedEarnings==1) equmTRANS(:)%taxrev = equmTRANS(:)%taxrev + (1.0-corptax)*equmTRANS(:)%profit-equmTRANS(:)%dividend 
	
		END IF
	ELSE IF(fsptransition==.true.) THEN	
	! 	!tax increase guess
	! 	ltaxincr = 0.003
	! 	equmTRANS(:)%labtax = equmINITSS%labtax
	! 	DO it = 1,Ttransition
	! 		IF(cumdeltatrans(it)>=fspointer%labtaxstart .and. cumdeltatrans(it)<fspointer%labtaxend) THEN
	! 			equmTRANS(it)%labtax = equmINITSS%labtax +ltaxincr
	! 		END  IF
	! 	END DO
		
	END IF

	!household bonds
	IF(RetainedEarningsInBondMarket==1) equmTRANS(:)%bond = -equmTRANS(:)%worldbond - equmTRANS(:)%govbond - equmTRANS(:)%fundbond - equmTRANS(:)%intfirmbond
	IF(RetainedEarningsInBondMarket==0) equmTRANS(:)%bond = -equmTRANS(:)%worldbond - equmTRANS(:)%govbond - equmTRANS(:)%fundbond
	
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

END SUBROUTINE IterateTransitionStickyRb