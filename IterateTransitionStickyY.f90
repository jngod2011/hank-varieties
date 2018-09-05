SUBROUTINE IterateTransitionStickyY

USE Parameters
USE Globals
USE Procedures

IMPLICIT NONE

INTEGER 	:: it,ii,infix,insmooth
REAL(8) 	:: lldY,lldK,ldiffY,ldiffK,ldiffT,lstepT,ltaxincr,ltaxincr_implied,ltaxincr1,lpvgovbc,lpvtaxableinc,lpvprofit,lpvfirmrb,loutputprodFINALSS
REAL(8), DIMENSION(Ttransition) :: lcapital,loutput,lfirmdiscount,lgovbond,loutputprod,lexcessdemand,lbondmarket,lKYratio,lKYstat1,lKYstat2,lKYstat3
REAL(8), DIMENSION(Ttransition) :: lcapital1,lcapital_xold,lcapital_fold,loutput1,loutput_xold,loutput_fold

iteratingtransition = .true.


IF(fsptransition==.false.) THEN
	IF(Display>=1 .and. stickytransition==.true.) write(*,*)' Solving for sticky price transition without ZLB'
	IF(Display>=1 .and. zlbtransition==.true.) write(*,*)' Solving for sticky price transition with ZLB'
ELSE IF (fsptransition==.true.) THEN
	IF(Display>=1 .and. stickytransition==.true.) write(*,*)' Solving for sticky price transition without ZLB (with FSP)'
	IF(Display>=1 .and. zlbtransition==.true.) write(*,*)' Solving for sticky price transition with ZLB (with FSP)'
END IF


!guess capital demand and aggregate output
IF(fsptransition==.false. .and. UseFlexTransitionAsGuess==0) THEN

	!initial capital pins down labor and KN ratio in first period of transition
	equmTRANS(1)%capital = equmINITSS%capital
	!construct sequence of guesses of capital: assume log linear in capital (constant if a temporary transition)
	lldK = (log(equmFINALSS%capital)  - log(equmTRANS(1)%capital)) / real(Ttransition)
	DO it = 2,Ttransition
		equmTRANS(it)%capital = equmTRANS(1)%capital  *exp(lldK*it)
	END DO

	!construct sequence of guesses of output: assume log linear in output) (constant if a temporary transition)
	lldY = (log(equmFINALSS%output)  - log(equmINITSS%output)) / real(Ttransition+1)
	DO it = 1,Ttransition
		equmTRANS(it)%output = equmINITSS%output  *exp(lldY*it)
	END DO

ELSE IF(fsptransition==.false. .and. UseFlexTransitionAsGuess==1 ) THEN
	equmTRANS%capital = irfpointer%equmFLEX%capital
	equmTRANS%output = irfpointer%equmFLEX%output

ELSE IF	(fsptransition==.true.) THEN
	IF(stickytransition==.true.) THEN
		equmTRANS%capital = irfpointer%equmSTICKY%capital
		equmTRANS%output = irfpointer%equmSTICKY%output
	ELSE IF(zlbtransition==.true.) THEN
		equmTRANS%capital = irfpointer%equmZLB%capital
		equmTRANS%output = irfpointer%equmZLB%output
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

!labor from production function
loutputprodFINALSS  = equmFINALSS%output - rhousing*housefrac*(1.0-equmFINALSS%fundlev)*equmFINALSS%capital/(1.0-housefrac)	
loutputprod  = equmTRANS(:)%output - rhousing*housefrac*(1.0-equmTRANS(:)%fundlev)*equmTRANS(:)%capital/(1.0-housefrac)	
equmTRANS(:)%labor = (loutputprod/(equmTRANS(:)%capital**alpha)) ** (1.0/(1.0-alpha))
equmTRANS(:)%KNratio = equmTRANS(:)%capital/equmTRANS(:)%labor
equmTRANS(:)%KYratio = (equmTRANS(:)%KNratio**(1.0-alpha)) / equmTRANS(:)%tfp

!wage from GHH
equmTRANS(:)%netwage = chi * (equmTRANS(:)%labor**(1.0/frisch))
equmTRANS(:)%wage = equmTRANS(:)%netwage/(1.0-equmTRANS(:)%labtax)

!marginal costs and rental rate
equmTRANS(:)%mc = (equmTRANS(:)%KNratio ** (-alpha))*equmTRANS(:)%wage /(equmTRANS(:)%tfp*(1.0-alpha))
equmTRANS(:)%gap = equmTRANS(:)%elast*equmTRANS(:)%mc / (equmTRANS(:)%elast-1.0) - 1.0
equmTRANS(:)%rcapital = equmTRANS(:)%mc * alpha / equmTRANS(:)%KYratio

!solve phillips curve for inflation
IF (FirmDiscountRate==1) lfirmdiscount = equmTRANS(:)%rho
IF (FirmDiscountRate==2) lfirmdiscount = equmINITSS%rb

it = Ttransition
equmTRANS(it)%pi = deltatransvec(it)*(equmTRANS(it)%elast-1.0)*equmTRANS(it)%gap/theta &
					+ equmFINALSS%pi * (1.0+loutputprodFINALSS-loutputprod(it) -lfirmdiscount(it)*deltatransvec(it))
DO it = Ttransition-1,1,-1
	equmTRANS(it)%pi = deltatransvec(it)*(equmTRANS(it)%elast-1.0)*equmTRANS(it)%gap/theta &
						+ equmTRANS(it+1)%pi * (1.0+loutputprod(it+1)-loutputprod(it) -lfirmdiscount(it)*deltatransvec(it))
END DO	


!real and nominal interest rates
IF(AssumePathRealRate==0) THEN
	IF(zlbtransition==.false.) equmTRANS(:)%rb = (phitaylor-1.0)*equmTRANS(:)%pi + equmINITSS%rnom + equmTRANS(:)%mpshock
	IF(zlbtransition==.true.) equmTRANS(:)%rb = max(phitaylor*equmTRANS(:)%pi + equmINITSS%rnom + equmTRANS(:)%mpshock,0.0_8) - equmTRANS(:)%pi
	
END IF

! DO it = 1,Ttransition
! 	write(*,*) 'it, rb, pi: ',equmTRANS(it)%rb, equmTRANS(it)%pi
! 	write(*,*) 'it, mc, capital: ',equmTRANS(it)%mc, equmTRANS(it)%capital
! 	write(*,*) 'it, output,gap: ',equmTRANS(it)%output,equmTRANS(it)%gap
! 	write(*,*) 'it, wage,labor: ',equmTRANS(it)%wage,equmTRANS(it)%labor
! 	write(*,*) 'it, netwage,labor: ',equmTRANS(it)%netwage,equmTRANS(it)%labor
! END DO
! STOP

equmTRANS(:)%rnom = equmTRANS(:)%rb + equmTRANS(:)%pi
equmTRANS(:)%rborr = equmTRANS(:)%rb + equmTRANS(:)%borrwedge

!fund and world bond holdings
equmTRANS(:)%fundbond = -equmTRANS(:)%capital*equmTRANS(:)%fundlev

! IF(SolveWorldBondImplicitMethod==0)THEN
	equmTRANS(1)%worldbond = equmINITSS%worldbond
	DO it = 1,Ttransition-1
		CALL WorldBondFunction2(equmTRANS(it)%rb,equmTRANS(it+1)%worldbond,equmINITSS%worldbond,equmINITSS%rb,bondelast)
		equmTRANS(it+1)%worldbond = equmTRANS(it)%worldbond + deltatransvec(it)*bondadjust*(equmTRANS(it+1)%worldbond - equmTRANS(it)%worldbond)
	END DO	
! ELSEIF(SolveWorldBondImplicitMethod==1)THEN
! 	equmTRANS(1)%worldbond = equmINITSS%worldbond
! 	DO it = 1,Ttransition-1
! 		CALL WorldBondFunction2(equmTRANS(it)%rb,equmTRANS(it+1)%worldbond,equmINITSS%worldbond,equmINITSS%rb,bondelast)
! 		equmTRANS(it+1)%worldbond = (equmTRANS(it)%worldbond + deltatransvec(it)*bondadjust*equmTRANS(it+1)%worldbond)/(1.0+deltatransvec(it))
! 	END DO
!
! END IF
	
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


ii = 1	 
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
	ldiffK= maxval(abs(lcapital/equmTRANS(:)%capital - 1.0))
! 	lcapital = equmTRANS(:)%capital
	
	
	lbondmarket = statsTRANS(:)%Eb + equmTRANS(:)%worldbond + equmTRANS(:)%fundbond + equmTRANS(:)%govbond + equmTRANS(:)%intfirmbond
	DO it = 1,Ttransition-1
  		lexcessdemand(it) = equmTRANS(it)%rb*lbondmarket(it) - (lbondmarket(it+1)-lbondmarket(it))/deltatransvec(it)
! 		lexcessdemand(it) = equmTRANS(it)%rb*lbondmarket(it)
	END DO
! 	lexcessdemand(Ttransition) = equmTRANS(Ttransition)%rb*lbondmarket(Ttransition) + lbondmarket(Ttransition)/deltatransvec(Ttransition)
	lexcessdemand(Ttransition) = equmTRANS(Ttransition)%rb*lbondmarket(Ttransition)
	loutputprod = equmTRANS(:)%tfp * (lcapital**alpha) * (equmTRANS(:)%labor**(1.0-alpha))
	
	ldiffY = maxval(abs((loutputprod+lexcessdemand)/loutputprod - 1.0))
!
! 	laggdemand = statsTRANS(:)%Ec + ((1.0-housefrac)/(1-equmTRANS(:)%fundlev))* deprec*statsTRANS(:)%Ea &
! 		 			+ statsTRANS(:)%Erent + statsTRANS(:)%Eadjcost + equmTRANS%govexp &
! 					-statsTRANS%EbN*(equmTRANS%rborr-equmTRANS%rb)
! 					 + (equmTRANS%profit-equmTRANS%dividend) + INTEREST PAYMENTS
! 					 +equmTRANS%priceadjust
	!firm retained earnings, net exports,hh investment (need change in Ea, price adjustment costs
	
	
	IF(fsptransition==.true.) ldiffT = lpvgovbc
	IF (Display>=1) write(*,"(A,I,A)") '  Transition iter ',ii, ':'
	IF (Display>=1 .and. fsptransition==.false.) write(*,"(A,E10.3,A,E10.3,A,E10.3)") '   K err',ldiffK, ',  Y err',ldiffY

!  	write(*,*) 'lcapital, ',lcapital
! 	write(*,*) ' '
!  	write(*,*) 'lbondmarket, ',lbondmarket
! 	write(*,*) ' '
!  	write(*,*) 'lexcessdemand, ',lexcessdemand
! 	write(*,*) ' '
! 	write(*,*) 'loutputprod, ',loutputprod
! 	DO it = 1,Ttransition
! 		write(*,*) 'it, hhbond, ',statsTRANS(it)%Eb
! 		write(*,*) 'it, worldbond, ',equmTRANS(it)%worldbond
! 		write(*,*) 'it, fundbond, ',equmTRANS(it)%fundbond
! 		write(*,*) 'it, govbond, ',equmTRANS(it)%govbond
! 		write(*,*) 'it, intfirmbond, ',equmTRANS(it)%intfirmbond
! 		write(*,*) ' '
!
!
! 	END DO


	!update capital and output
	IF (ii<maxitertranssticky .and. max(ldiffK,ldiffY)>toltransition ) THEN
		CALL PartialUpdate(Ttransition-1,stepstickytransK,equmTRANS(2:Ttransition)%capital,lcapital(2:Ttransition),lcapital1(2:Ttransition))
		lcapital1(1)  = equmTRANS(1)%capital
		lcapital_xold = equmTRANS(:)%capital
		lcapital_fold = lcapital
 		equmTRANS(:)%capital = lcapital1

		loutput = loutputprod + lexcessdemand + rhousing*housefrac*(1.0-equmTRANS(:)%fundlev)*lcapital/(1.0-housefrac)

		lKYstat1 = equmTRANS(:)%output/equmTRANS(:)%labor
		lKYstat2 = loutput/equmTRANS(:)%labor
		CALL PartialUpdate(Ttransition,stepstickytransY,lKYstat1,lKYstat2,lKYstat3)
		equmTRANS(:)%labor = equmTRANS(:)%capital*(lKYstat3**(-1.0/alpha))
		equmTRANS(:)%output = (equmTRANS(:)%capital**alpha) * (equmTRANS(:)%labor**(1.0-alpha))

!  		lKYratio = 	lcapital/(loutputprod + lexcessdemand)
! 		equmTRANS(:)%output = lcapital1/lKYratio + rhousing*housefrac*(1.0-equmTRANS(:)%fundlev)*lcapital1/(1.0-housefrac)
!
!  		CALL PartialUpdate(Ttransition,stepstickytransY,equmTRANS(:)%output,loutput,loutput1)
! 		loutput_xold = equmTRANS(:)%output
! 		loutput_fold = loutput
!  		equmTRANS(:)%output = loutput1

		!impose it settles down
		IF(ImposeGuessSettlesDown==1) THEN
			infix = 2
			equmTRANS(Ttransition-infix+1:Ttransition)%output = equmFINALSS%output
			equmTRANS(Ttransition-infix+1:Ttransition)%capital = equmFINALSS%capital
			insmooth = Ttransition-10
 			CALL MASmooth(insmooth-infix+1,equmTRANS(Ttransition-insmooth:Ttransition-infix)%output, deltatransvec(Ttransition-insmooth:Ttransition-infix), 3, equmTRANS(Ttransition-insmooth:Ttransition-infix)%output)
 			CALL MASmooth(insmooth-infix+1,equmTRANS(Ttransition-insmooth:Ttransition-infix)%capital, deltatransvec(Ttransition-insmooth:Ttransition-infix), 3, equmTRANS(Ttransition-insmooth:Ttransition-infix)%capital)
			
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
		equmTRANS(:)%output = loutput
		equmTRANS(:)%capital = lcapital
		
	END IF
	

	!labor from production function
! 	loutputprodFINALSS  = equmFINALSS%output - rhousing*housefrac*(1.0-equmFINALSS%fundlev)*equmFINALSS%capital/(1.0-housefrac)
! 	loutputprod  = equmTRANS(:)%output - rhousing*housefrac*(1.0-equmTRANS(:)%fundlev)*equmTRANS(:)%capital/(1.0-housefrac)
! 	equmTRANS(:)%labor = (loutputprod/(equmTRANS(:)%capital**alpha)) ** (1.0/(1.0-alpha))
	equmTRANS(:)%KNratio = equmTRANS(:)%capital/equmTRANS(:)%labor
	equmTRANS(:)%KYratio = (equmTRANS(:)%KNratio**(1.0-alpha)) / equmTRANS(:)%tfp

	!wage from GHH
	equmTRANS(:)%netwage = chi * (equmTRANS(:)%labor**(1.0/frisch))
	equmTRANS(:)%wage = equmTRANS(:)%netwage/(1.0-equmTRANS(:)%labtax)

	!marginal costs and rental rate
	equmTRANS(:)%mc = (equmTRANS(:)%KNratio ** (-alpha))*equmTRANS(:)%wage /(equmTRANS(:)%tfp*(1.0-alpha))
	equmTRANS(:)%gap = equmTRANS(:)%elast*equmTRANS(:)%mc / (equmTRANS(:)%elast-1.0) - 1.0
	equmTRANS(:)%rcapital = equmTRANS(:)%mc * alpha / equmTRANS(:)%KYratio

	!solve phillips curve for inflation
	IF (FirmDiscountRate==1) lfirmdiscount = equmTRANS(:)%rho
	IF (FirmDiscountRate==2) lfirmdiscount = equmINITSS%rb

	it = Ttransition
	equmTRANS(it)%pi = deltatransvec(it)*(equmTRANS(it)%elast-1.0)*equmTRANS(it)%gap/theta &
						+ equmFINALSS%pi * (1.0+loutputprodFINALSS-loutputprod(it) -lfirmdiscount(it)*deltatransvec(it))
	DO it = Ttransition-1,1,-1
		equmTRANS(it)%pi = deltatransvec(it)*(equmTRANS(it)%elast-1.0)*equmTRANS(it)%gap/theta &
							+ equmTRANS(it+1)%pi * (1.0+loutputprod(it+1)-loutputprod(it) -lfirmdiscount(it)*deltatransvec(it))
	END DO	


	!real and nominal interest rates
	IF(AssumePathRealRate==0) THEN
		IF(zlbtransition==.false.) equmTRANS(:)%rb = (phitaylor-1.0)*equmTRANS(:)%pi + equmINITSS%rnom + equmTRANS(:)%mpshock
		IF(zlbtransition==.true.) equmTRANS(:)%rb = max(phitaylor*equmTRANS(:)%pi + equmINITSS%rnom + equmTRANS(:)%mpshock,0.0_8) - equmTRANS(:)%pi
	
	END IF
	

! 	DO it = 1,Ttransition
! 		write(*,*) 'it, rb, pi: ',equmTRANS(it)%rb, equmTRANS(it)%pi
! 		write(*,*) 'it, mc, capital: ',equmTRANS(it)%mc, equmTRANS(it)%capital
! 		write(*,*) 'it, output,gap: ',equmTRANS(it)%output,equmTRANS(it)%gap
! 		write(*,*) 'it, wage,labor: ',equmTRANS(it)%wage,equmTRANS(it)%labor
! 		write(*,*) 'it, netwage,outputprod: ',equmTRANS(it)%netwage,loutputprod(it)
! 		write(*,*) ' '
! 	END DO
! 	PAUSE
	
	equmTRANS(:)%rnom = equmTRANS(:)%rb + equmTRANS(:)%pi
	equmTRANS(:)%rborr = equmTRANS(:)%rb + equmTRANS(:)%borrwedge

	!fund and world bond holdings
	equmTRANS(:)%fundbond = -equmTRANS(:)%capital*equmTRANS(:)%fundlev

	equmTRANS(1)%worldbond = equmINITSS%worldbond
	DO it = 1,Ttransition-1
		CALL WorldBondFunction2(equmTRANS(it)%rb,equmTRANS(it+1)%worldbond,equmINITSS%worldbond,equmINITSS%rb,bondelast)
		equmTRANS(it+1)%worldbond = equmTRANS(it)%worldbond + deltatransvec(it)*bondadjust*(equmTRANS(it+1)%worldbond - equmTRANS(it)%worldbond)
	END DO	

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

END SUBROUTINE IterateTransitionStickyY