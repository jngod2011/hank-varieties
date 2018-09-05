SUBROUTINE IterateTransitionFlex

USE Parameters
USE Globals
USE Procedures

IMPLICIT NONE

INTEGER 	:: it,ii
REAL(8) 	:: lldB,lldK,ldiffB,ldiffK,ldiffT,lstepT,lqa,lqb,lqc,lminmargcost,ltaxincr,ltaxincr_implied,ltaxincr1,lpvgovbc,lpvtaxableinc
REAL(8), DIMENSION(Ttransition) :: lcapital,lbond,lfirmdiscount,lgovbond
REAL(8), DIMENSION(Ttransition) :: lcapital1,lcapital_xold,lcapital_fold,lbond1,lbond_xold,lbond_fold

iteratingtransition = .true.

lminmargcost = 0.01

IF(fsptransition==.false.) write(*,*)' Solving for flexible price transition'
IF (fsptransition==.true.) write(*,*)' Solving for flexible price transition (with FSP)'


!guess capital and bond demands
IF((fsptransition==.false.) .and. (UseFlexTransitionAsGuess==0 .or. flextransition==.true.)) THEN

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

ELSE IF	(fsptransition==.true.) THEN
	equmTRANS%capital = irfpointer%equmFLEX%capital
	equmTRANS%bond = irfpointer%equmFLEX%bond
END IF

!government expenditure path (in baseline these are implied by tax revenues)
IF	(fsptransition==.true.) THEN
	equmTRANS%govexp = irfpointer%equmFLEX%govexp		
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

! write(*,*)'equmTRANS%rb',equmTRANS(:)%rb

equmTRANS(:)%rborr = equmTRANS(:)%rb + equmTRANS(:)%borrwedge


!inflation and nominal interest rates
equmTRANS(:)%pi = (equmTRANS(:)%rb - equmINITSS%rnom - equmTRANS(:)%mpshock) / (phitaylor-1.0) !taylor rule
equmTRANS(:)%rnom = equmTRANS(:)%rb + equmTRANS(:)%pi !fisher equn

equmTRANS(:)%mc = equmINITSS%mc
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
equmTRANS(:)%priceadjust = 0.0
 
!solve backward for investment
it = Ttransition
equmTRANS(it)%investment = (equmFINALSS%capital-equmTRANS(it)%capital)/deltatransvec(it) + deprec*equmTRANS(it)%capital
DO it = Ttransition-1,1,-1
	equmTRANS(it)%investment = (equmTRANS(it+1)%capital-equmTRANS(it)%capital)/deltatransvec(it) + deprec*equmTRANS(it)%capital
END DO

IF(FixProfitsOutOfSteadyState==0) THEN
	equmTRANS(:)%ra = (equmTRANS(:)%rcapital - deprec + (1.0-equmTRANS(:)%mc)/equmTRANS(:)%KYratio - equmTRANS(:)%priceadjust/equmTRANS(:)%capital - equmTRANS(:)%fundlev*equmTRANS(:)%rb) / (1.0-equmTRANS(:)%fundlev)
ELSE IF(FixProfitsOutOfSteadyState==1) THEN
	equmTRANS(:)%ra = (equmTRANS(:)%rcapital - deprec + (1.0-equmINITSS%mc)*equmINITSS%capital/(equmINITSS%KYratio*equmTRANS(:)%capital) - equmTRANS(:)%fundlev*equmTRANS(:)%rb) / (1.0-equmTRANS(:)%fundlev)
END IF

ii = 1
ldiffK = 1.0
ldiffB = 1.0

DO WHILE (ii<=maxitertransflex .and. max(ldiffK,ldiffB)>toltransition )
	!solve for transtion
	CALL Transition
	
	!compute implied govt bond demand
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
! 	IF (Display>=1 .and. fsptransition==.true.) write(*,"(A,E10.3,A,E10.3,A,E10.3)") '   K err',ldiffK, ',  B err',ldiffB, ',  PV gbc',ldiffT
	
	!update capital, labor and bond supply
	IF (ii<maxitertransflex .and. max(ldiffK,ldiffB)>toltransition ) THEN
		CALL PartialUpdate(Ttransition-1,stepflextransK,equmTRANS(2:Ttransition)%capital,lcapital(2:Ttransition),lcapital1(2:Ttransition))
		CALL PartialUpdate(Ttransition-1,stepflextransB,equmTRANS(2:Ttransition)%bond,lbond(2:Ttransition),lbond1(2:Ttransition))
				
		lcapital1(1)  = equmTRANS(1)%capital
		lbond1(1)  = equmTRANS(1)%bond
		
		lcapital_xold = equmTRANS(:)%capital
		lcapital_fold = lcapital
		equmTRANS(:)%capital = lcapital1

! write(*,*)'equmTRANS%bond',equmTRANS(:)%bond

		lbond_xold = equmTRANS(:)%bond
		lbond_fold = lbond
		equmTRANS(:)%bond = lbond1
! write(*,*)'equmTRANS%bond',equmTRANS(:)%bond

		!impose it settles down
! 		equmTRANS(Ttransition)%bond = equmFINALSS%bond
! 		equmTRANS(Ttransition)%capital = equmFINALSS%capital
		
		!update labor tax rate
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
		equmTRANS(:)%bond = lbond
		equmTRANS(:)%capital = lcapital
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
	equmTRANS(:)%pi = (equmTRANS(:)%rb - equmINITSS%rnom - equmTRANS(:)%mpshock) / (phitaylor-1.0) !taylor rule
	equmTRANS(:)%rnom = equmTRANS(:)%rb + equmTRANS(:)%pi !fisher equn

	equmTRANS(:)%mc = equmINITSS%mc
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
	equmTRANS(:)%priceadjust = 0.0


	!solve backward for investment
	it = Ttransition
	equmTRANS(it)%investment = (equmFINALSS%capital-equmTRANS(it)%capital)/deltatransvec(it) + deprec*equmTRANS(it)%capital
	DO it = Ttransition-1,1,-1
		equmTRANS(it)%investment = (equmTRANS(it+1)%capital-equmTRANS(it)%capital)/deltatransvec(it) + deprec*equmTRANS(it)%capital
	END DO

	IF(FixProfitsOutOfSteadyState==0) THEN
		equmTRANS(:)%ra = (equmTRANS(:)%rcapital - deprec + (1.0-equmTRANS(:)%mc)/equmTRANS(:)%KYratio - equmTRANS(:)%priceadjust/equmTRANS(:)%capital - equmTRANS(:)%fundlev*equmTRANS(:)%rb) / (1.0-equmTRANS(:)%fundlev)
	ELSE IF(FixProfitsOutOfSteadyState==1) THEN
		equmTRANS(:)%ra = (equmTRANS(:)%rcapital - deprec + (1.0-equmINITSS%mc)*equmINITSS%capital/(equmINITSS%KYratio*equmTRANS(:)%capital) - equmTRANS(:)%fundlev*equmTRANS(:)%rb) / (1.0-equmTRANS(:)%fundlev)
	END IF
	
 	ii = ii+1	
END DO

! write(*,*)'equmTRANS%bond',equmTRANS(:)%bond



IF(fsptransition==.false.) THEN 
	irfpointer%equmFLEX = equmTRANS
	irfpointer%statsFLEX = statsTRANS
	irfpointer%solnFLEX = solnTRANS
ELSE IF(fsptransition==.true.) THEN 
	irfpointer_fs%equmFLEX = equmTRANS
	irfpointer_fs%statsFLEX = statsTRANS
	irfpointer_fs%solnFLEX = solnTRANS
END IF

END SUBROUTINE IterateTransitionFlex