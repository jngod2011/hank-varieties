SUBROUTINE IterateTransition

USE Parameters
USE Globals
USE Procedures

IMPLICIT NONE

INTEGER 	:: it,ii,itfg
REAL(8) 	:: lldK,lldB,lldLN,lldLY,ldiffB,ldiffK,ldiffL,stepK,stepB,stepL
REAL(8) 	:: lminprice_W,lpvgovbc,lpvlumpincr,linitlumpincr
REAL(8), DIMENSION(Ttransition) :: lcapital,lcapital1,lbond,lbond1,lrb,lrb1,llabor_Y,llabor_Y1,llabor_N,llabor_N1
REAL(8), DIMENSION(Ttransition) :: lfirmdiscount,lworldbond,lrgov,ltransfer,lKYY,lKNN

iteratingtransition = .true.
IF(flextransition==.true.) THEN
	stepK = stepflextransK
	stepB = stepflextransB
	stepL = stepflextransL
ELSE IF((stickytransition==.true.) .or. (zlbtransition==.true.)) THEN
	stepK = stepstickytransK
	stepB = stepstickytransB
	stepL = stepstickytransL
END IF
lminprice_W = 0.01

IF(Display>=1)  write(*,*)' ';
IF(Display>=1 .and. flextransition==.true.) write(*,*)' Solving for flexible price transition'
IF(Display>=1 .and. stickytransition==.true.) write(*,*)' Solving for sticky price transition without ZLB'
IF(Display>=1 .and. zlbtransition==.true.) write(*,*)' Solving for sticky price transition with ZLB'

!forward guidance time
itfg = MINLOC(cumdeltatrans, 1, MASK = cumdeltatrans>=ForwardGuideShockQtrs)

!guess capital demand
IF((flextransition==.true.) .or. (UseFlexTransitionAsGuess==0)) THEN
	!construct sequence of guesses of capital: assume log linear in capital (constant if a temporary transition)
	equmTRANS(1)%capital = equmINITSS%capital
	lldK = (log(equmFINALSS%capital)  - log(equmTRANS(1)%capital)) / real(Ttransition)
	DO it = 2,Ttransition
		equmTRANS(it)%capital = equmTRANS(1)%capital  *exp(lldK*it)
	END DO
ELSE
	equmTRANS%capital = irfpointer%equmFLEX%capital		
END IF

!bond guess
IF (flextransition==.true. .and. UpdateFlexUsingBond==1) THEN !guess constant bond
	equmTRANS(1)%bond = equmINITSS%bond
	lldB = (log(equmFINALSS%bond)  - log(equmTRANS(1)%bond)) / real(Ttransition)
	DO it = 2,Ttransition
		equmTRANS(it)%bond = equmTRANS(1)%bond  *exp(lldB*it)
	END DO
	equmTRANS(:)%rb = equmINITSS%rb	!needed for govt budget constraint

ELSE IF	(monetaryshock==.false.) THEN !guess constant real rate
	equmTRANS(:)%rb = equmINITSS%rb		 
	
ELSE IF (monetaryshock==.true.) THEN !guess zero inflation
	equmTRANS(:)%pi = equmINITSS%pi

	IF(forwardguide==.false.) THEN
		IF(zlbtransition==.false.) equmTRANS(:)%rnom = equmINITSS%rnom +phitaylor*equmTRANS(:)%pi + equmTRANS(:)%mpshock
		IF(zlbtransition==.true.) equmTRANS(:)%rnom = max(equmINITSS%rnom +phitaylor*equmTRANS(:)%pi + equmTRANS(:)%mpshock, 0.0)
				
	ELSE IF(forwardguide==.true. ) THEN
		IF(zlbtransition==.false.) equmTRANS(1:itfg-1)%rnom = equmINITSS%rnom +phifg*equmTRANS(1:itfg-1)%pi + equmTRANS(1:itfg-1)%mpshock
		IF(zlbtransition==.true.) equmTRANS(1:itfg-1)%rnom = max(equmINITSS%rnom +phifg*equmTRANS(1:itfg-1)%pi + equmTRANS(1:itfg-1)%mpshock, 0.0)

		IF(zlbtransition==.false.) equmTRANS(itfg:Ttransition)%rnom = equmINITSS%rnom +phitaylor*equmTRANS(itfg:Ttransition)%pi + equmTRANS(itfg:Ttransition)%mpshock
		IF(zlbtransition==.true.) equmTRANS(itfg:Ttransition)%rnom = max(equmINITSS%rnom +phitaylor*equmTRANS(itfg:Ttransition)%pi + equmTRANS(itfg:Ttransition)%mpshock, 0.0)			
	END IF	
	equmTRANS(:)%rb = equmTRANS(:)%rnom - equmTRANS(:)%pi 
END IF

!guess labor
IF((flextransition==.true.) .or. (UseFlexTransitionAsGuess==0)) THEN
	lldLY = (log(equmFINALSS%labor_Y)  - log(equmINITSS%labor_Y)) / real(Ttransition+1)
	lldLN = (log(equmFINALSS%labor_N)  - log(equmINITSS%labor_N)) / real(Ttransition+1)
	DO it = 1,Ttransition
		equmTRANS(it)%labor_Y = equmINITSS%labor_Y  *exp(lldLY*it)
		equmTRANS(it)%labor_N = equmINITSS%labor_N  *exp(lldLN*it)
	END DO
ELSE
	equmTRANS%labor_Y = irfpointer%equmFLEX%labor_Y
	equmTRANS%labor_N = irfpointer%equmFLEX%labor_N
END IF

!
! !other eqm variables
! equmTRANS(:)%gap = equmTRANS(:)%elast*equmTRANS(:)%mc / (equmTRANS(:)%elast-1.0) - 1.0
! equmTRANS(:)%tfpadj = (equmTRANS(:)%tfp**((1.0+utilelast)/utilelastalpha)) * ((equmTRANS(:)%mc*alpha/equmINITSS%rcapital)**(alpha*utilelast/utilelastalpha))
! equmTRANS(:)%KNratio = equmTRANS(:)%capital/equmTRANS(:)%labor
! equmTRANS(:)%wage = equmTRANS(:)%mc*(1.0-alpha)* equmTRANS(:)%tfpadj * (equmTRANS(:)%KNratio**(alpha/utilelastalpha))
! equmTRANS(:)%netwage = (1.0-equmTRANS(:)%labtax)*equmTRANS(:)%wage
! equmTRANS(:)%caputil = ((equmTRANS(:)%mc*alpha*equmTRANS(:)%tfp/equmINITSS%rcapital) * equmTRANS(:)%KNratio**(alpha-1.0)) ** (utilelast/utilelastalpha)
! equmTRANS(:)%output = equmTRANS(:)%tfpadj * (equmTRANS(:)%capital**(alpha/utilelastalpha)) * (equmTRANS(:)%labor**((1.0-alpha)*(1.0+utilelast)/utilelastalpha))
! equmTRANS(:)%KYratio = (equmTRANS(:)%KNratio**(1.0-alpha)) / (equmTRANS(:)%tfp* (equmTRANS(:)%caputil**alpha))
! equmTRANS(:)%rcapital = ((equmINITSS%rcapital**utilelast) * equmTRANS(:)%mc * alpha / equmTRANS(:)%KYratio ) ** (1.0/(1.0+utilelast))
! IF (flextransition == .true.) equmTRANS(:)%priceadjust = 0.0
! IF (flextransition == .false.) equmTRANS(:)%priceadjust = (priceadjcost/2.0)*(equmTRANS(:)%pi**2)*equmTRANS(:)%capital/equmTRANS(:)%KYratio
! equmTRANS(:)%profit = (1.0-equmTRANS(:)%mc)*equmTRANS(:)%capital/equmTRANS(:)%KYratio - equmTRANS(:)%priceadjust
!
! equmTRANS(:)%deprec = equmINITSS%deprec + (utilelast*equmINITSS%rcapital/(1.0+ utilelast)) * ((equmTRANS(:)%rcapital/equmINITSS%rcapital)**(1.0+utilelast) -1.0)



ii = 1	 
ldiffK = 1.0
ldiffB = 1.0
ldiffL = 1.0
nblviolated = .false.
DO WHILE (ii<=maxitertranssticky .and. max(ldiffK,ldiffB,ldiffL)>toltransition .and. nblviolated == .false.)
		
	!inflation and nominal interest 
	IF (flextransition == .true.) THEN
		equmTRANS(:)%pi = (equmTRANS(:)%rb - equmINITSS%rnom - equmTRANS(:)%mpshock) / (phitaylor-1.0) !taylor rule
		equmTRANS(:)%rnom = equmTRANS(:)%rb + equmTRANS(:)%pi !fisher equn
	ELSE IF (flextransition == .false.) THEN
		IF(zlbtransition==.false.) THEN
			IF(BackwardTermInTaylorRule==0) THEN
				IF(forwardguide==.false.) THEN
					equmTRANS(:)%pi = (equmTRANS(:)%rb - equmINITSS%rnom - equmTRANS(:)%mpshock) / (phitaylor-1.0) !taylor rule
					equmTRANS(:)%rnom = equmTRANS(:)%rb + equmTRANS(:)%pi !fisher equn
			
				ELSE IF(forwardguide==.true.) THEN

					equmTRANS(1:itfg-1)%pi = (equmTRANS(1:itfg-1)%rb - equmINITSS%rnom - equmTRANS(1:itfg-1)%mpshock) / (phifg-1.0) !taylor rule
					equmTRANS(1:itfg-1)%rnom = equmTRANS(1:itfg-1)%rb + equmTRANS(1:itfg-1)%pi !fisher equn
			
					equmTRANS(itfg:Ttransition)%pi = (equmTRANS(itfg:Ttransition)%rb - equmINITSS%rnom - equmTRANS(itfg:Ttransition)%mpshock) / (phitaylor-1.0) !taylor rule
					equmTRANS(itfg:Ttransition)%rnom = equmTRANS(itfg:Ttransition)%rb + equmTRANS(itfg:Ttransition)%pi !fisher equn
				
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
	
	!price level 
	equmTRANS(1)%pricelev = 1.0
	DO it = 1,Ttransition-1
		equmTRANS(it+1)%pricelev = equmTRANS(it)%pricelev / (1.0 - deltatransvec(it)*equmTRANS(it)%pi)
	END DO
	
	!wholesale price
	IF (flextransition == .true.) THEN
		equmTRANS(:)%price_W = 1.0 -1.0/equmTRANS(:)%elast

	ELSE IF (flextransition == .false.) THEN

		!solve phillips curve backwards for marginal costs
		IF (FirmDiscountRate==1) lfirmdiscount = equmTRANS(:)%rho
		IF (FirmDiscountRate==2) lfirmdiscount = equmINITSS%rb
		IF (FirmDiscountRate==3) lfirmdiscount = equmINITSS%ra
		IF (FirmDiscountRate==4) lfirmdiscount = equmTRANS(:)%rb
		IF (FirmDiscountRate==5) lfirmdiscount = equmTRANS(:)%ra

		!final period of transition
		it = Ttransition
		equmTRANS(it)%price_W = (lfirmdiscount(it) 	- (equmFINALSS%totoutput-equmTRANS(it)%totoutput)/(equmTRANS(it)%totoutput*deltatransvec(it)) &
												 ) * equmFINALSS%pi * priceadjcost * equmTRANS(it)%varieties/ equmTRANS(it)%elast &
												+ (equmTRANS(it)%elast-1.0)/equmTRANS(it)%elast - ((equmFINALSS%pi-equmTRANS(it)%pi)/deltatransvec(it)) * priceadjcost * equmTRANS(it)%varieties / equmTRANS(it)%elast
		equmTRANS(it)%price_W = max(lminprice_W,equmTRANS(it)%price_W)

		!solve backwards
		DO it = Ttransition-1,1,-1
			equmTRANS(it)%price_W = (lfirmdiscount(it) 	- (equmTRANS(it+1)%totoutput-equmTRANS(it)%totoutput)/(equmTRANS(it)%totoutput*deltatransvec(it)) &
													 ) * equmTRANS(it+1)%pi * priceadjcost * equmTRANS(it)%varieties/ equmTRANS(it)%elast &
													+ (equmTRANS(it)%elast-1.0)/equmTRANS(it)%elast - ((equmTRANS(it+1)%pi-equmTRANS(it)%pi)/deltatransvec(it)) * priceadjcost * equmTRANS(it)%varieties / equmTRANS(it)%elast
			equmTRANS(it)%price_W = max(lminprice_W,equmTRANS(it)%price_W)

		END DO
	END IF
	
	!quantities
	equmTRANS(:)%K_totoutput_ratio = equmTRANS(:)%capital/equmTRANS(:)%totoutput
	
	equmTRANS(:)%rcapital = (equmTRANS(:)%price_W*alpha_Y*drs_Y + (1.0-equmTRANS(:)%price_W)*alpha_N*drs_N)/equmTRANS(:)%K_totoutput_ratio
	lKYY = equmTRANS(:)%price_W*alpha_Y*drs_Y/equmTRANS(:)%rcapital
	equmTRANS(:)%output = (equmTRANS(:)%tfp_Y**(1.0/(1.0-alpha_Y*drs_Y))) *((lKYY**alpha_Y)*(equmTRANS(:)%labor_Y**(1.0-alpha_Y)))**(drs_Y/(1.0-alpha_Y*drs_Y))
	lKNN = (1.0-equmTRANS(:)%price_W)*alpha_N*drs_N*equmTRANS(:)%output/equmTRANS(:)%rcapital
	equmTRANS(:)%varieties = (equmTRANS(:)%tfp_N**(1.0/(1.0-alpha_N*drs_N))) *((lKNN**alpha_N)*(equmTRANS(:)%labor_N**(1.0-alpha_N)))**(drs_N/(1.0-alpha_N*drs_N))
	equmTRANS(:)%totoutput = equmTRANS(:)%output*equmTRANS(:)%varieties

	!profits
	equmTRANS(:)%grossprofit_W = equmTRANS(:)%price_W * equmTRANS(:)%output*(1.0-drs_Y)
	equmTRANS(:)%netprofit_W = equmTRANS(:)%varieties * equmTRANS(:)%grossprofit_W
	equmTRANS(:)%grossprofit_R = (1.0-equmTRANS(:)%price_W)*equmTRANS(:)%output
	equmTRANS(:)%netprofit_R = equmTRANS(:)%varieties * (1.0-drs_N) * equmTRANS(:)%grossprofit_R
	equmTRANS(:)%profit = equmTRANS(:)%netprofit_R + equmTRANS(:)%netprofit_W
	equmTRANS(:)%dividend_A = profdistfracA * equmTRANS(:)%profit*(1.0-corptax)
	equmTRANS(:)%dividend_B = profdistfracB * equmTRANS(:)%profit*(1.0-corptax)
	
	
	!wages
	equmTRANS(:)%wage_Y = equmTRANS(:)%price_W*(1.0-alpha_Y)*drs_Y * equmTRANS(:)%output/ equmTRANS(:)%labor_Y
	equmTRANS(:)%mc_Y = (1.0/tfp_Y)*((equmTRANS(:)%rcapital/alpha_Y)**alpha_Y) *((equmTRANS(:)%wage_Y/(1.0-alpha_Y))**(1.0-alpha_Y))
	equmTRANS(:)%wage_N = equmTRANS(:)%grossprofit_R * (1.0-alpha_N)*drs_N* equmTRANS(:)%varieties / equmTRANS(:)%labor_N
	equmTRANS(:)%mc_N = (1.0/equmTRANS(:)%tfp_N)*((equmTRANS(:)%rcapital/alpha_N)**alpha_N) *((equmTRANS(:)%wage_N/(1.0-alpha_N))**(1.0-alpha_N))

	equmTRANS(:)%ra = equmTRANS(:)%rcapital - deprec
	
	!value of equity
	it = Ttransition
	equmTRANS(it)%equity_A = (equmFINALSS%equity_A + equmTRANS(it)%dividend_A*deltatransvec(it)) / (1.0+deltatransvec(it)*equmTRANS(it)%ra)		
	equmTRANS(it)%equity_B = (equmFINALSS%equity_B + equmTRANS(it)%dividend_B*deltatransvec(it)) / (1.0+deltatransvec(it)*equmTRANS(it)%rb)		
	DO it = Ttransition-1,1,-1
		equmTRANS(it)%equity_A = (equmTRANS(it+1)%equity_A + equmTRANS(it)%dividend_A*deltatransvec(it)) / (1.0+deltatransvec(it)*equmTRANS(it)%ra)		
		equmTRANS(it)%equity_B = (equmTRANS(it+1)%equity_B + equmTRANS(it)%dividend_B*deltatransvec(it)) / (1.0+deltatransvec(it)*equmTRANS(it)%rb)		
	END DO
	equmTRANS(:)%assetdrop_A = (equmTRANS(1)%capital + equmTRANS(1)%equity_A) / (equmINITSS%capital + equmINITSS%equity_A)
	equmTRANS(:)%assetdrop_B = equmTRANS(1)%bond / equmINITSS%bond

	!investment
	equmTRANS(:)%caputil = equmINITSS%caputil
	
	it = Ttransition
	equmTRANS(it)%investment = (equmFINALSS%capital-equmTRANS(it)%capital)/deltatransvec(it) + deprec*equmTRANS(it)%capital
	DO it = Ttransition-1,1,-1
		equmTRANS(it)%investment = (equmTRANS(it+1)%capital-equmTRANS(it)%capital)/deltatransvec(it) + deprec*equmTRANS(it)%capital
	END DO
	
	!government budget constraint,expenditures and tax rates
	IF (AdjGovBudgetConstraint==1) THEN !adjust spending
		equmTRANS(:)%govbond = equmINITSS%govbond * (equmTRANS(:)%pricelev**-fixnomgovdebt)
		equmTRANS(:)%labtax = equmINITSS%labtax
		equmTRANS(:)%lumptransfer = equmINITSS%lumptransfer		
		equmTRANS(:)%lumptransfer = equmTRANS(:)%lumptransfer * equmTRANS(:)%transfershock	
	
		equmTRANS(:)%taxrev = equmTRANS(:)%labtax*(equmTRANS(:)%wage_Y*equmTRANS(:)%labor_Y +equmTRANS(:)%wage_N*equmTRANS(:)%labor_N) - equmTRANS(:)%lumptransfer + corptax*equmTRANS(:)%profit
		IF(TaxHHProfitIncome == 1) equmTRANS(:)%taxrev = equmTRANS(:)%taxrev + equmTRANS(:)%labtax*profdistfracW*equmTRANS(:)%profit*(1.0-corptax)
	
		equmTRANS(:)%govexp = equmTRANS(:)%taxrev + (equmTRANS(:)%rb + fixnomgovdebt*equmTRANS(:)%pi) * equmTRANS(:)%govbond

	ELSE IF(AdjGovBudgetConstraint==2) THEN  !adjust lump sum taxes
		equmTRANS(:)%govbond = equmINITSS%govbond * (equmTRANS(:)%pricelev**-fixnomgovdebt)
		equmTRANS(:)%govexp = equmINITSS%govexp * equmTRANS(:)%govshock
		equmTRANS(:)%labtax = equmINITSS%labtax
		equmTRANS(:)%taxrev = equmTRANS(:)%govexp - (equmTRANS(:)%rb + fixnomgovdebt*equmTRANS(:)%pi) * equmTRANS(:)%govbond
		equmTRANS(:)%lumptransfer = equmTRANS(:)%labtax*(equmTRANS(:)%wage_Y*equmTRANS(:)%labor_Y +equmTRANS(:)%wage_N*equmTRANS(:)%labor_N) + corptax*equmTRANS(:)%profit +  (equmTRANS(:)%rb + fixnomgovdebt*equmTRANS(:)%pi) * equmTRANS(:)%govbond - equmTRANS(:)%govexp
		IF(TaxHHProfitIncome == 1) equmTRANS(:)%lumptransfer = equmTRANS(:)%lumptransfer + equmTRANS(:)%labtax*profdistfracW*equmTRANS(:)%profit*(1.0-corptax)

			
	ELSE IF(AdjGovBudgetConstraint==3) THEN !adjust debt
		IF(GovExpConstantFracOutput==0) equmTRANS(:)%govexp = equmINITSS%govexp * equmTRANS(:)%govshock
		IF(GovExpConstantFracOutput==1) equmTRANS(:)%govexp = equmTRANS(:)%output*equmINITSS%govexp/equmINITSS%output

		equmTRANS(:)%lumptransfer = equmINITSS%lumptransfer	 * equmTRANS(:)%transfershock	
	 		
		equmTRANS(:)%taxrev = equmTRANS(:)%labtax*(equmTRANS(:)%wage_Y*equmTRANS(:)%labor_Y +equmTRANS(:)%wage_N*equmTRANS(:)%labor_N) - equmTRANS(:)%lumptransfer + corptax*equmTRANS(:)%profit
		IF(TaxHHProfitIncome == 1) equmTRANS(:)%taxrev = equmTRANS(:)%taxrev + equmTRANS(:)%labtax*profdistfracW*equmTRANS(:)%profit*(1.0-corptax)

		!compute required increase in lumptransfer
		lrgov = equmTRANS(:)%rb
		lpvgovbc = equmFINALSS%govbond * (equmTRANS(Ttransition)%pricelev**-fixnomgovdebt)
		lpvlumpincr = 0.0
		DO it = Ttransition,1,-1
			lpvgovbc = (lpvgovbc + deltatransvec(it)*(equmTRANS(it)%govexp - equmTRANS(it)%taxrev))/(1.0+deltatransvec(it)*lrgov(it))
			IF(cumdeltatrans(it)>=taxincrstart) lpvlumpincr = (lpvlumpincr + deltatransvec(it))/(1.0+deltatransvec(it)*(lrgov(it)+taxincrdecay))
			IF(cumdeltatrans(it)<taxincrstart) lpvlumpincr = lpvlumpincr/(1.0+deltatransvec(it)*lrgov(it))
		END DO	

		IF (GovBCFDScheme==1) THEN
			linitlumpincr = (equmINITSS%govbond-lpvgovbc) / lpvlumpincr
			DO it = 1,Ttransition
				IF(cumdeltatrans(it)>=taxincrstart) equmTRANS(it)%lumptransfer = equmTRANS(it)%lumptransfer + linitlumpincr*exp(-taxincrdecay*(cumdeltatrans(it)-taxincrstart))
			END DO
		ELSE IF (GovBCFDScheme==2) THEN
			linitlumpincr = (equmINITSS%govbond-lpvgovbc) / lpvlumpincr
			ltransfer(1) = linitlumpincr
			DO it = 1,Ttransition-1
				ltransfer(it+1) = ltransfer(it) / (1.0 + deltatransvec(it)*taxincrdecay)
			END DO
			equmTRANS(:)%lumptransfer =  equmTRANS(:)%lumptransfer + ltransfer
		END IF
		
		equmTRANS(:)%taxrev = equmTRANS(:)%labtax*(equmTRANS(:)%wage_Y*equmTRANS(:)%labor_Y +equmTRANS(:)%wage_N*equmTRANS(:)%labor_N) - equmTRANS(:)%lumptransfer + corptax*equmTRANS(:)%profit
		IF(TaxHHProfitIncome == 1) equmTRANS(:)%taxrev = equmTRANS(:)%taxrev + equmTRANS(:)%labtax*profdistfracW*equmTRANS(:)%profit*(1.0-corptax)
		
		IF (GovBCFDScheme==1) THEN
			equmTRANS(Ttransition)%govbond = equmFINALSS%govbond * (equmTRANS(Ttransition)%pricelev**-fixnomgovdebt)
			DO it = Ttransition-1,2,-1
				equmTRANS(it)%govbond = (equmTRANS(it+1)%govbond - deltatransvec(it)*(equmTRANS(it)%taxrev-equmTRANS(it)%govexp)) / (1.0+deltatransvec(it)*lrgov(it))
			END DO
			equmTRANS(1)%govbond = equmINITSS%govbond
		ELSE IF (GovBCFDScheme==2) THEN
			equmTRANS(1)%govbond = equmINITSS%govbond
			DO it = 1,Ttransition-1
				equmTRANS(it+1)%govbond = (1.0+deltatransvec(it)*lrgov(it))*equmTRANS(it)%govbond + deltatransvec(it)*(equmTRANS(it)%taxrev-equmTRANS(it)%govexp) 
			END DO
		END IF

		equmTRANS(:)%lumptransfer = equmTRANS(:)%lumptransfer + (equmTRANS(:)%rb-lrgov(:))*equmTRANS(:)%govbond
		equmTRANS(:)%taxrev = equmTRANS(:)%labtax*(equmTRANS(:)%wage_Y*equmTRANS(:)%labor_Y +equmTRANS(:)%wage_N*equmTRANS(:)%labor_N) - equmTRANS(:)%lumptransfer + corptax*equmTRANS(:)%profit
		IF(TaxHHProfitIncome == 1) equmTRANS(:)%taxrev = equmTRANS(:)%taxrev + equmTRANS(:)%labtax*profdistfracW*equmTRANS(:)%profit*(1.0-corptax)
	
		
	
	
	ELSE IF(AdjGovBudgetConstraint==4) THEN  !adjust proportional tax rate
		equmTRANS(:)%govbond = equmINITSS%govbond * (equmTRANS(:)%pricelev**-fixnomgovdebt)
		equmTRANS(:)%govexp = equmINITSS%govexp * equmTRANS(:)%govshock
		equmTRANS(:)%lumptransfer = equmINITSS%lumptransfer * equmTRANS(:)%transfershock
		equmTRANS(:)%taxrev = equmTRANS(:)%govexp - (equmTRANS(:)%rb + fixnomgovdebt*equmTRANS(:)%pi)*equmTRANS(:)%govbond

		IF(TaxHHProfitIncome == 0) equmTRANS(:)%labtax  = (equmTRANS(:)%lumptransfer - corptax*equmTRANS(:)%profit - (equmTRANS(:)%rb + fixnomgovdebt*equmTRANS(:)%pi)*equmTRANS(:)%govbond + equmTRANS(:)%govexp) / (equmTRANS(:)%wage_Y*equmTRANS(:)%labor_Y + equmTRANS(:)%wage_N*equmTRANS(:)%labor_N)
		IF(TaxHHProfitIncome == 1) equmTRANS(:)%labtax  = (equmTRANS(:)%lumptransfer - corptax*equmTRANS(:)%profit - (equmTRANS(:)%rb + fixnomgovdebt*equmTRANS(:)%pi)*equmTRANS(:)%govbond + equmTRANS(:)%govexp) / (equmTRANS(:)%wage_Y*equmTRANS(:)%labor_Y + equmTRANS(:)%wage_N*equmTRANS(:)%labor_N + profdistfracW*equmTRANS(:)%profit*(1.0-corptax))

	END IF

	!household bond  or interest rate update  (because of govt bond)	
	IF (flextransition==.true. .and. UpdateFlexUsingBond==1) THEN 
		
		!world bond
		equmTRANS(:)%worldbond = -equmTRANS(:)%bond - equmTRANS(:)%govbond

		!interest rate
		it = Ttransition
		CALL WorldBondInverse2( (equmFINALSS%worldbond-equmTRANS(it)%worldbond)/(bondadjust*deltatransvec(it)) + equmTRANS(it)%worldbond, equmTRANS(it)%rb,equmINITSS%worldbond,equmINITSS%rb,bondelast)
		DO it = Ttransition-1,1,-1
			CALL WorldBondInverse2( (equmTRANS(it+1)%worldbond-equmTRANS(it)%worldbond)/(bondadjust*deltatransvec(it)) + equmTRANS(it)%worldbond, equmTRANS(it)%rb,equmINITSS%worldbond,equmINITSS%rb,bondelast)
		END DO	
		
	ELSE	
		!world bond
		equmTRANS(1)%worldbond = equmINITSS%worldbond
		DO it = 1,Ttransition-1
			CALL WorldBondFunction2( equmTRANS(it)%rb,equmTRANS(it+1)%worldbond,equmINITSS%worldbond,equmINITSS%rb,bondelast)
			equmTRANS(it+1)%worldbond = equmTRANS(it)%worldbond + bondadjust*deltatransvec(it)*(equmTRANS(it+1)%worldbond-equmTRANS(it)%worldbond)
		END DO
		
		!household bond
		equmTRANS(:)%bond = -equmTRANS(:)%worldbond - equmTRANS(:)%govbond
	END IF
	
	!borrowing rate
	IF(FixBorrowRateTransition==0) equmTRANS(:)%rborr = equmTRANS(:)%rb + equmTRANS(:)%borrwedge
	IF(FixBorrowRateTransition==1) equmTRANS(:)%rborr = equmINITSS%rb + equmTRANS(:)%borrwedge !allows for shock to borrowing wedge


	!CHECK
! 	write(*,*) ''
! 	write(*,*) ' ra ',equmINITSS%ra,equmTRANS(:)%ra
! 	write(*,*) ' rb ',equmINITSS%rb,equmTRANS(:)%rb
! 	write(*,*) ' price_W',equmINITSS%price_W,equmTRANS(:)%price_W
! 	write(*,*) ' dividend_A',equmINITSS%equity_A,equmTRANS(:)%dividend_A
! 	write(*,*) ' equity_A',equmINITSS%equity_A,equmTRANS(:)%equity_A
! 	write(*,*) ' profit',equmINITSS%profit,equmTRANS(:)%profit
! 	write(*,*) ' netprofit_R',equmINITSS%netprofit_R,equmTRANS(:)%netprofit_R
! 	write(*,*) ' netprofit_W',equmINITSS%netprofit_W,equmTRANS(:)%netprofit_W
! 	write(*,*) ' wage_Y',equmINITSS%wage_Y,equmTRANS(:)%wage_Y
! 	write(*,*) ' wage_N',equmINITSS%wage_N,equmTRANS(:)%wage_N
! 	write(*,*) ' pi',equmINITSS%pi,equmTRANS(:)%pi
! 	write(*,*) ' bond',equmINITSS%bond,equmTRANS(:)%bond
! 	write(*,*) ' capital',equmINITSS%capital,equmTRANS(:)%capital
! 	write(*,*) ' capital_N',equmINITSS%capital_N,equmTRANS(:)%capital_N
! 	write(*,*) ' capital_Y',equmINITSS%capital_N,equmTRANS(:)%capital_Y
! 	write(*,*) ' output',equmINITSS%output,equmTRANS(:)%output
! 	write(*,*) ' lump',equmINITSS%lumptransfer,equmTRANS(:)%lumptransfer
! 	write(*,*) ' labor_Y',equmINITSS%labor_Y,equmTRANS(:)%labor_Y
! 	write(*,*) ' rcapital',equmINITSS%rcapital,equmTRANS(:)%rcapital
! 	write(*,*) ' '
!
! 	STOP

	!transition
	CALL Transition
	
	!implied capital, bond, labor
	lcapital = statsTRANS(:)%Ea - equmTRANS(:)%equity_A		
	lbond = statsTRANS(:)%Eb
	llabor_Y = statsTRANS(:)%Elabor_Y/equmTRANS(:)%varieties
	llabor_N = statsTRANS(:)%Elabor_N

	
	!check convergence: needs to be fixed: 
	IF(ConvergenceRelToOutput==0) THEN
		ldiffK= maxval(abs(lcapital/equmTRANS(:)%capital - 1.0))
		ldiffB= maxval(abs(lbond/equmTRANS(:)%bond - 1.0))
	ELSEIF(ConvergenceRelToOutput==1) THEN
		ldiffK= maxval(abs(lcapital-equmTRANS(:)%capital)/equmINITSS%output)
		ldiffB= maxval(abs(lbond-equmTRANS(:)%bond)/equmINITSS%output)
	END IF
	ldiffL= max(maxval(abs(llabor_Y/equmTRANS(:)%labor_Y - 1.0)),maxval(abs(llabor_N/equmTRANS(:)%labor_N - 1.0)) )

	
	IF (Display>=1) THEN
		write(*,"(A,I,A)") '  Transition iter ',ii, ':'
		write(*,"(A,E10.3,A,E10.3,A,E10.3)") '   K err',ldiffK, ',  B err',ldiffB, ',  L err',ldiffL
		write(*,"(A,E14.7,A,E14.7)") '   rb , t=1',equmTRANS(1)%rb,'rb , t=2',equmTRANS(2)%rb
		write(*,"(A,E14.7,A,E14.7)") '   bond, t=2',lbond(2), ',  target',equmTRANS(2)%bond
		write(*,"(A,E14.7,A,E14.7)") '   capital, t=2',lcapital(2), ',  target',equmTRANS(2)%capital
		write(*,"(A,E14.7,A,E14.7)") '   labor_Y, t=1',llabor_Y(1), ',  target',equmTRANS(1)%labor_Y
		write(*,"(A,E14.7,A,E14.7)") '   labor_N, t=1',llabor_N(1), ',  target',equmTRANS(1)%labor_N
	END IF

! 	write(*,*) 'lcapital',lcapital
! 	write(*,*) 'equmTRANS(:)%capital', equmTRANS(:)%capital
! 	write(*,*) 'lbond',lbond
! 	write(*,*) 'equmTRANS(:)%bond', equmTRANS(:)%bond
!
	!CHECK
!
! 	write(*,*) '**************************'
! 	write(*,*) ' ra ',equmINITSS%ra,equmTRANS(:)%ra
! 	write(*,*) ' rb ',equmINITSS%rb,equmTRANS(:)%rb
! 	write(*,*) ' price_W',equmINITSS%price_W,equmTRANS(:)%price_W
! 	write(*,*) ' dividend_A',equmINITSS%dividend_A,equmTRANS(:)%dividend_A
! 	write(*,*) ' equity_A',equmINITSS%equity_A,equmTRANS(:)%equity_A
! 	write(*,*) ' profit',equmINITSS%profit,equmTRANS(:)%profit
! 	write(*,*) ' netprofit_R',equmINITSS%netprofit_R,equmTRANS(:)%netprofit_R
! 	write(*,*) ' netprofit_W',equmINITSS%netprofit_W,equmTRANS(:)%netprofit_W
! 	write(*,*) ' wage_Y',equmINITSS%wage_Y,equmTRANS(:)%wage_Y
! 	write(*,*) ' wage_N',equmINITSS%wage_N,equmTRANS(:)%wage_N
! 	write(*,*) ' pi',equmINITSS%pi,equmTRANS(:)%pi
! 	write(*,*) ' bond',equmINITSS%bond,equmTRANS(:)%bond
! 	write(*,*) ' capital',equmINITSS%capital,equmTRANS(:)%capital
! 	write(*,*) ' capital_N',equmINITSS%capital_N,equmTRANS(:)%capital_N
! 	write(*,*) ' capital_Y',equmINITSS%capital_N,equmTRANS(:)%capital_Y
! 	write(*,*) ' output',equmINITSS%output,equmTRANS(:)%output
! 	write(*,*) ' lump',equmINITSS%lumptransfer,equmTRANS(:)%lumptransfer
! 	write(*,*) ' labor_Y',equmINITSS%labor_Y,equmTRANS(:)%labor_Y
! 	write(*,*) ' rcapital',equmINITSS%rcapital,equmTRANS(:)%rcapital
! 	write(*,*) ' hour',equmINITSS%capital_N,equmTRANS(:)%capital_N
! 	write(*,*) ' labor_Y',equmINITSS%labor_Y,llabor_Y
! 	write(*,*) ' labor_N',equmINITSS%labor_Y,llabor_N
! 	write(*,*) ' assetdrop_A',equmTRANS(:)%assetdrop_A
! 	write(*,*) ' assetdrop_B',equmTRANS(:)%assetdrop_B
!
! 	write(*,*) ' '
! 	write(*,*) '**************************'
! 	write(*,*) ' '
! 	write(*,*) ' V INITSS',solnINITSS%V(1,:,5)
! 	write(*,*) ' V FINALSS',solnFINALSS%V(1,:,5)
! 	write(*,*) ' V TRANS T',solnTRANS(Ttransition)%V(1,:,5)

! STOP

	!updates
	IF (ii<maxitertranssticky .and. max(ldiffK,ldiffB,ldiffL)>toltransition ) THEN

		!update capital
		CALL PartialUpdate(Ttransition-1,stepK,equmTRANS(2:Ttransition)%capital,lcapital(2:Ttransition),lcapital1(2:Ttransition))
		lcapital1(1)  = equmTRANS(1)%capital
		equmTRANS(:)%capital = lcapital1
	
		!update rb using bond
		lworldbond = -lbond - equmTRANS(:)%govbond
		it = Ttransition
		CALL WorldBondInverse2( (equmFINALSS%worldbond-lworldbond(it))/(bondadjust*deltatransvec(it)) + lworldbond(it) ,lrb(it),equmINITSS%worldbond,equmINITSS%rb,bondelast)
		DO it = Ttransition-1,1,-1
			CALL WorldBondInverse2( (lworldbond(it+1)-lworldbond(it))/(bondadjust*deltatransvec(it)) + lworldbond(it) ,lrb(it),equmINITSS%worldbond,equmINITSS%rb,bondelast)
		END DO		
		CALL PartialUpdate(Ttransition,stepB,equmTRANS(:)%rb,lrb,lrb1)
		equmTRANS(:)%rb = lrb1

		!update labor
		CALL PartialUpdate(Ttransition,stepL,equmTRANS(:)%labor_Y,llabor_Y,llabor_Y1)
		CALL PartialUpdate(Ttransition,stepL,equmTRANS(:)%labor_N,llabor_N,llabor_N1)
		equmTRANS(:)%labor_Y = llabor_Y1
		equmTRANS(:)%labor_N = llabor_N1

	END IF		

	


! !
! !
! !
! ! 	!labor
! ! 	equmTRANS(:)%labor = statsTRANS(:)%Elabor
! !
! ! 	!other eqm variables
! ! 	equmTRANS(:)%gap = equmTRANS(:)%elast*equmTRANS(:)%mc / (equmTRANS(:)%elast-1.0) - 1.0
! ! 	equmTRANS(:)%tfpadj = (equmTRANS(:)%tfp**((1.0+utilelast)/utilelastalpha)) * ((equmTRANS(:)%mc*alpha/equmINITSS%rcapital)**(alpha*utilelast/utilelastalpha))
! ! 	equmTRANS(:)%KNratio = equmTRANS(:)%capital/equmTRANS(:)%labor
! ! 	equmTRANS(:)%wage = equmTRANS(:)%mc*(1.0-alpha)* equmTRANS(:)%tfpadj * (equmTRANS(:)%KNratio**(alpha/utilelastalpha))
! ! 	equmTRANS(:)%netwage = (1.0-equmTRANS(:)%labtax)*equmTRANS(:)%wage
! ! 	equmTRANS(:)%caputil = ((equmTRANS(:)%mc*alpha*equmTRANS(:)%tfp/equmINITSS%rcapital) * equmTRANS(:)%KNratio**(alpha-1.0)) ** (utilelast/utilelastalpha)
! ! 	equmTRANS(:)%output = equmTRANS(:)%tfpadj * (equmTRANS(:)%capital**(alpha/utilelastalpha)) * (equmTRANS(:)%labor**((1.0-alpha)*(1.0+utilelast)/utilelastalpha))
! ! 	equmTRANS(:)%KYratio = (equmTRANS(:)%KNratio**(1.0-alpha)) / (equmTRANS(:)%tfp* (equmTRANS(:)%caputil**alpha))
! ! 	equmTRANS(:)%rcapital = ((equmINITSS%rcapital**utilelast) * equmTRANS(:)%mc * alpha / equmTRANS(:)%KYratio ) ** (1.0/(1.0+utilelast))
! ! 	IF (flextransition == .true.) equmTRANS(:)%priceadjust = 0.0
! ! 	IF (flextransition == .false.) equmTRANS(:)%priceadjust = (priceadjcost/2.0)*(equmTRANS(:)%pi**2)*equmTRANS(:)%capital/equmTRANS(:)%KYratio
! ! 	equmTRANS(:)%profit = (1.0-equmTRANS(:)%mc)*equmTRANS(:)%capital/equmTRANS(:)%KYratio - equmTRANS(:)%priceadjust
! !
! ! 	equmTRANS(:)%deprec = equmINITSS%deprec + (utilelast*equmINITSS%rcapital/(1.0+ utilelast)) * ((equmTRANS(:)%rcapital/equmINITSS%rcapital)**(1.0+utilelast) -1.0)
! !
!
!
! 	!dividends and illiquid return
!
! 	equmTRANS(:)%ra = (equmTRANS(:)%rcapital*equmTRANS(:)%caputil - equmTRANS(:)%deprec + equmTRANS(:)%divrate - equmTRANS(:)%fundlev*equmTRANS(:)%rb) / (1.0-equmTRANS(:)%fundlev)
!
! 	!value of equity component of investmemt fund
! 	it = Ttransition
! 	equmTRANS(it)%equity = (equmFINALSS%equity + equmTRANS(it)%dividend*deltatransvec(it)) / (1.0+deltatransvec(it)*equmTRANS(it)%ra)
! 	DO it = Ttransition-1,1,-1
! 		equmTRANS(it)%equity = (equmTRANS(it+1)%equity + equmTRANS(it)%dividend*deltatransvec(it)) / (1.0+deltatransvec(it)*equmTRANS(it)%ra)
! 	END DO
! 	equmTRANS(:)%illassetdrop = ((1.0-equmTRANS(1)%fundlev)*equmTRANS(1)%capital + equmTRANS(1)%equity) / ((1.0-equmINITSS%fundlev)*equmINITSS%capital + equmINITSS%equity)
!

	
	ii = ii+1	
END DO

!run distribution stats with full
iteratingtransition = .false.
CALL Transition


IF(stickytransition==.true.) THEN
	irfpointer%equmSTICKY = equmTRANS
	irfpointer%statsSTICKY = statsTRANS
	irfpointer%solnSTICKY = solnTRANS
ELSE IF(flextransition==.true.) THEN
	irfpointer%equmFLEX = equmTRANS
	irfpointer%statsFLEX = statsTRANS
	irfpointer%solnFLEX = solnTRANS	
ELSE IF	(zlbtransition==.true.) THEN
	irfpointer%equmZLB = equmTRANS
	irfpointer%statsZLB = statsTRANS
	irfpointer%solnZLB = solnTRANS
END IF


END SUBROUTINE IterateTransition
