SUBROUTINE IterateTransition

USE Parameters
USE Globals
USE Procedures

IMPLICIT NONE

INTEGER 	:: it,ii,itfg
REAL(8) 	:: lldK,lldB,ldiffB,ldiffK,lminmargcost,lpvgovbc,lpvlumpincr,linitlumpincr,stepK,stepB
REAL(8), DIMENSION(Ttransition) :: lcapital,lcapital1,lbond,lbond1,lfirmdiscount,lrb,lrb1,lworldbond,lrgov,ltransfer

iteratingtransition = .true.
IF(flextransition==.true.) THEN
	stepK = stepflextransK
	stepB = stepflextransB
ELSE IF((stickytransition==.true.) .or. (zlbtransition==.true.)) THEN
	stepK = stepstickytransK
	stepB = stepstickytransB
END IF
lminmargcost = 0.01

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

!fund bond
equmTRANS(:)%fundbond = -equmTRANS(:)%capital*equmTRANS(:)%fundlev


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

!marginal costs
IF (flextransition == .true.) THEN
	equmTRANS(:)%mc = (equmTRANS(:)%elast-1.0)/equmTRANS(:)%elast

ELSE IF (flextransition == .false.) THEN

	!solve phillips curve backwards for marginal costs
	IF (FirmDiscountRate==1) lfirmdiscount = equmTRANS(:)%rho
	IF (FirmDiscountRate==2) lfirmdiscount = equmINITSS%rb
	IF (FirmDiscountRate==3) lfirmdiscount = equmINITSS%ra
	IF (FirmDiscountRate==4) lfirmdiscount = equmTRANS(:)%rb
	IF (FirmDiscountRate==5) lfirmdiscount = equmTRANS(:)%ra

	!final period of transition
	it = Ttransition
	equmTRANS(it)%mc = (lfirmdiscount(it) 	- (equmFINALSS%tfp-equmTRANS(it)%tfp)/(equmTRANS(it)%tfp*deltatransvec(it)) &
											- alpha*(equmFINALSS%capital-equmTRANS(it)%capital)/(equmTRANS(it)%capital*deltatransvec(it)) &
											- alpha*(equmFINALSS%caputil-equmTRANS(it)%caputil)/(equmTRANS(it)%caputil*deltatransvec(it)) &
											- (1.0-alpha)*(equmFINALSS%labor-equmTRANS(it)%labor)/(equmTRANS(it)%labor*deltatransvec(it)) ) *equmFINALSS%pi * theta/ equmTRANS(it)%elast &
											+ (equmTRANS(it)%elast-1.0)/equmTRANS(it)%elast - ((equmFINALSS%pi-equmTRANS(it)%pi)/deltatransvec(it)) * theta/ equmTRANS(it)%elast
	equmTRANS(it)%mc = max(lminmargcost,equmTRANS(it)%mc)

	!solve backwards
	DO it = Ttransition-1,1,-1
		equmTRANS(it)%mc = (lfirmdiscount(it) 	- (equmTRANS(it+1)%tfp-equmTRANS(it)%tfp)/(equmTRANS(it)%tfp*deltatransvec(it)) &
												- alpha*(equmTRANS(it+1)%capital-equmTRANS(it)%capital)/(equmTRANS(it)%capital*deltatransvec(it)) &
												- alpha*(equmTRANS(it+1)%caputil-equmTRANS(it)%caputil)/(equmTRANS(it)%caputil*deltatransvec(it)) &
												- (1.0-alpha)*(equmTRANS(it+1)%labor-equmTRANS(it)%labor)/(equmTRANS(it)%labor*deltatransvec(it)) ) *equmTRANS(it+1)%pi * theta/ equmTRANS(it)%elast &
												+ (equmTRANS(it)%elast-1.0)/equmTRANS(it)%elast - ((equmTRANS(it+1)%pi-equmTRANS(it)%pi)/deltatransvec(it)) * theta/ equmTRANS(it)%elast
		equmTRANS(it)%mc = max(lminmargcost,equmTRANS(it)%mc)

	END DO
END IF

!labor
equmTRANS(:)%labor = equmINITSS%labor

!other eqm variables
equmTRANS(:)%gap = equmTRANS(:)%elast*equmTRANS(:)%mc / (equmTRANS(:)%elast-1.0) - 1.0
equmTRANS(:)%tfpadj = (equmTRANS(:)%tfp**((1.0+utilelast)/utilelastalpha)) * ((equmTRANS(:)%mc*alpha/equmINITSS%rcapital)**(alpha*utilelast/utilelastalpha))
equmTRANS(:)%KNratio = equmTRANS(:)%capital/equmTRANS(:)%labor
equmTRANS(:)%wage = equmTRANS(:)%mc*(1.0-alpha)* equmTRANS(:)%tfpadj * (equmTRANS(:)%KNratio**(alpha/utilelastalpha))
equmTRANS(:)%netwage = (1.0-equmTRANS(:)%labtax)*equmTRANS(:)%wage
equmTRANS(:)%caputil = ((equmTRANS(:)%mc*alpha*equmTRANS(:)%tfp/equmINITSS%rcapital) * equmTRANS(:)%KNratio**(alpha-1.0)) ** (utilelast/utilelastalpha)
equmTRANS(:)%output = equmTRANS(:)%tfpadj * (equmTRANS(:)%capital**(alpha/utilelastalpha)) * (equmTRANS(:)%labor**((1.0-alpha)*(1.0+utilelast)/utilelastalpha))
equmTRANS(:)%KYratio = (equmTRANS(:)%KNratio**(1.0-alpha)) / (equmTRANS(:)%tfp* (equmTRANS(:)%caputil**alpha))
equmTRANS(:)%rcapital = ((equmINITSS%rcapital**utilelast) * equmTRANS(:)%mc * alpha / equmTRANS(:)%KYratio ) ** (1.0/(1.0+utilelast))
IF (flextransition == .true.) equmTRANS(:)%priceadjust = 0.0
IF (flextransition == .false.) equmTRANS(:)%priceadjust = (theta/2.0)*(equmTRANS(:)%pi**2)*equmTRANS(:)%capital/equmTRANS(:)%KYratio
equmTRANS(:)%profit = (1.0-equmTRANS(:)%mc)*equmTRANS(:)%capital/equmTRANS(:)%KYratio - equmTRANS(:)%priceadjust

equmTRANS(:)%deprec = equmINITSS%deprec + (utilelast*equmINITSS%rcapital/(1.0+ utilelast)) * ((equmTRANS(:)%rcapital/equmINITSS%rcapital)**(1.0+utilelast) -1.0)

!solve backward for investment
it = Ttransition
equmTRANS(it)%investment = (equmFINALSS%capital-equmTRANS(it)%capital)/deltatransvec(it) + equmTRANS(it)%deprec*equmTRANS(it)%capital
DO it = Ttransition-1,1,-1
	equmTRANS(it)%investment = (equmTRANS(it+1)%capital-equmTRANS(it)%capital)/deltatransvec(it) + equmTRANS(it)%deprec*equmTRANS(it)%capital
END DO

!dividends and illiquid return
IF(FixProfitsOutOfSteadyState==0) equmTRANS(:)%dividend 	= equmTRANS(:)%profit*(1.0-corptax)
IF(FixProfitsOutOfSteadyState==1) equmTRANS(:)%dividend 	= equmINITSS%profit*(1.0-corptax)
IF(DistributeProfitsInProportion==1) equmTRANS(:)%dividend  = profdistfrac*equmTRANS(:)%dividend 

IF(DividendFundLumpSum==1) equmTRANS(:)%divrate = 0.0
IF(DividendFundLumpSum==0) equmTRANS(:)%divrate = equmTRANS(:)%dividend/equmTRANS(:)%capital


equmTRANS(:)%ra = (equmTRANS(:)%rcapital*equmTRANS(:)%caputil - equmTRANS(:)%deprec + equmTRANS(:)%divrate - equmTRANS(:)%fundlev*equmTRANS(:)%rb) / (1.0-equmTRANS(:)%fundlev)

!value of equity component of investmemt fund
IF(DividendFundLumpSum==0) THEN
	equmTRANS(:)%equity = 0.0
	equmTRANS(:)%illassetdrop = 1.0	
ELSE IF(DividendFundLumpSum==1) THEN
	it = Ttransition
	equmTRANS(it)%equity = (equmFINALSS%equity + equmTRANS(it)%dividend*deltatransvec(it)) / (1.0+deltatransvec(it)*equmTRANS(it)%ra)		
	DO it = Ttransition-1,1,-1
		equmTRANS(it)%equity = (equmTRANS(it+1)%equity + equmTRANS(it)%dividend*deltatransvec(it)) / (1.0+deltatransvec(it)*equmTRANS(it)%ra)		
	END DO
	equmTRANS(:)%illassetdrop = ((1.0-equmTRANS(1)%fundlev)*equmTRANS(1)%capital + equmTRANS(1)%equity) / ((1.0-equmINITSS%fundlev)*equmINITSS%capital + equmINITSS%equity)
END IF

	
!government budget constraint,expenditures and tax rates
IF (AdjGovBudgetConstraint==1) THEN !adjust spending
	equmTRANS(:)%govbond = equmINITSS%govbond * (equmTRANS(:)%pricelev**-fixnomgovdebt)
	equmTRANS(:)%labtax = equmINITSS%labtax
	IF(RebateCorpTaxLumpSum==0) equmTRANS(:)%lumptransfer = equmINITSS%lumptransfer		
	IF(RebateCorpTaxLumpSum==1) THEN
		equmTRANS(:)%lumptransfer = equmINITSS%lumptransfer + corptax*(equmTRANS(:)%profit - equmINITSS%profit)
	END IF
	IF(RebateCorpTaxLumpSum==2) THEN
		equmTRANS(:)%lumptransfer = equmINITSS%lumptransfer + (1.0-labtax)*corptax*(equmTRANS(:)%profit - equmINITSS%profit)
	END IF 
	equmTRANS(:)%lumptransfer = equmTRANS(:)%lumptransfer * equmTRANS(:)%transfershock	
	
	equmTRANS(:)%taxrev = equmTRANS(:)%labtax*equmTRANS(:)%wage*equmTRANS(:)%labor - equmTRANS(:)%lumptransfer + corptax*equmTRANS(:)%profit
	IF(DistributeProfitsInProportion == 1 .and. TaxHHProfitIncome == 1) equmTRANS(:)%taxrev = equmTRANS(:)%taxrev + equmTRANS(:)%labtax*(1.0-profdistfrac)*equmTRANS(:)%profit*(1.0-corptax)
	
	equmTRANS(:)%govexp = equmTRANS(:)%taxrev + (equmTRANS(:)%rb + fixnomgovdebt*equmTRANS(:)%pi) * equmTRANS(:)%govbond

ELSE IF(AdjGovBudgetConstraint==2) THEN  !adjust lump sum taxes
	equmTRANS(:)%govbond = equmINITSS%govbond * (equmTRANS(:)%pricelev**-fixnomgovdebt)
	equmTRANS(:)%govexp = equmINITSS%govexp * equmTRANS(:)%govshock
	equmTRANS(:)%labtax = equmINITSS%labtax
	equmTRANS(:)%taxrev = equmTRANS(:)%govexp - (equmTRANS(:)%rb + fixnomgovdebt*equmTRANS(:)%pi) * equmTRANS(:)%govbond
	equmTRANS(:)%lumptransfer = equmTRANS(:)%labtax*equmTRANS(:)%wage*equmTRANS(:)%labor + corptax*equmTRANS(:)%profit +  (equmTRANS(:)%rb + fixnomgovdebt*equmTRANS(:)%pi) * equmTRANS(:)%govbond - equmTRANS(:)%govexp
	IF(DistributeProfitsInProportion == 1 .and. TaxHHProfitIncome == 1) equmTRANS(:)%lumptransfer = equmTRANS(:)%lumptransfer + equmTRANS(:)%labtax*(1.0-profdistfrac)*equmTRANS(:)%profit*(1.0-corptax)

			
ELSE IF(AdjGovBudgetConstraint==3) THEN !adjust debt
	IF(GovExpConstantFracOutput==0) equmTRANS(:)%govexp = equmINITSS%govexp * equmTRANS(:)%govshock
	IF(GovExpConstantFracOutput==1) equmTRANS(:)%govexp = equmTRANS(:)%output*equmINITSS%govexp/equmINITSS%output

	IF(RebateCorpTaxLumpSum==0) equmTRANS(:)%lumptransfer = equmINITSS%lumptransfer		
	IF(RebateCorpTaxLumpSum==1) THEN
		equmTRANS(:)%lumptransfer = equmINITSS%lumptransfer + corptax*(equmTRANS(:)%profit - equmINITSS%profit)
	END IF
	IF(RebateCorpTaxLumpSum==2) THEN
		equmTRANS(:)%lumptransfer = equmINITSS%lumptransfer + (1.0-labtax)*corptax*(equmTRANS(:)%profit - equmINITSS%profit)
	END IF
	equmTRANS(:)%lumptransfer = equmTRANS(:)%lumptransfer * equmTRANS(:)%transfershock	
	 		
	equmTRANS(:)%taxrev = equmTRANS(:)%labtax*equmTRANS(:)%wage*equmTRANS(:)%labor - equmTRANS(:)%lumptransfer + corptax*equmTRANS(:)%profit
	IF(DistributeProfitsInProportion == 1 .and. TaxHHProfitIncome == 1) equmTRANS(:)%taxrev = equmTRANS(:)%taxrev + equmTRANS(:)%labtax*(1.0-profdistfrac)*equmTRANS(:)%profit*(1.0-corptax)

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
	
	equmTRANS(:)%taxrev = equmTRANS(:)%labtax*equmTRANS(:)%wage*equmTRANS(:)%labor - equmTRANS(:)%lumptransfer + corptax*equmTRANS(:)%profit
	IF(DistributeProfitsInProportion == 1 .and. TaxHHProfitIncome == 1) equmTRANS(:)%taxrev = equmTRANS(:)%taxrev + equmTRANS(:)%labtax*(1.0-profdistfrac)*equmTRANS(:)%profit*(1.0-corptax)
	
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
	equmTRANS(:)%taxrev = equmTRANS(:)%labtax*equmTRANS(:)%wage*equmTRANS(:)%labor - equmTRANS(:)%lumptransfer + corptax*equmTRANS(:)%profit
	IF(DistributeProfitsInProportion == 1 .and. TaxHHProfitIncome == 1) equmTRANS(:)%taxrev = equmTRANS(:)%taxrev + equmTRANS(:)%labtax*(1.0-profdistfrac)*equmTRANS(:)%profit*(1.0-corptax)
	
	
	
	
ELSE IF(AdjGovBudgetConstraint==4) THEN  !adjust proportional tax rate
	equmTRANS(:)%govbond = equmINITSS%govbond * (equmTRANS(:)%pricelev**-fixnomgovdebt)
	equmTRANS(:)%govexp = equmINITSS%govexp * equmTRANS(:)%govshock
	equmTRANS(:)%lumptransfer = equmINITSS%lumptransfer * equmTRANS(:)%transfershock
	equmTRANS(:)%taxrev = equmTRANS(:)%govexp - (equmTRANS(:)%rb + fixnomgovdebt*equmTRANS(:)%pi)*equmTRANS(:)%govbond

	IF(DistributeProfitsInProportion == 0 .or. TaxHHProfitIncome == 0) equmTRANS(:)%labtax  = (equmTRANS(:)%lumptransfer - corptax*equmTRANS(:)%profit - (equmTRANS(:)%rb + fixnomgovdebt*equmTRANS(:)%pi)*equmTRANS(:)%govbond + equmTRANS(:)%govexp) / (equmTRANS(:)%wage*equmTRANS(:)%labor)
	IF(DistributeProfitsInProportion == 1 .and. TaxHHProfitIncome == 1) equmTRANS(:)%labtax  = (equmTRANS(:)%lumptransfer - corptax*equmTRANS(:)%profit - (equmTRANS(:)%rb + fixnomgovdebt*equmTRANS(:)%pi)*equmTRANS(:)%govbond + equmTRANS(:)%govexp) / (equmTRANS(:)%wage*equmTRANS(:)%labor + (1.0-profdistfrac)*equmTRANS(:)%profit*(1.0-corptax))

END IF

!household bond  or interest rate guess
IF (flextransition==.true. .and. UpdateFlexUsingBond==1) THEN !interest rate implied by constant household bond (because of govt bond)
 
	!world bond
	equmTRANS(:)%worldbond = -equmTRANS(:)%bond - equmTRANS(:)%govbond - equmTRANS(:)%fundbond

	!interest rate
	it = Ttransition
	CALL WorldBondInverse2( (equmFINALSS%worldbond-equmTRANS(it)%worldbond)/(bondadjust*deltatransvec(it)) + equmTRANS(it)%worldbond, equmTRANS(it)%rb,equmINITSS%worldbond,equmINITSS%rb,bondelast)
	DO it = Ttransition-1,1,-1
		CALL WorldBondInverse2( (equmTRANS(it+1)%worldbond-equmTRANS(it)%worldbond)/(bondadjust*deltatransvec(it)) + equmTRANS(it)%worldbond, equmTRANS(it)%rb,equmINITSS%worldbond,equmINITSS%rb,bondelast)
	END DO	
	
ELSE !household bond implied by constant interest rate
	
	!world bond from real rate
	equmTRANS(1)%worldbond = equmINITSS%worldbond
	DO it = 1,Ttransition-1
		CALL WorldBondFunction2( equmTRANS(it)%rb,equmTRANS(it+1)%worldbond,equmINITSS%worldbond,equmINITSS%rb,bondelast)
		equmTRANS(it+1)%worldbond = equmTRANS(it)%worldbond + bondadjust*deltatransvec(it)*(equmTRANS(it+1)%worldbond-equmTRANS(it)%worldbond)
	END DO
	
	!household bond
	equmTRANS(:)%bond = -equmTRANS(:)%worldbond - equmTRANS(:)%govbond - equmTRANS(:)%fundbond
END IF

!borrowing rate
IF(FixBorrowRateTransition==0) equmTRANS(:)%rborr = equmTRANS(:)%rb + equmTRANS(:)%borrwedge
IF(FixBorrowRateTransition==1) equmTRANS(:)%rborr = equmINITSS%rb + equmTRANS(:)%borrwedge !allows for shock to borrowing wedge


ii = 1	 
ldiffK = 1.0
ldiffB = 1.0
nblviolated = .false.
DO WHILE (ii<=maxitertranssticky .and. max(ldiffK,ldiffB)>toltransition .and. nblviolated == .false.)
	!solve for transtion
	CALL Transition
	
	!computed implied equilibrium quantities
	lbond = statsTRANS(:)%Eb
	lcapital = (statsTRANS(:)%Ea - equmTRANS(:)%equity)/ (1.0 - equmTRANS(:)%fundlev)
	IF(ConvergenceRelToOutput==0) THEN
		ldiffK= maxval(abs(lcapital/equmTRANS(:)%capital - 1.0))
		ldiffB= maxval(abs(lbond/equmTRANS(:)%bond - 1.0))
	ELSEIF(ConvergenceRelToOutput==1) THEN
		ldiffK= maxval(abs(lcapital-equmTRANS(:)%capital)/equmINITSS%output)
		ldiffB= maxval(abs(lbond-equmTRANS(:)%bond)/equmINITSS%output)
	END IF
	IF (Display>=1) write(*,"(A,I,A)") '  Transition iter ',ii, ':'
	IF (Display>=1) write(*,"(A,E10.3,A,E10.3,A,E10.3)") '   K err',ldiffK, ',  B err',ldiffB
	IF (Display>=1) write(*,"(A,E14.7,A,E14.7)") '   rb , t=1',equmTRANS(1)%rb,'rb , t=2',equmTRANS(2)%rb
	IF (Display>=1) write(*,"(A,E14.7,A,E14.7)") '   hh bond, t=2',lbond(2), ',  target',equmTRANS(2)%bond
	IF (Display>=1) write(*,"(A,E14.7,A,E14.7)") '   capital, t=2',lcapital(2), ',  target',equmTRANS(2)%capital
	
	IF (ii<maxitertranssticky .and. max(ldiffK,ldiffB)>toltransition ) THEN

		!update capital
		CALL PartialUpdate(Ttransition-1,stepK,equmTRANS(2:Ttransition)%capital,lcapital(2:Ttransition),lcapital1(2:Ttransition))
		lcapital1(1)  = equmTRANS(1)%capital
		equmTRANS(:)%capital = lcapital1
		
		!fund bond
		equmTRANS(:)%fundbond = -equmTRANS(:)%capital*equmTRANS(:)%fundlev
		
		
		IF((flextransition==.true.) .and. (UpdateFlexUsingBond==1)) THEN !update using B

			CALL PartialUpdate(Ttransition-1,stepB,equmTRANS(2:Ttransition)%bond,lbond(2:Ttransition),lbond1(2:Ttransition))
			lbond1(1)  = equmTRANS(1)%bond
			equmTRANS(:)%bond = lbond1

		ELSE  !update using rb

			!interest rate implied by world bond
			lworldbond = -lbond - equmTRANS(:)%govbond - equmTRANS(:)%fundbond

			it = Ttransition
			CALL WorldBondInverse2( (equmFINALSS%worldbond-lworldbond(it))/(bondadjust*deltatransvec(it)) + lworldbond(it) ,lrb(it),equmINITSS%worldbond,equmINITSS%rb,bondelast)
			DO it = Ttransition-1,1,-1
				CALL WorldBondInverse2( (lworldbond(it+1)-lworldbond(it))/(bondadjust*deltatransvec(it)) + lworldbond(it) ,lrb(it),equmINITSS%worldbond,equmINITSS%rb,bondelast)
			END DO
		
			IF((flextransition==.true.) .and. (UpdateFlexUsingBond==1)) THEN !update using B
				equmTRANS(:)%rb = lrb
			ELSE	!update using rb
				CALL PartialUpdate(Ttransition,stepB,equmTRANS(:)%rb,lrb,lrb1)
				equmTRANS(:)%rb = lrb1			
			END IF

		END IF
			
	ElSE
! 		!run distribution stats with full
		iteratingtransition = .false.
		CALL Transition
		equmTRANS(:)%capital = lcapital
		equmTRANS(:)%bond = lbond
!   		equmTRANS(:)%rb = lrb

	END IF
	
	

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
	
	!marginal costs
	IF (flextransition == .true.) THEN
		equmTRANS(:)%mc = (equmTRANS(:)%elast-1.0)/equmTRANS(:)%elast

	ELSE IF (flextransition == .false.) THEN

		!solve phillips curve backwards for marginal costs
		IF (FirmDiscountRate==1) lfirmdiscount = equmTRANS(:)%rho
		IF (FirmDiscountRate==2) lfirmdiscount = equmINITSS%rb
		IF (FirmDiscountRate==3) lfirmdiscount = equmINITSS%ra
		IF (FirmDiscountRate==4) lfirmdiscount = equmTRANS(:)%rb
		IF (FirmDiscountRate==5) lfirmdiscount = equmTRANS(:)%ra

		!final period of transition
		it = Ttransition
		equmTRANS(it)%mc = (lfirmdiscount(it) 	- (equmFINALSS%tfp-equmTRANS(it)%tfp)/(equmTRANS(it)%tfp*deltatransvec(it)) &
												- alpha*(equmFINALSS%capital-equmTRANS(it)%capital)/(equmTRANS(it)%capital*deltatransvec(it)) &
												- alpha*(equmFINALSS%caputil-equmTRANS(it)%caputil)/(equmTRANS(it)%caputil*deltatransvec(it)) &
												- (1.0-alpha)*(equmFINALSS%labor-equmTRANS(it)%labor)/(equmTRANS(it)%labor*deltatransvec(it)) ) *equmFINALSS%pi * theta/ equmTRANS(it)%elast &
												+ (equmTRANS(it)%elast-1.0)/equmTRANS(it)%elast - ((equmFINALSS%pi-equmTRANS(it)%pi)/deltatransvec(it)) * theta/ equmTRANS(it)%elast
		equmTRANS(it)%mc = max(lminmargcost,equmTRANS(it)%mc)

		!solve backwards
		DO it = Ttransition-1,1,-1
			equmTRANS(it)%mc = (lfirmdiscount(it) 	- (equmTRANS(it+1)%tfp-equmTRANS(it)%tfp)/(equmTRANS(it)%tfp*deltatransvec(it)) &
													- alpha*(equmTRANS(it+1)%capital-equmTRANS(it)%capital)/(equmTRANS(it)%capital*deltatransvec(it)) &
													- alpha*(equmTRANS(it+1)%caputil-equmTRANS(it)%caputil)/(equmTRANS(it)%caputil*deltatransvec(it)) &
													- (1.0-alpha)*(equmTRANS(it+1)%labor-equmTRANS(it)%labor)/(equmTRANS(it)%labor*deltatransvec(it)) ) *equmTRANS(it+1)%pi * theta/ equmTRANS(it)%elast &
													+ (equmTRANS(it)%elast-1.0)/equmTRANS(it)%elast - ((equmTRANS(it+1)%pi-equmTRANS(it)%pi)/deltatransvec(it)) * theta/ equmTRANS(it)%elast
			equmTRANS(it)%mc = max(lminmargcost,equmTRANS(it)%mc)

		END DO
	END IF

	!labor
	equmTRANS(:)%labor = statsTRANS(:)%Elabor

	!other eqm variables
	equmTRANS(:)%gap = equmTRANS(:)%elast*equmTRANS(:)%mc / (equmTRANS(:)%elast-1.0) - 1.0
	equmTRANS(:)%tfpadj = (equmTRANS(:)%tfp**((1.0+utilelast)/utilelastalpha)) * ((equmTRANS(:)%mc*alpha/equmINITSS%rcapital)**(alpha*utilelast/utilelastalpha))
	equmTRANS(:)%KNratio = equmTRANS(:)%capital/equmTRANS(:)%labor
	equmTRANS(:)%wage = equmTRANS(:)%mc*(1.0-alpha)* equmTRANS(:)%tfpadj * (equmTRANS(:)%KNratio**(alpha/utilelastalpha))
	equmTRANS(:)%netwage = (1.0-equmTRANS(:)%labtax)*equmTRANS(:)%wage
	equmTRANS(:)%caputil = ((equmTRANS(:)%mc*alpha*equmTRANS(:)%tfp/equmINITSS%rcapital) * equmTRANS(:)%KNratio**(alpha-1.0)) ** (utilelast/utilelastalpha)
	equmTRANS(:)%output = equmTRANS(:)%tfpadj * (equmTRANS(:)%capital**(alpha/utilelastalpha)) * (equmTRANS(:)%labor**((1.0-alpha)*(1.0+utilelast)/utilelastalpha))
	equmTRANS(:)%KYratio = (equmTRANS(:)%KNratio**(1.0-alpha)) / (equmTRANS(:)%tfp* (equmTRANS(:)%caputil**alpha))
	equmTRANS(:)%rcapital = ((equmINITSS%rcapital**utilelast) * equmTRANS(:)%mc * alpha / equmTRANS(:)%KYratio ) ** (1.0/(1.0+utilelast))
	IF (flextransition == .true.) equmTRANS(:)%priceadjust = 0.0
	IF (flextransition == .false.) equmTRANS(:)%priceadjust = (theta/2.0)*(equmTRANS(:)%pi**2)*equmTRANS(:)%capital/equmTRANS(:)%KYratio
	equmTRANS(:)%profit = (1.0-equmTRANS(:)%mc)*equmTRANS(:)%capital/equmTRANS(:)%KYratio - equmTRANS(:)%priceadjust

	equmTRANS(:)%deprec = equmINITSS%deprec + (utilelast*equmINITSS%rcapital/(1.0+ utilelast)) * ((equmTRANS(:)%rcapital/equmINITSS%rcapital)**(1.0+utilelast) -1.0)

	!solve backward for investment
	it = Ttransition
	equmTRANS(it)%investment = (equmFINALSS%capital-equmTRANS(it)%capital)/deltatransvec(it) + equmTRANS(it)%deprec*equmTRANS(it)%capital
	DO it = Ttransition-1,1,-1
		equmTRANS(it)%investment = (equmTRANS(it+1)%capital-equmTRANS(it)%capital)/deltatransvec(it) + equmTRANS(it)%deprec*equmTRANS(it)%capital
	END DO

	!dividends and illiquid return
	IF(FixProfitsOutOfSteadyState==0) equmTRANS(:)%dividend 	= equmTRANS(:)%profit*(1.0-corptax)
	IF(FixProfitsOutOfSteadyState==1) equmTRANS(:)%dividend 	= equmINITSS%profit*(1.0-corptax)
	IF(DistributeProfitsInProportion==1) equmTRANS(:)%dividend  = profdistfrac*equmTRANS(:)%dividend

	IF(DividendFundLumpSum==1) equmTRANS(:)%divrate = 0.0
	IF(DividendFundLumpSum==0) equmTRANS(:)%divrate = equmTRANS(:)%dividend/equmTRANS(:)%capital

	equmTRANS(:)%ra = (equmTRANS(:)%rcapital*equmTRANS(:)%caputil - equmTRANS(:)%deprec + equmTRANS(:)%divrate - equmTRANS(:)%fundlev*equmTRANS(:)%rb) / (1.0-equmTRANS(:)%fundlev)

	!value of equity component of investmemt fund
	IF(DividendFundLumpSum==0) THEN
		equmTRANS(:)%equity = 0.0
		equmTRANS(:)%illassetdrop = 1.0	
	ELSE IF(DividendFundLumpSum==1) THEN
		it = Ttransition
		equmTRANS(it)%equity = (equmFINALSS%equity + equmTRANS(it)%dividend*deltatransvec(it)) / (1.0+deltatransvec(it)*equmTRANS(it)%ra)		
		DO it = Ttransition-1,1,-1
			equmTRANS(it)%equity = (equmTRANS(it+1)%equity + equmTRANS(it)%dividend*deltatransvec(it)) / (1.0+deltatransvec(it)*equmTRANS(it)%ra)		
		END DO
		equmTRANS(:)%illassetdrop = ((1.0-equmTRANS(1)%fundlev)*equmTRANS(1)%capital + equmTRANS(1)%equity) / ((1.0-equmINITSS%fundlev)*equmINITSS%capital + equmINITSS%equity)
	END IF

	!government budget constraint,expenditures and tax rates
	IF (AdjGovBudgetConstraint==1) THEN !adjust spending
		equmTRANS(:)%govbond = equmINITSS%govbond * (equmTRANS(:)%pricelev**-fixnomgovdebt)
		equmTRANS(:)%labtax = equmINITSS%labtax
		IF(RebateCorpTaxLumpSum==0) equmTRANS(:)%lumptransfer = equmINITSS%lumptransfer		
		IF(RebateCorpTaxLumpSum==1) THEN
			equmTRANS(:)%lumptransfer = equmINITSS%lumptransfer + corptax*(equmTRANS(:)%profit - equmINITSS%profit)
		END IF
		IF(RebateCorpTaxLumpSum==2) THEN
			equmTRANS(:)%lumptransfer = equmINITSS%lumptransfer + (1.0-labtax)*corptax*(equmTRANS(:)%profit - equmINITSS%profit)
		END IF 
		equmTRANS(:)%lumptransfer = equmTRANS(:)%lumptransfer * equmTRANS(:)%transfershock	
	
		equmTRANS(:)%taxrev = equmTRANS(:)%labtax*equmTRANS(:)%wage*equmTRANS(:)%labor - equmTRANS(:)%lumptransfer + corptax*equmTRANS(:)%profit
		IF(DistributeProfitsInProportion == 1 .and. TaxHHProfitIncome == 1) equmTRANS(:)%taxrev = equmTRANS(:)%taxrev + equmTRANS(:)%labtax*(1.0-profdistfrac)*equmTRANS(:)%profit*(1.0-corptax)
	
		equmTRANS(:)%govexp = equmTRANS(:)%taxrev + (equmTRANS(:)%rb + fixnomgovdebt*equmTRANS(:)%pi) * equmTRANS(:)%govbond

	ELSE IF(AdjGovBudgetConstraint==2) THEN  !adjust lump sum taxes
		equmTRANS(:)%govbond = equmINITSS%govbond * (equmTRANS(:)%pricelev**-fixnomgovdebt)
		equmTRANS(:)%govexp = equmINITSS%govexp * equmTRANS(:)%govshock
		equmTRANS(:)%labtax = equmINITSS%labtax
		equmTRANS(:)%taxrev = equmTRANS(:)%govexp - (equmTRANS(:)%rb + fixnomgovdebt*equmTRANS(:)%pi) * equmTRANS(:)%govbond
		equmTRANS(:)%lumptransfer = equmTRANS(:)%labtax*equmTRANS(:)%wage*equmTRANS(:)%labor + corptax*equmTRANS(:)%profit +  (equmTRANS(:)%rb + fixnomgovdebt*equmTRANS(:)%pi) * equmTRANS(:)%govbond - equmTRANS(:)%govexp
		IF(DistributeProfitsInProportion == 1 .and. TaxHHProfitIncome == 1) equmTRANS(:)%lumptransfer = equmTRANS(:)%lumptransfer + equmTRANS(:)%labtax*(1.0-profdistfrac)*equmTRANS(:)%profit*(1.0-corptax)

			
	ELSE IF(AdjGovBudgetConstraint==3) THEN !adjust debt
		IF(GovExpConstantFracOutput==0) equmTRANS(:)%govexp = equmINITSS%govexp * equmTRANS(:)%govshock
		IF(GovExpConstantFracOutput==1) equmTRANS(:)%govexp = equmTRANS(:)%output*equmINITSS%govexp/equmINITSS%output

		IF(RebateCorpTaxLumpSum==0) equmTRANS(:)%lumptransfer = equmINITSS%lumptransfer		
		IF(RebateCorpTaxLumpSum==1) THEN
			equmTRANS(:)%lumptransfer = equmINITSS%lumptransfer + corptax*(equmTRANS(:)%profit - equmINITSS%profit)
		END IF
		IF(RebateCorpTaxLumpSum==2) THEN
			equmTRANS(:)%lumptransfer = equmINITSS%lumptransfer + (1.0-labtax)*corptax*(equmTRANS(:)%profit - equmINITSS%profit)
		END IF
		equmTRANS(:)%lumptransfer = equmTRANS(:)%lumptransfer * equmTRANS(:)%transfershock	
	 		
		equmTRANS(:)%taxrev = equmTRANS(:)%labtax*equmTRANS(:)%wage*equmTRANS(:)%labor - equmTRANS(:)%lumptransfer + corptax*equmTRANS(:)%profit
		IF(DistributeProfitsInProportion == 1 .and. TaxHHProfitIncome == 1) equmTRANS(:)%taxrev = equmTRANS(:)%taxrev + equmTRANS(:)%labtax*(1.0-profdistfrac)*equmTRANS(:)%profit*(1.0-corptax)

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
		
		equmTRANS(:)%taxrev = equmTRANS(:)%labtax*equmTRANS(:)%wage*equmTRANS(:)%labor - equmTRANS(:)%lumptransfer + corptax*equmTRANS(:)%profit
		IF(DistributeProfitsInProportion == 1 .and. TaxHHProfitIncome == 1) equmTRANS(:)%taxrev = equmTRANS(:)%taxrev + equmTRANS(:)%labtax*(1.0-profdistfrac)*equmTRANS(:)%profit*(1.0-corptax)
		
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
		equmTRANS(:)%taxrev = equmTRANS(:)%labtax*equmTRANS(:)%wage*equmTRANS(:)%labor - equmTRANS(:)%lumptransfer + corptax*equmTRANS(:)%profit
		IF(DistributeProfitsInProportion == 1 .and. TaxHHProfitIncome == 1) equmTRANS(:)%taxrev = equmTRANS(:)%taxrev + equmTRANS(:)%labtax*(1.0-profdistfrac)*equmTRANS(:)%profit*(1.0-corptax)
	
	
	
	
	ELSE IF(AdjGovBudgetConstraint==4) THEN  !adjust proportional tax rate
		equmTRANS(:)%govbond = equmINITSS%govbond * (equmTRANS(:)%pricelev**-fixnomgovdebt)
		equmTRANS(:)%govexp = equmINITSS%govexp * equmTRANS(:)%govshock
		equmTRANS(:)%lumptransfer = equmINITSS%lumptransfer * equmTRANS(:)%transfershock
		equmTRANS(:)%taxrev = equmTRANS(:)%govexp - (equmTRANS(:)%rb + fixnomgovdebt*equmTRANS(:)%pi)*equmTRANS(:)%govbond

		IF(DistributeProfitsInProportion == 0 .or. TaxHHProfitIncome == 0) equmTRANS(:)%labtax  = (equmTRANS(:)%lumptransfer - corptax*equmTRANS(:)%profit - (equmTRANS(:)%rb + fixnomgovdebt*equmTRANS(:)%pi)*equmTRANS(:)%govbond + equmTRANS(:)%govexp) / (equmTRANS(:)%wage*equmTRANS(:)%labor)
		IF(DistributeProfitsInProportion == 1 .and. TaxHHProfitIncome == 1) equmTRANS(:)%labtax  = (equmTRANS(:)%lumptransfer - corptax*equmTRANS(:)%profit - (equmTRANS(:)%rb + fixnomgovdebt*equmTRANS(:)%pi)*equmTRANS(:)%govbond + equmTRANS(:)%govexp) / (equmTRANS(:)%wage*equmTRANS(:)%labor + (1.0-profdistfrac)*equmTRANS(:)%profit*(1.0-corptax))

	END IF


	!household bond  or interest rate update  (because of govt bond)	
	IF (flextransition==.true. .and. UpdateFlexUsingBond==1) THEN 
		
		!world bond
		equmTRANS(:)%worldbond = -equmTRANS(:)%bond - equmTRANS(:)%govbond - equmTRANS(:)%fundbond

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
		equmTRANS(:)%bond = -equmTRANS(:)%worldbond - equmTRANS(:)%govbond - equmTRANS(:)%fundbond
	END IF
	
	!borrowing rate
	IF(FixBorrowRateTransition==0) equmTRANS(:)%rborr = equmTRANS(:)%rb + equmTRANS(:)%borrwedge
	IF(FixBorrowRateTransition==1) equmTRANS(:)%rborr = equmINITSS%rb + equmTRANS(:)%borrwedge !allows for shock to borrowing wedge
	
	
	ii = ii+1	
END DO

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
