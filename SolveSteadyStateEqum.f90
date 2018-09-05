SUBROUTINE SolveSteadyStateEqum

USE Parameters
USE Globals
USE Procedures

IMPLICIT NONE

REAL(8) 	:: ldiffKN, lKNratio,lKYratio,lwage,lrcapital,lstepKN

OPEN(3, FILE = trim(OutputDir) // 'SteadyStateEqumIteration.txt', STATUS = 'replace'); CLOSE(3)

converged = .false.
neqmiter = 1
IF(calibrating == .false.) lstepKN = stepequmss
IF(calibrating == .true. ) lstepKN = 0.2*stepequmss
ldiffKN = 1.0

DO WHILE (neqmiter<=maxiterequmss .and. ldiffKN>tolequmss )

	IF (Display>=2) THEN
		WRITE(*,*) '*******************************************'
		WRITE(*,*) ' ITERATION : 			',neqmiter
		WRITE(*,*) ' r guess: 			',rcapital
		WRITE(*,*) '  implied ra: 		',ra
		WRITE(*,*) '  implied borr rate: ',rborr
		WRITE(*,*) '  implied wage: ',wage
		WRITE(*,*) '  implied KY ratio: ',KYratio
		WRITE(*,*) '  implied KN firm: 	',KNratio
	END IF

	IF(initialSS == .true.) CALL Grids
	CALL IterateBellman
	CALL StationaryDistribution
	CALL DistributionStatistics



	IF(DividendFundLumpSum==0) capital = (1.0-housefrac)*Ea/(1.0-fundlev)
	IF(DividendFundLumpSum==1) capital = (1.0-housefrac)*Ea/(1.0-fundlev + (1.0-mc-operatecost)/(ra*KYratio))
	bond = Eb
	labor = (netwage/chi)**frisch
	investment = deprec*capital
	lKNratio = capital / labor
	lKYratio = (lKNratio**(1.0-alpha)) / tfp
	lwage = mc*(1.0-alpha)* tfp * (lKNratio**alpha)
	lrcapital = mc*alpha/lKYratio

	ldiffKN= abs(lKNratio/KNratio - 1.0)
	IF (Display>=1) write(*,"(A,I2,A,E11.4)") ' Steady state equm iter ',neqmiter, ', K/N error',lKNratio/KNratio - 1.0

	OPEN(3, FILE = trim(OutputDir) // 'SteadyStateEqumIteration.txt', ACCESS = 'append')
	WRITE(3,*) '*******************************************'
	WRITE(3,*) ' ITERATION : 			',neqmiter
	WRITE(3,*) ' r guess: 	',rcapital
	WRITE(3,*) '  implied ra: 	',ra
	WRITE(3,*) '  implied borr rate: ',rborr	
	WRITE(3,*) '  implied wage: ',wage
	WRITE(3,*) '  implied KY ratio: ',KYratio
	WRITE(3,*) '  actual KY ratio: ',lKYratio
	WRITE(3,*) '  implied KN firm: 	',KNratio
	WRITE(3,*) '  implied KN hh: 	',lKNratio
	WRITE(3,*) '  relative error: ',lKNratio/KNratio - 1.0
	CLOSE(3)
		
	!update KN ratio
	IF (neqmiter<=maxiterequmss .and. ldiffKN>tolequmss ) THEN
		KNratio = (1.0-lstepKN)*KNratio +lstepKN*lKNratio
	ELSE
		KNratio = lKNratio
	END IF

	KYratio = (KNratio**(1.0-alpha)) / tfp
	rcapital = mc * alpha / KYratio
	wage = mc*(1.0-alpha)* tfp * (KNratio**alpha)
	netwage = (1.0-labtax)*wage
	rborr = rb + borrwedge
	labor = (netwage/chi)**frisch
	profit = (1.0-mc-operatecost)*capital/KYratio - priceadjust
	dividend = profit*(1.0-corptax)
	IF(DividendFundLumpSum==0) divrate = dividend/capital
	IF(DividendFundLumpSum==1) divrate = 0.0
	ra = (rcapital - deprec + divrate - fundlev*rb)/(1.0-fundlev)
	IF(DividendFundLumpSum==0) equity = 0.0
	IF(DividendFundLumpSum==1)equity = profit/ra
	output = tfp*(capital**alpha)*(labor**(1.0-alpha)) + rhousing*housefrac*((1.0-fundlev)*capital+equity)/(1.0-housefrac)
	
	taxrev = labtax*wage*labor - lumptransfer + corptax*profit
	fundbond = -capital*fundlev
	IF(GovBondResidualZeroWorld==0) THEN
		govbond = -ssdebttogdp*output
		govexp = taxrev + rb*govbond 
		worldbond = -bond-govbond-fundbond
	ELSE IF(GovBondResidualZeroWorld==1) THEN
		worldbond = 0.0
		govbond = -bond-worldbond-fundbond
		govexp = taxrev + rb*govbond		
	END IF
	bondelast = bondelastrelgdp*output	
	priceadjust = 0.0
	intfirmbond = 0.0
	
	
	directdepmax = directdepmaxfrac*output
	directdepmin = directdepminfrac*output
	
	neqmiter = neqmiter+1 

END DO

END SUBROUTINE SolveSteadyStateEqum
