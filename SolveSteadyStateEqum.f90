SUBROUTINE SolveSteadyStateEqum

USE Parameters
USE Globals
USE Procedures

IMPLICIT NONE

! REAL(8) 	:: lwage,lrcapital,ldiffprof,lprofit
REAL(8) 	:: ldiff,llabor_Y,llabor_N,lK_totoutput_ratio,lstep,lstep_lab

OPEN(3, FILE = trim(OutputDir) // 'SteadyStateEqumIteration.txt', STATUS = 'replace'); CLOSE(3)

converged = .false.
neqmiter = 1
ldiff = 1.0
lstep 	= stepequmss
lstep_lab = 0.5

IF(initialSS == .true. .and. CalibrateDiscountRate==0) THEN !initialize grids and prices at target
	 CALL Grids
	 CALL InitialPrices
END IF


IF(OneAssetNoCapital==0) THEN
	DO WHILE (neqmiter<=maxiterequmss .and. ldiff>tolequmss )




		CALL IterateBellman
		CALL StationaryDistribution
		CALL DistributionStatistics
	
		!implied capital-output ratio
		CALL KNYfromANY(Ea,lK_totoutput_ratio,price_W)

		!implied total labor
		llabor_Y = Elabor_Y/varieties
		llabor_N = Elabor_N
	
		!check convergence
		ldiff = abs(lK_totoutput_ratio/K_totoutput_ratio - 1.0) + abs(llabor_Y/labor_Y - 1.0) + abs(llabor_N/labor_N - 1.0)
		IF (Display>=1) THEN
			write(*,"(9999(G20.6,:,','))") 'Steady state iter ',neqmiter
			write(*,"(9999(G14.6,:,','))") ' K-NY ratio:',K_totoutput_ratio,', labor_Y:',labor_Y,', labor_N:',labor_N
			write(*,"(9999(G14.6,:,','))") ' K-NY ratio:',lK_totoutput_ratio,', labor_Y:',llabor_Y,', labor_N:',llabor_N
			write(*,*) ' Err: ',ldiff
			write(*,*) ' '
		END IF
	
		!update
		K_totoutput_ratio = lstep*lK_totoutput_ratio + (1.0-lstep)*K_totoutput_ratio
		labor_Y = lstep_lab*llabor_Y + (1.0-lstep_lab)*labor_Y
		labor_N = lstep_lab*llabor_N + (1.0-lstep_lab)*labor_N
	
		!implied prices and equm objects
		capital = totoutput * K_totoutput_ratio
		rcapital = (price_W*alpha_Y*drs_Y + (grossprofit_R/output)*alpha_N*drs_N)/K_totoutput_ratio

		capital_Y = price_W*alpha_Y*drs_Y*output/rcapital
		tfp_Y = output / (((capital_Y**alpha_Y)*(labor_Y**(1.0-alpha_Y))) ** drs_Y)
		wage_Y = price_W*(1.0-alpha_Y)*drs_Y*output/labor_Y
		mc_Y = (1.0/tfp_Y)*((rcapital/alpha_Y)**alpha_Y) *((wage_Y/(1.0-alpha_Y))**(1.0-alpha_Y))

		capital_N = grossprofit_R*alpha_N*drs_N*varieties/rcapital
		tfp_N = varieties / (((capital_N**alpha_N)*(labor_N**(1.0-alpha_N))) ** drs_N)
		wage_N = grossprofit_R*(1.0-alpha_N)*drs_N*varieties/labor_N
		mc_N = (1.0/tfp_N)*((rcapital/alpha_N)**alpha_N) *((wage_N/(1.0-alpha_N))**(1.0-alpha_N))

		ra = rcapital - deprec
		dividend_A = profdistfracA*profit*(1.0-corptax)
		equity_A = dividend_A/ra
		dividend_B = profdistfracB*profit*(1.0-corptax)
		equity_B = dividend_B/rb

		netwagegrid = (1.0-labtax) * yprodgrid * ( wage_N*yoccgrid + wage_Y*(1.0-yoccgrid) )
		neqmiter = neqmiter+1 

	END DO

! ELSE IF(OneAssetNoCapital==1) THEN
!
!
! 	ldiffprof=1
! 	DO WHILE (neqmiter<=maxiterequmss .and. ldiffprof>tolequmss )
! 		IF(initialSS == .true.)CALL Grids
! 		CALL IterateBellman
! 		CALL StationaryDistribution
! 		CALL DistributionStatistics
! 		labor = Elabor
! 		lprofit = (1.0-mc)*tfp*labor
!
! 		ldiffprof= abs(lprofit/profit - 1.0)
! 		IF (Display>=1) write(*,"(A,I2,A,E11.4)") ' Steady state equm iter ',neqmiter, ', Profit error',lprofit/profit - 1.0
!
! 		profit = lprofit
! 		neqmiter = neqmiter+1
!
! 	END DO
!
! 	capital = 0.0
! 	equity = 0.0
! 	KYratio = 0.0
! 	KNratio = 0.0
END IF

bond = Eb
investment = deprec*capital
priceadjust = 0.0
bondelast = bondelastrelgdp*output

taxrev = labtax*wage_Y*labor_Y*varieties + labtax*wage_N*labor_N - lumptransfer + corptax*profit
IF(TaxHHProfitIncome == 1) taxrev = taxrev + labtax*profdistfracW*profit*(1.0-corptax)

caputil 	= 1.0
assetdrop_A = 1.0
assetdrop_B = 1.0

IF(GovBondResidualZeroWorld==0) THEN
	govbond = -ssdebttogdp*output
	govexp = taxrev + rb*govbond 
	worldbond = equity_B - bond - govbond
ELSE IF(GovBondResidualZeroWorld==1) THEN
	worldbond = 0.0
	govbond = equity_B - bond - worldbond
	govexp = taxrev + rb*govbond		
END IF


END SUBROUTINE SolveSteadyStateEqum

