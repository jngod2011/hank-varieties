REAL(8) FUNCTION FnDiscountRate(lrhoT)

USE Parameters
USE Globals
USE Procedures

IMPLICIT NONE

REAL(8), INTENT(IN) :: lrhoT
REAL(8)				:: lstep_lab,llabor_Y,llabor_N

lstep_lab = 0.5

rho = -log(logistic(lrhoT))

CALL IterateBellman
CALL StationaryDistribution
CALL DistributionStatistics


!implied capital-output ratio
CALL KNYfromANY(Ea,K_totoutput_ratio,price_W)

!implied labor
llabor_Y = Elabor_Y/varieties
llabor_N = Elabor_N

!check convergence
FnDiscountRate = K_totoutput_ratio/target_K_totoutput_ratio - 1.0

IF (Display>=1) THEN
	write(*,"(9999(G20.6,:,','))") 'Rho calib iter: ',nrhoiter,' rho:',rho
	write(*,"(9999(G14.6,:,','))") ' K-NY target :',target_K_totoutput_ratio,' labor_Y:',labor_Y,' labor_N:',labor_N
	write(*,"(9999(G14.6,:,','))") ' K-NY current:',K_totoutput_ratio,' impl labor_Y:',llabor_Y,' impl labor_N:',llabor_N
	write(*,*) ' Err: ',FnDiscountRate
	write(*,*) ' '
END IF

IF (nrhoiter<=2) THEN
	labor_Y = hourtarget*meanlabeff*DOT_PRODUCT(1.0-occgrid,occdist)/varieties
	labor_N = hourtarget*meanlabeff*DOT_PRODUCT(occgrid,occdist)
ELSE
	labor_Y = lstep_lab*llabor_Y + (1.0-lstep_lab)*labor_Y
	labor_N = lstep_lab*llabor_N + (1.0-lstep_lab)*labor_N
END IF

!update other prices
capital = totoutput * target_K_totoutput_ratio
rcapital = (price_W*alpha_Y*drs_Y + (1.0-price_W)*alpha_N*drs_N)/target_K_totoutput_ratio

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

! IF(exploring==.true.) write(4,"(A,I2,A,E11.4,A,E11.4)") '   Rho iter ',neqmiter, ', rho ',rho, ', K/N err',FnDiscountRate


!
! IF(CalibrateDiscountRate==1) OPEN(3, FILE = trim(OutputDir) // 'DiscountRateCalibration.txt', ACCESS = 'append')
! IF(CalibrateRhoAtInitialGuess==1) OPEN(3, FILE = trim(OutputDir) // 'DiscountRateAtInitialGuess.txt', ACCESS = 'append')
! WRITE(3,*) '*******************************************'
! WRITE(3,*) ' ITERATION : 			',neqmiter
! WRITE(3,*) ' rho guess: 			',rho
! WRITE(3,*) '  target KY ratio: ',KYratio
! WRITE(3,*) '  implied KY ratio: ',lKYratio
! WRITE(3,*) '  target KN ratio: ',KNratio
! WRITE(3,*) '  implied KN ratio: ',lKNratio
! WRITE(3,*) '  relative error: ',FnDiscountRate
! CLOSE(3)

nrhoiter = nrhoiter+1
END FUNCTION FnDiscountRate
