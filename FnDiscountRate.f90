REAL(8) FUNCTION FnDiscountRate(lrhoT)

USE Parameters
USE Globals
USE Procedures

IMPLICIT NONE

REAL(8), INTENT(IN) :: lrhoT
REAL(8)				:: lstep_lab,llabor_Y,llabor_N,lK_totoutput_ratio

lstep_lab = 1.0

rho = -log(logistic(lrhoT))

CALL Grids
CALL InitialPrices
CALL IterateBellman
CALL StationaryDistribution
CALL DistributionStatistics


!implied capital-output ratio
CALL KNYfromANY(Ea,lK_totoutput_ratio,price_W)

!update labor
llabor_Y = Elabor_Y/varieties
llabor_N = Elabor_N
labor_Y = lstep_lab*llabor_Y + (1.0-lstep_lab)*labor_Y
labor_N = lstep_lab*llabor_N + (1.0-lstep_lab)*labor_N

FnDiscountRate = lK_totoutput_ratio/target_K_totoutput_ratio - 1.0
IF (Display>=1) THEN
	write(*,*) 'Rho calirbation iter ',neqmiter
	write(*,"(9999(G14.6,:,','))") ' K-NY ratio:',target_K_totoutput_ratio,', labor_Y:',labor_Y,', labor_N:',labor_N
	write(*,"(9999(G14.6,:,','))") ' K-NY ratio:',lK_totoutput_ratio,', labor_Y:',llabor_Y,', labor_N:',llabor_N
	write(*,*) ' Err: ',FnDiscountRate
	write(*,*) ' '
END IF


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

neqmiter = neqmiter+1
END FUNCTION FnDiscountRate
