REAL(8) FUNCTION FnDiscountRate(lrhoT)

USE Parameters
USE Globals
USE Procedures

IMPLICIT NONE

REAL(8), INTENT(IN) :: lrhoT
REAL(8)				:: lKNratio,lKYratio,lcapital

rho = -log(logistic(lrhoT))

IF (Display>=2) THEN
	WRITE(*,*) '*******************************************'
	WRITE(*,*) ' ITERATION : 			',neqmiter
	WRITE(*,*) ' rho guess: 			',rho
	WRITE(*,*) '  target KY ratio: ',KYratio
	WRITE(*,*) '  target KN ratio: ',KNratio
END IF

CALL Grids
CALL IterateBellman
CALL StationaryDistribution
CALL DistributionStatistics

IF(DividendFundLumpSum==0) lcapital = (1.0-housefrac)*Ea/(1.0-fundlev)
IF(DividendFundLumpSum==1) lcapital = (1.0-housefrac)*Ea/(1.0-fundlev + (1.0-mc-operatecost)/(ra*KYratio))

lKNratio = lcapital / labor
lKYratio = (lKNratio**(1.0-alpha)) / tfp

IF(Display>=1) write(*,"(A,I2,A,E11.4,A,E11.4)") '  Rho iter ',neqmiter, ', rho ',rho, ', K/N err',lKNratio/KNratio - 1.0
IF(exploring==.true.) write(4,"(A,I2,A,E11.4,A,E11.4)") '   Rho iter ',neqmiter, ', rho ',rho, ', K/N err',lKNratio/KNratio - 1.0

FnDiscountRate = lKNratio/KNratio - 1.0


IF(CalibrateDiscountRate==1) OPEN(3, FILE = trim(OutputDir) // 'DiscountRateCalibration.txt', ACCESS = 'append')
IF(CalibrateRhoAtInitialGuess==1) OPEN(3, FILE = trim(OutputDir) // 'DiscountRateAtInitialGuess.txt', ACCESS = 'append')
WRITE(3,*) '*******************************************'
WRITE(3,*) ' ITERATION : 			',neqmiter
WRITE(3,*) ' rho guess: 			',rho
WRITE(3,*) '  target KY ratio: ',KYratio
WRITE(3,*) '  implied KY ratio: ',lKYratio
WRITE(3,*) '  target KN ratio: ',KNratio
WRITE(3,*) '  implied KN ratio: ',lKNratio
WRITE(3,*) '  relative error: ',FnDiscountRate
CLOSE(3)

neqmiter = neqmiter+1
    
END FUNCTION FnDiscountRate
