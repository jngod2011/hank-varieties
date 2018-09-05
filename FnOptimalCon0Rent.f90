REAL(8) FUNCTION  FnOptimalCon0Rent(lc)

USE Parameters
USE Globals
USE Procedures

IMPLICIT NONE

REAL(8), INTENT(IN) :: lc

gc = lc
FnOptimalCon0Rent =  -utilfn1inv(gVb) + ((gc**nondurwgt)*((housefrac*flowamnt*gill)**(1.0-nondurwgt))-glabdisutil) * utilfn1inv(nondurwgt*((gc/(housefrac*flowamnt*gill))**(nondurwgt-1.0)))

END FUNCTION FnOptimalCon0Rent