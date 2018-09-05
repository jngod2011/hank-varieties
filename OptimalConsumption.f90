SUBROUTINE OptimalConsumption(lVb,lc,lp,ls,lHc)

USE Parameters
USE Globals
USE Procedures

IMPLICIT NONE

REAL(8),INTENT(IN)	:: lVb
REAL(8),INTENT(OUT)	:: lc,lp,ls,lHc

IF(nondurwgt==1.0) THEN
	lp = 0.0
	lc = utilfn1inv(lVb) + glabdisutil
ELSE
	lp = max( utilfn1inv(lVb/nondurwgttilde)*(1.0-nondurwgt)/nondurwgttilde - housefrac*flowamnt*gill, 0.0_8)
	IF(lp>0.0) lc = (lp+housefrac*flowamnt*gill)*nondurwgt/(1.0-nondurwgt) + glabdisutil
	IF(lp==0) lc = (lVb / (nondurwgt * (housefrac*flowamnt*gill)**((1.0-gam)*(1.0-nondurwgt))) ) ** (1.0/(nondurwgt-1.0-nondurwgt*gam)) + glabdisutil
END IF
ls = gbdrift-lc-lp
lHc = utilfn(((lc-glabdisutil)**nondurwgt)*((lp+housefrac*flowamnt*gill)**(1.0-nondurwgt))) + lVb*ls

END SUBROUTINE OptimalConsumption