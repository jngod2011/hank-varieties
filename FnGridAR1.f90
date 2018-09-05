REAL(8) FUNCTION  FnGridAR1(lx)

USE Parameters
USE Globals
USE Procedures

IMPLICIT NONE

REAL(8), INTENT(IN) :: lx
REAL(8)	:: lydist(ngpy),lydist1(ngpy),lSDlogy,lwidth,lmean,ltemp1,ltemp2,ltemp3,ltemp4,ldiff,lvarlog
INTEGER :: iy

lSDlogy = sqrt(ysig**2/(1.0-yrho**2))

CALL LinearSpacedGrid (ngpy,-lx*lSDlogy,lx*lSDlogy,logygrid)
lwidth = logygrid(2)-logygrid(1)

!transitions
ytrans = 0.0
DO iy = 1,ngpy
	lmean = yrho*logygrid(iy)
	IF(iy==1) CALL cumnor ((logygrid(1)+0.5*lwidth-lmean)/ysig, ytrans(iy,1), ytrans(iy,2))    
	IF(iy>1 .and. iy<ngpy) THEN
		CALL cumnor ((logygrid(iy)+0.5*lwidth-lmean)/ysig, ltemp1, ltemp2)
        CALL cumnor ((logygrid(iy-1)+0.5*lwidth-lmean)/ysig, ltemp3, ltemp4)
		ytrans(iy,iy) = ltemp1 - ltemp3	
		ytrans(iy,iy-1) = ltemp3
		ytrans(iy,iy+1) = 1.0-ltemp1
	END IF
	IF(iy==ngpy) CALL cumnor ((logygrid(ngpy)-0.5*lwidth-lmean)/ysig, ytrans(iy,ngpy-1),ytrans(iy,ngpy))

	!adjust for poisson component
	ytrans(iy,:) = ylambda*ytrans(iy,:) 
	ytrans(iy,iy) = ytrans(iy,iy) + 1.0-ylambda 

END DO

lydist = 1.0/real(ngpy)
ldiff = 1.0
DO WHILE (ldiff>1.0e-10) 
	lydist1 = MATMUL(lydist,ytrans)
	ldiff = maxval(abs(lydist-lydist1))
	lydist = lydist1
END DO
ydist = lydist1

!adjust so mean in levels is 1
ygrid = exp(logygrid)
lmean = DOT_PRODUCT(ygrid,ydist)
ygrid = ygrid/lmean
logygrid = log(ygrid)

lvarlog = DOT_PRODUCT(logygrid**2,ydist) - DOT_PRODUCT(logygrid,ydist)**2

!moment
FnGridAR1  = (lvarlog - ylogvar)**2
IF(Display>=2) write(*,*) 'AR1 grid:',lx,FnGridAR1

END FUNCTION  FnGridAR1