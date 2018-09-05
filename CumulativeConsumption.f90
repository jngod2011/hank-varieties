SUBROUTINE CumulativeConsumption

USE Parameters
USE Globals
USE Procedures
USE mUMFPACK

IMPLICIT NONE

INTEGER			:: iy,iaby,it,iw(naby),nsteps
REAL(8) 		:: ldelta
REAL(8), DIMENSION(nab) 	:: ldiag,lvec
REAL(8), DIMENSION(nab,ngpy) 	:: lcvec,lcumvec,lcumvec1
TYPE(tCSR_di), DIMENSION(ngpy) 	:: AUMF     !umfpack type
TYPE(CSR)						:: lA

IF (Display>=1) write(*,*) "Computing cumulative consumption for MPC"

!$OMP PARALLEL DO
DO iaby = 1,naby
	lcvec(abfromab(afromaby(iaby),bfromaby(iaby)),yfromaby(iaby)) = c(afromaby(iaby),bfromaby(iaby),yfromaby(iaby))
END DO
!$OMP END PARALLEL DO

ldelta = 0.01 !0.1

DO iy = 1,ngpy
	lA = AUCSR(iy)
	lA%val = -ldelta*lA%val
	ldiag = 1.0 - ldelta*ymarkovdiag(iy,iy)	
	CALL apldia (nab, 0, lA%val, lA%col, lA%row, ldiag, lA%val, lA%col, lA%row, iw )
	AUMF(iy) = tCSR_di(lA%row-1,lA%col-1,lA%val) 
END DO

!one quarter cumulative expected consumption
nsteps = int(1/ldelta)
lcumvec = 0.0
DO it = 1,nsteps

	!$OMP PARALLEL DO PRIVATE(lvec)
	DO iy = 1,ngpy	
		lvec = lcumvec(:,iy) + ldelta*lcvec(:,iy) + ldelta*MATMUL(lcumvec(:,:),ymarkovoff(:,iy))
		lcumvec1(:,iy) = AUMF(iy) .umfpack. lvec
	END DO
	!$OMP END PARALLEL DO 

	lcumvec = lcumvec1
END DO

!$OMP PARALLEL DO
DO iaby = 1,naby
	ccum1(afromaby(iaby),bfromaby(iaby),yfromaby(iaby)) = lcumvec(abfromab(afromaby(iaby),bfromaby(iaby)),yfromaby(iaby))
END DO
!$OMP END PARALLEL DO

!four quarter cumulative expected consumption
nsteps = int(4/ldelta)
lcumvec = 0.0
DO it = 1,nsteps

	!$OMP PARALLEL DO PRIVATE(lvec)
	DO iy = 1,ngpy	
		lvec = lcumvec(:,iy) + ldelta*lcvec(:,iy) + ldelta*MATMUL(lcumvec(:,:),ymarkovoff(:,iy))
		lcumvec1(:,iy) = AUMF(iy) .umfpack. lvec
	END DO
	!$OMP END PARALLEL DO 

	lcumvec = lcumvec1
END DO


!$OMP PARALLEL DO
DO iaby = 1,naby
	ccum4(afromaby(iaby),bfromaby(iaby),yfromaby(iaby)) = lcumvec(abfromab(afromaby(iaby),bfromaby(iaby)),yfromaby(iaby))
END DO
!$OMP END PARALLEL DO


!two quarter cumulative expected consumption
nsteps = int(2/ldelta)
lcumvec = 0.0
DO it = 1,nsteps

	!$OMP PARALLEL DO PRIVATE(lvec)
	DO iy = 1,ngpy	
		lvec = lcumvec(:,iy) + ldelta*lcvec(:,iy) + ldelta*MATMUL(lcumvec(:,:),ymarkovoff(:,iy))
		lcumvec1(:,iy) = AUMF(iy) .umfpack. lvec
	END DO
	!$OMP END PARALLEL DO 

	lcumvec = lcumvec1
END DO


!$OMP PARALLEL DO
DO iaby = 1,naby
	ccum2(afromaby(iaby),bfromaby(iaby),yfromaby(iaby)) = lcumvec(abfromab(afromaby(iaby),bfromaby(iaby)),yfromaby(iaby))
END DO
!$OMP END PARALLEL DO

END SUBROUTINE CumulativeConsumption