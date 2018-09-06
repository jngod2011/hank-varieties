SUBROUTINE Grids

USE Parameters
USE Globals
USE Procedures

IMPLICIT NONE

INTEGER			:: ia,ib,iy,iaby,iab,io,ip,iy2,io2,ip2
REAL(8)         :: lmean
REAL(8)         :: leye(ngpprod,ngpprod)

REAL(8), EXTERNAL :: golden

!identity matrix
leye = 0.0; 
DO ip = 1,ngpprod
	leye(ip,ip) = 1.0
END DO

!create occ, prod indices
iy = 0
DO io = 1,ngpocc
	DO ip = 1,ngpprod
		iy = iy+1
		occfromy(iy) = io
		prodfromy(iy) = ip
		yfromoccprod(io,ip) = iy
	END DO
END DO	

!productivity process
IF(ngpprod==1) THEN
	logprodgrid = 0.0
	prodgrid = 1.0
	prodmarkov = 0.0
	prodtrans = 1.0
	proddist = 1.0

ELSE IF (ReadEarningsProcess==1) THEN
	OPEN(1, FILE = trim(EarningsProcessDir) // '/ygrid_combined.txt');READ(1,*) logprodgrid;CLOSE(1)
	OPEN(1, FILE = trim(EarningsProcessDir) // '/ydist_combined.txt');READ(1,*) proddist;CLOSE(1)
	OPEN(1, FILE = trim(EarningsProcessDir) // '/ymarkov_combined.txt');READ(1,*) prodmarkov;CLOSE(1)
	IF (AdjustProdGridFrisch==1) logprodgrid = logprodgrid/ (1.0+adjfricshgridfrac*frisch)
	prodgrid = exp(logprodgrid)
	prodmarkov = TRANSPOSE(prodmarkov) !since fortran reads in column major order	
	!fix up rounding in markov matrix
	DO ip = 1,ngpprod
		prodmarkov(ip,ip) = prodmarkov(ip,ip) - SUM(prodmarkov(ip,:))
	END DO
	
ELSE IF (TwoPointWageProcess==1) THEN
	prodgrid(1) = 0.8
	prodgrid(2) = 1.2
	logprodgrid = log(prodgrid)
	
	prodtrans(1,1) = 1.0-0.06667
	prodtrans(1,2) = 0.06667
	prodtrans(2,1) = 0.06667
	prodtrans(2,2) = 1.0-0.06667
	proddist(1) = 0.5
	proddist(2) = 0.5
	prodmarkov = (prodtrans-leye)/1.0 !assumes ytrans is quarterly	
END IF

!adjust mean productivity
lmean = DOT_PRODUCT(proddist,prodgrid)
prodgrid = meanlabeff*prodgrid/lmean


!occupation types
IF(ngpocc==1) THEN
	occgrid = 0.5
	occdist = 1.0
ELSE !equally spaced [0,1]
	occgrid(1) = 0.0
	occgrid(ngpocc) = 1.0
	IF(ngpocc>=2) THEN
		DO io = 2,ngpocc-1
			occgrid(io) = occgrid(io-1) + 1.0/real(ngpocc-1)
		END DO	
	END IF
	occdist = 1.0/ngpocc
END IF	


!construct combined markov process and indicators
ymarkov = 0
iy = 0
DO iy = 1,ngpy
	io = occfromy(iy)
	ip = prodfromy(iy)
	yprodgrid(iy) = prodgrid(ip)
	yoccgrid(iy) = occgrid(io)
	ydist(iy) = proddist(ip)*occdist(io)
	DO iy2 = 1,ngpy
		io2 = occfromy(iy2)
		ip2 = prodfromy(iy2)
		IF (io==io2) THEN
			ymarkov(iy,iy2) = prodmarkov(ip,ip2)
		END IF
	END DO
END DO		

	
ymarkovdiag = 0.0
DO iy = 1,ngpy
	ymarkovdiag(iy,iy) = ymarkov(iy,iy)
END DO
ymarkovoff = ymarkov-ymarkovdiag


!profit distribution shares: compute in advance so can hold fix when labor productivity distribution moves
profsharegrid = yprodgrid/meanlabeff

!out of steady-state adjustments
prodgridscale = 1.0
prodmarkovscale = 1.0

CALL PowerSpacedGrid (ngpa,agridparam,0.0_8,amax,agrid)

!with low gridparam points get bunched close to zero, so evenly space first 8 points;
IF(ngpa>10) THEN
	DO ia = 1,10-1
		agrid(ia) = (ia-1)*agrid(10)/(10.0-1.0)
	END DO
END IF
dagrid = agrid(2:ngpa)-agrid(1:ngpa-1)

!liquid grid
IF(Borrowing==0) THEN
	CALL PowerSpacedGrid (ngpb,bgridparam,0.0_8,bmax,bgrid)
ELSE IF(Borrowing==1) THEN
	CALL PowerSpacedGrid (ngpbPOS,bgridparam,0.0_8,bmax,bgrid(ngpbNEG+1:ngpb))
	nbl = -lumptransfer/(rborr+PerfectAnnuityMarkets*deathrate)
 	abl = max(nbl + cmin,blim)
	IF (Display>=2) write(*,*) 'natural borrowing limit = ',nbl
	IF (Display>=2) write(*,*) 'actual borrowing limit = ',abl
	CALL PowerSpacedGrid (ngpbNEG/2+1,bgridparamNEG,abl,(abl+bgrid(ngpbNEG+1))/2.0,bgrid(1:ngpbNEG/2+1))
	DO ib = ngpbNEG/2+2,ngpbNEG
		bgrid(ib) = bgrid(ngpbNEG+1) -(bgrid(ngpbNEG+2-ib)-bgrid(1))
	END DO
END IF
dbgrid = bgrid(2:ngpb)-bgrid(1:ngpb-1)


adelta(1) = 0.5*dagrid(1)
adelta(2:ngpa-1) = 0.5*(dagrid(1:ngpa-2)+dagrid(2:ngpa-1))
adelta(ngpa) = 0.5*dagrid(ngpa-1)

bdelta(1) = 0.5*dbgrid(1)
bdelta(2:ngpb-1) = 0.5*(dbgrid(1:ngpb-2)+dbgrid(2:ngpb-1))
bdelta(ngpb) = 0.5*dbgrid(ngpb-1)


!combined grids for openmp loops
iab = 1
DO ia = 1,ngpa
DO ib = 1,ngpb
	afromab(iab) = ia
	bfromab(iab) = ib
	abfromab(ia,ib) = iab
	abdelta(iab) = adelta(ia)*bdelta(ib)
	iab = iab+1
END DO
END DO

iaby = 1
DO ia = 1,ngpa
DO ib = 1,ngpb
DO iy = 1,ngpy
	afromaby(iaby) = ia
	bfromaby(iaby) = ib
	yfromaby(iaby) = iy
	abyfromaby(ia,ib,iy) = iaby
	abydelta(iaby) = adelta(ia)*bdelta(ib)
	abfromaby(iaby) = abfromab(ia,ib)
	iaby = iaby+1
END DO
END DO
END DO



END SUBROUTINE Grids

