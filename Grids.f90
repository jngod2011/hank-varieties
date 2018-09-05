SUBROUTINE Grids

USE Parameters
USE Globals
USE Procedures

IMPLICIT NONE

INTEGER			:: ia,ib,iy,iaby,iab
REAL(8)         :: lElogy,lVlogy,lSDlogy,lmean,lwidth,ltemp1,ltemp2,ltemp3,ltemp4,lwidthSD,lf(ngpa),lg(ngpa)
REAL(8)         :: lm,lxa,lxb,lxc,lfa,lfb,lfc,lfval,leye(ngpy,ngpy),lmeangrosslabinc,lydistcum(ngpy)

REAL(8), EXTERNAL :: FnGridAR1,golden

!identity matrix
leye = 0.0; DO iy = 1,ngpy
	leye(iy,iy) = 1.0
END DO

!productivity process
IF(ngpy==1) THEN
	logygrid = 0.0
	ygrid = 1.0
	ymarkov = 0.0
	ytrans = 1.0
	ydist = 1.0

ELSE IF (ReadEarningsProcess==1) THEN
	OPEN(1, FILE = trim(EarningsProcessDir) // '/ygrid_combined.txt');READ(1,*) logygrid;CLOSE(1)
	OPEN(1, FILE = trim(EarningsProcessDir) // '/ydist_combined.txt');READ(1,*) ydist;CLOSE(1)
	OPEN(1, FILE = trim(EarningsProcessDir) // '/ymarkov_combined.txt');READ(1,*) ymarkov;CLOSE(1)
	ymarkov = TRANSPOSE(ymarkov) !since fortran reads in column major order
	
	!adjust for mean 1
	ygrid = exp(logygrid)
	lmean = DOT_PRODUCT(ydist,ygrid)
	ygrid = ygrid/lmean

	!fix up rounding in markov matrix
	DO iy = 1,ngpy
		ymarkov(iy,iy) = ymarkov(iy,iy) - SUM(ymarkov(iy,:))
	END DO
	
ELSE IF (TwoPointWageProcess==1) THEN
	ygrid(1) = 0.8
	ygrid(2) = 1.2
	ytrans(1,1) = 1.0-0.06667
	ytrans(1,2) = 0.06667
	ytrans(2,1) = 0.06667
	ytrans(2,2) = 1.0-0.06667
	ymarkov = (ytrans-leye)/1.0 !assumes ytrans is quarterly
	
ELSE IF (IncludeUnemploymentState==1 ) THEN! 1 is unemployment, 2:ngpy are wage shocks
	logygrid(1) = log(uiben)
	lElogy = -(ysig**2)/(2*(1.0-yrho))
	lVlogy = ysig**2/(1.0-yrho**2)
	lSDlogy = sqrt(lVlogy)
	IF(ngpy==4) lwidthSD = 1.0
	IF(ngpy>=6) lwidthSD = 1.5
	
	CALL LinearSpacedGrid (ngpy-1,lElogy-lwidthSD*lSDlogy,lElogy+lwidthSD*lSDlogy,logygrid(2:ngpy))
	lwidth = logygrid(3)-logygrid(2)
	ygrid = exp(logygrid)
	
	!transitions if employed
	ytrans = 0.0
	DO iy = 2,ngpy
		lmean = yrho*logygrid(iy)-(ysig**2)/2

		IF(iy==2) CALL cumnor ((logygrid(2)+0.5*lwidth-lmean)/ysig, ytrans(iy,2), ytrans(iy,3))    
		IF(iy>2 .and. iy<ngpy) THEN
			CALL cumnor ((logygrid(iy)+0.5*lwidth-lmean)/ysig, ltemp1, ltemp2)
	        CALL cumnor ((logygrid(iy-1)+0.5*lwidth-lmean)/ysig, ltemp3, ltemp4)
			ytrans(iy,iy) = ltemp1 - ltemp3	
			ytrans(iy,iy-1) = ltemp3
			ytrans(iy,iy+1) = 1.0-ltemp1
		END IF
		IF(iy==ngpy) CALL cumnor ((logygrid(ngpy)-0.5*lwidth-lmean)/ysig, ytrans(iy,ngpy-1),ytrans(iy,ngpy))

		!adjust for poisson component
		ytrans(iy,2:ngpy) = ylambda*ytrans(iy,2:ngpy) 
		ytrans(iy,iy) = ytrans(iy,iy) + 1.0-ylambda 

		ytrans(iy,2:ngpy) = (1.0-jobloss)*ytrans(iy,2:ngpy)

		ytrans(iy,1) = jobloss
	END DO

	!transitions out of unemployment
	ytrans(1,1) = 1.0-jobfind
	ytrans(1,2) = jobfind
	ytrans(1,3:ngpy) = 0.0
	
	ymarkov = (ytrans-leye)/1.0 !assumes ytrans is quarterly
	
ELSE IF (IncludeUnemploymentState==0 ) THEN !use poisson AR(1)
	lxa = 1.0e-8
 	lxb = 1.0
 	CALL MNBRAK(lxa,lxb,lxc,lfa,lfb,lfc,FnGridAR1)
 	lfval = GOLDEN(lxa,lxb,lxc,FnGridAR1,1.0e-2_8,lm)
 	lfval = FnGridAR1(lm)
	IF(Display>=2) THEN
		write(*,*) 'ygrid: ',ygrid
		write(*,*) 'ydist: ',ydist
	END IF
	ymarkov = (ytrans-leye)/1.0 !assumes ytrans is quarterly
END IF

ymarkovdiag = 0.0
DO iy = 1,ngpy
	ymarkovdiag(iy,iy) = ymarkov(iy,iy)
END DO
ymarkovoff = ymarkov-ymarkovdiag

!labor transfer and taxes
DO iy = 1,ngpy
	lydistcum(iy) = SUM(ydist(1:iy))
END DO
IF(ngpy>1) THEN
	CALL LinInterp1 (ngpy,lydistcum,ygrid,lumptransferpc,lumptransfercutoff)
	lumptransfer = lumptransfercutoff * labtax * wage * ((netwage/chi)**frisch) 
ELSE IF (ngpy==1) THEN
	lumptransfer = 0.0
END IF
!hours worked: note these will change along transition and calibration
IF (LaborSupply==0) THEN
	lgrid = 1.0
	labdisutilgrid = 0.0
END IF
IF (LaborSupply==1) THEN
	IF(ScaleGHHIdiosyncratic==0) THEN	
		lgrid = (netwage*ygrid/chi)**frisch
		labdisutilgrid = chi*(lgrid**(1+1.0/frisch))/(1+1.0/frisch)
	ELSE IF(ScaleGHHIdiosyncratic==1) THEN
		lgrid = (netwage/chi)**frisch
		labdisutilgrid = chi*ygrid*(lgrid**(1+1.0/frisch))/(1+1.0/frisch)
	END IF
END IF	

grosslabincgrid = wage*lgrid*ygrid
lmeangrosslabinc = DOT_PRODUCT(grosslabincgrid,ydist)
netlabincgrid = netwage*lgrid*ygrid + lumptransfer
netinc2illgrid = directdepfrac*(min(max(netlabincgrid,directdepmin),directdepmax)-directdepmin)
netinc2liqgrid = netlabincgrid - netinc2illgrid

aendogmax = amax

CALL PowerSpacedGrid (ngpa,agridparam,0.0_8,min(amax*lmeangrosslabinc,aendogmax),agrid)
!with low gridparam points get bunched close to zero, so evenly space first 8 points;
DO ia = 1,10-1
	agrid(ia) = (ia-1)*agrid(10)/(10.0-1.0)
END DO
dagrid = agrid(2:ngpa)-agrid(1:ngpa-1)

!liquid grid
IF(Borrowing==0) THEN
	CALL PowerSpacedGrid (ngpb,bgridparam,0.0_8,bmax*lmeangrosslabinc,bgrid)
ELSE IF(Borrowing==1) THEN
	CALL PowerSpacedGrid (ngpbPOS,bgridparam,0.0_8,bmax*lmeangrosslabinc,bgrid(ngpbNEG+1:ngpb))
	nbl = maxval(labdisutilgrid-netinc2liqgrid)/(rborr+PerfectAnnuityMarkets*deathrate)
 	abl = max(nbl+ cmin,blim*lmeangrosslabinc)
	IF (Display>=2) write(*,*) 'natural borrowing limit = ',nbl
	IF (Display>=2) write(*,*) 'actual borrowing limit = ',abl
	CALL PowerSpacedGrid (ngpbNEG/2+1,bgridparamNEG,abl,(abl+bgrid(ngpbNEG+1))/2.0,bgrid(1:ngpbNEG/2+1))
	DO ib = ngpbNEG/2+2,ngpbNEG
		bgrid(ib) = bgrid(ngpbNEG+1) -(bgrid(ngpbNEG+2-ib)-bgrid(1))
	END DO
END IF
dbgrid = bgrid(2:ngpb)-bgrid(1:ngpb-1)

!areas for density approximation
! adelta(1:ngpa-1) = dagrid
! adelta(ngpa) = dagrid(ngpa-1)
! bdelta(1:ngpb-1) = dbgrid
! bdelta(ngpb) = dbgrid(ngpb-1)

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

!
! lf = exp(-agrid+110.0)
! lf = lf/sum(lf*adelta)
! write(*,*) 'sum lf',sum(lf*adelta),'mean lf',sum(lf*adelta*agrid)
!
! CALL AdjustDistProportionately(agrid,adelta,lf,0.8_8,lg)
!
! write(*,*) 'sum lf',sum(lf*adelta),'sum lg',sum(lg*adelta)
! write(*,*) 'mean lf',sum(lf*adelta*agrid),'mean lg',sum(lg*adelta*agrid)
! write(*,*) 'lf',lf
! write(*,*) 'lg',lg
! write(*,*) 'lf*adelta',lf*adelta
! write(*,*) 'lg*adelta',lg*adelta
! STOP

END SUBROUTINE Grids

