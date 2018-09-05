SUBROUTINE IterateBellman

USE Parameters
USE Globals
USE Procedures
USE mUMFPACK

IMPLICIT NONE

INTEGER			:: iaby,ii
REAL(8)			:: lVdiff,ltau,ltau0,ladrift

delta = deltass
	
!set drifts
ltau = 15.0
ltau0 = (ra+PerfectAnnuityMarkets*deathrate)*(1.0-housefrac)*(agrid(ngpa)*0.999)**(1.0-ltau)
!ltau0 = 0.0
! adrift = (ra+PerfectAnnuityMarkets*deathrate)*(1.0-housefrac)*agrid - ltau0*(agrid**ltau)
adrift = ra*(1.0-housefrac)*agrid + PerfectAnnuityMarkets*deathrate*agrid - ltau0*(agrid**ltau)
! adrift = (ra+PerfectAnnuityMarkets*deathrate)*(1.0-housefrac)*agrid
! adrift(ngpa) = -1.0
bdrift = MERGE((rb+PerfectAnnuityMarkets*deathrate)*bgrid,(rborr+PerfectAnnuityMarkets*deathrate)*bgrid,bgrid>0.0)
fspamount = 0.0
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
netlabincgrid = netwage*lgrid*ygrid + lumptransfer
netinc2illgrid = directdepfrac*(min(max(netlabincgrid,directdepmin),directdepmax)-directdepmin)
netinc2liqgrid = netlabincgrid - netinc2illgrid

IF (Borrowing==1 .and. bgrid(1) < maxval(labdisutilgrid-netinc2liqgrid)/(rborr+PerfectAnnuityMarkets*deathrate)) THEN
	write(*,*) 'Warning: natural borrowing limit violated'
END IF


!Initial Guess
IF	(	(EquilibriumR==0 .and. calibrating==.false.) &
 	.or. (EquilibriumR==1 .and. neqmiter<=3 .and. calibrating==.false.) &
 	.or. (CalibrateDiscountRate==1 .and. neqmiter<=3 .and. calibrating==.false.) &
 	.or. (CalibrateRhoAtInitialGuess==1 .and. neqmiter<=3 .and. calibrating==.false.) &
	.or. (calibrating==.true. .and. ImposeEqumInCalibration==1  .and. neqmiter==1 ) &
	.or. (calibrating==.true. .and. ImposeEqumInCalibration==0 ))  THEN
	
	!$OMP PARALLEL DO PRIVATE (ladrift)
	DO iaby = 1,naby
		glabdisutil = labdisutilgrid(yfromaby(iaby))
!   		ladrift = (ra+PerfectAnnuityMarkets*deathrate)*(1.0-housefrac)*agrid(afromaby(iaby))
  		ladrift = (ra*(1.0-housefrac) + PerfectAnnuityMarkets*deathrate)*agrid(afromaby(iaby))
! 		ladrift = max(adrift(afromaby(iaby)),0.0)
! 		V(afromaby(iaby),bfromaby(iaby),yfromaby(iaby)) = (utilfn(netlabincgrid(yfromaby(iaby)) + ladrift + bdrift(bfromaby(iaby)) - glabdisutil) + utilfnflow(housefrac*flowamnt*agrid(afromaby(iaby))) ) / (rho)
! 		V(afromaby(iaby),bfromaby(iaby),yfromaby(iaby)) = (utilfn(netlabincgrid(yfromaby(iaby)) + ladrift + (rb+PerfectAnnuityMarkets*deathrate)*bgrid(bfromaby(iaby)) - glabdisutil) + utilfnflow(housefrac*flowamnt*agrid(afromaby(iaby))) ) / (rho+deathrate)
! 		V(afromaby(iaby),bfromaby(iaby),yfromaby(iaby)) = (utilfn(netlabincgrid(yfromaby(iaby)) + ladrift - glabdisutil) + utilfnflow(housefrac*flowamnt*agrid(afromaby(iaby))) ) / (rho)
! 		V(afromaby(iaby),bfromaby(iaby),yfromaby(iaby)) = utilfn(netlabincgrid(yfromaby(iaby))- glabdisutil) / (rho)
! 		V(afromaby(iaby),bfromaby(iaby),yfromaby(iaby)) = (utilfn(netlabincgrid(yfromaby(iaby)) + (rb+PerfectAnnuityMarkets*deathrate)*bgrid(bfromaby(iaby)) - glabdisutil) + utilfnflow(housefrac*flowamnt*agrid(afromaby(iaby))) ) / (rho)
! 		V(afromaby(iaby),bfromaby(iaby),yfromaby(iaby)) = (utilfn(netlabincgrid(yfromaby(iaby)) + (rb+PerfectAnnuityMarkets*deathrate)*bgrid(bfromaby(iaby)) - glabdisutil) ) / (rho+deathrate)
! 		V(afromaby(iaby),bfromaby(iaby),yfromaby(iaby)) = (utilfn(netlabincgrid(yfromaby(iaby)) + (rb+PerfectAnnuityMarkets*deathrate)*bgrid(bfromaby(iaby)) - glabdisutil) + utilfnflow(housefrac*flowamnt*agrid(afromaby(iaby))) ) / (rho+deathrate)
! 		V(afromaby(iaby),bfromaby(iaby),yfromaby(iaby)) = (utilfn(netinc2liqgrid(yfromaby(iaby)) + (rb+PerfectAnnuityMarkets*deathrate)*bgrid(bfromaby(iaby)) - glabdisutil) + utilfnflow(housefrac*flowamnt*agrid(afromaby(iaby))) ) / (rho+deathrate)		
		V(afromaby(iaby),bfromaby(iaby),yfromaby(iaby)) = (utilfn(netinc2liqgrid(yfromaby(iaby)) + (rb+PerfectAnnuityMarkets*deathrate)*bgrid(bfromaby(iaby)) + housefrac*flowamnt*agrid(afromaby(iaby)) - glabdisutil)) / (rho+deathrate)
	END DO
	!$OMP END PARALLEL DO
	
END IF



ii = 1
lVdiff = 1.0
DO WHILE (ii<=maxiter .and. lVdiff>Vtol)

	CALL HJBUpdate

	!check for convergence
	lVdiff = maxval(abs(Vnew-V))
	IF(Display>=2) write(*,*) "Iteration: ",ii," max V change: ",lVdiff
	
	!update
	V = Vnew
	ii = ii+1
	
END DO


END SUBROUTINE IterateBellman