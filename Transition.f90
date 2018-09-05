SUBROUTINE Transition

USE Parameters
USE Globals
USE Procedures

IMPLICIT NONE

INTEGER						:: it,iw(nab),iaby,ia,ib,iy,iab
REAL(8), DIMENSION(nab) 	:: ldavec,ldbvec,ldiag
REAL(8) 					:: ltau,ltau0
TYPE(tCSR_di),DIMENSION(ngpy) 	 :: AUMF     !umfpack type
LOGICAL 					:: nblviolated
REAL(8) 					:: lbgrid(ngpb,Ttransition),lbmin,lbmax,holdbgrid(ngpb)
REAL(8), DIMENSION(nab,ngpy) 	:: lgmat,lgmat1
REAL(8), DIMENSION(ngpa,ngpb,ngpy) 	:: lgjoint
REAL(8), DIMENSION(ngpy,ngpy) 	:: leye,lmat

!identity matrix
leye = 0.0; DO iy = 1,ngpy
	leye(iy,iy) = 1.0
END DO

!$OMP PARALLEL DO
DO iab = 1,nab
	ldavec(iab) = adelta(afromab(iab))
	ldbvec(iab) = bdelta(bfromab(iab))
END DO
!$OMP END PARALLEL DO

holdbgrid = bgrid

!solve backward
DO it = Ttransition,1,-1
	IF(Display>=2) WRITE(*,*) "  Solving transition backward: ",it
	
	
	IF(it==Ttransition) V = solnFINALSS%V
	IF(it<Ttransition) V = solnTRANS(it+1)%V

	!set up fiscal stimulus
	fspamount = 0.0
	IF (fsptransition== .true.) THEN
	 	IF(cumdeltatrans(it)>=fspointer%fspstart .and. cumdeltatrans(it)<fspointer%fspend) fspamount = fspointer%fspamount
	END IF
	
	!set drifts and globals
	rho = equmTRANS(it)%rho
	ra = equmTRANS(it)%ra
	rborr = equmTRANS(it)%rborr
	borrwedge = equmTRANS(it)%borrwedge
	wage = equmTRANS(it)%wage
	labor = equmTRANS(it)%labor
	netwage = equmTRANS(it)%netwage
	labtax = equmTRANS(it)%labtax
	lumptransfer = equmTRANS(it)%lumptransfer
	rb = equmTRANS(it)%rb
	tfp = equmTRANS(it)%tfp
	kappafc_w = equmTRANS(it)%kappafc_w
	mpshock = equmTRANS(it)%mpshock
	ysig = equmTRANS(it)%ysig
	prefshock = equmTRANS(it)%prefshock
	fundlev = equmTRANS(it)%fundlev
	fundbond = equmTRANS(it)%fundbond
	worldbond = equmTRANS(it)%worldbond
	elast = equmTRANS(it)%elast
	gam = equmTRANS(it)%gam
	
! 	IF(it==1) THEN
! 		write(*,*) ' ra ',equmINITSS%ra,ra
! 		write(*,*) ' rb ',equmINITSS%rb,rb
! 		write(*,*) ' wage ',equmINITSS%wage,wage
! 		write(*,*) ' equity',equmINITSS%equity,equmTRANS(1)%equity
! 		write(*,*) ' profit',equmINITSS%profit,equmTRANS(1)%profit
! 		write(*,*) ' illequitydrop',illequitydrop
! 	END IF
	
	ltau = 15.0
	ltau0 = (ra+PerfectAnnuityMarkets*deathrate)*(1.0-housefrac)*(agrid(ngpa)*0.999)**(1.0-ltau)
! 	adrift = (ra+PerfectAnnuityMarkets*deathrate)*(1.0-housefrac)*agrid - ltau0*(agrid**ltau)
	adrift = (ra*(1.0-housefrac)+PerfectAnnuityMarkets*deathrate)*agrid - ltau0*(agrid**ltau)
	bdrift = MERGE((rb+PerfectAnnuityMarkets*deathrate)*bgrid,(rborr+PerfectAnnuityMarkets*deathrate)*bgrid,bgrid>0.0)
	IF (LaborSupply==0) THEN
		lgrid = 1.0
		labdisutilgrid = 0.0
	END IF
	IF (LaborSupply==1) THEN
		IF(ScaleGHHIdiosyncratic==0) THEN
			lgrid = (netwage*ygrid/chi)**frisch
			labdisutilgrid = chi*(lgrid**(1+1.0/frisch))/(1+1.0/frisch)
		ELSE IF(ScaleGHHIdiosyncratic==1) THEN
			lgrid = labor
! 			labdisutilgrid = chi*ygrid*(lgrid**(1+1.0/frisch))/(1+1.0/frisch)
			IF((testingGHH == .false.) .and. (NoChangeLabDisutility==0)) labdisutilgrid = chi*ygrid*(labor**(1+1.0/frisch))/(1+1.0/frisch)
			IF((testingGHH == .true.) .or. (NoChangeLabDisutility==1)) labdisutilgrid = chi*ygrid*(equmINITSS%labor**(1+1.0/frisch))/(1+1.0/frisch)
			
		END IF
	END IF	
	grosslabincgrid = wage*lgrid*ygrid
	netlabincgrid = netwage*lgrid*ygrid + lumptransfer
	netinc2illgrid = directdepfrac*(min(max(netlabincgrid,directdepmin),directdepmax)-directdepmin)
	netinc2liqgrid = netlabincgrid - netinc2illgrid
	
	nblviolated = .false.
	lbgrid(:,it) = bgrid
	IF (Borrowing==1 .and. rborr+PerfectAnnuityMarkets*deathrate>0.0 .and. bgrid(1) < maxval(labdisutilgrid-netinc2liqgrid)/(rborr+PerfectAnnuityMarkets*deathrate)) THEN
		IF(Display>=1) write(*,*) '  Warning: natural borrowing limit violated in transition'
		IF(Display>=1) write(*,*) '  Steady State ABL: ',abl, 'Current NBL: ',maxval(labdisutilgrid-netinc2liqgrid)/(rborr+PerfectAnnuityMarkets*deathrate)
		IF(Display>=1) THEN
			write(*,*) ' pi ',equmTRANS(it)%pi
			write(*,*) ' rborr ',rborr
			write(*,*) ' rb ',rb
			write(*,*) ' wage ',wage
			write(*,*) ' mc ',equmTRANS(it)%mc
			write(*,*) ' KYratio ',equmTRANS(it)%KYratio
			write(*,*) ' mpshock ',equmTRANS(it)%mpshock
		END IF
		STOP
		nblviolated = .true.
		lbmin = maxval(labdisutilgrid-netinc2liqgrid)/(rborr+PerfectAnnuityMarkets*deathrate) + cmin
		
		
		CALL PowerSpacedGrid (ngpbNEG/2+1,bgridparamNEG,lbmin,(lbmin+lbgrid(ngpbNEG+1,it))/2.0,lbgrid(1:ngpbNEG/2+1,it))
		DO ib = ngpbNEG/2+2,ngpbNEG
			lbgrid(ib,it) = lbgrid(ngpbNEG+1,it) -(lbgrid(ngpbNEG+2-ib,it)-lbgrid(1,it))
		END DO
		
		bdrift = MERGE((rb+PerfectAnnuityMarkets*deathrate)*lbgrid(:,it),(rborr+PerfectAnnuityMarkets*deathrate)*lbgrid(:,it),lbgrid(:,it)>0.0)
		
	END IF
	IF (Borrowing==1 .and. rb+PerfectAnnuityMarkets*deathrate<0.0 .and. bgrid(ngpb) > minval(netinc2liqgrid-labdisutilgrid)/(-rb-PerfectAnnuityMarkets*deathrate)) THEN
		IF(Display>=1) write(*,*) '  Warning: upper bound on b violated in transition with rb<0'
		IF(Display>=1) write(*,*) '  Steady State upperbound: ',bgrid(ngpb), 'Current upper bound: ',minval(netinc2liqgrid-labdisutilgrid)/(-rb-PerfectAnnuityMarkets*deathrate)

! 		nblviolated = .true.
! 		lbmax = minval(netinc2liqgrid-labdisutilgrid)/(-rb-PerfectAnnuityMarkets*deathrate) - cmin		
! 		CALL PowerSpacedGrid (ngpbPOS,bgridparam,0.0_8,lbmax,lbgrid(ngpbNEG+1:ngpb,it))		
! 		bdrift = MERGE((rb+PerfectAnnuityMarkets*deathrate)*lbgrid(:,it),(rborr+PerfectAnnuityMarkets*deathrate)*lbgrid(:,it),lbgrid(:,it)>0.0)
	END IF

	delta = deltatransvec(it)
	CALL HJBUpdate
	
! 	
! 	!store value functions and A matrix
	solnTRANS(it)%V = Vnew
	solnTRANS(it)%u = u 
	solnTRANS(it)%c = c 
	solnTRANS(it)%p = p
	solnTRANS(it)%d = d 
	solnTRANS(it)%bdot = bdot
! 	
! 	$OMP PARALLEL DO
	DO iy = 1,ngpy
		IF (ALLOCATED(solnTRANS(it)%A(iy)%val)) DEALLOCATE(solnTRANS(it)%A(iy)%val) 
		IF (ALLOCATED(solnTRANS(it)%A(iy)%row)) DEALLOCATE(solnTRANS(it)%A(iy)%row) 
		IF (ALLOCATED(solnTRANS(it)%A(iy)%col)) DEALLOCATE(solnTRANS(it)%A(iy)%col) 
		ALLOCATE(solnTRANS(it)%A(iy)%val(ACSR(iy)%nz),solnTRANS(it)%A(iy)%row(ACSR(iy)%n+1),solnTRANS(it)%A(iy)%col(ACSR(iy)%nz))

		IF (ALLOCATED(solnTRANS(it)%AU(iy)%val)) DEALLOCATE(solnTRANS(it)%AU(iy)%val) 
		IF (ALLOCATED(solnTRANS(it)%AU(iy)%row)) DEALLOCATE(solnTRANS(it)%AU(iy)%row) 
		IF (ALLOCATED(solnTRANS(it)%AU(iy)%col)) DEALLOCATE(solnTRANS(it)%AU(iy)%col) 
 		ALLOCATE(solnTRANS(it)%AU(iy)%val(AUCSR(iy)%nz),solnTRANS(it)%AU(iy)%row(AUCSR(iy)%n+1),solnTRANS(it)%AU(iy)%col(AUCSR(iy)%nz))
		
	END DO
! 	!$OMP END PARALLEL DO
	solnTRANS(it)%A = ACSR
	solnTRANS(it)%AU = AUCSR
		
END DO

!simulate forward
DO it = 1,Ttransition
	IF(Display>=2) WRITE(*,*) "Iterating transition forward: ",it
	
	!set up fiscal stimulus
	IF (fsptransition== .true.) THEN
	 	IF(cumdeltatrans(it)>=fspointer%fspstart .and. cumdeltatrans(it)<fspointer%fspend) fspamount = fspointer%fspamount
	END IF	
	
	!$OMP PARALLEL DO PRIVATE(ldiag, iw)
	DO iy = 1,ngpy
	
		!construct B matrix
		solnTRANS(it)%B(iy)%nz = solnTRANS(it)%AU(iy)%nz
		solnTRANS(it)%B(iy)%n = solnTRANS(it)%AU(iy)%n
		IF (ALLOCATED(solnTRANS(it)%B(iy)%val)) DEALLOCATE(solnTRANS(it)%B(iy)%val) 
		IF (ALLOCATED(solnTRANS(it)%B(iy)%row)) DEALLOCATE(solnTRANS(it)%B(iy)%row) 
		IF (ALLOCATED(solnTRANS(it)%B(iy)%col)) DEALLOCATE(solnTRANS(it)%B(iy)%col) 
		ALLOCATE(solnTRANS(it)%B(iy)%val(solnTRANS(it)%B(iy)%nz),solnTRANS(it)%B(iy)%row(solnTRANS(it)%B(iy)%n+1),solnTRANS(it)%B(iy)%col(solnTRANS(it)%B(iy)%nz))

		!transpose
		solnTRANS(it)%B(iy) = solnTRANS(it)%AU(iy)
		CALL csrcsc (nab, 1, 1, solnTRANS(it)%AU(iy)%val, solnTRANS(it)%AU(iy)%col, solnTRANS(it)%AU(iy)%row, solnTRANS(it)%B(iy)%val, solnTRANS(it)%B(iy)%col, solnTRANS(it)%B(iy)%row )

		!adjust A' matrix for non-linearly spaced grids
		ldiag = 1.0/(ldavec*ldbvec)
		CALL diamua (nab, 1, solnTRANS(it)%B(iy)%val, solnTRANS(it)%B(iy)%col, solnTRANS(it)%B(iy)%row, ldiag,solnTRANS(it)%B(iy)%val, solnTRANS(it)%B(iy)%col, solnTRANS(it)%B(iy)%row )
		ldiag = ldavec*ldbvec
		CALL amudia (nab, 1, solnTRANS(it)%B(iy)%val, solnTRANS(it)%B(iy)%col, solnTRANS(it)%B(iy)%row, ldiag,solnTRANS(it)%B(iy)%val, solnTRANS(it)%B(iy)%col, solnTRANS(it)%B(iy)%row )


		solnTRANS(it)%B(iy)%val = -deltatransvec(it)*solnTRANS(it)%B(iy)%val
		ldiag = 1.0  - deltatransvec(it)*ymarkovdiag(iy,iy) + deltatransvec(it)*deathrate

		CALL apldia (nab, 0, solnTRANS(it)%B(iy)%val, solnTRANS(it)%B(iy)%col, solnTRANS(it)%B(iy)%row, ldiag, solnTRANS(it)%B(iy)%val, solnTRANS(it)%B(iy)%col, solnTRANS(it)%B(iy)%row, iw )
	
		AUMF(iy) = tCSR_di(solnTRANS(it)%B(iy)%row-1,solnTRANS(it)%B(iy)%col-1,solnTRANS(it)%B(iy)%val) 	
	END DO
	!$OMP END PARALLEL DO
	
	IF(it==1 .and. DividendFundLumpSum==0) lgmat = solnINITSS%gmat
	IF(it==1 .and. DividendFundLumpSum==1) THEN
		
		!$OMP PARALLEL DO PRIVATE(ib)
		DO iy = 1,ngpy
			DO ib = 1,ngpb
				CALL AdjustDistProportionately(agrid,adelta,solnINITSS%gjoint(:,ib,iy),illequitydrop,lgjoint(:,ib,iy))
			END DO
		END DO
		!$OMP END PARALLEL DO

		!$OMP PARALLEL DO PRIVATE(ia,ib,iy,iab)	
		DO iaby = 1,naby
			ia = afromaby(iaby)
			ib = bfromaby(iaby)
			iy = yfromaby(iaby)
			iab = 	abfromab(ia,ib)
			lgmat(iab,iy) = lgjoint(ia,ib,iy)
		END DO
		!$OMP END PARALLEL DO
		
	END IF
	
	IF(it>1) lgmat = solnTRANS(it-1)%gmat

	lmat = leye + deltatransvec(it)*TRANSPOSE(ymarkovoff) 
	
	!$OMP PARALLEL DO
	DO iy = 1,ngpy
	
		lgmat1(:,iy) = MATMUL(lgmat(:,:),lmat(iy,:))
		lgmat1(:,iy) = lgmat1(:,iy) + deltatransvec(it)*deathrate*MERGE(SUM(lgmat(:,iy)*ldavec*ldbvec),0.0_8,afromab==1 .and. bfromab==ngpbNEG+1 )/(ldavec*ldbvec)
		lgmat1(:,iy) = AUMF(iy) .umfpack. lgmat1(:,iy)
		
	END DO
	!$OMP END PARALLEL DO
	
	solnTRANS(it)%gmat = lgmat1
	
	!$OMP PARALLEL DO PRIVATE(ia,ib,iy,iab)
	DO iaby = 1,ngpa*ngpb*ngpy
		ia = afromaby(iaby)
		ib = bfromaby(iaby)
		iy = yfromaby(iaby)
		iab = 	abfromab(ia,ib)
		solnTRANS(it)%gvec(iaby) = solnTRANS(it)%gmat(iab,iy)
		solnTRANS(it)%gjoint(ia,ib,iy) = solnTRANS(it)%gmat(iab,iy)		
	END DO
	!$OMP END PARALLEL DO
	
	
	!marginal distributions
	!$OMP PARALLEL DO PRIVATE(ia,ib)
	DO iy = 1,ngpy
		DO ia = 1,ngpa
			solnTRANS(it)%gamarg(ia,iy) = SUM(solnTRANS(it)%gjoint(ia,:,iy)*bdelta)
		END DO

		DO ib = 1,ngpb
			solnTRANS(it)%gbmarg(ib,iy) = SUM(solnTRANS(it)%gjoint(:,ib,iy)*adelta)
		END DO
	END DO
	!$OMP END PARALLEL DO
	
	!set globals: apply current period policy function to last period distribution of state variables
 	IF(it==1  .and. DividendFundLumpSum==0) THEN
		gvec = solnINITSS%gvec
		gmat = solnINITSS%gmat
		gamarg = solnINITSS%gamarg
		gbmarg = solnINITSS%gbmarg
 	ELSE IF(it==1  .and. DividendFundLumpSum==1) THEN
		gmat = lgmat
		gjoint = lgjoint
		
		!$OMP PARALLEL DO PRIVATE(ia,ib,iy,iab)	
		DO iaby = 1,naby
			ia = afromaby(iaby)
			ib = bfromaby(iaby)
			iy = yfromaby(iaby)
			iab = 	abfromab(ia,ib)
			gvec(iaby) = gmat(iab,iy)
		END DO
		!$OMP END PARALLEL DO

		DO ia = 1,ngpa
		DO iy = 1,ngpy
			gamarg(ia,iy) = SUM(gjoint(ia,:,iy)*bdelta)
		END DO
		END DO

		DO ib = 1,ngpb
		DO iy = 1,ngpy
			gbmarg(ib,iy) = SUM(gjoint(:,ib,iy)*adelta)
		END DO
		END DO
		
	ELSE IF(it>1) THEN
		gvec = solnTRANS(it-1)%gvec
		gmat = solnTRANS(it-1)%gmat
		gamarg = solnTRANS(it-1)%gamarg
		gbmarg = solnTRANS(it-1)%gbmarg
	END IF	
		
	wage = equmTRANS(it)%wage
	labor = equmTRANS(it)%labor
	netwage = equmTRANS(it)%netwage
	lumptransfer = equmTRANS(it)%lumptransfer	
	kappafc_w = equmTRANS(it)%kappafc_w	
	rb = equmTRANS(it)%rb
	ra = equmTRANS(it)%ra
	rborr = equmTRANS(it)%rborr
	
	c = solnTRANS(it)%c
	p = solnTRANS(it)%p
	d = solnTRANS(it)%d
	IF (LaborSupply==0) THEN
		lgrid = 1.0
		labdisutilgrid = 0.0
	END IF
	IF (LaborSupply==1) THEN
		IF(ScaleGHHIdiosyncratic==0) THEN
			lgrid = (netwage*ygrid/chi)**frisch
			labdisutilgrid = chi*(lgrid**(1+1.0/frisch))/(1+1.0/frisch)
		ELSE IF(ScaleGHHIdiosyncratic==1) THEN
			lgrid = labor
! 			labdisutilgrid = chi*ygrid*(lgrid**(1+1.0/frisch))/(1+1.0/frisch)
			IF((testingGHH == .false.) .and. (NoChangeLabDisutility==0)) labdisutilgrid = chi*ygrid*(labor**(1+1.0/frisch))/(1+1.0/frisch)
			IF((testingGHH == .true.) .or. (NoChangeLabDisutility==1)) labdisutilgrid = chi*ygrid*(equmINITSS%labor**(1+1.0/frisch))/(1+1.0/frisch)
		END IF
	END IF
		
	grosslabincgrid = wage*lgrid*ygrid
	netlabincgrid = netwage*lgrid*ygrid + lumptransfer
	netinc2illgrid = directdepfrac*(min(max(netlabincgrid,directdepmin),directdepmax)-directdepmin)
	netinc2liqgrid = netlabincgrid - netinc2illgrid

	IF(it>1) bgrid = lbgrid(:,it-1)
	!PROBLEM HERE SINCE bdelta IS USED IN DISTRIBUTION STATISTICS

	CALL DistributionStatistics
	bgrid = holdbgrid
	statsTRANS(it) = DistributionStatsType(Ea,Eb,Ec,Erent,Ed,Ewage,Enetlabinc,Egrosslabinc,Einc,Ehours,Enw,FRACa0,FRACa0close,FRACb0,FRACb0close,FRACb0a0,FRACb0aP,FRACbN,FRACnw0,FRACnw0close,FRACb0a0close, &
										EbN,EbP,Eadjcost,Efsp,PERCa,PERCb,PERCnw,PERCc,PERCinc,GINIa,GINIb,GINInw,GINIc,GINIinc, &
										Ea_nwQ,Eb_nwQ,Ec_nwQ,Einc_nwQ,Ea_incQ,Eb_incQ,Ec_incQ,Einc_incQ,Ec_bN,Ec_b0close,Ec_b0far)

END DO

END SUBROUTINE Transition


