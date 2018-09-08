SUBROUTINE Transition

USE Parameters
USE Globals
USE Procedures

IMPLICIT NONE

INTEGER						:: it,iw(nab),iaby,ia,ib,iy,iab
REAL(8), DIMENSION(nab) 	:: ldavec,ldbvec,ldiag
REAL(8) 					:: ltau,ltau0
TYPE(tCSR_di),DIMENSION(ngpy) 	 :: AUMF     !umfpack type
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


!solve backward
DO it = Ttransition,1,-1
	IF(Display>=2) WRITE(*,*) "  Solving transition backward: ",it
	
	
	IF(it==Ttransition) V = solnFINALSS%V
	IF(it<Ttransition) V = solnTRANS(it+1)%V

	
	!set drifts and globals
	rho = equmTRANS(it)%rho
	ra = equmTRANS(it)%ra
	rb = equmTRANS(it)%rb
	rborr = equmTRANS(it)%rborr
	borrwedge = equmTRANS(it)%borrwedge
	wage_Y = equmTRANS(it)%wage_Y
	wage_N = equmTRANS(it)%wage_N
	labtax = equmTRANS(it)%labtax
	lumptransfer = equmTRANS(it)%lumptransfer
	tfp_Y = equmTRANS(it)%tfp_Y
	tfp_N = equmTRANS(it)%tfp_N
	kappa0_w = equmTRANS(it)%kappa0_w	
	IF (SymmetricAdjustmentCost) kappa0_d = kappa0_w
	kappa1_w = equmTRANS(it)%kappa1_w	
	IF (SymmetricAdjustmentCost) kappa1_d = kappa1_w
	mpshock = equmTRANS(it)%mpshock
	prefshock = equmTRANS(it)%prefshock
	worldbond = equmTRANS(it)%worldbond
	elast = equmTRANS(it)%elast
	gam = equmTRANS(it)%gam
	profit = equmTRANS(it)%profit
	finwedge = equmTRANS(it)%finwedge
	labwedge = equmTRANS(it)%labwedge
	yprodgrid = equmTRANS(it)%yprodgrid
	prodmarkovscale = equmTRANS(it)%prodmarkovscale
	netwagegrid = (1.0-labtax) * equmTRANS(it)%yprodgrid * ( equmTRANS(it)%wage_N*yoccgrid + equmTRANS(it)%wage_Y*(1.0-yoccgrid) )
	
	
	
    IF(it==1) THEN
  		write(*,*) ' '
   		write(*,*) ' ra ',equmINITSS%ra,ra
   		write(*,*) ' rb ',equmINITSS%rb,rb
   		write(*,*) ' equity',equmINITSS%equity_A,equmTRANS(1)%equity_A
   		write(*,*) ' profit',equmINITSS%profit,equmTRANS(1)%profit
   		write(*,*) ' wage_Y',equmINITSS%wage_Y,equmTRANS(1)%wage_Y
   		write(*,*) ' wage_N',equmINITSS%wage_N,equmTRANS(1)%wage_N
   		write(*,*) ' pi',equmINITSS%pi,equmTRANS(1)%pi
   		write(*,*) ' bond',equmINITSS%bond,equmTRANS(1:3)%bond
   		write(*,*) ' lump',equmINITSS%lumptransfer,equmTRANS(1:2)%lumptransfer
  		write(*,*) ' '
	END IF

	!set drifts
	ltau = 15.0
	ltau0 = (ra+PerfectAnnuityMarkets*deathrate)*(agrid(ngpa)*0.999)**(1.0-ltau)
	adrift = ra*agrid + PerfectAnnuityMarkets*deathrate*agrid - ltau0*(agrid**ltau)
	bdrift = MERGE((rb-finwedge+PerfectAnnuityMarkets*deathrate)*bgrid,(rborr+PerfectAnnuityMarkets*deathrate)*bgrid,bgrid>0.0)
	
	
	IF (Borrowing==1 .and. rborr+PerfectAnnuityMarkets*deathrate>0.0 .and. bgrid(1) < -lumptransfer/(rborr+PerfectAnnuityMarkets*deathrate)) THEN
		IF(Display>=1) THEN
			write(*,*) '  Warning: natural borrowing limit violated in transition'
			write(*,*) '  Steady State ABL: ',abl, 'Current NBL: ',-lumptransfer/(rborr+PerfectAnnuityMarkets*deathrate)
			write(*,*) rborr, lumptransfer, bgrid(1)
	  		write(*,*) ' it: ',it
	   		write(*,*) ' ra ',equmINITSS%ra,ra
	   		write(*,*) ' rb ',equmINITSS%rb,rb
	   		write(*,*) ' equity',equmINITSS%equity_A,equmTRANS(it)%equity_A
	   		write(*,*) ' profit',equmINITSS%profit,equmTRANS(it)%profit
	   		write(*,*) ' wage_Y',equmINITSS%wage_Y,wage_Y
	   		write(*,*) ' wage_Y',equmINITSS%wage_N,wage_N
	   		write(*,*) ' pi',equmINITSS%pi,equmTRANS(it)%pi
	   		write(*,*) ' bond',equmINITSS%bond,equmTRANS(1:3)%bond
	   		write(*,*) ' lump',equmINITSS%lumptransfer,lumptransfer
	  		write(*,*) ' '
			
		END IF
		nblviolated = .true.
		STOP
	END IF

	delta = deltatransvec(it)
	CALL HJBUpdate
 	
! 	!store value functions and A matrix
	solnTRANS(it)%V = Vnew
	solnTRANS(it)%u = u 
	solnTRANS(it)%c = c 
	solnTRANS(it)%s = s
	solnTRANS(it)%h = h
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
		ldiag = 1.0  - deltatransvec(it)*prodmarkovscale*ymarkovdiag(iy,iy) + deltatransvec(it)*deathrate

		CALL apldia (nab, 0, solnTRANS(it)%B(iy)%val, solnTRANS(it)%B(iy)%col, solnTRANS(it)%B(iy)%row, ldiag, solnTRANS(it)%B(iy)%val, solnTRANS(it)%B(iy)%col, solnTRANS(it)%B(iy)%row, iw )
	
		AUMF(iy) = tCSR_di(solnTRANS(it)%B(iy)%row-1,solnTRANS(it)%B(iy)%col-1,solnTRANS(it)%B(iy)%val) 	
	END DO
	!$OMP END PARALLEL DO
	
	!adjust for asset price drops on impact
	IF(it==1 .and. OneAssetNoCapital==1) lgmat = solnINITSS%gmat
	IF(it==1 .and. OneAssetNoCapital==0) THEN
		
		!$OMP PARALLEL DO PRIVATE(ib)
		DO iy = 1,ngpy
			DO ib = 1,ngpb
				CALL AdjustDistProportionately(agrid,adelta,solnINITSS%gjoint(:,ib,iy),equmTRANS(1)%assetdrop_A,lgjoint(:,ib,iy))
			END DO
		END DO
		!$OMP END PARALLEL DO

		!$OMP PARALLEL DO PRIVATE(ia)
		DO iy = 1,ngpy
			DO ia = 1,ngpa
				CALL AdjustDistProportionately(bgrid,bdelta,solnINITSS%gjoint(ia,:,iy),equmTRANS(1)%assetdrop_B,lgjoint(ia,:,iy))
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

	lmat = leye + deltatransvec(it)*TRANSPOSE(prodmarkovscale*ymarkovoff) 
	
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
	
	!set globals for distribution statistics: apply current period policy function to last period distribution of state variables
	IF(it==1 .and. OneAssetNoCapital==1)  THEN
		gvec = solnINITSS%gvec
		gmat = solnINITSS%gmat
		gamarg = solnINITSS%gamarg
		gbmarg = solnINITSS%gbmarg
 	ELSEIF(it==1 .and. OneAssetNoCapital==0) THEN
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
		
	wage_Y = equmTRANS(it)%wage_Y
	wage_N = equmTRANS(it)%wage_N
	lumptransfer = equmTRANS(it)%lumptransfer	
	labtax = equmTRANS(it)%labtax
	kappa0_w = equmTRANS(it)%kappa0_w	
	IF (SymmetricAdjustmentCost) kappa0_d = kappa0_w
	kappa1_w = equmTRANS(it)%kappa1_w	
	IF (SymmetricAdjustmentCost) kappa1_d = kappa1_w
	labtax = equmTRANS(it)%labtax
	rb = equmTRANS(it)%rb
	ra = equmTRANS(it)%ra
	rborr = equmTRANS(it)%rborr
	finwedge = equmTRANS(it)%finwedge
	profit = equmTRANS(it)%profit
	yprodgrid = equmTRANS(it)%yprodgrid
	
	c = solnTRANS(it)%c
	h = solnTRANS(it)%h 
	d = solnTRANS(it)%d


	CALL DistributionStatistics
	statsTRANS(it) = DistributionStatsType(Ea,Eb,Ec,Elabor,Elabor_N,Elabor_Y,Ed,Ewage,Enetlabinc,Egrosslabinc,Enetprofinc,Egrossprofinc,Einc,Ehours,Ehours_N,Ehours_Y,Enw, &
									FRACa0,FRACa0close,FRACb0,FRACb0close,FRACb0a0,FRACb0aP,FRACbN,FRACnw0,FRACnw0close,FRACb0a0close, &
 									EbN,EbP,Eadjcost,PERCa,PERCb,PERCnw,PERCc,PERCinc,GINIa,GINIb,GINInw,GINIc,GINIinc, &
 									Ea_nwQ,Eb_nwQ,Ec_nwQ,Einc_nwQ,Ea_incQ,Eb_incQ,Ec_incQ,Einc_incQ,Ec_bN,Ec_b0close,Ec_b0far,Ec_nwQ_add)

END DO

END SUBROUTINE Transition


