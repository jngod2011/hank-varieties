SUBROUTINE IRFSequence(lIRFDir)

USE Parameters
USE Globals
USE Procedures

IMPLICIT NONE

CHARACTER, INTENT(IN) :: lIRFDir*20
CHARACTER	:: lstring*80
INTEGER	:: it,ipe
REAL(8) :: lequity

irfpointer => irfstruct		

IF(SolveFlexPriceTransition==1) THEN
	flextransition = .true.
! 	CALL IterateTransitionFlex
	CALL IterateTransition
! 	IF (SaveCumPolicyFnsIRF==1) CALL CumulativeConsTransition
! 	IF(ComputeDiscountedMPC==1) CALL DiscountedMPCTransition
	flextransition = .false.				
END IF
IF(SolveStickyPriceTransition==1) THEN
	stickytransition = .true.
! 	IF(OneAssetNoCapital==0) CALL IterateTransitionStickyRb
	IF(OneAssetNoCapital==0) CALL IterateTransition
! 	IF(OneAssetNoCapital==1) CALL IterateTransOneAssetStickyRb
! 	IF (SaveCumPolicyFnsIRF==1) CALL CumulativeConsTransition
! 	IF(ComputeDiscountedMPC==1) CALL DiscountedMPCTransition
	stickytransition = .false.		
END IF	
IF(SolveZLBTransition==1) THEN
	zlbtransition = .true.
! 	IF(OneAssetNoCapital==0) CALL IterateTransitionStickyRb
	IF(OneAssetNoCapital==0) CALL IterateTransition
! 	IF(OneAssetNoCapital==1) CALL IterateTransOneAssetStickyRb
! 	IF (SaveCumPolicyFnsIRF==1) CALL CumulativeConsTransition
! 	IF(ComputeDiscountedMPC==1) CALL DiscountedMPCTransition
	zlbtransition = .false.
END IF

irfsave => irfstruct
OutputDirIRF = trim(OutputDir) // "IRF_" // trim(lIRFDir) // "/NOFS/"
CALL SaveIRFOutput

! IF (DoPriceExperiments==1) THEN
! 	irfpointer => irfpriceexp
!
! 	DO ipe = 1,16
! 	IF (WhichPriceExperiment(ipe)==1) THEN
! 		testingGHH = .false.
!
! 		IF(SolveFlexPriceTransition==1) THEN
! 			write(*,*) "Solving for flex price transition, price experiment ", ipe
!
! 			DO it = 1,Ttransition
! 				equmTRANS(it) = equmINITSS
! 			END DO
!
! 			SELECT CASE(ipe)
! 				CASE(1) !change wage and labor tax only
! 					equmTRANS%wage = irfstruct%equmFLEX%wage
! 					equmTRANS%netwage = irfstruct%equmFLEX%netwage
! 					equmTRANS%labtax = irfstruct%equmFLEX%labtax
!
! 				CASE(2) !only change profits
! 					equmTRANS%profit = irfstruct%equmFLEX%profit
!
! 				CASE(3) !only change profits and wage
! 					equmTRANS%profit = irfstruct%equmFLEX%profit
! 					equmTRANS%wage = irfstruct%equmFLEX%wage
! 					equmTRANS%netwage = irfstruct%equmFLEX%netwage
! 					equmTRANS%labtax = irfstruct%equmFLEX%labtax
!
! 				CASE(4) !only change rb (and rborr if the wedge is fixed out of steady state)
! 					equmTRANS%rb = irfstruct%equmFLEX%rb
! 					IF(FixBorrowRateTransition==0) equmTRANS%rborr = irfstruct%equmFLEX%rborr
!
! 				CASE(5) !only change ra only
! 					equmTRANS%ra = irfstruct%equmFLEX%ra
!
! 				CASE(6) !only change illiquid asset drop
! 					equmTRANS%illassetdrop  = irfstruct%equmFLEX%illassetdrop
!
! 				CASE(7) !only change transfers
! 					equmTRANS%lumptransfer = irfstruct%equmFLEX%lumptransfer
!
! 				CASE(8) !change all
! 					equmTRANS%wage = irfstruct%equmFLEX%wage
! 					equmTRANS%netwage = irfstruct%equmFLEX%netwage
! 					equmTRANS%labtax = irfstruct%equmFLEX%labtax
! 					equmTRANS%lumptransfer = irfstruct%equmFLEX%lumptransfer
! 					equmTRANS%rb = irfstruct%equmFLEX%rb
! 					equmTRANS%rborr = irfstruct%equmFLEX%rborr
! 					equmTRANS%ra = irfstruct%equmFLEX%ra
! 					equmTRANS%illassetdrop  = irfstruct%equmFLEX%illassetdrop
! 					equmTRANS%profit = irfstruct%equmFLEX%profit
! 					equmTRANS%labwedge = irfstruct%equmFLEX%labwedge
! 					equmTRANS%finwedge = irfstruct%equmFLEX%finwedge
! 					equmTRANS%prodgridscale = irfstruct%equmFLEX%prodgridscale
! 					DO it = 1,Ttransition
! 						equmTRANS(it)%ygrid = irfstruct%equmFLEX(it)%ygrid
! 					END DO
! 					equmTRANS%prodmarkovscale = irfstruct%equmFLEX%prodmarkovscale
! 					equmTRANS%kappa0_w = irfstruct%equmFLEX%kappa0_w
! 					equmTRANS%kappa1_w = irfstruct%equmFLEX%kappa1_w
! 					equmTRANS%prefshock = irfstruct%equmFLEX%prefshock
!
!  				CASE(9) !change lump transfer by direct effect from government interest payments
! 					equmTRANS%lumptransfer = equmINITSS%lumptransfer + (irfstruct%equmFLEX%rb - equmTRANS%rb)*equmINITSS%govbond
!
!  				CASE(10) !change lump transfer, not including direct effect from govt interest payments
! 					equmTRANS%lumptransfer = irfstruct%equmFLEX%lumptransfer - (irfstruct%equmFLEX%rb - equmTRANS%rb)*equmINITSS%govbond
!
! 				CASE(11) !change rb, rborr, and change ra by same amount as rb
! 					equmTRANS%rb = irfstruct%equmFLEX%rb
! 					equmTRANS%rborr = irfstruct%equmFLEX%rborr
! 					equmTRANS%ra = equmINITSS%ra + (irfstruct%equmFLEX%rb - equmINITSS%rb)
!
! 				CASE(12) !change ra, and change rb, rborr by same amount as ra
! 					equmTRANS%ra = irfstruct%equmFLEX%ra
! 					equmTRANS%rb = irfstruct%equmFLEX%rb + (irfstruct%equmFLEX%ra - equmINITSS%ra)
! 					equmTRANS%rborr = irfstruct%equmFLEX%rborr + (irfstruct%equmFLEX%ra - equmINITSS%ra)
!
! 				CASE(13) !change only proportional labor tax
! 					equmTRANS%labtax = irfstruct%equmFLEX%labtax
! 					equmTRANS%netwage = equmINITSS%wage*(1.0-irfstruct%equmFLEX%labtax)
!
! 				CASE(14) !change rb, rborr, and change ra by same amount as rb, and discount eqm profits at implied ra
! 					equmTRANS%rb = irfstruct%equmFLEX%rb
! 					equmTRANS%rborr = irfstruct%equmFLEX%rborr
! 					equmTRANS%ra = equmINITSS%ra + (irfstruct%equmFLEX%rb - equmINITSS%rb)
! 					IF(DividendFundLumpSum==0) THEN
! 						equmTRANS(:)%illassetdrop = 1.0
! 					ELSE IF(DividendFundLumpSum==1) THEN
! 						it = Ttransition
! 						lequity = (equmFINALSS%equity + irfstruct%equmFLEX(it)%dividend*deltatransvec(it)) / (1.0+deltatransvec(it)*equmTRANS(it)%ra)
! 						DO it = Ttransition-1,1,-1
! 							lequity = (lequity + irfstruct%equmFLEX(it)%dividend*deltatransvec(it)) / (1.0+deltatransvec(it)*equmTRANS(it)%ra)
! 						END DO
! 						equmTRANS(:)%illassetdrop = ((1.0-irfstruct%equmFLEX(1)%fundlev) * irfstruct%equmFLEX(1)%capital + lequity) / ((1.0-equmINITSS%fundlev)*equmINITSS%capital + equmINITSS%equity)
! 					END IF
!
! 				CASE(15) !change rb, rborr, and change ra by same amount as rb, and discount initimake al profits at implied ra
! 					equmTRANS%rb = irfstruct%equmFLEX%rb
! 					equmTRANS%rborr = irfstruct%equmFLEX%rborr
! 					equmTRANS%ra = equmINITSS%ra + (irfstruct%equmFLEX%rb - equmINITSS%rb)
! 					IF(DividendFundLumpSum==0) THEN
! 						equmTRANS(:)%illassetdrop = 1.0
! 					ELSE IF(DividendFundLumpSum==1) THEN
! 						it = Ttransition
! 						lequity = (equmFINALSS%equity + equmINITSS%dividend*deltatransvec(it)) / (1.0+deltatransvec(it)*equmTRANS(it)%ra)
! 						DO it = Ttransition-1,1,-1
! 							lequity = (lequity + equmINITSS%dividend*deltatransvec(it)) / (1.0+deltatransvec(it)*equmTRANS(it)%ra)
! 						END DO
! 						equmTRANS(:)%illassetdrop = ((1.0-irfstruct%equmFLEX(1)%fundlev) * irfstruct%equmFLEX(1)%capital + lequity) / ((1.0-equmINITSS%fundlev)*equmINITSS%capital + equmINITSS%equity)
! 					END IF
!
! 				CASE(16) !only change the direct shock
! 					equmTRANS%labwedge = irfstruct%equmFLEX%labwedge
! 					equmTRANS%finwedge = irfstruct%equmFLEX%finwedge
! 					IF(FixBorrowRateTransition==1) equmTRANS%rborr = irfstruct%equmFLEX%rborr
! 					equmTRANS%prodgridscale = irfstruct%equmFLEX%prodgridscale
! 					DO it = 1,Ttransition
! 						equmTRANS(it)%ygrid = irfstruct%equmFLEX(it)%ygrid
! 					END DO
! 					equmTRANS%prodmarkovscale = irfstruct%equmFLEX%prodmarkovscale
! 					equmTRANS%kappa0_w = irfstruct%equmFLEX%kappa0_w
! 					equmTRANS%kappa1_w = irfstruct%equmFLEX%kappa1_w
! 					equmTRANS%prefshock = irfstruct%equmFLEX%prefshock
!
! 			END SELECT
! 			CALL Transition
!
! 			IF (SaveCumPolicyFnsIRF==1) CALL CumulativeConsTransition
!
! 			irfpointer%equmFLEX = equmTRANS
! 			irfpointer%statsFLEX = statsTRANS
! 			irfpointer%solnFLEX = solnTRANS
! 			irfpointer%cumFLEX = CumulativePolicyType(ccum1,ccum2,ccum4,dcum1,dcum2,dcum4)
!
! 			IF(ComputeDiscountedMPC==1) CALL DiscountedMPCTransition
!
! 		END IF
!
! 		IF(SolveStickyPriceTransition==1) THEN
!
! 			write(*,*) "Solving for sticky price transition without ZLB, price experiment ", ipe
!
! 			DO it = 1,Ttransition
! 				equmTRANS(it) = equmINITSS
! 			END DO
!
! 			SELECT CASE(ipe)
! 				CASE(1) !change wage and labor tax only
! 					equmTRANS%wage = irfstruct%equmSTICKY%wage
! 					equmTRANS%netwage = irfstruct%equmSTICKY%netwage
! 					equmTRANS%labtax = irfstruct%equmSTICKY%labtax
!
! 				CASE(2) !only change profits
! 					equmTRANS%profit = irfstruct%equmSTICKY%profit
!
! 				CASE(3) !only change profits and wage
! 					equmTRANS%profit = irfstruct%equmSTICKY%profit
! 					equmTRANS%wage = irfstruct%equmSTICKY%wage
! 					equmTRANS%netwage = irfstruct%equmSTICKY%netwage
! 					equmTRANS%labtax = irfstruct%equmSTICKY%labtax
!
! 				CASE(4) !only change rb (and rborr if the wedge is fixed out of steady state)
! 					equmTRANS%rb = irfstruct%equmSTICKY%rb
! 					IF(FixBorrowRateTransition==0) equmTRANS%rborr = irfstruct%equmSTICKY%rborr
!
! 				CASE(5) !only change ra only
! 					equmTRANS%ra = irfstruct%equmSTICKY%ra
!
! 				CASE(6) !only change illiquid asset drop
! 					equmTRANS%illassetdrop  = irfstruct%equmSTICKY%illassetdrop
!
! 				CASE(7) !only change transfers
! 					equmTRANS%lumptransfer = irfstruct%equmSTICKY%lumptransfer
!
! 				CASE(8) !change all
! 					equmTRANS%wage = irfstruct%equmSTICKY%wage
! 					equmTRANS%netwage = irfstruct%equmSTICKY%netwage
! 					equmTRANS%labtax = irfstruct%equmSTICKY%labtax
! 					equmTRANS%lumptransfer = irfstruct%equmSTICKY%lumptransfer
! 					equmTRANS%rb = irfstruct%equmSTICKY%rb
! 					equmTRANS%rborr = irfstruct%equmSTICKY%rborr
! 					equmTRANS%ra = irfstruct%equmSTICKY%ra
! 					equmTRANS%illassetdrop  = irfstruct%equmSTICKY%illassetdrop
! 					equmTRANS%profit = irfstruct%equmSTICKY%profit
! 					equmTRANS%labwedge = irfstruct%equmSTICKY%labwedge
! 					equmTRANS%finwedge = irfstruct%equmSTICKY%finwedge
! 					equmTRANS%prodgridscale = irfstruct%equmSTICKY%prodgridscale
! 					DO it = 1,Ttransition
! 						equmTRANS(it)%ygrid = irfstruct%equmSTICKY(it)%ygrid
! 					END DO
! 					equmTRANS%prodmarkovscale = irfstruct%equmSTICKY%prodmarkovscale
! 					equmTRANS%kappa0_w = irfstruct%equmSTICKY%kappa0_w
! 					equmTRANS%kappa1_w = irfstruct%equmSTICKY%kappa1_w
! 					equmTRANS%prefshock = irfstruct%equmSTICKY%prefshock
!
!  				CASE(9) !change lump transfer by direct effect from government interest payments
! 					equmTRANS%lumptransfer = equmINITSS%lumptransfer + (irfstruct%equmSTICKY%rb - equmTRANS%rb)*equmINITSS%govbond
!
!  				CASE(10) !change lump transfer, not including direct effect from govt interest payments
! 					equmTRANS%lumptransfer = irfstruct%equmSTICKY%lumptransfer - (irfstruct%equmSTICKY%rb - equmTRANS%rb)*equmINITSS%govbond
!
! 				CASE(11) !change rb, rborr, and change ra by same amount as rb
! 					equmTRANS%rb = irfstruct%equmSTICKY%rb
! 					equmTRANS%rborr = irfstruct%equmSTICKY%rborr
! 					equmTRANS%ra = equmINITSS%ra + (irfstruct%equmSTICKY%rb - equmINITSS%rb)
!
! 				CASE(12) !change ra, and change rb, rborr by same amount as ra
! 					equmTRANS%ra = irfstruct%equmSTICKY%ra
! 					equmTRANS%rb = irfstruct%equmSTICKY%rb + (irfstruct%equmSTICKY%ra - equmINITSS%ra)
! 					equmTRANS%rborr = irfstruct%equmSTICKY%rborr + (irfstruct%equmSTICKY%ra - equmINITSS%ra)
!
! 				CASE(13) !change only proportional labor tax
! 					equmTRANS%labtax = irfstruct%equmSTICKY%labtax
! 					equmTRANS%netwage = equmINITSS%wage*(1.0-irfstruct%equmSTICKY%labtax)
!
! 				CASE(14) !change rb, rborr, and change ra by same amount as rb, and discount eqm profits at implied ra
! 					equmTRANS%rb = irfstruct%equmSTICKY%rb
! 					equmTRANS%rborr = irfstruct%equmSTICKY%rborr
! 					equmTRANS%ra = equmINITSS%ra + (irfstruct%equmSTICKY%rb - equmINITSS%rb)
! 					IF(DividendFundLumpSum==0) THEN
! 						equmTRANS(:)%illassetdrop = 1.0
! 					ELSE IF(DividendFundLumpSum==1) THEN
! 						it = Ttransition
! 						lequity = (equmFINALSS%equity + irfstruct%equmSTICKY(it)%dividend*deltatransvec(it)) / (1.0+deltatransvec(it)*equmTRANS(it)%ra)
! 						DO it = Ttransition-1,1,-1
! 							lequity = (lequity + irfstruct%equmSTICKY(it)%dividend*deltatransvec(it)) / (1.0+deltatransvec(it)*equmTRANS(it)%ra)
! 						END DO
! 						equmTRANS(:)%illassetdrop = ((1.0-irfstruct%equmSTICKY(1)%fundlev) * irfstruct%equmSTICKY(1)%capital + lequity) / ((1.0-equmINITSS%fundlev)*equmINITSS%capital + equmINITSS%equity)
! 					END IF
!
! 				CASE(15) !change rb, rborr, and change ra by same amount as rb, and discount initimake al profits at implied ra
! 					equmTRANS%rb = irfstruct%equmSTICKY%rb
! 					equmTRANS%rborr = irfstruct%equmSTICKY%rborr
! 					equmTRANS%ra = equmINITSS%ra + (irfstruct%equmSTICKY%rb - equmINITSS%rb)
! 					IF(DividendFundLumpSum==0) THEN
! 						equmTRANS(:)%illassetdrop = 1.0
! 					ELSE IF(DividendFundLumpSum==1) THEN
! 						it = Ttransition
! 						lequity = (equmFINALSS%equity + equmINITSS%dividend*deltatransvec(it)) / (1.0+deltatransvec(it)*equmTRANS(it)%ra)
! 						DO it = Ttransition-1,1,-1
! 							lequity = (lequity + equmINITSS%dividend*deltatransvec(it)) / (1.0+deltatransvec(it)*equmTRANS(it)%ra)
! 						END DO
! 						equmTRANS(:)%illassetdrop = ((1.0-irfstruct%equmSTICKY(1)%fundlev) * irfstruct%equmSTICKY(1)%capital + lequity) / ((1.0-equmINITSS%fundlev)*equmINITSS%capital + equmINITSS%equity)
! 					END IF
!
! 				CASE(16) !only change the direct shock
! 					equmTRANS%labwedge = irfstruct%equmSTICKY%labwedge
! 					equmTRANS%finwedge = irfstruct%equmSTICKY%finwedge
! 					IF(FixBorrowRateTransition==1) equmTRANS%rborr = irfstruct%equmSTICKY%rborr
! 					equmTRANS%prodgridscale = irfstruct%equmSTICKY%prodgridscale
! 					DO it = 1,Ttransition
! 						equmTRANS(it)%ygrid = irfstruct%equmSTICKY(it)%ygrid
! 					END DO
! 					equmTRANS%prodmarkovscale = irfstruct%equmSTICKY%prodmarkovscale
! 					equmTRANS%kappa0_w = irfstruct%equmSTICKY%kappa0_w
! 					equmTRANS%kappa1_w = irfstruct%equmSTICKY%kappa1_w
! 					equmTRANS%prefshock = irfstruct%equmSTICKY%prefshock
!
! 			END SELECT
! 			CALL Transition
!
! 			IF (SaveCumPolicyFnsIRF==1) CALL CumulativeConsTransition
!
! 			irfpointer%equmSTICKY = equmTRANS
! 			irfpointer%statsSTICKY = statsTRANS
! 			irfpointer%solnSTICKY = solnTRANS
! 			irfpointer%cumSTICKY = CumulativePolicyType(ccum1,ccum2,ccum4,dcum1,dcum2,dcum4)
!
! 			IF(ComputeDiscountedMPC==1) CALL DiscountedMPCTransition
!
! 		END IF
!
! 		IF(SolveZLBTransition==1) THEN
! 			write(*,*) "Solving for sticky price transition with ZLB, price experiment ", ipe
!
! 			DO it = 1,Ttransition
! 				equmTRANS(it) = equmINITSS
! 			END DO
!
! 			SELECT CASE(ipe)
! 				CASE(1) !change wage and labor tax only
! 					equmTRANS%wage = irfstruct%equmZLB%wage
! 					equmTRANS%netwage = irfstruct%equmZLB%netwage
! 					equmTRANS%labtax = irfstruct%equmZLB%labtax
!
! 				CASE(2) !only change profits
! 					equmTRANS%profit = irfstruct%equmZLB%profit
!
! 				CASE(3) !only change profits and wage
! 					equmTRANS%profit = irfstruct%equmZLB%profit
! 					equmTRANS%wage = irfstruct%equmZLB%wage
! 					equmTRANS%netwage = irfstruct%equmZLB%netwage
! 					equmTRANS%labtax = irfstruct%equmZLB%labtax
!
! 				CASE(4) !only change rb (and rborr if the wedge is fixed out of steady state)
! 					equmTRANS%rb = irfstruct%equmZLB%rb
! 					IF(FixBorrowRateTransition==0) equmTRANS%rborr = irfstruct%equmZLB%rborr
!
! 				CASE(5) !only change ra only
! 					equmTRANS%ra = irfstruct%equmZLB%ra
!
! 				CASE(6) !only change illiquid asset drop
! 					equmTRANS%illassetdrop  = irfstruct%equmZLB%illassetdrop
!
! 				CASE(7) !only change transfers
! 					equmTRANS%lumptransfer = irfstruct%equmZLB%lumptransfer
!
! 				CASE(8) !change all
! 					equmTRANS%wage = irfstruct%equmZLB%wage
! 					equmTRANS%netwage = irfstruct%equmZLB%netwage
! 					equmTRANS%labtax = irfstruct%equmZLB%labtax
! 					equmTRANS%lumptransfer = irfstruct%equmZLB%lumptransfer
! 					equmTRANS%rb = irfstruct%equmZLB%rb
! 					equmTRANS%rborr = irfstruct%equmZLB%rborr
! 					equmTRANS%ra = irfstruct%equmZLB%ra
! 					equmTRANS%illassetdrop  = irfstruct%equmZLB%illassetdrop
! 					equmTRANS%profit = irfstruct%equmZLB%profit
! 					equmTRANS%labwedge = irfstruct%equmZLB%labwedge
! 					equmTRANS%finwedge = irfstruct%equmZLB%finwedge
! 					equmTRANS%prodgridscale = irfstruct%equmZLB%prodgridscale
! 					DO it = 1,Ttransition
! 						equmTRANS(it)%ygrid = irfstruct%equmZLB(it)%ygrid
! 					END DO
! 					equmTRANS%prodmarkovscale = irfstruct%equmZLB%prodmarkovscale
! 					equmTRANS%kappa0_w = irfstruct%equmZLB%kappa0_w
! 					equmTRANS%kappa1_w = irfstruct%equmZLB%kappa1_w
! 					equmTRANS%prefshock = irfstruct%equmZLB%prefshock
!
!  				CASE(9) !change lump transfer by direct effect from government interest payments
! 					equmTRANS%lumptransfer = equmINITSS%lumptransfer + (irfstruct%equmZLB%rb - equmTRANS%rb)*equmINITSS%govbond
!
!  				CASE(10) !change lump transfer, not including direct effect from govt interest payments
! 					equmTRANS%lumptransfer = irfstruct%equmZLB%lumptransfer - (irfstruct%equmZLB%rb - equmTRANS%rb)*equmINITSS%govbond
!
! 				CASE(11) !change rb, rborr, and change ra by same amount as rb
! 					equmTRANS%rb = irfstruct%equmZLB%rb
! 					equmTRANS%rborr = irfstruct%equmZLB%rborr
! 					equmTRANS%ra = equmINITSS%ra + (irfstruct%equmZLB%rb - equmINITSS%rb)
!
! 				CASE(12) !change ra, and change rb, rborr by same amount as ra
! 					equmTRANS%ra = irfstruct%equmZLB%ra
! 					equmTRANS%rb = irfstruct%equmZLB%rb + (irfstruct%equmZLB%ra - equmINITSS%ra)
! 					equmTRANS%rborr = irfstruct%equmZLB%rborr + (irfstruct%equmZLB%ra - equmINITSS%ra)
!
! 				CASE(13) !change only proportional labor tax
! 					equmTRANS%labtax = irfstruct%equmZLB%labtax
! 					equmTRANS%netwage = equmINITSS%wage*(1.0-irfstruct%equmZLB%labtax)
!
! 				CASE(14) !change rb, rborr, and change ra by same amount as rb, and discount eqm profits at implied ra
! 					equmTRANS%rb = irfstruct%equmZLB%rb
! 					equmTRANS%rborr = irfstruct%equmZLB%rborr
! 					equmTRANS%ra = equmINITSS%ra + (irfstruct%equmZLB%rb - equmINITSS%rb)
! 					IF(DividendFundLumpSum==0) THEN
! 						equmTRANS(:)%illassetdrop = 1.0
! 					ELSE IF(DividendFundLumpSum==1) THEN
! 						it = Ttransition
! 						lequity = (equmFINALSS%equity + irfstruct%equmZLB(it)%dividend*deltatransvec(it)) / (1.0+deltatransvec(it)*equmTRANS(it)%ra)
! 						DO it = Ttransition-1,1,-1
! 							lequity = (lequity + irfstruct%equmZLB(it)%dividend*deltatransvec(it)) / (1.0+deltatransvec(it)*equmTRANS(it)%ra)
! 						END DO
! 						equmTRANS(:)%illassetdrop = ((1.0-irfstruct%equmZLB(1)%fundlev) * irfstruct%equmZLB(1)%capital + lequity) / ((1.0-equmINITSS%fundlev)*equmINITSS%capital + equmINITSS%equity)
! 					END IF
!
! 				CASE(15) !change rb, rborr, and change ra by same amount as rb, and discount initimake al profits at implied ra
! 					equmTRANS%rb = irfstruct%equmZLB%rb
! 					equmTRANS%rborr = irfstruct%equmZLB%rborr
! 					equmTRANS%ra = equmINITSS%ra + (irfstruct%equmZLB%rb - equmINITSS%rb)
! 					IF(DividendFundLumpSum==0) THEN
! 						equmTRANS(:)%illassetdrop = 1.0
! 					ELSE IF(DividendFundLumpSum==1) THEN
! 						it = Ttransition
! 						lequity = (equmFINALSS%equity + equmINITSS%dividend*deltatransvec(it)) / (1.0+deltatransvec(it)*equmTRANS(it)%ra)
! 						DO it = Ttransition-1,1,-1
! 							lequity = (lequity + equmINITSS%dividend*deltatransvec(it)) / (1.0+deltatransvec(it)*equmTRANS(it)%ra)
! 						END DO
! 						equmTRANS(:)%illassetdrop = ((1.0-irfstruct%equmZLB(1)%fundlev) * irfstruct%equmZLB(1)%capital + lequity) / ((1.0-equmINITSS%fundlev)*equmINITSS%capital + equmINITSS%equity)
! 					END IF
!
! 				CASE(16) !only change the direct shock
! 					equmTRANS%labwedge = irfstruct%equmZLB%labwedge
! 					equmTRANS%finwedge = irfstruct%equmZLB%finwedge
! 					IF(FixBorrowRateTransition==1) equmTRANS%rborr = irfstruct%equmZLB%rborr
! 					equmTRANS%prodgridscale = irfstruct%equmZLB%prodgridscale
! 					DO it = 1,Ttransition
! 						equmTRANS(it)%ygrid = irfstruct%equmZLB(it)%ygrid
! 					END DO
! 					equmTRANS%prodmarkovscale = irfstruct%equmZLB%prodmarkovscale
! 					equmTRANS%kappa0_w = irfstruct%equmZLB%kappa0_w
! 					equmTRANS%kappa1_w = irfstruct%equmZLB%kappa1_w
! 					equmTRANS%prefshock = irfstruct%equmZLB%prefshock
!
! 			END SELECT
! 			CALL Transition
!
! 			IF (SaveCumPolicyFnsIRF==1) CALL CumulativeConsTransition
!
! 			irfpointer%equmZLB = equmTRANS
! 			irfpointer%statsZLB = statsTRANS
! 			irfpointer%solnZLB = solnTRANS
! 			irfpointer%cumZLB = CumulativePolicyType(ccum1,ccum2,ccum4,dcum1,dcum2,dcum4)
!
! 			IF(ComputeDiscountedMPC==1) CALL DiscountedMPCTransition
!
! 		END IF
!
! 		irfsave => irfpriceexp
! 		IF(ipe<10) WRITE(UNIT=lstring, FMT='(I1)') ipe
! 		IF(ipe>=10) WRITE(UNIT=lstring, FMT='(I2)') ipe
! 		OutputDirIRF = trim(OutputDir) // "IRF_" // trim(lIRFDir) // "/PE"// trim(lstring) // "/"
! 		CALL SaveIRFOutput
!
! 	END IF
! 	END DO
!
! END IF
testingGHH = .false.


END SUBROUTINE IRFSequence