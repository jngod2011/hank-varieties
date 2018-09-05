SUBROUTINE IRFSequence(lIRFDir)

USE Parameters
USE Globals
USE Procedures

IMPLICIT NONE

CHARACTER, INTENT(IN) :: lIRFDir*20
CHARACTER	:: lstring*80
INTEGER	:: ifs,it,ipe

irfpointer => irfstruct(0)		

fsptransition = .false.
IF(SolveFlexPriceTransition==1) THEN
	flextransition = .true.
	CALL IterateTransitionFlex
	flextransition = .false.				
END IF
IF(SolveStickyPriceTransition==1) THEN
	stickytransition = .true.
	IF(StickyPriceAlgorithm==1) CALL IterateTransitionStickyB
	IF(StickyPriceAlgorithm==2) CALL IterateTransitionStickyPi
	IF(StickyPriceAlgorithm==3) CALL IterateTransitionStickyY
	IF(StickyPriceAlgorithm==4) CALL IterateTransitionStickyRb
	stickytransition = .false.		
END IF	
IF(SolveZLBTransition==1) THEN
	zlbtransition = .true.
	IF(ZLBAlgorithm==1) CALL IterateTransitionStickyB
	IF(ZLBAlgorithm==2) CALL IterateTransitionStickyPi
	IF(ZLBAlgorithm==3) CALL IterateTransitionStickyY
	IF(ZLBAlgorithm==4) CALL IterateTransitionStickyRb
	zlbtransition = .false.
END IF

irfsave => irfstruct(0)
OutputDirIRF = trim(OutputDir) // "IRF_" // trim(lIRFDir) // "/NOFS/"
CALL SaveIRFOutput

IF (DoPriceExperiments==1) THEN
	irfpointer => irfpriceexp
	
	DO ipe = 1,13
		testingGHH = .false.
				
		IF(SolveFlexPriceTransition==1) THEN
			write(*,*) "Solving for flex price transition, price experiment ", ipe
			
			DO it = 1,Ttransition
				equmTRANS(it) = equmINITSS
			END DO
		
			SELECT CASE(ipe)
				CASE(1) !change wage and lump
					equmTRANS%wage = irfstruct(0)%equmFLEX%wage
					equmTRANS%netwage = irfstruct(0)%equmFLEX%netwage
					equmTRANS%lumptransfer = irfstruct(0)%equmFLEX%lumptransfer
					equmTRANS%labor = irfstruct(0)%equmFLEX%labor
				
				CASE(2) !only change rb, rborr
					equmTRANS%rb = irfstruct(0)%equmFLEX%rb
					equmTRANS%rborr = irfstruct(0)%equmFLEX%rborr

				CASE(3) !only change ra
					equmTRANS%ra = irfstruct(0)%equmFLEX%ra

				CASE(4) !keep wage fixed, change ra, rb, rborr
					equmTRANS%ra = irfstruct(0)%equmFLEX%ra
					equmTRANS%rb = irfstruct(0)%equmFLEX%rb
					equmTRANS%rborr = irfstruct(0)%equmFLEX%rborr

				CASE(5) !keep ra fixed, change wage, rb, rborr
					equmTRANS%wage = irfstruct(0)%equmFLEX%wage
					equmTRANS%netwage = irfstruct(0)%equmFLEX%netwage
					equmTRANS%lumptransfer = irfstruct(0)%equmFLEX%lumptransfer
					equmTRANS%labor = irfstruct(0)%equmFLEX%labor
					equmTRANS%rb = irfstruct(0)%equmFLEX%rb
					equmTRANS%rborr = irfstruct(0)%equmFLEX%rborr

				CASE(6) !keep wage fixed, change rb, rborr, change ra by same amount as rb
					equmTRANS%rb = irfstruct(0)%equmFLEX%rb
					equmTRANS%rborr = irfstruct(0)%equmFLEX%rborr
					equmTRANS%ra = equmINITSS%ra + (irfstruct(0)%equmFLEX%rb - equmINITSS%rb)

				CASE(7) !keep wage fixed, change ra, change rb, rborr by same amount as ra
					equmTRANS%ra = irfstruct(0)%equmFLEX%ra
					equmTRANS%rb = irfstruct(0)%equmFLEX%rb + (irfstruct(0)%equmFLEX%ra - equmINITSS%ra)
					equmTRANS%rborr = irfstruct(0)%equmFLEX%rborr + (irfstruct(0)%equmFLEX%ra - equmINITSS%ra)

				CASE(8) !change all 3
					equmTRANS%wage = irfstruct(0)%equmFLEX%wage
					equmTRANS%netwage = irfstruct(0)%equmFLEX%netwage
					equmTRANS%lumptransfer = irfstruct(0)%equmFLEX%lumptransfer
					equmTRANS%labor = irfstruct(0)%equmFLEX%labor
					equmTRANS%rb = irfstruct(0)%equmFLEX%rb
					equmTRANS%rborr = irfstruct(0)%equmFLEX%rborr
					equmTRANS%ra = irfstruct(0)%equmFLEX%ra
				
 				CASE(9) !only change wage, not lump transfer
 					equmTRANS%wage = irfstruct(0)%equmFLEX%wage
 					equmTRANS%netwage = irfstruct(0)%equmFLEX%netwage
					equmTRANS%labor = irfstruct(0)%equmFLEX%labor

 				CASE(10) !only lump transfer, not including direct effect from govt interest payments
					equmTRANS%lumptransfer = irfstruct(0)%equmFLEX%lumptransfer - (irfstruct(0)%equmFLEX%rb - equmTRANS%rb)*equmINITSS%govbond

 				CASE(11) !change lump transfer by direct effect from government interest payments
					equmTRANS%lumptransfer = equmINITSS%lumptransfer + (irfstruct(0)%equmFLEX%rb - equmTRANS%rb)*equmINITSS%govbond

 				CASE(12) !only lump transfer by indirect effect from change in tax revenues
					equmTRANS%lumptransfer = equmINITSS%lumptransfer + (irfstruct(0)%equmFLEX%labtax*irfstruct(0)%equmFLEX%wage*irfstruct(0)%equmFLEX%labor &
											- equmINITSS%labtax*equmINITSS%wage*equmINITSS%labor)

				CASE(13) !only change wage, not lump transfer, but do not change disutility
					equmTRANS%wage = irfstruct(0)%equmFLEX%wage
					equmTRANS%netwage = irfstruct(0)%equmFLEX%netwage
					equmTRANS%labor = irfstruct(0)%equmFLEX%labor
					testingGHH = .true.

			END SELECT	
			CALL Transition

			irfpointer%equmFLEX = equmTRANS
			irfpointer%statsFLEX = statsTRANS
			irfpointer%solnFLEX = solnTRANS
				
		END IF
		
		IF(SolveStickyPriceTransition==1) THEN
			write(*,*) "Solving for sticky price transition without ZLB, price experiment ", ipe
		
			DO it = 1,Ttransition
				equmTRANS(it) = equmINITSS
			END DO

			SELECT CASE(ipe)
				CASE(1) !only change wage
					equmTRANS%wage = irfstruct(0)%equmSTICKY%wage
					equmTRANS%netwage = irfstruct(0)%equmSTICKY%netwage
					equmTRANS%lumptransfer = irfstruct(0)%equmSTICKY%lumptransfer
					equmTRANS%labor = irfstruct(0)%equmSTICKY%labor

				CASE(2) !only change rb, rborr
					equmTRANS%rb = irfstruct(0)%equmSTICKY%rb
					equmTRANS%rborr = irfstruct(0)%equmSTICKY%rborr

				CASE(3) !only change ra
					equmTRANS%ra = irfstruct(0)%equmSTICKY%ra
					
				CASE(4) !keep wage fixed, change ra, rb, rborr
					equmTRANS%ra = irfstruct(0)%equmSTICKY%ra
					equmTRANS%rb = irfstruct(0)%equmSTICKY%rb
					equmTRANS%rborr = irfstruct(0)%equmSTICKY%rborr

				CASE(5) !keep ra fixed, change wage, rb, rborr
					equmTRANS%wage = irfstruct(0)%equmSTICKY%wage
					equmTRANS%netwage = irfstruct(0)%equmSTICKY%netwage
					equmTRANS%lumptransfer = irfstruct(0)%equmSTICKY%lumptransfer
					equmTRANS%labor = irfstruct(0)%equmSTICKY%labor
					equmTRANS%rb = irfstruct(0)%equmSTICKY%rb
					equmTRANS%rborr = irfstruct(0)%equmSTICKY%rborr
					
				CASE(6) !keep wage fixed, change rb, rborr, change ra by same amount as rb
					equmTRANS%rb = irfstruct(0)%equmSTICKY%rb
					equmTRANS%rborr = irfstruct(0)%equmSTICKY%rborr
					equmTRANS%ra = equmINITSS%ra + (irfstruct(0)%equmSTICKY%rb - equmINITSS%rb)

				CASE(7) !keep wage fixed, change ra, change rb, rborr by same amount as ra
					equmTRANS%ra = irfstruct(0)%equmSTICKY%ra
					equmTRANS%rb = irfstruct(0)%equmSTICKY%rb + (irfstruct(0)%equmSTICKY%ra - equmINITSS%ra)
					equmTRANS%rborr = irfstruct(0)%equmSTICKY%rborr + (irfstruct(0)%equmSTICKY%ra - equmINITSS%ra)

				CASE(8) !change all 3
					equmTRANS%wage = irfstruct(0)%equmSTICKY%wage
					equmTRANS%netwage = irfstruct(0)%equmSTICKY%netwage
					equmTRANS%lumptransfer = irfstruct(0)%equmSTICKY%lumptransfer
					equmTRANS%labor = irfstruct(0)%equmSTICKY%labor
					equmTRANS%rb = irfstruct(0)%equmSTICKY%rb
					equmTRANS%rborr = irfstruct(0)%equmSTICKY%rborr
					equmTRANS%ra = irfstruct(0)%equmSTICKY%ra

 				CASE(9) !only change wage, not lump transfer
 					equmTRANS%wage = irfstruct(0)%equmSTICKY%wage
 					equmTRANS%netwage = irfstruct(0)%equmSTICKY%netwage
					equmTRANS%labor = irfstruct(0)%equmSTICKY%labor

 				CASE(10) !only lump transfer, not including direct effect from govt interest payments
					equmTRANS%lumptransfer = irfstruct(0)%equmSTICKY%lumptransfer - (irfstruct(0)%equmSTICKY%rb - equmTRANS%rb)*equmINITSS%govbond

 				CASE(11) !change lump transfer by direct effect from government interest payments
					equmTRANS%lumptransfer = equmINITSS%lumptransfer + (irfstruct(0)%equmSTICKY%rb - equmTRANS%rb)*equmINITSS%govbond
				
 				CASE(12) !only lump transfer by indirect effect from change in tax revenues
					equmTRANS%lumptransfer = equmINITSS%lumptransfer + (irfstruct(0)%equmSTICKY%labtax*irfstruct(0)%equmSTICKY%wage*irfstruct(0)%equmSTICKY%labor &
												- equmINITSS%labtax*equmINITSS%wage*equmINITSS%labor)
				
				CASE(13) !only change wage, not lump transfer, but do not change disutility
					equmTRANS%wage = irfstruct(0)%equmSTICKY%wage
					equmTRANS%netwage = irfstruct(0)%equmSTICKY%netwage
					equmTRANS%labor = irfstruct(0)%equmSTICKY%labor
					testingGHH = .true.
								
				
			END SELECT	
			CALL Transition

			irfpointer%equmSTICKY = equmTRANS
			irfpointer%statsSTICKY = statsTRANS
			irfpointer%solnSTICKY = solnTRANS

		END IF

		IF(SolveZLBTransition==1) THEN
			DO it = 1,Ttransition
				equmTRANS(it) = equmINITSS
			END DO

			SELECT CASE(ipe)
				CASE(1) !only change wage
					equmTRANS%wage = irfstruct(0)%equmZLB%wage
					equmTRANS%netwage = irfstruct(0)%equmZLB%netwage
					equmTRANS%lumptransfer = irfstruct(0)%equmZLB%lumptransfer
					equmTRANS%labor = irfstruct(0)%equmZLB%labor
		
				CASE(2) !only change rb, rborr
					equmTRANS%rb = irfstruct(0)%equmZLB%rb
					equmTRANS%rborr = irfstruct(0)%equmZLB%rborr

				CASE(3) !only change ra
					equmTRANS%ra = irfstruct(0)%equmZLB%ra
					
				CASE(4) !keep wage fixed, change ra, rb, rborr
					equmTRANS%ra = irfstruct(0)%equmZLB%ra
					equmTRANS%rb = irfstruct(0)%equmZLB%rb
					equmTRANS%rborr = irfstruct(0)%equmZLB%rborr

				CASE(5) !keep ra fixed, change wage, rb, rborr
					equmTRANS%wage = irfstruct(0)%equmZLB%wage
					equmTRANS%netwage = irfstruct(0)%equmZLB%netwage
					equmTRANS%lumptransfer = irfstruct(0)%equmZLB%lumptransfer
					equmTRANS%labor = irfstruct(0)%equmZLB%labor
					equmTRANS%rb = irfstruct(0)%equmZLB%rb
					equmTRANS%rborr = irfstruct(0)%equmZLB%rborr
					
				CASE(6) !keep wage fixed, change rb, rborr, change ra by same amount as rb
					equmTRANS%rb = irfstruct(0)%equmZLB%rb
					equmTRANS%rborr = irfstruct(0)%equmZLB%rborr
					equmTRANS%ra = equmINITSS%ra + (irfstruct(0)%equmZLB%rb - equmINITSS%rb)

				CASE(7) !keep wage fixed, change ra, change rb, rborr by same amount as ra
					equmTRANS%ra = irfstruct(0)%equmZLB%ra
					equmTRANS%rb = irfstruct(0)%equmZLB%rb + (irfstruct(0)%equmZLB%ra - equmINITSS%ra)
					equmTRANS%rborr = irfstruct(0)%equmZLB%rborr + (irfstruct(0)%equmZLB%ra - equmINITSS%ra)	

				CASE(8) !change all 3
					equmTRANS%wage = irfstruct(0)%equmZLB%wage
					equmTRANS%netwage = irfstruct(0)%equmZLB%netwage
					equmTRANS%lumptransfer = irfstruct(0)%equmZLB%lumptransfer
					equmTRANS%labor = irfstruct(0)%equmZLB%labor
					equmTRANS%rb = irfstruct(0)%equmZLB%rb
					equmTRANS%rborr = irfstruct(0)%equmZLB%rborr
					equmTRANS%ra = irfstruct(0)%equmZLB%ra

 				CASE(9) !only change wage, not lump transfer
 					equmTRANS%wage = irfstruct(0)%equmZLB%wage
 					equmTRANS%netwage = irfstruct(0)%equmZLB%netwage
					equmTRANS%labor = irfstruct(0)%equmZLB%labor

 				CASE(10) !only lump transfer, not including direct effect from govt interest payments
					equmTRANS%lumptransfer = irfstruct(0)%equmZLB%lumptransfer - (irfstruct(0)%equmZLB%rb - equmTRANS%rb)*equmINITSS%govbond

 				CASE(11) !change lump transfer by direct effect from government interest payments
					equmTRANS%lumptransfer = equmINITSS%lumptransfer + (irfstruct(0)%equmZLB%rb - equmTRANS%rb)*equmINITSS%govbond

				CASE(12) !only lump transfer by indirect effect from change in tax revenues
					equmTRANS%lumptransfer = equmINITSS%lumptransfer + (irfstruct(0)%equmZLB%labtax*irfstruct(0)%equmZLB%wage*irfstruct(0)%equmZLB%labor &
											- equmINITSS%labtax*equmINITSS%wage*equmINITSS%labor)

				CASE(13) !only change wage, not lump transfer, but do not change disutility
					equmTRANS%wage = irfstruct(0)%equmZLB%wage
					equmTRANS%netwage = irfstruct(0)%equmZLB%netwage
					equmTRANS%labor = irfstruct(0)%equmZLB%labor
					testingGHH = .true.
									
			END SELECT	
			CALL Transition

			irfpointer%equmZLB = equmTRANS
			irfpointer%statsZLB = statsTRANS
			irfpointer%solnZLB = solnTRANS
		
		END IF
				
		irfsave => irfpriceexp
		IF(ipe<10) WRITE(UNIT=lstring, FMT='(I1)') ipe
		IF(ipe>=10) WRITE(UNIT=lstring, FMT='(I2)') ipe
		OutputDirIRF = trim(OutputDir) // "IRF_" // trim(lIRFDir) // "/PE"// trim(lstring) // "/"		
		CALL SaveIRFOutput
		
	END DO
END IF
testingGHH = .false.


IF (DoFiscalStimulus==1) THEN
	fsptransition = .true.
	DO ifs = 1,nfs
		irfpointer_fs => irfstruct(ifs)	
		fspointer => fsconfig(ifs)

		IF(SolveFlexPriceTransition==1) THEN
			flextransition = .true.
			CALL IterateTransitionFlex
			flextransition = .false.				
		END IF
		IF(SolveStickyPriceTransition==1) THEN
			stickytransition = .true.
			IF(StickyPriceAlgorithm==1) CALL IterateTransitionStickyB
			IF(StickyPriceAlgorithm==2) CALL IterateTransitionStickyPi
			stickytransition = .false.		
		END IF	
		IF(SolveZLBTransition==1) THEN
			zlbtransition = .true.
			IF(ZLBAlgorithm==1) CALL IterateTransitionStickyB
			IF(ZLBAlgorithm==2) CALL IterateTransitionStickyPi
			zlbtransition = .false.
		END IF

		irfsave => irfstruct(ifs)
		WRITE(UNIT=lstring, FMT='(I1)') ifs
		OutputDirIRF = trim(OutputDir) // "IRF_" // trim(lIRFDir) // "/FS"// trim(lstring) // "/"		
		CALL SaveIRFOutput
	END DO
	fsptransition = .false.
END IF

END SUBROUTINE IRFSequence