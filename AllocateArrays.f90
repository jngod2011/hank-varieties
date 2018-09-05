SUBROUTINE AllocateArrays

USE Parameters
USE Globals
USE Procedures

IMPLICIT NONE

INTEGER	:: it,ifs

write(*,*) 'Allocating arrays for transition'

DO it = 1,Ttransition
 	CALL AllocateSolutionType(solnTRANS(it))
END DO

DO it = 1,Ttransition
	ifs = 0
	IF(SolveFlexPriceTransition==1) CALL AllocateSolutionType(irfstruct(ifs)%solnFLEX(it))
	IF(SolveStickyPriceTransition==1) CALL AllocateSolutionType(irfstruct(ifs)%solnSTICKY(it))
	IF(SolveZLBTransition==1) CALL AllocateSolutionType(irfstruct(ifs)%solnZLB(it))

	IF (DoPriceExperiments==1) THEN
		IF(SolveFlexPriceTransition==1) CALL AllocateSolutionType(irfpriceexp%solnFLEX(it))
		IF(SolveStickyPriceTransition==1) CALL AllocateSolutionType(irfpriceexp%solnSTICKY(it))
		IF(SolveZLBTransition==1) CALL AllocateSolutionType(irfpriceexp%solnZLB(it))
	END IF
	
	IF (DoFiscalStimulus==1) THEN
		DO ifs = 1,nfs
			IF(SolveFlexPriceTransition==1) CALL AllocateSolutionType(irfstruct(ifs)%solnFLEX(it))
			IF(SolveStickyPriceTransition==1) CALL AllocateSolutionType(irfstruct(ifs)%solnSTICKY(it))
			IF(SolveZLBTransition==1) CALL AllocateSolutionType(irfstruct(ifs)%solnZLB(it))
		END DO
	END IF
	
END DO	


write(*,*) 'Finished allocating arrays'

END SUBROUTINE