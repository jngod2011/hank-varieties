SUBROUTINE AllocateArrays

USE Parameters
USE Globals
USE Procedures

IMPLICIT NONE

INTEGER	:: it

write(*,*) 'Allocating arrays for transition'

DO it = 1,Ttransition
 	CALL AllocateSolutionType(solnTRANS(it))
END DO

DO it = 1,Ttransition
	IF(SolveFlexPriceTransition==1) CALL AllocateSolutionType(irfstruct%solnFLEX(it))
	IF(SolveStickyPriceTransition==1) CALL AllocateSolutionType(irfstruct%solnSTICKY(it))
	IF(SolveZLBTransition==1) CALL AllocateSolutionType(irfstruct%solnZLB(it))

	IF (DoPriceExperiments==1) THEN
		IF(SolveFlexPriceTransition==1) CALL AllocateSolutionType(irfpriceexp%solnFLEX(it))
		IF(SolveStickyPriceTransition==1) CALL AllocateSolutionType(irfpriceexp%solnSTICKY(it))
		IF(SolveZLBTransition==1) CALL AllocateSolutionType(irfpriceexp%solnZLB(it))
	END IF
	
END DO	


IF (SaveCumPolicyFnsIRF==1)	THEN
	IF(SolveFlexPriceTransition==1) CALL AllocateCumulativePolicyType(irfstruct%cumFLEX)
	IF(SolveStickyPriceTransition==1) CALL AllocateCumulativePolicyType(irfstruct%cumSTICKY)
	IF(SolveZLBTransition==1) CALL AllocateCumulativePolicyType(irfstruct%cumZLB)

	IF (DoPriceExperiments==1) THEN
		IF(SolveFlexPriceTransition==1) CALL AllocateCumulativePolicyType(irfpriceexp%cumFLEX)
		IF(SolveStickyPriceTransition==1) CALL AllocateCumulativePolicyType(irfpriceexp%cumSTICKY)
		IF(SolveZLBTransition==1) CALL AllocateCumulativePolicyType(irfpriceexp%cumZLB)
	END IF
END IF	


write(*,*) 'Finished allocating arrays'

END SUBROUTINE