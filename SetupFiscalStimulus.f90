SUBROUTINE SetupFiscalStimulus

USE Parameters
USE Globals

IMPLICIT NONE
INTEGER	:: ifs

!allocate fiscal stimulus
DO ifs = 1,nfs
	ALLOCATE(fsconfig(ifs)%fspamount(ngpa,ngpb,ngpy))
END DO

!configuration 1: lump sum transfers
fsconfig(1)%fspamount 	= 0.05
fsconfig(1)%fspstart 	= 0.0	!quarters into the transtion. if less than deltatransmin then it will be on impact. for now must be zero
fsconfig(1)%fspend 		= 1.0	!quarters into the transtion. if less than deltatransmin then it will be on impact. for now must be zero
fsconfig(1)%labtaxstart = 40	!quarters
fsconfig(1)%labtaxend 	= 60 !quarters
fsconfig(1)%govamount 	= 0.0
fsconfig(1)%govstart 	= 0.0
fsconfig(1)%govend 		= 0.0

!configuration 2: government spending
fsconfig(2)%fspamount 	= 0.0
fsconfig(2)%fspstart 	= 0.0	!quarters into the transtion. if less than deltatransmin then it will be on impact. for now must be zero
fsconfig(2)%fspend 		= 0.0	!quarters into the transtion. if less than deltatransmin then it will be on impact. for now must be zero
fsconfig(2)%labtaxstart = 40	!quarters
fsconfig(2)%labtaxend 	= 60 !quarters
fsconfig(2)%govamount 	= 0.05
fsconfig(2)%govstart 	= 0.0
fsconfig(2)%govend 		= 1.0

END SUBROUTINE