SUBROUTINE ImpulseResponses

USE Parameters
USE Globals
USE Procedures

IMPLICIT NONE
INTEGER		:: it,it1
CHARACTER	:: IRFDir*20

!allocate arrays
CALL AllocateArrays

flextransition = .false.
zlbtransition = .false.
stickytransition = .false.
fsptransition = .false.

!set up deltatransvec
CALL PowerSpacedGrid (Ttransition,deltatransparam,deltatransmin,deltatransmax,deltatransvec)

nendtrans = min(50,Ttransition)
cumdeltatrans(1) = deltatransvec(1)
DO it = 2,Ttransition
	cumdeltatrans(it) = cumdeltatrans(it-1) + deltatransvec(it)
END DO
OPEN(3, FILE = trim(OutputDir) // 'deltatransvec.txt', STATUS = 'replace'); CALL WriteMatrix(3,Ttransition,1,deltatransvec)

IF(DoFiscalStimulus==1) CALL SetupFiscalStimulus

!TFP shock
IF(IncludeTFPShock==1) THEN
	IF(Display>=1) write(*,*)'Solving for TFP shock IRF'	
	equmTRANS(:) = equmINITSS	
	
	equmTRANS(1)%tfp = equmINITSS%tfp * exp(TFPShockSize)
	DO it = 2,Ttransition-nendtrans
		equmTRANS(it)%tfp = equmTRANS(it-1)%tfp ** (TFPShockPers**deltatransvec(it-1))
	END DO	
	equmTRANS(Ttransition-nendtrans+1:Ttransition)%tfp = equmINITSS%tfp

	IRFDir = "TFP"
	CALL IRFSequence(IRFDir)
END IF

!News shock: this needs to be fixed up because of the way that deltatransvec works - it=1 is NOT 1 quarter in (see forward guidance shock)
IF (IncludeNewsShock==1) THEN

	IF(Display>=1) write(*,*)'Solving for news shock IRF'
	
	equmTRANS(:) = equmINITSS	
	equmTRANS(1:NewsShockQtrs)%tfp = equmINITSS%tfp
	equmTRANS(NewsShockQtrs+1)%tfp = equmINITSS%tfp * exp(NewsShockSize)
	DO it = NewsShockQtrs+2,Ttransition-nendtrans
		equmTRANS(it)%tfp = equmTRANS(it-1)%tfp ** (NewsShockPers**deltatransvec(it-1))
	END DO	
	equmTRANS(Ttransition-nendtrans+1:Ttransition)%tfp = equmINITSS%tfp

	IRFDir = "News"
	CALL IRFSequence(IRFDir)

END IF

!Kappa fc shock
IF(IncludeKappafcShock==1) THEN
	IF(Display>=1) write(*,*)'Solving for fixed cost shock IRF'	
	
	equmTRANS(:) = equmINITSS		
	equmTRANS(1)%kappafc_w = equmINITSS%kappafc_w * exp(KappafcShockSize)
	DO it = 2,Ttransition-nendtrans
		equmTRANS(it)%kappafc_w =(equmINITSS%kappafc_w **(1.0-KappafcShockPers**deltatransvec(it-1))) * (equmTRANS(it-1)%kappafc_w ** (KappafcShockPers**deltatransvec(it-1)))
		
	END DO	
	equmTRANS(Ttransition-nendtrans+1:Ttransition)%kappafc_w = equmINITSS%kappafc_w

	IRFDir = "Kappafc"
	CALL IRFSequence(IRFDir)

END IF

!Monetary policy shock
IF(IncludeMonetaryShock==1) THEN
	IF(Display>=1) write(*,*)'Solving for monetary policy shock IRF'	
	irfpointer => irfstruct(0)

	equmTRANS(:) = equmINITSS		
	equmTRANS(1)%mpshock = equmINITSS%mpshock + MonetaryShockSize
	DO it = 2,Ttransition-nendtrans
		equmTRANS(it)%mpshock =equmINITSS%mpshock *(1.0-MonetaryShockPers**deltatransvec(it-1)) + equmTRANS(it-1)%mpshock * (MonetaryShockPers**deltatransvec(it-1))

	END DO	
	equmTRANS(Ttransition-nendtrans+1:Ttransition)%mpshock = equmINITSS%mpshock

	IRFDir = "Monetary"
	CALL IRFSequence(IRFDir)

END IF

!Forward Guidance shock
IF(IncludeForwardGuideShock==1) THEN
	forwardguide = .true.
	IF(Display>=1) write(*,*)'Solving for forward guidance shock IRF'	
	irfpointer => irfstruct(0)

	equmTRANS(:) = equmINITSS
	it1 = MINLOC(cumdeltatrans, 1, MASK = cumdeltatrans>=ForwardGuideShockQtrs)
	equmTRANS(1:it1-1)%mpshock = equmINITSS%mpshock
	equmTRANS(it1)%mpshock = equmINITSS%mpshock + ForwardGuideShockSize
	DO it = it1+1,Ttransition-nendtrans
		equmTRANS(it)%mpshock =equmINITSS%mpshock *(1.0-MonetaryShockPers**deltatransvec(it-1)) + equmTRANS(it-1)%mpshock * (MonetaryShockPers**deltatransvec(it-1))

	END DO	
	equmTRANS(Ttransition-nendtrans+1:Ttransition)%mpshock = equmINITSS%mpshock

	IRFDir = "ForwardGuide"
	CALL IRFSequence(IRFDir)
	forwardguide = .false.

END IF



! StDevYShockSize

!Preference shock
IF(IncludePrefShock==1) THEN
	IF(Display>=1) write(*,*)'Solving for Preference shock IRF'	
	
	equmTRANS(:) = equmINITSS	
	equmTRANS(1)%prefshock = equmINITSS%prefshock * exp(PrefShockSize)
	DO it = 2,Ttransition-nendtrans
		equmTRANS(it)%prefshock = equmTRANS(it-1)%prefshock ** (PrefShockPers**deltatransvec(it-1))
	END DO	
	equmTRANS(Ttransition-nendtrans+1:Ttransition)%prefshock = equmINITSS%prefshock
	
	IRFDir = "Pref"
	CALL IRFSequence(IRFDir)

	
END IF

!Fund Leverage shock
IF(IncludeFundLevShock==1) THEN
	IF(Display>=1) write(*,*)'Solving for fund leverage shock IRF'	

	equmTRANS(:) = equmINITSS		
	equmTRANS(1)%fundlev = equmINITSS%fundlev + FundLevShockSize
	DO it = 2,Ttransition-nendtrans
		equmTRANS(it)%fundlev =equmINITSS%fundlev *(1.0-FundLevShockPers**deltatransvec(it-1)) + equmTRANS(it-1)%fundlev * (FundLevShockPers**deltatransvec(it-1))

	END DO	
	equmTRANS(Ttransition-nendtrans+1:Ttransition)%fundlev = equmINITSS%fundlev

	IRFDir = "FundLev"
	CALL IRFSequence(IRFDir)

END IF

!Risk Aversion shock
IF(IncludeRiskAversionShock==1) THEN
	IF(Display>=1) write(*,*)'Solving for risk aversion shock IRF'	

	equmTRANS(:) = equmINITSS		
	equmTRANS(1)%gam = equmINITSS%gam + RiskAversionShockSize
	DO it = 2,Ttransition-nendtrans
		equmTRANS(it)%gam =equmINITSS%gam *(1.0-RiskAversionShockPers**deltatransvec(it-1)) + equmTRANS(it-1)%gam * (RiskAversionShockPers**deltatransvec(it-1))

	END DO	
	equmTRANS(Ttransition-nendtrans+1:Ttransition)%gam = equmINITSS%gam

	IRFDir = "RiskAversion"
	CALL IRFSequence(IRFDir)

END IF

!Markup shock
IF(IncludeMarkupShock==1) THEN
	IF(Display>=1) write(*,*)'Solving for markup shock IRF'	

	equmTRANS(:) = equmINITSS		
	equmTRANS(1)%elast = 1.0/equmINITSS%elast + MarkupShockSize
	equmTRANS(1)%elast = 1.0/equmTRANS(1)%elast
	DO it = 2,Ttransition-nendtrans
		equmTRANS(it)%elast = (1.0/equmINITSS%elast) *(1.0-MarkupShockPers**deltatransvec(it-1)) + (1.0/equmTRANS(it-1)%elast) * (MarkupShockPers**deltatransvec(it-1))
		equmTRANS(it)%elast = 1.0/equmTRANS(it)%elast
	END DO	
	equmTRANS(Ttransition-nendtrans+1:Ttransition)%elast = equmINITSS%elast

	IRFDir = "Markup"
	CALL IRFSequence(IRFDir)

END IF

!Borrowing wedge shock
IF(IncludeBorrWedgeShock==1) THEN
	IF(Display>=1) write(*,*)'Solving for borrrowing wedge shock IRF'
	equmTRANS(:) = equmINITSS	
	
	equmTRANS(1)%borrwedge = equmINITSS%borrwedge + BorrWedgeShockSize
	DO it = 2,Ttransition-nendtrans
		equmTRANS(it)%borrwedge =equmINITSS%borrwedge *(1.0-BorrWedgeShockPers**deltatransvec(it-1)) + equmTRANS(it-1)%borrwedge * (BorrWedgeShockPers**deltatransvec(it-1))
	END DO	
	equmTRANS(Ttransition-nendtrans+1:Ttransition)%borrwedge = equmINITSS%borrwedge

	IRFDir = "BorrWedge"
	CALL IRFSequence(IRFDir)

END IF


END SUBROUTINE ImpulseResponses