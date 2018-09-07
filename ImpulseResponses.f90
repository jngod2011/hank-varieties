SUBROUTINE ImpulseResponses

USE Parameters
USE Globals
USE Procedures

IMPLICIT NONE
INTEGER		:: it,it1
CHARACTER	:: IRFDir*20
REAL(8) 	:: linitchi,lminygrid
REAL(8), DIMENSION(Ttransition) 	:: lchi,lmpshock1,lmpshock2

!allocate arrays
CALL AllocateArrays

flextransition = .false.
zlbtransition = .false.
stickytransition = .false.
monetaryshock = .false.
forwardguide = .false.
prodispshock = .false.

! nendtrans = min(50,Ttransition)
nendtrans = min(10,Ttransition)

!set up deltatransvec
IF (TransitionTimeStepType==1) THEN
	CALL PowerSpacedGrid (Ttransition,deltatransparam,deltatransmin,deltatransmax,deltatransvec)

	cumdeltatrans(1) = deltatransvec(1)
	DO it = 2,Ttransition
		cumdeltatrans(it) = cumdeltatrans(it-1) + deltatransvec(it)
	END DO
ELSEIF (TransitionTimeStepType==2) THEN
	cumdeltatrans(1) = 1.0/real(Ttransition+1) 
	DO it = 2,Ttransition
		cumdeltatrans(it) = cumdeltatrans(it-1) + 1.0/real(Ttransition+1) 
	END DO
	
	cumdeltatrans = (cumdeltatrans/(1.0-cumdeltatrans)) / deltatransnu
	deltatransvec(1) = cumdeltatrans(1)
	deltatransvec(2:Ttransition) = cumdeltatrans(2:Ttransition) - cumdeltatrans(1:Ttransition-1)

END IF
OPEN(3, FILE = trim(OutputDir) // 'deltatransvec.txt', STATUS = 'replace'); CALL WriteMatrix(3,Ttransition,1,deltatransvec)

!TFP shock
IF(IncludeTFPShock==1) THEN
	IF(Display>=1) write(*,*)'Solving for TFP shock IRF'	
	equmTRANS(:) = equmINITSS	
	
	equmTRANS(1)%tfp_Y = equmINITSS%tfp_Y * exp(TFPShockSize)
	DO it = 2,Ttransition-nendtrans
		equmTRANS(it)%tfp_Y = equmTRANS(it-1)%tfp_Y ** (TFPShockPers**deltatransvec(it-1))
	END DO	
	equmTRANS(Ttransition-nendtrans+1:Ttransition)%tfp_Y = equmINITSS%tfp_Y

	IRFDir = "TFP"
	CALL IRFSequence(IRFDir)
END IF

!News shock
IF (IncludeNewsShock==1) THEN

	IF(Display>=1) write(*,*)'Solving for news shock IRF'
	
	equmTRANS(:) = equmINITSS	
	it1 = MINLOC(cumdeltatrans, 1, MASK = cumdeltatrans>=NewsShockQtrs)
	equmTRANS(1:it1-1)%tfp_Y = equmINITSS%tfp_Y
	equmTRANS(it1)%tfp_Y = equmINITSS%tfp_Y * exp(NewsShockSize)
	DO it = it1+1,Ttransition-nendtrans
		equmTRANS(it)%tfp_Y = equmTRANS(it-1)%tfp_Y ** (NewsShockPers**deltatransvec(it-1))
	END DO	
	equmTRANS(Ttransition-nendtrans+1:Ttransition)%tfp_Y = equmINITSS%tfp_Y

	IRFDir = "News"
	CALL IRFSequence(IRFDir)
	
END IF

!Kappa 0 shock
IF(IncludeKappa0Shock==1) THEN
	IF(Display>=1) write(*,*)'Solving for fixed cost shock IRF'	
	
	equmTRANS(:) = equmINITSS		
	equmTRANS(1)%kappa0_w = equmINITSS%kappa0_w  + Kappa0ShockSize
	DO it = 2,Ttransition-nendtrans
		equmTRANS(it)%kappa0_w = equmINITSS%kappa0_w *(1.0-Kappa0ShockPers**deltatransvec(it-1)) + equmTRANS(it-1)%kappa0_w * (Kappa0ShockPers**deltatransvec(it-1))
	END DO	
	equmTRANS(Ttransition-nendtrans+1:Ttransition)%kappa0_w = equmINITSS%kappa0_w

	IRFDir = "Kappa0"
	CALL IRFSequence(IRFDir)

END IF

!Kappa 1 shock
IF(IncludeKappa1Shock==1) THEN
	IF(Display>=1) write(*,*)'Solving for Adjustment cost shock IRF'	
	
	equmTRANS(:) = equmINITSS		
	linitchi = (equmINITSS%kappa1_w**(-kappa2_w)) / (1.0+kappa2_w)
	lchi(1) = linitchi * exp(Kappa1ShockSize)
	DO it = 2,Ttransition-nendtrans
! 		lchi(it) = linitchi *(1.0-Kappa1ShockPers**deltatransvec(it-1)) + lchi(it-1) * (Kappa1ShockPers**deltatransvec(it-1))
		lchi(it) = lchi(it-1) ** (Kappa1ShockPers**deltatransvec(it-1))		
	END DO	
	lchi(Ttransition-nendtrans+1:Ttransition) = linitchi
	equmTRANS(:)%kappa1_w = ((1.0+kappa2_w)*lchi(:)) ** (-1.0/kappa2_w)
	
	IRFDir = "Kappa1"
	CALL IRFSequence(IRFDir)

END IF

!Monetary policy shock
IF(IncludeMonetaryShock==1) THEN
	monetaryshock =.true.
	IF(Display>=1) write(*,*)'Solving for monetary policy shock IRF'	

	equmTRANS(:) = equmINITSS		
	equmTRANS(1)%mpshock = equmINITSS%mpshock + MonetaryShockSize
	DO it = 2,Ttransition-nendtrans
		equmTRANS(it)%mpshock =equmINITSS%mpshock *(1.0-MonetaryShockPers**deltatransvec(it-1)) + equmTRANS(it-1)%mpshock * (MonetaryShockPers**deltatransvec(it-1))

	END DO	
	equmTRANS(Ttransition-nendtrans+1:Ttransition)%mpshock = equmINITSS%mpshock

	IRFDir = "Monetary"
	CALL IRFSequence(IRFDir)
	monetaryshock =.false.

END IF

!Forward Guidance shock
IF(IncludeForwardGuideShock==1) THEN
	forwardguide = .true.
	monetaryshock= .true.
	IF(Display>=1) write(*,*)'Solving for forward guidance shock IRF'	

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
	monetaryshock =.false.

END IF

!Taylor rule path shock
IF(IncludeTaylorPathShock==1) THEN
	monetaryshock =.true.
	IF(Display>=1) write(*,*)'Solving for taylor rule path shock IRF'	

	equmTRANS(:) = equmINITSS		

	lmpshock1(1) = TaylorPathShockSize1
	lmpshock2(1) = TaylorPathShockSize2
	DO it = 2,Ttransition-nendtrans
		lmpshock1(it) = lmpshock1(it-1) * (TaylorPathShockPers1**deltatransvec(it-1))
		lmpshock2(it) = lmpshock2(it-1) * (TaylorPathShockPers2**deltatransvec(it-1))

	END DO	
	
	lmpshock1(Ttransition-nendtrans+1:Ttransition) = 0.0
	lmpshock2(Ttransition-nendtrans+1:Ttransition) = 0.0

	equmTRANS(:)%mpshock = equmINITSS%mpshock + lmpshock1 + lmpshock2


	IRFDir = "TaylorPath"
	CALL IRFSequence(IRFDir)
	monetaryshock =.false.

END IF

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

!Government Expenditure Shock
IF(IncludeGovExpShock==1) THEN
	IF(Display>=1) write(*,*)'Solving for govt expenditure shock IRF'	
	equmTRANS(:) = equmINITSS	
	
	equmTRANS(1)%govshock = equmINITSS%govshock * exp(GovExpShockSize)
	DO it = 2,Ttransition-nendtrans
		equmTRANS(it)%govshock = equmTRANS(it-1)%govshock ** (GovExpShockPers**deltatransvec(it-1))
	END DO	
	equmTRANS(Ttransition-nendtrans+1:Ttransition)%govshock = equmINITSS%govshock

	IRFDir = "GovExp"
	CALL IRFSequence(IRFDir)
END IF

!Lump Transfer Shock
IF(IncludeTransferShock==1) THEN
	IF(Display>=1) write(*,*)'Solving for lump transfer shock IRF'	
	equmTRANS(:) = equmINITSS	
	
	equmTRANS(1)%transfershock = equmINITSS%transfershock * exp(TransferShockSize)
	DO it = 2,Ttransition-nendtrans
		equmTRANS(it)%transfershock = equmTRANS(it-1)%transfershock ** (TransferShockPers**deltatransvec(it-1))
	END DO	
	equmTRANS(Ttransition-nendtrans+1:Ttransition)%transfershock = equmINITSS%transfershock

	IRFDir = "Transfer"
	CALL IRFSequence(IRFDir)
END IF

!Financial wedge shock
IF(IncludeFinWedgeShock==1) THEN
	IF(Display>=1) write(*,*)'Solving for Financial wedge shock IRF'	

	equmTRANS(:) = equmINITSS		
	equmTRANS(1)%finwedge = equmINITSS%finwedge + FinWedgeShockSize
	DO it = 2,Ttransition-nendtrans
		equmTRANS(it)%finwedge = equmINITSS%finwedge *(1.0-FinWedgeShockPers**deltatransvec(it-1)) + equmTRANS(it-1)%finwedge * (FinWedgeShockPers**deltatransvec(it-1))
	END DO	
	equmTRANS(Ttransition-nendtrans+1:Ttransition)%finwedge = equmINITSS%finwedge

	IRFDir = "FinWedge"
	CALL IRFSequence(IRFDir)

END IF

!Labor wedge shock
IF(IncludeLabWedgeShock==1) THEN
	IF(Display>=1) write(*,*)'Solving for Labor wedge shock IRF'	
	equmTRANS(:) = equmINITSS	
	
	equmTRANS(1)%labwedge = equmINITSS%labwedge * exp(LabWedgeShockSize)
	DO it = 2,Ttransition-nendtrans
		equmTRANS(it)%labwedge = equmTRANS(it-1)%labwedge ** (LabWedgeShockPers**deltatransvec(it-1))
	END DO	
	equmTRANS(Ttransition-nendtrans+1:Ttransition)%labwedge = equmINITSS%labwedge

	IRFDir = "LabWedge"
	CALL IRFSequence(IRFDir)
END IF

!Productivity dispersion shock
IF(IncludeProdDispShock==1) THEN
	prodispshock = .true.
	IF(Display>=1) write(*,*)'Solving for producitivity dispersion shock IRF'	
	equmTRANS(:) = equmINITSS	
	
	equmTRANS(1)%prodgridscale = equmINITSS%prodgridscale * exp(ProdDispShockSize)
	DO it = 2,Ttransition-nendtrans
		equmTRANS(it)%prodgridscale = equmTRANS(it-1)%prodgridscale ** (ProdDispShockPers**deltatransvec(it-1))
	END DO	
	equmTRANS(Ttransition-nendtrans+1:Ttransition)%prodgridscale = equmINITSS%prodgridscale

	lminygrid = 1.0e-8_8
	equmTRANS(1)%yprodgrid(:) = equmTRANS(1)%prodgridscale * equmINITSS%yprodgrid(:) - meanlabeff*(equmTRANS(1)%prodgridscale-1)
	equmTRANS(1)%yprodgrid(:) = max(equmTRANS(1)%yprodgrid(:),lminygrid)
	DO it = 2,Ttransition
		equmTRANS(it)%yprodgrid(:) = equmTRANS(it)%prodgridscale * equmINITSS%yprodgrid(:) - meanlabeff*(equmTRANS(it)%prodgridscale-1)
		equmTRANS(it)%yprodgrid(:) = max(equmTRANS(it)%yprodgrid(:),lminygrid)
	END DO	


	IRFDir = "ProdDisp"
	CALL IRFSequence(IRFDir)
	prodispshock = .false.
END IF

!Productivity persistence shock
IF(IncludeProdPersShock==1) THEN
	IF(Display>=1) write(*,*)'Solving for producitivity persistence shock IRF'	
	equmTRANS(:) = equmINITSS	
	
	equmTRANS(1)%prodmarkovscale = equmINITSS%prodmarkovscale * exp(ProdPersShockSize)
	DO it = 2,Ttransition-nendtrans
		equmTRANS(it)%prodmarkovscale = equmTRANS(it-1)%prodmarkovscale ** (ProdPersShockPers**deltatransvec(it-1))
	END DO	
	equmTRANS(Ttransition-nendtrans+1:Ttransition)%prodmarkovscale = equmINITSS%prodmarkovscale

	IRFDir = "ProdPers"
	CALL IRFSequence(IRFDir)
END IF

END SUBROUTINE ImpulseResponses