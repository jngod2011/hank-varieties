SUBROUTINE MomentConditions(lnpar,ly,lf,lfvec)

USE Parameters
USE Globals
USE Procedures

IMPLICIT NONE

INTEGER, INTENT(IN)     					:: lnpar
REAL(8), DIMENSION(lnpar), INTENT(IN)  	:: ly   !length lnpar
REAL(8), INTENT(OUT) 					:: lf,lfvec(nmoments)
INTEGER     	:: ip,iflag
REAL(8)			:: lx(lnpar),loutput,lKNratio,lub
REAL(8), EXTERNAL	:: FnCapitalEquity

lx = ly/paramscale
! write(*,*) 'lx is: ',lx

write(4,*) '*********************************'
write(4,*) 'EVALUATION NUMBER: ',objeval

!Extract parameters
ip = 0
IF(EstimateKappafc==1) THEN
    ip=ip+1
	kappafc_w = exp(lx(ip))
	write(4,*) ' kappafc_w guess: ',kappafc_w
END IF
IF(EstimateKappa0==1) THEN
    ip=ip+1
	kappa0_w = logistic(lx(ip))
	write(4,*) ' kappa0_w guess: ',kappa0_w
END IF
IF(EstimateKappa1==1) THEN
    ip=ip+1
	kappa1_w = exp(lx(ip))
	write(4,*) ' kappa1_w guess: ',kappa1_w
END IF
IF(EstimateKappa2==1) THEN
    ip=ip+1
	kappa2_w = kappa2min + exp(lx(ip))
	write(4,*) ' kappa2_w guess: ',kappa2_w
END IF
IF(EstimateKappa3==1) THEN
    ip=ip+1
	kappa3 = exp(lx(ip))
	write(4,*) ' kappa3 guess: ',kappa3
END IF
IF(EstimateKappa4==1) THEN
    ip=ip+1
	kappa4_w = exp(lx(ip))
	write(4,*) ' kappa4_w guess: ',kappa4_w
END IF
IF(EstimateRho==1) THEN
    ip=ip+1
    rho = -log(logistic(lx(ip)))
	write(4,*) ' rho: ',rho
END IF
IF(EstimateBorrWedge==1) THEN
    ip=ip+1
	borrwedge = borrwedgemax*logistic(lx(ip))	
	write(4,*) ' borrwedge guess: ',borrwedge
END IF
IF(EstimateGamma==1) THEN
    ip=ip+1
	gam = exp(lx(ip))
	write(4,*) ' gam guess: ',gam
END IF

IF(PinKappa1ByKappa02==1)THEN
	kappa1_w = ((1.0-kappa0_w)*(1.0+kappa2_w))**(-1.0/kappa2_w)
END IF	

IF(SymmetricAdjustmentCost==1) THEN
	kappafc_d = kappafc_w
	kappa0_d = kappa0_w
	kappa1_d = kappa1_w
	kappa2_d = kappa2_w
	kappa4_d = kappa4_w
ELSE IF(NoDepositCost==1) THEN
	kappafc_d = 0.0
	kappa0_d = 0.0
	kappa1_d = 100.0	!this is the optimal depost rate, i.e. higher for lower cost
	kappa2_d = 1.0
	kappa4_d = 0.0
END IF


IF(ImposeEqumInCalibration==0) THEN
	CALL Grids
	CALL InitialPrices
	CALL IterateBellman
	CALL StationaryDistribution
ELSE IF(ImposeEqumInCalibration==1) THEN
	initialSS = .true.
	CALL SolveSteadyStateEqum
END IF

CALL DistributionStatistics

!implied capital-output ratio
la = -deprec 
lb = Ea*deprec + (1.0-1.0/elast)*alpha_Y*drs_Y + (1.0/elast)*alpha_N*drs_N + ((1.0-1.0/elast)*(1.0-drs_Y) + (1.0/elast)*(1.0-drs_N)) *(1.0-corptax)*profdistfracA
lc = - Ea* ( (1.0-1.0/elast)*alpha_Y*drs_Y + (1.0/elast)*alpha_N*drs_N )
lK_totoutput_ratio = (-lb+sqrt(lb**2-4*la*lc)) / (2*la)

!implied total labor
llabor_Y = Elabor_Y/varieties
llabor_N = Elabor_N

!recall normalize steady-state output to 1
PROBLEM HERE BECAUSE IT ONLY WORKS WITH RHO CALIBRATION BECAUSE WE UPDATE LABOR.

!model implied moments
IF(MatchRelativeToTargetOutput==0) THEN
	loutput = tfp*(capital**alpha)*(labor**(1.0-alpha)) 
	loutput = max(loutput, 0.01)
	modelMeanIll = Ea/loutput
	modelMeanLiq = Eb/loutput
	modelMedianIll = PERCa(6)/loutput
	modelP75Ill = PERCa(7)/loutput
	modelMedianLiq = PERCb(6)/loutput
ELSE !use for no labor supply
! 	loutput = KYratio/max(capital,1.0e-8)
	loutput = 1.0
	modelMeanIll = Ea/loutput
	modelMeanLiq = Eb/loutput
	modelMedianIll = PERCa(6)/loutput
	modelP75Ill = PERCa(7)/loutput
	modelMedianLiq = PERCb(6)/loutput
END IF
modelFracIll0 = FRACa0close 
modelFracLiq0 = FRACb0close
modelFracIll0Liq0 = FRACb0a0close 
modelFracLiqNeg = FRACbN

write(*,*) ra,Ea,Eb,labor,capital, equity, loutput

!moments to match: need to include weights
ip = 0
IF (MatchMeanIll == 1) THEN
	ip = ip+1
! 	lfvec(ip) = (modelMeanIll - targetMeanIll)
	lfvec(ip) = (modelMeanIll/targetMeanIll - 1.0)
	write(4,*) 'MeanIll, target: ',targetMeanIll, ' model: ',modelMeanIll
	
! 	IF(UseDiagWeight==1) diagweight = 1.0/(stdevElogh**2.0)
END IF
IF (MatchKYratio == 1 .or. ImposeEqumInCalibration==1) THEN
	ip = ip+1
! 	lfvec(ip) = (modelKYratio - targetKYratio)
	lfvec(ip) = (modelKYratio/targetKYratio - 1.0)
	write(4,*) 'KYratio, target: ',targetKYratio, ' model: ',modelKYratio
	
	!for now just upweight this by 10
	lfvec(ip) = lfvec(ip)*10.0
! 	IF(UseDiagWeight==1) diagweight = 1.0/(stdevElogh**2.0)
END IF
IF (MatchMedianIll == 1) THEN
	ip = ip+1
! 	lfvec(ip) = (modelMedianIll - targetMedianIll)
	lfvec(ip) = (modelMedianIll/targetMedianIll - 1.0)
	write(4,*) 'MedianIll, target: ',targetMedianIll, ' model: ',modelMedianIll
END IF
IF (MatchP75Ill == 1) THEN
	ip = ip+1
! 	lfvec(ip) = (modelMedianIll - targetMedianIll)
	lfvec(ip) = (modelP75Ill/targetP75Ill - 1.0)
	write(4,*) 'P75Ill, target: ',targetP75Ill, ' model: ',modelP75Ill
END IF
IF (MatchFracIll0 == 1) THEN
	ip = ip+1
! 	lfvec(ip) = (modelFracIll0 - targetFracIll0)
	lfvec(ip) = (modelFracIll0/targetFracIll0 - 1.0)
	write(4,*) 'FracIll0, target: ',targetFracIll0, ' model: ',modelFracIll0
END IF
IF (MatchMeanLiq == 1) THEN
	ip = ip+1
! 	lfvec(ip) = (modelMeanLiq - targetMeanLiq)
	lfvec(ip) = (modelMeanLiq/targetMeanLiq -1.0)
	write(4,*) 'MeanLiq, target: ',targetMeanLiq, ' model: ',modelMeanLiq
! 	IF(UseDiagWeight==1) diagweight = 1.0/
!upweight by 2                                                                                                                                                                                                                                                                 
       lfvec(ip) = lfvec(ip)*2.0
END IF
IF (MatchMedianLiq == 1) THEN
	ip = ip+1
! 	lfvec(ip) = (modelMedianLiq - targetMedianLiq)
	lfvec(ip) = (modelMedianLiq/targetMedianLiq - 1.0)
	write(4,*) 'MedianLiq, target: ',targetMedianLiq, ' model: ',modelMedianLiq
! 	IF(UseDiagWeight==1) diagweight = 1.0/(stdevElogh**2.0)
END IF
IF (MatchFracLiq0 == 1) THEN
	ip = ip+1
! 	lfvec(ip) = (modelFracLiq0 - targetFracLiq0)
	lfvec(ip) = (modelFracLiq0/targetFracLiq0 - 1.0)
	write(4,*) 'FracLiq0, target: ',targetFracLiq0, ' model: ',modelFracLiq0
! 	IF(UseDiagWeight==1) diagweight = 1.0/(stdevElogh**2.0)

	!for now just upweight this by 2
	lfvec(ip) = lfvec(ip)*2.0

END IF
IF (MatchFracLiqNeg == 1) THEN
	ip = ip+1
! 	lfvec(ip) = (modelFracLiqNeg - targetFracLiqNeg)
	lfvec(ip) = (modelFracLiqNeg/targetFracLiqNeg-1.0)
	write(4,*) 'FracLiqNeg, target: ',targetFracLiqNeg, ' model: ',modelFracLiqNeg
END IF
IF (MatchFracIll0Liq0 == 1) THEN
	ip = ip+1
	lfvec(ip) = (modelFracIll0Liq0/targetFracIll0Liq0 - 1.0)
	write(4,*) 'FracIll0Liq0, target: ',targetFracIll0Liq0, ' model: ',modelFracIll0Liq0
	
    lfvec(ip) = lfvec(ip)
	
END IF



!construct objective function manually as well:
lf = sum((lfvec**2.0)*diagweight)/real(nmoments)
write(4,*) 'objective fun: ',lf
write(4,*) ' '

objeval = objeval+1

END SUBROUTINE MomentConditions