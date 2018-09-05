SUBROUTINE Exploration

USE Parameters
USE Globals
USE Procedures

IMPLICIT NONE

INTEGER					:: ip,is,iflag,ismin(1)
REAL(8)					:: lrhoL,lrhoU
REAL(8)			:: lcapital,loutput,lKNratio
REAL(8), EXTERNAL 		:: FnDiscountRate

exploring = .true.

!set up parameters
nparam = 0
IF (EstimateKappafc == 1) 		nparam = nparam+1
IF (EstimateKappa0 == 1) 		nparam = nparam+1
IF (EstimateKappa1 == 1) 		nparam = nparam+1
IF (EstimateKappa2 == 1) 		nparam = nparam+1
IF (EstimateKappa3 == 1) 		nparam = nparam+1
IF (EstimateRho == 1 .and. CalibrateRhoInExploration==0) nparam = nparam+1
IF (EstimateBorrWedge == 1) 	nparam = nparam+1
IF (EstimateGamma == 1) 	nparam = nparam+1

IF (ALLOCATED(paramout)) DEALLOCATE(paramout)
IF (ALLOCATED(paramub)) DEALLOCATE(paramub)
IF (ALLOCATED(paramlb)) DEALLOCATE(paramlb)
ALLOCATE(paramout(nparam),paramub(nparam),paramlb(nparam))

!set up moments
nmoments = 0
IF (MatchMeanIll == 1)		nmoments = nmoments+1
IF (MatchKYratio == 1 .and. CalibrateRhoInExploration==0) nmoments = nmoments+1
IF (MatchMedianIll == 1)	nmoments = nmoments+1
IF (MatchP75Ill == 1)	nmoments = nmoments+1
IF (MatchFracIll0 == 1)		nmoments = nmoments+1
IF (MatchMeanLiq == 1)		nmoments = nmoments+1
IF (MatchMedianLiq == 1)	nmoments = nmoments+1
IF (MatchFracLiq0 == 1)		nmoments = nmoments+1
IF (MatchFracLiqNeg == 1)	nmoments = nmoments+1
IF (MatchFracIll0Liq0 == 1)	nmoments = nmoments+1

!get sobol sequence
nsobol = 200 !500
nsobolskip = floor(log(real(nsobol))/log(2.0)) !largest power of 2 smaller than nsobol
IF (ALLOCATED(sobolseq)) DEALLOCATE(sobolseq)
IF (ALLOCATED(paramsobol)) DEALLOCATE(paramsobol)
IF (ALLOCATED(sobolmom)) DEALLOCATE(sobolmom)
IF (ALLOCATED(sobolobj)) DEALLOCATE(sobolobj)

ALLOCATE(sobolseq(nparam,nsobol),paramsobol(nparam,nsobol),sobolmom(nmoments,nsobol),sobolobj(nsobol))
CALL i8_sobol_generate (nparam,nsobol,nsobolskip,sobolseq )

!set up bounds
ip = 0
IF(EstimateKappafc==1) THEN
    ip=ip+1
	paramlb(ip) = 0.0
	paramub(ip) = 0.5
END IF
IF(EstimateKappa0==1) THEN
    ip=ip+1
	paramlb(ip) = 0.0
	paramub(ip) = 0.95
END IF
IF(EstimateKappa1==1) THEN
    ip=ip+1
    paramlb(ip) = 0.0
    paramub(ip) = 3.0
END IF
IF(EstimateKappa2==1) THEN
    ip=ip+1
    paramlb(ip) = kappa2min
	paramub(ip) = 15.0
END IF
IF(EstimateKappa3==1) THEN
    ip=ip+1
	paramlb(ip) = log(0.01)
    paramub(ip) = log(10.0)
END IF
IF(EstimateKappa4==1) THEN
    ip=ip+1
	paramlb(ip) = 1.0e-5
    paramub(ip) = 0.5
END IF
IF(EstimateRho==1 .and. CalibrateRhoInExploration==0) THEN
    ip=ip+1
    paramlb(ip) = 0.001
    paramub(ip) = 0.1
END IF
IF(EstimateBorrWedge==1) THEN
    ip=ip+1
    paramlb(ip) = ra-rb
    paramub(ip) = borrwedgemax
END IF
IF(EstimateGamma==1) THEN
    ip=ip+1
	paramlb(ip) = 1.0
    paramub(ip) = 3.0
END IF

!create sobol paramter points
!$OMP PARALLEL DO PRIVATE(ip)
DO is = 1,nsobol
	DO ip = 1,nparam
		paramsobol(ip,is) = paramlb(ip) + sobolseq(ip,is)*(paramub(ip)-paramlb(ip))
	END DO
END DO
!$OMP END PARALLEL DO
OPEN(3, FILE = trim(OutputDir) // 'sobolseq.txt', STATUS = 'replace'); CALL WriteMatrixLong(3,nparam,nsobol,sobolseq)
OPEN(3, FILE = trim(OutputDir) // 'paramsobol.txt', STATUS = 'replace'); CALL WriteMatrixLong(3,nparam,nsobol,paramsobol)


!evaluate model at each point in sobol sequence
OPEN(4, FILE = trim(OutputDir) // 'sobol_iterations.txt', STATUS = 'replace')
DO is = 1,nsobol	
	
	IF(Display>=1) WRITE(*,*)"Exploring parameter space, sobol sequence point ",is
	
	write(4,*) '*********************************'
	write(4,*) 'EVALUATION NUMBER: ',is

	!extract parameters
	ip = 0
	IF(EstimateKappafc==1) THEN
	    ip=ip+1
		kappafc_w = paramsobol(ip,is)
		write(4,*) ' kappafc_w guess: ',kappafc_w		
	END IF
	IF(EstimateKappa0==1) THEN
	    ip=ip+1
		kappa0_w = paramsobol(ip,is)
		write(4,*) ' kappa0_w guess: ',kappa0_w
	END IF
	IF(EstimateKappa1==1) THEN
	    ip=ip+1
		kappa1_w = paramsobol(ip,is)
		write(4,*) ' kappa1_w guess: ',kappa1_w
	END IF
	IF(EstimateKappa2==1) THEN
	    ip=ip+1
		kappa2_w = paramsobol(ip,is)
		write(4,*) ' kappa2_w guess: ',kappa2_w
	END IF
	IF(EstimateKappa3==1) THEN
	    ip=ip+1
		kappa3 = exp(paramsobol(ip,is))
		write(4,*) ' kappa3 guess: ',kappa3
	END IF
	IF(EstimateKappa4==1) THEN
	    ip=ip+1
		kappa4_w = paramsobol(ip,is)
		write(4,*) ' kappa4_w guess: ',kappa4_w
	END IF
	IF(EstimateRho==1 .and. CalibrateRhoInExploration==0) THEN
	    ip=ip+1
		rho = paramsobol(ip,is)
		write(4,*) ' rho guess: ',rho
	END IF
	IF(EstimateBorrWedge==1) THEN
	    ip=ip+1
		borrwedge = paramsobol(ip,is)
		write(4,*) ' borrwedge guess: ',borrwedge	
	END IF
	IF(EstimateGamma==1) THEN
	    ip=ip+1
		gam = paramsobol(ip,is)
		write(4,*) ' gam guess: ',gam
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

	KYratio = targetKYratio
	KNratio = (tfp*KYratio)**(1.0/(1.0-alpha))
	rcapital = mc*alpha/KYratio
	wage = mc*(1.0-alpha)*tfp*(KNratio**alpha)
	netwage = (1.0-labtax)*wage
	rborr = rb + borrwedge
	labor = (netwage/chi)**frisch
	capital = KNratio*labor
	profit = (1.0-mc-operatecost)*capital/KYratio - priceadjust
	dividend = profit*(1.0-corptax)
	IF(DividendFundLumpSum==0) divrate = dividend/capital
	IF(DividendFundLumpSum==1) divrate = 0.0
	ra = (rcapital - deprec + divrate - fundlev*rb)/(1.0-fundlev)
	IF(DividendFundLumpSum==0) equity = 0.0
	IF(DividendFundLumpSum==1)equity = profit/ra
	output = tfp*(capital**alpha)*(labor**(1.0-alpha)) + rhousing*housefrac*((1.0-fundlev)*capital+equity)/(1.0-housefrac)
	fundbond = -capital*fundlev
	govbond = -ssdebttogdp*output
	priceadjust = 0.0
	intfirmbond = 0.0

	directdepmax = directdepmaxfrac*output
	directdepmin = directdepminfrac*output

	IF (CalibrateRhoInExploration==0) THEN
		CALL Grids
		CALL IterateBellman
		CALL StationaryDistribution	
		CALL DistributionStatistics
	ELSE IF (CalibrateRhoInExploration==1) THEN
		converged = .false.
		neqmiter = 1
		OPEN(3, FILE = trim(OutputDir) // 'DiscountRateCalibration.txt', STATUS = 'replace'); CLOSE(3)
		lrhoL = invlogistic(exp(-0.03_8))
		lrhoU = invlogistic(exp(-0.012_8))
	 	CALL rtflsp(FnDiscountRate,lrhoL,lrhoU,1.0e-8_8,tolrho,iflag,maxiterrho)
		converged = .true.
		write(4,*) ' calibrated rho: ',rho
	END IF

	!model implied moments
	IF(DividendFundLumpSum==0) lcapital = (1.0-housefrac)*Ea/(1.0-fundlev)
	IF(DividendFundLumpSum==1) lcapital = (1.0-housefrac)*Ea/(1.0-fundlev + (1.0-mc-operatecost)/(ra*KYratio))
	
	loutput = tfp*(lcapital**alpha)*(labor**(1.0-alpha)) + rhousing*housefrac*((1.0-fundlev)*lcapital+equity)/(1.0-housefrac)
		
	lKNratio = lcapital/labor
	modelKYratio = (lKNratio**(1.0-alpha)) / tfp	
	modelMeanIll = Ea/loutput
	modelMedianIll = PERCa(6)/loutput
	modelP75Ill = PERCa(7)/loutput
	modelMeanLiq = Eb/loutput
	modelMedianLiq = PERCb(6)/loutput
	modelFracIll0 = FRACa0close 
	modelFracLiq0 = FRACb0close
	modelFracIll0Liq0 = FRACb0a0close 
	modelFracLiqNeg = FRACbN

	!moments to match: need to include weights
	ip = 0
	IF (MatchMeanIll == 1) THEN
		ip = ip+1
		sobolmom(ip,is) = (modelMeanIll/targetMeanIll - 1.0)
		write(4,*) 'MeanIll, target: ',targetMeanIll, ' model: ',modelMeanIll
	END IF
	IF (MatchKYratio == 1 .and. CalibrateRhoInExploration==0) THEN
		ip = ip+1
		sobolmom(ip,is) = (modelKYratio/targetKYratio - 1.0)
		write(4,*) 'KYratio, target: ',targetKYratio, ' model: ',modelKYratio
	END IF
	IF (MatchMedianIll == 1) THEN
		ip = ip+1
		sobolmom(ip,is) = (modelMedianIll/targetMedianIll - 1.0)
		write(4,*) 'MedianIll, target: ',targetMedianIll, ' model: ',modelMedianIll
	END IF
	IF (MatchP75Ill == 1) THEN
		ip = ip+1
		sobolmom(ip,is) = (modelP75Ill/targetP75Ill - 1.0)
		write(4,*) 'P75Ill, target: ',targetP75Ill, ' model: ',modelP75Ill
	END IF
	IF (MatchFracIll0 == 1) THEN
		ip = ip+1
		sobolmom(ip,is) = (modelFracIll0/targetFracIll0 - 1.0)
		write(4,*) 'FracIll0, target: ',targetFracIll0, ' model: ',modelFracIll0
	END IF
	IF (MatchMeanLiq == 1) THEN
		ip = ip+1
		sobolmom(ip,is) = (modelMeanLiq/targetMeanLiq -1.0)
		write(4,*) 'MeanLiq, target: ',targetMeanLiq, ' model: ',modelMeanLiq
	END IF
	IF (MatchMedianLiq == 1) THEN
		ip = ip+1
		sobolmom(ip,is) = (modelMedianLiq/targetMedianLiq - 1.0)
		write(4,*) 'MedianLiq, target: ',targetMedianLiq, ' model: ',modelMedianLiq
	END IF
	IF (MatchFracLiq0 == 1) THEN
		ip = ip+1
		sobolmom(ip,is) = (modelFracLiq0/targetFracLiq0 - 1.0)
		write(4,*) 'FracLiq0, target: ',targetFracLiq0, ' model: ',modelFracLiq0
	END IF
	IF (MatchFracLiqNeg == 1) THEN
		ip = ip+1
		sobolmom(ip,is) = (modelFracLiqNeg/targetFracLiqNeg-1.0)
		write(4,*) 'FracLiqNeg, target: ',targetFracLiqNeg, ' model: ',modelFracLiqNeg
	END IF
	IF (MatchFracIll0Liq0 == 1) THEN
		ip = ip+1
		sobolmom(ip,is) = (modelFracIll0Liq0/targetFracIll0Liq0 - 1.0)
		write(4,*) 'FracFracIll0Liq0, target: ',targetFracIll0Liq0, ' model: ',modelFracIll0Liq0
	END IF


	IF (CalibrateRhoInExploration==1) write(4,*) 'KYratio, target: ',targetKYratio, ' model: ',modelKYratio

	!construct objective function manually as well:
	sobolobj(is) = sum((sobolmom(:,is)**2.0))/real(nmoments)
	write(4,*) 'objective fun (unweighted): ',sobolobj(is)
	write(4,*) ' '

END DO


OPEN(6, FILE = trim(OutputDir) // 'sobolmom.txt', STATUS = 'replace'); CALL WriteMatrixLong(6,nmoments,nsobol,sobolmom)
OPEN(6, FILE = trim(OutputDir) // 'sobolobj.txt', STATUS = 'replace'); CALL WriteMatrixLong(6,nsobol,1,sobolobj)


!extract parameters at minimum objective
write(4,*) '*********************************'
ismin = MINLOC(sobolobj)
is =ismin(1)
write(4,*) 'MINIMUM EVALUATION NUMBER: ',is
write(4,*) 'MINIMUM OBJECTIVE: ',MINVAL(sobolobj)


ip = 0
IF(EstimateKappafc==1) THEN
    ip=ip+1
	kappafc_w = paramsobol(ip,is)
	write(4,*) ' kappafc_w guess: ',kappafc_w		
END IF
IF(EstimateKappa0==1) THEN
    ip=ip+1
	kappa0_w = paramsobol(ip,is)
	write(4,*) ' kappa0_w guess: ',kappa0_w
END IF
IF(EstimateKappa1==1) THEN
    ip=ip+1
	kappa1_w = paramsobol(ip,is)
	write(4,*) ' kappa1_w guess: ',kappa1_w
END IF
IF(EstimateKappa2==1) THEN
    ip=ip+1
	kappa2_w = paramsobol(ip,is)
	write(4,*) ' kappa2_w guess: ',kappa2_w
END IF
IF(EstimateKappa3==1) THEN
    ip=ip+1
	kappa3 = exp(paramsobol(ip,is))
	write(4,*) ' kappa3 guess: ',kappa3
END IF
IF(EstimateKappa4==1) THEN
    ip=ip+1
	kappa4_w = paramsobol(ip,is)
	write(4,*) ' kappa4_w guess: ',kappa4_w
END IF
IF(EstimateRho==1 .and. CalibrateRhoInExploration==0) THEN
    ip=ip+1
	rho = paramsobol(ip,is)
	write(4,*) ' rho guess: ',rho
END IF
IF(EstimateBorrWedge==1) THEN
    ip=ip+1
	borrwedge = paramsobol(ip,is)
	write(4,*) ' borrwedge guess: ',borrwedge	
END IF
IF(EstimateGamma==1) THEN
    ip=ip+1
	gam = paramsobol(ip,is)
	write(4,*) ' gam guess: ',gam
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

ip = 0
IF (MatchMeanIll == 1) THEN
	ip = ip+1
	write(4,*) 'MeanIll, target: ',targetMeanIll, ' model: ',targetMeanIll*(1.0+sobolmom(ip,is))
END IF
IF (MatchKYratio == 1 .and. CalibrateRhoInExploration==0) THEN
	ip = ip+1
	write(4,*) 'KYratio, target: ',targetKYratio, ' model: ',targetKYratio*(1.0+sobolmom(ip,is))
END IF
IF (MatchMedianIll == 1) THEN
	ip = ip+1
	write(4,*) 'MedianIll, target: ',targetMedianIll, ' model: ',targetMedianIll*(1.0+sobolmom(ip,is))
END IF
IF (MatchP75Ill == 1) THEN
	ip = ip+1
	write(4,*) 'P75Ill, target: ',targetP75Ill, ' model: ',targetP75Ill*(1.0+sobolmom(ip,is))
END IF
IF (MatchFracIll0 == 1) THEN
	ip = ip+1
	write(4,*) 'FracIll0, target: ',targetFracIll0, ' model: ',targetFracIll0*(1.0+sobolmom(ip,is))
END IF
IF (MatchMeanLiq == 1) THEN
	ip = ip+1
	write(4,*) 'MeanLiq, target: ',targetMeanLiq, ' model: ',targetMeanLiq*(1.0+sobolmom(ip,is))
END IF
IF (MatchMedianLiq == 1) THEN
	ip = ip+1
	write(4,*) 'MedianLiq, target: ',targetMedianLiq, ' model: ',targetMedianLiq*(1.0+sobolmom(ip,is))
END IF
IF (MatchFracLiq0 == 1) THEN
	ip = ip+1
	write(4,*) 'FracLiq0, target: ',targetFracLiq0, ' model: ',targetFracLiq0*(1.0+sobolmom(ip,is))
END IF
IF (MatchFracLiqNeg == 1) THEN
	ip = ip+1
	write(4,*) 'FracLiqNeg, target: ',targetFracLiqNeg, ' model: ',targetFracLiqNeg*(1.0+sobolmom(ip,is))
END IF
IF (MatchFracIll0Liq0 == 1) THEN
	ip = ip+1
	write(4,*) 'FracIll0Liq0, target: ',targetFracIll0Liq0, ' model: ',targetFracIll0Liq0*(1.0+sobolmom(ip,is))
END IF

CLOSE(4)

exploring = .false.

KYratio = targetKYratio
KNratio = (tfp*KYratio)**(1.0/(1.0-alpha))
rcapital = mc*alpha/KYratio
wage = mc*(1.0-alpha)*tfp*(KNratio**alpha)
netwage = (1.0-labtax)*wage
labor = (netwage/chi)**frisch
capital = KNratio*labor
profit = (1.0-mc-operatecost)*capital/KYratio - priceadjust
dividend = profit*(1.0-corptax)
IF(DividendFundLumpSum==0) divrate = dividend/capital
IF(DividendFundLumpSum==1) divrate = 0.0
ra = (rcapital - deprec + divrate - fundlev*rb)/(1.0-fundlev)
IF(DividendFundLumpSum==0) equity = 0.0
IF(DividendFundLumpSum==1)equity = profit/ra
output = tfp*(capital**alpha)*(labor**(1.0-alpha)) + rhousing*housefrac*((1.0-fundlev)*capital+equity)/(1.0-housefrac)
rborr = rb + borrwedge
fundbond = -capital*fundlev
taxrev = labtax*wage*labor - lumptransfer + corptax*profit
govbond = -ssdebttogdp*output
govexp = taxrev + rb*govbond
directdepmax = directdepmaxfrac*output
directdepmin = directdepminfrac*output
priceadjust = 0.0

END SUBROUTINE Exploration

