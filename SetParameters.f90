SUBROUTINE SetParameters

USE Parameters
USE Globals
USE Procedures

IMPLICIT NONE

INTEGER 	:: iy
REAL 		:: la,lb,lc,lelastoc

!OUTPUT DIR
OutputDir =			"/Volumes/FILES/Large/ContinuousTimeAdjustment/check_christian3/"
!OutputDir =			"Output/"
EarningsProcessDir	= "/Volumes/FILES/Git/hank-main/earnings_input/2point_3_11"
! EarningsProcessDir	= "/Volumes/FILES/Git/hank-main/earnings_input/2point_3_3"
!EarningsProcessDir	= "/Volumes/FILES/Projects/MollContinuousTime/Fortran/DiscretizeEarnings2/2point_3_11_test"
!EarningsProcessDir	= "/home/gwk210/hank-main/earnings_input/2point_3_11"

CALL system ("mkdir -p " // trim(OutputDir))

!OPTIONS

!display options
Display              		= 2
ReportNonMonotonicity   	= 0

!run options
CalibrateDiscountRate	= 0 !at most 1 of EquilibriumR and CalibrateDiscountRate
EquilibriumR		 	= 0
ComputeCumulativeMPC 	= 1
DoImpulseResponses 		= 0!1
DoFiscalStimulus 		= 0
DoPriceExperiments		= 1
SaveTransitionPolicyFns = 0	!warning: will create a very large number of text files 
SaveTime1PolicyFns 		= 1

!labor risk options
ReadEarningsProcess 		= 1
LaborSupply					= 0
IncludeUnemploymentState	= 0
ScaleGHHIdiosyncratic		= 1
PerfectAnnuityMarkets		= 1
GovBondResidualZeroWorld    = 1	!imposes closed economy and imputes residual bond holdings to govt

!calibration options
CalibrateCostFunction		= 0
CalibrateRhoAtInitialGuess  = 1
SobolSequenceExploration	= 0
CalibrateRhoInExploration	= 0
CalibrateCloseToSobol 		= 0
MatchRelativeToTargetOutput = 1
ImposeEqumInCalibration 	= 0

!adjustment cost options
SymmetricAdjustmentCost 	= 1 !at most 1 of SymmetricAdjustmentCost and NoDepositCost
NoDepositCost 			= 0
AdjustCostFnRatio		= 1
Kappa3LinearOrMax 		= 2	!1 for linear, 2 for max (only valid if AdjustCostFnRatio==1)
ExponAdjustConstFn 		= 0
FixedValueForD 			= 0 !overides other adjustment cost options, must set dexog below

!transition computation options
PermanentShock 			= 0
SolveFlexPriceTransition	= 0
SolveStickyPriceTransition	= 1
SolveZLBTransition		= 0
StickyPriceAlgorithm		= 4 !1 for iterate on B, 2 for iterate on Pi, 3 for iterate on Y, 4 for iterate on Rb
ZLBAlgorithm			= 4 !1 for iterate on B, 2 for iterate on Pi, 3 for iterate on Y, 4 for iterate on Rb
UpdateUsingR			= 1 !only relevant if Algorithm==1
ConvergenceRelToOutput 		= 1
UpdateRbFromMarketClearing      = 0
UseFlexTransitionAsGuess	= 0
GuessZeroInflation 		= 1
ImposeGuessSettlesDown 		= 0
AssumePathRealRate 		= 0 !not active yet
FirmDiscountRate		= 5	!1 for rho, 2 for rb initial steady state, 3 for ra initial steady state, 4 for rb transition, 5 for ra transition
BackwardTermInTaylorRule 	= 0
NoChangeLabDisutility 		= 0 !hold labor disutility constant at steady state level when solving HJB (does not affect labor supply decision)
bondelastrelgdp 		= 1.0 !1.0 !bigger for smaller interest rate movements, closer to zero for larger interest rate movements. relative to steady state gdp
bondadjust 			= 0.1 !more responsive interest rate when closer to zero
wageflex			= 1.0 !1.0 for full flexibility, 0.0 for fixed real wage

!dividend options
FixProfitsOutOfSteadyState 	= 0
DividendSmoothing		= 0  !1 for smoothing problem, 2 for exponential decay, 3 for exponential decay with interest, 4 for alternate interest rates
RetainedEarningsInBondMarket    = 0
RetainedEarningsAsCapital	= 0
TaxRetainedEarnings 		= 0
firmgamma 			= 2.0	!only relevent if DividendSmoothing==1
divsmooth 			= 0.01 	!closer to zero for more closer to constant permanent dividends (smoother)
DividendFundLumpSum 		= 0 !1

!government bc options
AdjGovBudgetConstraint 		= 2 !1 for adjust spending, 2 for adjust lump sum taxes, 3 for let debt adjust (choose options below for financing)
GovExpConstantFracOutput 	= 0 !only active if AdjGovBudgetConstraint==3
taxincrstart 		= 1 !1 !quarters after shock that fiscal policy adjusts
taxincrdecay 		= 0.02 !0.1 !decay rate for tax increase higher for faster decay


!SHOCKS
IncludeTFPShock		= 0
IncludeBorrWedgeShock	= 0
IncludeNewsShock	= 0
IncludeKappafcShock	= 0
IncludeMonetaryShock	= 1
IncludeForwardGuideShock= 0!1
IncludeStDevYShock	= 0
IncludePrefShock	= 0!1
IncludeFundLevShock	= 0!1
IncludeRiskAversionShock= 0
IncludeMarkupShock	= 0
IncludeGovExpShock	= 0
IncludeTransferShock	= 0

TFPShockSize 		= -0.05 !log change
TFPShockPers 		= 0.7 !quarterly

BorrWedgeShockSize 		= 0.01 !percentage points
BorrWedgeShockPers 		= 0.5 !quarterly

NewsShockSize 			= -0.05 !log change
NewsShockPers 			= 0.5 !quarterly
NewsShockQtrs 			= 4 !number of quarters in advance

KappafcShockSize 		= 0.5 !log change
KappafcShockPers 		= 0.5 !quarterly

MonetaryShockSize 		= -0.0025 !percentage points
MonetaryShockPers 		= 0.5 !quarterly

ForwardGuideShockSize 	= -0.0025 !percentage points
ForwardGuideShockPers 	= 0.5 !quarterly

StDevYShockSize 		= 0.5 !log change
StDevYShockPers 		= 0.5 !quarterly

PrefShockSize 			= -0.6 !log change
PrefShockPers 			= 0.5 !quarterly

FundLevShockSize 		= -0.01 !percentage points
FundLevShockPers 		= 0.7 !quarterly

RiskAversionShockSize	= 0.1  !actual units
RiskAversionShockPers	= 0.5 !quarterly

MarkupShockSize 		= 0.05 !percentage points change in markup (converted later into elastici)
MarkupShockPers 		= 0.5 !quarterly

!FORWARD GUIDANCE OPTIONS
ForwardGuideShockQtrs 		= 9 !number of quarters in advance
ForwardGuideFixNomPreShock 	= 1 !if 1, then nominal rate is fixed before shock, and follows taylor rule after


!FISAL STIMULUS


!CALIBRATION OPTIONS
EstimateKappafc		= 0
EstimateKappa0		= 1
EstimateKappa1		= 0
EstimateKappa2		= 1
EstimateKappa3		= 0
EstimateKappa4		= 0
EstimateRho			= 1
EstimateBorrWedge	= 1
EstimateGamma		= 0

PinKappa1ByKappa02 	= 1

MatchMeanIll		= 0
MatchKYratio		= 1
MatchMedianIll		= 0
MatchP75Ill		= 0
MatchFracIll0		= 0
MatchMeanLiq		= 1
MatchMedianLiq		= 0
MatchFracLiq0		= 1
MatchFracLiqNeg		= 1
MatchFracIll0Liq0 	= 0

defnbclose			= 0.0 !0.0289 !0.01		!fraction of quarterly gross labor income (69100/4)
defnaclose 			= 0.0 !0.0289 !0.01

ndfls		= 1 

!SOLUTION PARAMETERS
maxiter 		= 500
Vtol			= 1.0e-8
maxiterKFE		= 2000
KFEtol			= 1.0e-8 !1.0e-7
deltass  		= 1.0e6
deltakfe 		= 1.0e6 !1.0 !1.0e4 
dVamin 			= 0.0
dVbmin 			= 1.0e-8

tolequmss		= 1.0e-8
stepequmss		= 0.1 !0.05
maxiterequmss	= 30 !20
maxiterrho 		= 30 !50
tolrho			= 1.0e-6

toltransition	= 1.0e-6
deltatransmin	= 1.0 !1.0e6 !1.0 !0.25
deltatransmax	= 40.0 !1.0e6 !40.0
deltatransparam	= 0.35 !0.3
maxitertransflex	= 2 !20 !50 !100 !300
maxitertranssticky	= 4000 !3000 !2000
stepflextransK  = 0.1 !0.05
stepflextransB  = 0.1 !0.051
stepstickytransK  = 0.01
stepstickytransB  = 0.001
stepstickytransRb  = 0.005
stepstickytransY  = 0.05

rho		=  0.019 !0.0145268586 !0.011860851 !0.011857549(fund lev) !0.011860851(base)
deathrate	= 0.0 !1.0/(4.0*45.0) !poisson death rate
gam		= 1.0	!risk aversion

flowgam		= 1.0 !10.0 !1.5 !risk aversion for flow consumption from illiquid assets
flowcon		= 0.05 !same order as ra (marginal utility of housing at a=0)

nondurwgt   = 1.0 !0.85
nondurwgttilde = (nondurwgt**nondurwgt) * ((1.0-nondurwgt)**(1.0-nondurwgt))

IF(LaborSupply==1) THEN
	frisch 		= 0.0 !0.5 	!frisch elasticity labor supply
	chi		= (1.0/3.0)**(-(1.0+frisch)/frisch) !sets labor income = 1 when hours = 1/3
ELSE IF(LaborSupply==0) THEN
	frisch = 0.0
	chi = 1.0
END IF

rb		= 0.02/4.0 !liquid return
prefshock	= 1.0

borrwedge 	 = 0.019875 !0.01611345 !quarterly wedge between rb and rborr: intermediation cost  
borrwedgemax     = 0.10 !0.125
blim 		 = -1.0 	!borrowing limit multiple of quarterly gross labor income

rborr = rb + borrwedge

directdepfrac 		= 0.0 !fraction of labor income deposited into illquid account
directdepmaxfrac 	= 5.0 	!as fraction of target output
directdepminfrac 	= 0.0  !as fraction of target output

!withdrawal fixed costs
kappafc_w	= 0.0 !10000.0 * 2.407e-5
kappa0_w	= 0.070046 !0.075355
kappa2_w	= 0.971050 !0.7359
kappa3		= 10000.0 * 2.407e-5
kappa4_w 	= 0.0

IF(PinKappa1ByKappa02==0) kappa1_w	= 0.535791 !0.525717
IF(PinKappa1ByKappa02==1) kappa1_w	= ((1.0-kappa0_w)*(1.0+kappa2_w))**(-1.0/kappa2_w)

dexog 	= 1.0	!only relevant if FixedValueForD==1
kappa2min   = 0.6 !0.5	!to make sure there is enough curvature for calibration

!deposit fixed costs
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
ELSE
	kappafc_d = 0.2
	kappa0_d = 0.0
	kappa1_d = 0.6
	kappa2_d = 1.0
	kappa4_d = 0.0
END IF	
	

theta 		= 100.0 !price adjustment cost
phitaylor 	= 1.25  ! inflation coefficient in taylor rule, 0.0 for fixed nominal interest rate at ss.
taylorpers 	= 0.1	!strong persistence is 0, low persistence is 1. only active if BackwardTermInTaylorRule=1
pi 		= 0.0		!inflation rate: this is steady state target	
rnom = rb + pi		!nominal interest rate (fisher equation): this will be constant in taylor rule
mpshock 	= 0.0		!shock to nominal interest rate in taylor rule
IF(BackwardTermInTaylorRule==1) phitaylor = phitaylor/taylorpers

!currently not in play:
yrho = 0.97		!ar1 parameter 
ylogvar = 0.4 !0.75 !cross-sectional variation of earnings
ysig =  sqrt(ylogvar*(1.0-yrho**2)) !st dev innovations 
ylambda = 1.0 !0.25 !1.0 !prob of a quarterly shock
jobloss = 0.04
jobfind = 0.666
uiben = 0.25

tfp 	= 1.0
elast 	= 1.0e8 !10.0 !elasticity of DS aggregator
gap 	= 0.0 !steady state output gap
mc = (elast-1.0)/elast

alpha 		= 0.4 !0.33 
alphatilde 	= (alpha**alpha) * ((1.0-alpha)**(1.0-alpha))
deprec 		= 0.08/4.0 !0.10/4.0 !0.0868/4.0 !0.10/4.0	!depreciation rate
operatecost = 0.0

housefrac = 0.0 !0.27  !fraction of illiquid assets that is non-productive

rhousing 	= 0.145 / 4.0
housemaint 	= 0.04 / 4.0
housedeprec = 0.09 / 4.0
flowamnt = rhousing - housedeprec - housemaint !quarterly flow consumption as fraction of non-productive illiquid asset

fundlev 	= 0.0 !0.11 !0.0
investadjcost 	= 0.0
capitaladjcost  = 0.0
utilelast 	= 5.0 !2.0
utilelastalpha  = 1.0 + utilelast-alpha*utilelast

!government
labtax 		= 0.25
corptax 	= 0.0 !0.15
ssdebttogdp 	= 0.26*4
lumptransferpc	= 0.4 !fraction of population that recevies more transfers than pays taxes

!calibration targets (relative to output, except KYratio)
targetMeanIll 		= 2.92 * 4.0
targetMeanLiq  		= 0.26 * 4.0
targetMedianIll 	= 0.21 * targetMeanIll
targetP75Ill 		= 0.71 * targetMeanIll
targetMedianLiq 	= 0.085 * targetMeanLiq
targetFracIll0 		= 0.115 + 0.12 !0.115
targetFracLiq0 		= 0.30 !0.36 !0.35
targetFracLiqNEG	= 0.15 !0.10
targetFracIll0Liq0 	= 0.10 !0.12 !0.115

IF (DividendFundLumpSum ==0) THEN
	targetKYtotalratio 	= targetMeanIll*(1.0-housefrac)/(1.0 - fundlev)
! 	targetKYratio 		= targetKYtotalratio / (1- flowamnt*(1.0 - fundlev)*housefrac*targetKYtotalratio/(1.0-housefrac))
	targetKYratio 		= targetKYtotalratio / (1- rhousing*(1.0 - fundlev)*housefrac*targetKYtotalratio/(1.0-housefrac))

ELSE IF (DividendFundLumpSum ==1) THEN
	lelastoc 	= elast/(1.0-operatecost*elast)
	
	la = -(deprec+rb*fundlev) * (1.0-rhousing*housefrac*targetMeanIll)
	lb = ((elast-1.0)/elast)*alpha*(1.0-rhousing*housefrac*targetMeanIll) &
	        + ((1.0-housefrac)/(1.0-fundlev))*targetMeanIll*(deprec+rb*fundlev) &
	        + (1.0/lelastoc)*(1.0-rhousing*housefrac*targetMeanIll)
	lc = -((1-housefrac)/(1-fundlev)) * targetMeanIll *((elast-1)/elast) *alpha
	
	targetKYratio = (-lb+sqrt(lb**2-4*la*lc)) / (2*la)
	
END IF

!if solving for equilibrium, these are guesses
KYratio = targetKYratio  !ratio to productive output
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
fundbond = -capital*fundlev
govbond = -ssdebttogdp*output
bondelast = bondelastrelgdp*output
priceadjust = 0.0

intfirmbond = 0.0
caputil 		= 1.0
tfpadj = ((tfp**(1.0+utilelast)) * (mc*alpha/rcapital)**(alpha*utilelast))**(1.0/utilelastalpha)

directdepmax = directdepmaxfrac*output
directdepmin = directdepminfrac*output

!allocate large arrays
ALLOCATE(V(ngpa,ngpb,ngpy),Vnew(ngpa,ngpb,ngpy))
ALLOCATE(u(ngpa,ngpb,ngpy),c(ngpa,ngpb,ngpy),p(ngpa,ngpb,ngpy),d(ngpa,ngpb,ngpy),H(ngpa,ngpb,ngpy),s(ngpa,ngpb,ngpy),bdot(ngpa,ngpb,ngpy),fspamount(ngpa,ngpb,ngpy))
ALLOCATE(ccum1(ngpa,ngpb,ngpy),ccum4(ngpa,ngpb,ngpy),ccum2(ngpa,ngpb,ngpy))
ALLOCATE(gjoint(ngpa,ngpb,ngpy),gamarg(ngpa,ngpy),gbmarg(ngpb,ngpy),gvec(naby),gmat(nab,ngpy),gabmarg(ngpa,ngpb),gabcum(ngpa,ngpb))

!allocate solution types
CALL AllocateSolutionType(solnINITSS)
CALL AllocateSolutionType(solnFINALSS)

!allocate COO matrices
DO iy = 1,ngpy
	ALLOCATE(ACOO(iy)%val(5*nab),ACOO(iy)%row(5*nab),ACOO(iy)%col(5*nab))
	ALLOCATE(AUCOO(iy)%val(3*nab),AUCOO(iy)%row(3*nab),AUCOO(iy)%col(3*nab))
END DO

exploring = .false.
calibrating = .false.
iteratingtransition = .false.
testingGHH = .false.
forwardguide = .false.

END SUBROUTINE SetParameters
