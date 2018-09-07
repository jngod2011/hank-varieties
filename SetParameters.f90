SUBROUTINE SetParameters

USE Parameters
USE Globals
USE Procedures

IMPLICIT NONE

INTEGER 	:: iy

!OUTPUT DIR
OutputDir =			"/Volumes/FILES/Large/HANKVarieties/test/"
! OutputDir ="Output/"
EarningsProcessDir	= "/Volumes/FILES/Git/hank-main/earnings_input/2point_3_11"
! EarningsProcessDir	= "/home/gkaplan/git/hank-main/earnings_input/2point_3_11"

!INPUT / OUTPUT DIR FOR FED IN PRICES
FeedPriceFile = "/Volumes/FILES/Dropbox/Projects/JEPMacroHeterogeneity/HANK_shocks/FeedPriceDirectories.txt" !file with paths on separate lines
! FeedPriceFile = "/Volumes/FILES/Dropbox/Projects/JEPMacroHeterogeneity/HANK_shocks/FeedPriceDirectoriesHANK.txt" !file with paths on separate lines
nfeed = 57


CALL system ("mkdir -p " // trim(OutputDir))

!OPTIONS

!display options
Display              		= 1
ReportNonMonotonicity   	= 0

!run options
CalibrateDiscountRate	= 1
EquilibriumR		 	= 1
ComputeCumulativeMPC 	= 0
DoImpulseResponses 		= 1
DoFeedInPrices 			= 0
SaveTransitionPolicyFns = 0	!warning: will create a very large number of text files 
SaveTime1PolicyFns 		= 0
SaveCumPolicyFnsIRF 	= 0
ComputeDiscountedMPC 	= 0
 
!labor risk options
ReadEarningsProcess = 0
NoLaborSupply		= 0	!only one of the labor supply options
LaborSupplySep		= 1	
LaborSupplyGHH 		= 0
ScaleDisutilityIdio	= 0
ImposeMaxHours 		= 1
PerfectAnnuityMarkets	= 1
GovBondResidualZeroWorld= 1	!imposes closed economy and imputes residual bond holdings to govt
AdjustProdGridFrisch 	= 1
adjfricshgridfrac = 0.85 !fraction of Frisch to adjust by

!calibration options
CalibrateCostFunction		= 0
CalibrateRhoAtInitialGuess  = 0!1
SobolSequenceExploration	= 0
CalibrateRhoInExploration	= 0
CalibrateCloseToSobol 		= 0
MatchRelativeToTargetOutput = 0
ImposeEqumInCalibration 	= 0

!adjustment cost options
OneAssetNoCapital		= 0
SymmetricAdjustmentCost	= 1 !at most 1 of SymmetricAdjustmentCost and NoDepositCost
NoDepositCost 			= 0
AdjustCostFnRatio		= 1
Kappa3LinearOrMax 		= 2	!1 for linear, 2 for max (only valid if AdjustCostFnRatio==1)
ExponAdjustConstFn 		= 0
FixedValueForD 			= 0 !overides other adjustment cost options, must set dexog below
PinKappa1ByKappa02 		= 1!0 !imposes that coefficient on convex term is equal to 1

!transition computation options
PermanentShock 				= 0
SolveFlexPriceTransition	= 0
SolveStickyPriceTransition	= 1
SolveZLBTransition			= 0
ConvergenceRelToOutput 		= 1
UpdateFlexUsingBond 		= 0
OppositeWorldBondFunction 	= 0 !Inverts slope of world bond demand for sticky prices. Use with fixed nominal rate (phitaylor=0)
UseFlexTransitionAsGuess	= 0
FirmDiscountRate			= 5	!1 for rho, 2 for rb initial steady state, 3 for ra initial steady state, 4 for rb transition, 5 for ra transition
BackwardTermInTaylorRule 	= 0
NoChangeLabDisutility 		= 0 !hold labor disutility constant at steady state level when solving HJB (does not affect labor supply decision)
bondelastrelgdp 			= 1.0 !1.0 !bigger for smaller interest rate movements (larger world bond movements), closer to zero for larger interest rate movements (smaller world bond movements). relative to steady state gdp
bondadjust 					= 0.1 !more responsive interest rate on impact when closer to zero
FixBorrowRateTransition 	= 1

!profit distribution options: fractions must sum to 1.0
TaxHHProfitIncome	= 1 !taxes labor subsidy component of profit income at labor tax rate
profdistfracA		= 0.33	!fraction of profits to illiquid equity (set to alpha)
profdistfracB		= 0.0	!fraction of profits to liquid equity
profdistfracW		= 0.67	!fraction of profits to labor subsidy (set to 1-alpha)
profdistfracL 		= 0.0 !0.0 	!fraction of profits lump-sum to households


!government bc options
AdjGovBudgetConstraint 		= 3 !1 for adjust spending, 2 for adjust lump sum taxes, 3 for let debt adjust (choose options below for financing), 4 for adjust proportional tax
GovExpConstantFracOutput 	= 0 !only active if AdjGovBudgetConstraint==3
fixnomgovdebt 				= 0.0 !0 for fixed real government debt, 1 for for fixed nominal government debt, or in between
taxincrstart 		= 0.0 !1.0 !quarters after shock that fiscal policy adjusts
taxincrdecay 		= 0.02 !0.1 !decay rate for tax increase higher for faster decay
GovBCFDScheme 		= 2 !1 for old, 2 for new. 2 requires taxincrstart = 0.0, ony affects B adjusts

!DECOMPOSITION
DoPriceExperiments		= 0
WhichPriceExperiment(1) = 1 !change wage and labor tax only
WhichPriceExperiment(2) = 1 !only change profits
WhichPriceExperiment(3) = 1 !only change profits and wage
WhichPriceExperiment(4) = 1 !only change rb (and rborr if the wedge is fixed out of steady state)
WhichPriceExperiment(5) = 1 !only change ra only
WhichPriceExperiment(6) = 1 !only change illiquid asset drop
WhichPriceExperiment(7) = 1 !only change transfers
WhichPriceExperiment(8) = 1 !change all
WhichPriceExperiment(9) = 0 !change lump transfer by direct effect from government interest payments
WhichPriceExperiment(10) = 0 !change lump transfer, not including direct effect from govt interest payments
WhichPriceExperiment(11) = 0 !change rb, rborr, and change ra by same amount as rb
WhichPriceExperiment(12) = 0 !change ra, and change rb, rborr by same amount as ra
WhichPriceExperiment(13) = 0 !change only proportional labor tax
WhichPriceExperiment(14) = 0 !change rb, rborr, and change ra by same amount as rb, and discount eqm profits at implied ra
WhichPriceExperiment(15) = 0 !change rb, rborr, and change ra by same amount as rb, and discount initimake al profits at implied ra
WhichPriceExperiment(16) = 1 !only change the direct shock

!SHOCKS
IncludeTFPShock=0
IncludeBorrWedgeShock=0
IncludeNewsShock=0
IncludeKappa0Shock=0
IncludeKappa1Shock=0
IncludeMonetaryShock=1
IncludeForwardGuideShock=0
IncludePrefShock=0
IncludeFundLevShock=0
IncludeRiskAversionShock=0
IncludeMarkupShock=0
IncludeGovExpShock=0
IncludeTransferShock=0
IncludeFinWedgeShock=0
IncludeLabWedgeShock=0
IncludeProdDispShock=0
IncludeProdPersShock=0
IncludeTaylorPathShock=0

!SHOCK SIZES AND PERSISTENCE
TFPShockSize 		= -0.02 !log change
TFPShockPers 		= exp(-0.3) !quarterly

BorrWedgeShockSize 		= 0.025 !percentage points
BorrWedgeShockPers 		= exp(-0.15) !quarterly

NewsShockSize 			= -0.01 !log change
NewsShockPers 			= exp(-0.3) !quarterly
NewsShockQtrs 			= 8 !number of quarters in advance

Kappa0ShockSize 		= 0.01 !percentage points
Kappa0ShockPers 		= exp(-0.5) !quarterly

Kappa1ShockSize 		= 0.05 !percentage points
Kappa1ShockPers 		= exp(-0.3) !quarterly

MonetaryShockSize 		= 0.0 !0.005 !percentage points
MonetaryShockPers 		= exp(-0.5) !0.5 !quarterly

mpshockimpact 			= 0.005
TaylorPathShockSize1 		= -10.0 !determines duration it is undone
TaylorPathShockPers1 		= exp(-0.1) !0.5 !quarterly
TaylorPathShockSize2        = -TaylorPathShockSize1 + mpshockimpact
TaylorPathShockPers2        = (TaylorPathShockPers1) **( (TaylorPathShockSize1-mpshockimpact)/TaylorPathShockSize1) !ensures integrates to zero  

ForwardGuideShockSize 	= -0.0025 !percentage points
ForwardGuideShockPers 	= exp(-0.5) !quarterly
ForwardGuideShockQtrs 		= 8 !number of quarters in advance (set phifg below)

PrefShockSize 			= -0.0075 !log change
PrefShockPers 			= exp(-0.3) !quarterly

FundLevShockSize 		= -0.01 !percentage points
FundLevShockPers 		= exp(-0.5) !quarterly

RiskAversionShockSize	= 0.1  !actual units
RiskAversionShockPers	= exp(-0.5) !quarterly

MarkupShockSize 		= 0.05 !percentage points change in markup (converted later into elasticity)
MarkupShockPers 		= exp(-0.3) !quarterly

GovExpShockSize 		= 0.05 !log change
GovExpShockPers 		= exp(-0.3) !quarterly

TransferShockSize 		= 0.40 !log change
TransferShockPers 		= exp(-0.8) !quarterly

FinWedgeShockSize 		= -0.0025 !percentage points
FinWedgeShockPers 		= exp(-0.3) !quarterly

LabWedgeShockSize 		= -0.05 !log change
LabWedgeShockPers 		= exp(-0.5) !quarterly

ProdDispShockSize 		= 0.20 !log change ( in cross-sectional standard deviation)
ProdDispShockPers 		= exp(-0.05) !quarterly
ProdDispScaleDisutilty 	= 1

ProdPersShockSize 		= 0.10 !log change (in transition matrix)
ProdPersShockPers 		= exp(-0.5) !quarterly


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


MatchMeanIll		= 0
MatchKYratio		= 1
MatchMedianIll		= 0
MatchP75Ill			= 0
MatchFracIll0		= 0
MatchMeanLiq		= 1
MatchMedianLiq		= 0
MatchFracLiq0		= 1
MatchFracLiqNeg		= 1
MatchFracIll0Liq0 	= 1!0

defnbclose			= 0.0 !500.0 /(69100.0/4.0)		!fraction of quarterly av gross labor income (69100/4)
defnaclose 			= 0.0 !500.0 /(69100.0/4.0)

ndfls		= 2

!SOLUTION PARAMETERS
maxiter 		= 500
Vtol			= 1.0e-8
maxiterKFE		= 2000
KFEtol			= 1.0e-12 !1.0e-8 !1.0e-7
deltass  		= 1.0e6
deltakfe 		= 1.0e6 !1.0 !1.0e4 
dVamin 			= 1.0e-8 !0.0
dVbmin 			= 1.0e-8

tolequmss		= 1.0e-6
stepequmss		= 0.05
maxiterequmss	= 40 !20
maxiterrho 		= 50 !30 !50
tolrho			= 1.0e-5 !1.0e-8

toltransition	= 1.0e-7
TransitionTimeStepType = 2
deltatransmin	= 1.0/3.0 !1.0 !use with TransitionTimeStepType==1
deltatransmax	= 40.0 !use with TransitionTimeStepType==1
deltatransparam	= 0.35 !use with TransitionTimeStepType==1
deltatransnu 	= 0.05 !use with TransitionTimeStepType==2
maxitertransflex	= 2000 !500 !200 !300
maxitertranssticky	= 2000 !500 !200 !500
stepflextransK  = 0.05
stepflextransB  = 0.001 !make smaller if update using B
stepflextransL  = 0.2
stepstickytransK  = 0.05 !0.01
stepstickytransB  = 0.001 !0.001 !0.00005
stepstickytransL  = 0.2

deltacumcon = 0.01 !deltatransmin !set to a low number like 0.01 for accurate steady state MPCs, and to deltratransmin for IRF consistency

!discount rates
rho		=  0.0125

!preferences
deathrate	= 1.0/(4.0*45.0) !poisson death rate
gam			= 1.0	!risk aversion
prefshock	= 1.0
labwedge 	= 1.0

BondsInUtility 	= 0
bondgam 		= 1.0
bondprefweight 	= 0.01

!liquid assets
rb			 = 0.02/4.0 !liquid return
borrwedge 	 = 0.015 !0.0148846 !0.019663 ! !quarterly wedge between rb and rborr: intermediation cost  
rborr = rb + borrwedge
borrwedgemax = 0.09
blim 		 = -1.0 	!borrowing limit multiple of quarterly output
finwedge	= 0.0 !financial intermediation wedge: households receive rb - finwedge


!withdrawal costs
kappafc_w	= 0.0
kappa0_w	= 0.0 !0.04383 !0.02306 !0.04383
kappa2_w	= 0.40176 !1.71078 !0.40176
kappa3		= 0.03 *2.92/4.0 !0.02 !0.03 *2.92/4.0 !as a fraction of average illiquid assets, approx $10,000
kappa4_w 	= 0.0

IF(PinKappa1ByKappa02==0) kappa1_w	=  0.04336 !0.48236
IF(PinKappa1ByKappa02==1) kappa1_w	= ((1.0-kappa0_w)*(1.0+kappa2_w))**(-1.0/kappa2_w)

dexog 	= 1.0	!only relevant if FixedValueForD==1
kappa2min   = 0.05 !0.5	!to make sure there is enough curvature for calibration

!deposit costs
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
	

priceadjcost 		= 100.0 !price adjustment cost
phitaylor 	= 1.25  ! inflation coefficient in taylor rule, 0.0 for fixed nominal interest rate at ss.
taylorpers 	= 0.1	!strong persistence is 0, low persistence is 1. only active if BackwardTermInTaylorRule=1
IF(BackwardTermInTaylorRule==1) phitaylor = phitaylor/taylorpers

phifg 		= 1.2 !1.1 !0.9 !with forward guidance use phifg instead of phitaylor pre-shock, and phitaylor after (not equal 1)

pi 		= 0.0		!inflation rate: this is steady state target	
pricelev = 1.0 	!steady stat price level if zero inflation 
rnom = rb + pi		!nominal interest rate (fisher equation): this will be constant in taylor rule
mpshock 	= 0.0		!shock to nominal interest rate in taylor rule

elast 	= 10.0 !elasticity of DS aggregator
gap 	= 0.0 !steady state output gap

meanlabeff = 1.0
frisch 		= 1.0 	!frisch elasticity labor supply

!production parameters
drs_Y 		= 0.6
alpha_Y		= 0.33
tfp_Y 		= 1.0

drs_N 		= 0.95
alpha_N		= 0.33
tfp_N 		= 1.0

deprec 		= 0.07/4.0 !depreciation rate

utilelast 	= 0.0
utilelastalpha  = 1.0 + utilelast-alpha_Y*utilelast

!government
labtax 			= 0.30
lumptransferpc 	= 0.05 !40000*labtax/(115000.0*2.92*4.0)
lumptransfer	= lumptransferpc*1.0 !initialize assuming totoutput=1
corptax 		= 0.0
ssdebttogdp 	= 0.26*4 !if foreign sector assumed to hold residual bonds or if equilibrium in liquid assets
govshock 		= 1.0
transfershock 	= 1.0

!calibration targets (relative to quarterly output)
targetMeanIll 		= 2.92 * 4.0 !3.0* 4.0
targetMeanLiq  		= 0.2* 4.0 !0.26 * 4.0 !0.2 * 4.0 
targetMedianIll 	= 0.21 * targetMeanIll
targetP75Ill 		= 0.71 * targetMeanIll
targetMedianLiq 	= 0.085 * targetMeanLiq
targetFracIll0 		= 0.115 + 0.12 !0.115
targetFracLiq0 		= 0.35 !0.30 !0.36 !0.35
targetFracLiqNEG	= 0.15 !0.10
targetFracIll0Liq0 	= 0.10 !0.12 !0.115




!allocate large arrays
ALLOCATE(V(ngpa,ngpb,ngpy),Vnew(ngpa,ngpb,ngpy))
ALLOCATE(u(ngpa,ngpb,ngpy),c(ngpa,ngpb,ngpy),h(ngpa,ngpb,ngpy),d(ngpa,ngpb,ngpy),s(ngpa,ngpb,ngpy),bdot(ngpa,ngpb,ngpy))
ALLOCATE(ccum1(ngpa,ngpb,ngpy),ccum2(ngpa,ngpb,ngpy),ccum4(ngpa,ngpb,ngpy))
ALLOCATE(dcum1(ngpa,ngpb,ngpy),dcum2(ngpa,ngpb,ngpy),dcum4(ngpa,ngpb,ngpy))
ALLOCATE(gjoint(ngpa,ngpb,ngpy),gamarg(ngpa,ngpy),gbmarg(ngpb,ngpy),gvec(naby),gmat(nab,ngpy),gabmarg(ngpa,ngpb),gabcum(ngpa,ngpb))
ALLOCATE(mpc(ngpa,ngpb,ngpy),subeff1ass(ngpa,ngpb,ngpy),subeff2ass(ngpa,ngpb,ngpy),wealtheff1ass(ngpa,ngpb,ngpy),wealtheff2ass(ngpa,ngpb,ngpy))

!allocate solution types
CALL AllocateSolutionType(solnINITSS)
CALL AllocateSolutionType(solnFINALSS)

!allocate cumulative policy types
CALL AllocateCumulativePolicyType(cumINITSS)

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
monetaryshock = .false.
prodispshock = .false.

END SUBROUTINE SetParameters
