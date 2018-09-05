MODULE Globals
USE Parameters

IMPLICIT NONE

!GLOBALS FOR DIRECTORIES
character(len=100)  OutputDir,OutputDirIRF,EarningsProcessDir
character(len=100)  InputParamFile

!OPTIONS GLOBALS
integer	:: Display, ReportNonMonotonicity, LaborSupply, EquilibriumR, CalibrateCostFunction, DoImpulseResponses,CalibrateDiscountRate
integer	:: PermanentShock, SolveFlexPriceTransition, SolveStickyPriceTransition, SolveZLBTransition, FirmDiscountRate, ImposeEqumInCalibration,UpdateRbFromMarketClearing
integer :: IncludeUnemploymentState, ComputeCumulativeMPC, ScaleGHHIdiosyncratic,SaveTransitionPolicyFns,SaveTime1PolicyFns,UseFlexTransitionAsGuess,DoPriceExperiments
integer :: BackwardTermInTaylorRule,FixProfitsOutOfSteadyState,DividendSmoothing,AssumePathRealRate,AdjGovBudgetConstraint,UpdateUsingR,GuessZeroInflation,ConvergenceRelToOutput
integer :: RetainedEarningsInBondMarket,GovExpConstantFracOutput,TaxRetainedEarnings,GovBondResidualZeroWorld,RetainedEarningsAsCapital,ExponAdjustConstFn,Kappa3LinearOrMax,PinKappa1ByKappa02
integer :: DoFiscalStimulus,ImposeGuessSettlesDown,ReadEarningsProcess,StickyPriceAlgorithm,ZLBAlgorithm,PerfectAnnuityMarkets,CalibrateRhoAtInitialGuess,FixedValueForD,NoChangeLabDisutility
integer :: SymmetricAdjustmentCost,NoDepositCost,AdjustCostFnRatio,SobolSequenceExploration,CalibrateRhoInExploration,MatchRelativeToTargetOutput,CalibrateCloseToSobol,DividendFundLumpSum
logical :: flextransition,stickytransition,zlbtransition,fsptransition,testingGHH
integer :: maxiter,maxitertransflex,maxitertranssticky,maxiterKFE,maxiterequmss,maxiterrho
real(8) :: Vtol,KFEtol,deltass,deltakfe,nendtrans,deltatransvec(Ttransition),cumdeltatrans(Ttransition),toltransition,deltatransparam,deltatransmin,deltatransmax
real(8) :: tolequmss,stepequmss,tolrho,stepflextransK,stepflextransB,stepstickytransK,stepstickytransB,stepstickytransY,stepstickytransRb,dVamin,dVbmin

!CALIBRATION OPTIONS
integer :: EstimateKappafc, EstimateKappa0, EstimateKappa1, EstimateKappa2, EstimateKappa3, EstimateKappa4, EstimateRho, EstimateBorrWedge,EstimateGamma
integer :: MatchMeanIll, MatchMedianIll, MatchP75Ill, MatchFracIll0, MatchMeanLiq, MatchFracIll0Liq0, MatchMedianLiq, MatchFracLiq0, MatchFracLiqNeg,MatchKYratio
logical :: calibrating,iteratingtransition,exploring
real(8) :: defnaclose,defnbclose

!SHOCK GLOBALS
integer :: IncludeTFPShock,IncludeBorrWedgeShock,IncludeNewsShock,NewsShockQtrs,IncludeKappafcShock,IncludeMonetaryShock
integer :: IncludeForwardGuideShock,ForwardGuideShockQtrs,IncludeStDevYShock,IncludePrefShock
integer :: IncludeFundLevShock,IncludeRiskAversionShock,IncludeMarkupShock,IncludeGovExpShock,IncludeTransferShock
integer :: ForwardGuideFixNomPreShock
real(8) :: TFPShockSize,BorrWedgeShockSize,NewsShockSize,KappafcShockSize,MonetaryShockSize,ForwardGuideShockSize,StDevYShockSize,PrefShockSize,FundLevShockSize,MarkupShockSize,RiskAversionShockSize
real(8) :: TFPShockPers,BorrWedgeShockPers,NewsShockPers,KappafcShockPers,MonetaryShockPers,ForwardGuideShockPers,StDevYShockPers,PrefShockPers,FundLevShockPers,MarkupShockPers,RiskAversionShockPers 
logical :: forwardguide

!GRIDS GLOBALS
real(8), dimension(ngpy) 		:: ygrid,netlabincgrid,grosslabincgrid,lgrid,logygrid,labdisutilgrid,netinc2illgrid,netinc2liqgrid		!wages,hours,earnings
real(8), dimension(ngpa)  		:: agrid,adelta,adrift        !illiquid asset
real(8), dimension(ngpb)  		:: bgrid,bdelta,bdrift        !liquid asset
real(8), dimension(nab)  		:: abdelta
real(8), dimension(naby)  		:: abydelta
integer, dimension(naby)  		:: afromaby,bfromaby,yfromaby,abfromaby
integer, dimension(nab)  		:: afromab,bfromab
integer, dimension(ngpa,ngpb,ngpy)  :: abyfromaby
integer, dimension(ngpa,ngpb)  		:: abfromab
real(8), dimension(ngpb-1)  		:: dbgrid        !liquid asset spacing
real(8), dimension(ngpa-1)  		:: dagrid        !illiquid asset spacing

!GLOBALS FOR INCOME RISK
real(8), dimension(ngpy)		    :: ydist
real(8), dimension(ngpy,ngpy)		:: ytrans,ymarkov,ymarkovdiag,ymarkovoff
real(8) 							:: yrho,ysig,ylogvar,ylambda,jobloss,jobfind

!GLOBALS FOR VALUE FUNCTIONS AND DECISION
real(8), dimension(:,:,:), allocatable		:: V,Vnew,u,gjoint,c,p,d,H,s,ccum1,ccum4,ccum2,bdot,fspamount
real(8), dimension(:,:), allocatable		:: gamarg,gbmarg,gmat,gabcum,gabmarg
real(8), dimension(:), allocatable			:: gvec

!ITERATION GLOBALS
real(8) :: delta

!PARAMETER GLOBALS
real(8)     :: rho,gam,utilcost,chi,frisch,blim,nbl,abl,uiben,prefshock,fundlev,fundbond,rhousing,housemaint,deathrate,flowcon,directdepfrac,directdepmax,directdepmin,directdepmaxfrac,directdepminfrac
real(8)     :: elast,alpha,deprec,alphatilde,theta,phitaylor,bondelast,borrwedge,mpshock,capitaladjcost,investadjcost,bondadjust,bondelastrelgdp,taylorpers,firmgamma,divsmooth
real(8) 	:: housefrac,flowamnt,flowgam,kappa0_d,kappa1_d,kappa2_d,kappafc_d,kappa0_w,kappa1_w,kappa2_w,kappafc_w,kappa3,kappa4_d,kappa4_w
real(8) 	:: dmin,aendogmax,nondurwgt,nondurwgttilde,taxincrstart,taxincrdecay,housedeprec,dexog,wageflex,utilelast,utilelastalpha,operatecost

!EQUILIBRIUM GLOBALS
real(8)     :: ra,rborr,rcapital,wage,netwage,KYratio,KNratio,mc,rb,tfp,pi,rnom,gap,bond,capital,labor,output,investment,govexp,taxrev,govbond,worldbond,profit,dividend,divrate,priceadjust,intfirmbond,equity
real(8)     :: labtax,lumptransfer,lumptransferpc,lumptransfercutoff,ssdebttogdp,corptax,illequitydrop,caputil,tfpadj
integer 	:: neqmiter
logical		:: converged,initialSS

!STATISTICS GLOBALS
real(8) 	:: Ea,Eb,Ec,Erent,Ed,Ewage,Enetlabinc,Egrosslabinc,Einc,Ehours,Enw,EbN,EbP,Eadjcost,Efsp
real(8) 	:: FRACa0,FRACa0close,FRACb0,FRACb0close,FRACb0a0,FRACb0aP,FRACbN,FRACnw0,FRACnw0close,FRACb0a0close
real(8) 	:: PERCa(11),PERCb(11),PERCnw(11),PERCc(11),PERCinc(11)
real(8) 	:: GINIa,GINIb,GINInw,GINIc,GINIinc
real(8) 	:: Ea_nwQ(4),Eb_nwQ(4),Ec_nwQ(4),Einc_nwQ(4),Ea_incQ(4),Eb_incQ(4),Ec_incQ(4),Einc_incQ(4)
real(8)		:: Ec_bN,Ec_b0close,Ec_b0far

!CALIBRATION GLOBALS
integer					:: nparam,nmoments,objeval,ndfls
real(8), allocatable	:: paramguess(:),paramout(:),paramscale(:),paramlb(:),paramub(:),diagweight(:),nobsvec(:)
real(8)		:: targetMeanIll,targetMeanLiq,targetMedianIll,targetP75Ill,targetMedianLiq,targetFracIll0,targetFracLiq0,targetFracIll0Liq0,targetFracLiqNEG,targetKYratio,targetKYtotalratio
real(8)		:: modelMeanIll,modelMeanLiq,modelMedianIll,modelP75Ill,modelMedianLiq,modelFracIll0,modelFracLiq0,modelFracIll0Liq0,modelFracLiqNEG,modelKYratio
real(8) 	:: kappa2min,borrwedgemax

!SOBOL EXPLORATION GLOBALS
integer			:: nsobol,nsobolskip
real(8), allocatable	:: sobolseq(:,:),paramsobol(:,:),sobolmom(:,:),sobolobj(:)


!SPARSE MATRIX TYPES
type COO
	integer	:: nz
	real(8), dimension(:), allocatable :: val
	integer, dimension(:), allocatable :: row, col
end type

type CSR
	integer	:: n,nz
	real(8), allocatable :: val(:)
	integer, allocatable :: row(:), col(:)
end type

!SPARSE MATRICES
type(COO), dimension(ngpy)	:: ACOO,AUCOO
type(CSR), dimension(ngpy)	:: ACSR,BCSR,AUCSR

!SOLUTION TYPE
type SolutionType
	real(8), dimension(:,:,:), allocatable	:: V,c,p,d,u,gjoint,bdot
	real(8), dimension(:), allocatable 		:: gvec	
	real(8), dimension(:,:), allocatable 	:: gamarg,gbmarg,gmat
	type(CSR), dimension(:), allocatable 	:: A,B,AU
end type

type(SolutionType)		:: solnINITSS,solnFINALSS,solnTRANS(Ttransition)


!EQUILIBRIUM TYPE
type EquilibriumType
	real(8)		:: ra,rborr,rcapital,wage,netwage,KYratio,KNratio,mc,rb,tfp,pi,rnom,gap,bond,capital,labor,output,investment,govexp,taxrev,govbond,worldbond,labtax,borrwedge,rho,kappafc_w,mpshock,ysig,prefshock,&
					priceadjust,fundlev,elast,gam,fundbond,profit,dividend,divrate,intfirmbond,lumptransfer,equity,caputil,deprec,tfpadj
end type

type(EquilibriumType)	:: equmINITSS,equmFINALSS,equmTRANS(Ttransition)


!DISTRIBUTION STATISTICS TYPE
type DistributionStatsType
	real(8) 	:: Ea,Eb,Ec,Erent,Ed,Ewage,Enetlabinc,Egrosslabinc,Einc,Ehours,Enw,FRACa0,FRACa0close,FRACb0,FRACb0close,FRACb0a0,FRACb0aP,FRACbN,FRACnw0,FRACnw0close,FRACb0a0close, &
					EbN,EbP,Eadjcost,Efsp,PERCa(11),PERCb(11),PERCnw(11),PERCc(11),PERCinc(11),GINIa,GINIb,GINInw,GINIc,GINIinc, &
					Ea_nwQ(4),Eb_nwQ(4),Ec_nwQ(4),Einc_nwQ(4),Ea_incQ(4),Eb_incQ(4),Ec_incQ(4),Einc_incQ(4),Ec_bN,Ec_b0close,Ec_b0far	
end type

type(DistributionStatsType)	::statsINITSS,statsFINALSS,statsTRANS(Ttransition)


!IMPULSE RESPONSE FUNCTION TYPE
type ImpulseResponseType
	type(SolutionType)			:: solnFLEX(Ttransition),solnSTICKY(Ttransition),solnZLB(Ttransition)
	type(EquilibriumType)		:: equmFLEX(Ttransition),equmSTICKY(Ttransition),equmZLB(Ttransition)
	type(DistributionStatsType)	:: statsFLEX(Ttransition),statsSTICKY(Ttransition),statsZLB(Ttransition)
end type	

type(ImpulseResponseType), dimension(0:nfs), target :: irfstruct
type(ImpulseResponseType), target					:: irfpriceexp
type(ImpulseResponseType), pointer 					:: irfsave,irfpointer,irfpointer_fs 

!FISCAL STIMULUS TYPE
type FiscalStimulusType
	real(8), dimension(:,:,:), allocatable	:: fspamount
	real(8)		:: govamount
	real(8)		:: fspstart,fspend,govstart,govend,labtaxincr,labtaxstart,labtaxend
end type

type(FiscalStimulusType), dimension(nfs), target 	::	fsconfig
type(FiscalStimulusType), pointer 					::	fspointer

!THREADPRIVATE GLOBALS
real(8)							:: gbdrift,glabdisutil,gill,gVb,gc
!$OMP THREADPRIVATE(gbdrift,glabdisutil,gill,gVb,gc)

END MODULE Globals