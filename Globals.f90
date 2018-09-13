MODULE Globals
USE Parameters

IMPLICIT NONE

!GLOBALS FOR DIRECTORIES
character(len=300)  OutputDir,OutputDirIRF,EarningsProcessDir,FeedPriceDir
character(len=300)  InputParamFile,FeedPriceFile

!OPTIONS GLOBALS
integer	:: Display, ReportNonMonotonicity, NoLaborSupply,LaborSupplySep,LaborSupplyGHH,EquilibriumR, CalibrateCostFunction, DoImpulseResponses,CalibrateDiscountRate,SaveCumPolicyFnsIRF,ComputeDiscountedMPC,DoFeedInPrices
integer	:: PermanentShock, SolveFlexPriceTransition, SolveStickyPriceTransition, SolveZLBTransition, FirmDiscountRate, ImposeEqumInCalibration,UpdateFlexUsingBond, OppositeWorldBondFunction
integer :: ComputeCumulativeMPC, ScaleDisutilityIdio,SaveTransitionPolicyFns,SaveTime1PolicyFns,UseFlexTransitionAsGuess,DoPriceExperiments,TaxHHProfitIncome
integer :: BackwardTermInTaylorRule,FixProfitsOutOfSteadyState,AdjGovBudgetConstraint,ConvergenceRelToOutput,AdjustProdGridFrisch,FixBorrowRateTransition,WhichPriceExperiment(16)
integer :: GovExpConstantFracOutput,GovBondResidualZeroWorld,ExponAdjustConstFn,Kappa3LinearOrMax,PinKappa1ByKappa02,ImposeMaxHours,OneAssetNoCapital,GovBCFDScheme
integer :: ReadEarningsProcess,StickyPriceAlgorithm,ZLBAlgorithm,PerfectAnnuityMarkets,CalibrateRhoAtInitialGuess,FixedValueForD,NoChangeLabDisutility,TransitionTimeStepType
integer :: SymmetricAdjustmentCost,NoDepositCost,AdjustCostFnRatio,SobolSequenceExploration,CalibrateRhoInExploration,MatchRelativeToTargetOutput,CalibrateCloseToSobol,BondsInUtility
logical :: flextransition,stickytransition,zlbtransition,testingGHH
integer :: maxiter,maxitertransflex,maxitertranssticky,maxiterKFE,maxiterequmss,maxiterrho,nfeed
real(8) :: Vtol,KFEtol,deltass,deltakfe,deltacumcon,nendtrans,deltatransvec(Ttransition),cumdeltatrans(Ttransition),toltransition,deltatransparam,deltatransmin,deltatransmax,adjfricshgridfrac,deltatransnu
real(8) :: tolequmss,stepequmss,tolrho,stepflextransK,stepflextransB,stepflextransL,stepstickytransK,stepstickytransB,stepstickytransL,dVamin,dVbmin

!CALIBRATION OPTIONS
integer :: EstimateKappafc, EstimateKappa0, EstimateKappa1, EstimateKappa2, EstimateKappa3, EstimateKappa4, EstimateRho, EstimateBorrWedge,EstimateGamma
integer :: MatchMeanIll, MatchMedianIll, MatchP75Ill, MatchFracIll0, MatchMeanLiq, MatchFracIll0Liq0, MatchMedianLiq, MatchFracLiq0, MatchFracLiqNeg,MatchKYratio
logical :: calibrating,iteratingtransition,exploring
real(8) :: defnaclose,defnbclose

!SHOCK GLOBALS
integer :: IncludeTFPShock,IncludeBorrWedgeShock,IncludeNewsShock,NewsShockQtrs,IncludeKappa0Shock,IncludeMonetaryShock,IncludeKappa1Shock
integer :: IncludeForwardGuideShock,ForwardGuideShockQtrs,IncludePrefShock,IncludeFinWedgeShock,IncludeLabWedgeShock,IncludeProdDispShock,IncludeProdPersShock
integer :: IncludeFundLevShock,IncludeRiskAversionShock,IncludeMarkupShock,IncludeGovExpShock,IncludeTransferShock,IncludeTaylorPathShock
real(8) :: TFPShockSize,BorrWedgeShockSize,NewsShockSize,Kappa0ShockSize,MonetaryShockSize,ForwardGuideShockSize,PrefShockSize,Kappa1ShockSize,TaylorPathShockSize1,TaylorPathShockSize2
real(8) :: FundLevShockSize,MarkupShockSize,RiskAversionShockSize,GovExpShockSize,TransferShockSize,FinWedgeShockSize,LabWedgeShockSize,ProdDispShockSize,ProdPersShockSize
real(8) :: TFPShockPers,BorrWedgeShockPers,NewsShockPers,Kappa0ShockPers,MonetaryShockPers,ForwardGuideShockPers,PrefShockPers,Kappa1ShockPers,TaylorPathShockPers1,TaylorPathShockPers2
real(8) :: FundLevShockPers,MarkupShockPers,RiskAversionShockPers,GovExpShockPers,TransferShockPers,FinWedgeShockPers,LabWedgeShockPers,ProdDispShockPers,ProdPersShockPers
integer :: ProdDispScaleDisutilty
logical :: forwardguide,monetaryshock,prodispshock,nblviolated


!GRIDS GLOBALS
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


!GLOBALS FOR PRODUCTIVITY AND OCCOUPATION PROCESS
real(8), dimension(ngpy) 		:: yprodgrid,yoccgrid,profsharegrid,ydist,netwagegrid		!individual productivity
real(8), dimension(ngpprod) 		:: prodgrid,logprodgrid,proddist
real(8), dimension(ngpocc) 		:: occgrid,occdist

integer, dimension(ngpy)	  		:: occfromy,prodfromy
integer, dimension(ngpocc,ngpprod)	  	:: yfromoccprod

real(8), dimension(ngpy,ngpy)		:: ytrans,ymarkov,ymarkovdiag,ymarkovoff
real(8), dimension(ngpprod,ngpprod)		:: prodtrans,prodmarkov

!GLOBALS FOR VALUE FUNCTIONS AND DECISION
real(8), dimension(:,:,:), allocatable		:: V,Vnew,u,gjoint,c,h,d,s,ccum1,ccum4,ccum2,dcum1,dcum4,dcum2,bdot,mpc,subeff1ass,subeff2ass,wealtheff1ass,wealtheff2ass
real(8), dimension(:,:), allocatable		:: gamarg,gbmarg,gmat,gabcum,gabmarg
real(8), dimension(:), allocatable			:: gvec

!ITERATION GLOBALS
real(8) :: delta

!PARAMETER GLOBALS
real(8)     :: rho,gam,utilcost,chi,frisch,blim,nbl,abl,prefshock,deathrate,meanlabeff,hourtarget
real(8)     :: profdistfracA,profdistfracB,profdistfracW,profdistfracL
real(8)     :: elast,alpha_Y,alpha_N,drs_Y,drs_N,deprec,priceadjcost,phitaylor,phifg,bondelast,borrwedge,mpshock,bondadjust,bondelastrelgdp,taylorpers
real(8) 	:: kappa0_d,kappa1_d,kappa2_d,kappafc_d,kappa0_w,kappa1_w,kappa2_w,kappafc_w,kappa3,kappa4_d,kappa4_w
real(8) 	:: dmin,taxincrstart,taxincrdecay,housedeprec,dexog,utilelast,utilelastalpha,fixnomgovdebt
real(8) 	:: bondprefweight,bondgam,mpshockimpact

!EQUILIBRIUM GLOBALS
real(8)     :: ra,rborr,rcapital,rb,pi,rnom,gap,bond,investment,govexp,taxrev,govbond,worldbond,profit,priceadjust
real(8)     :: totoutput,varieties,output,capital,K_totoutput_ratio,equity_A,equity_B,dividend_A,dividend_B
real(8)     :: capital_Y,labor_Y,wage_Y,mc_Y,tfp_Y,capital_N,labor_N,wage_N,mc_N,tfp_N,price_W,grossprofit_W,netprofit_W,grossprofit_R,netprofit_R
real(8)     :: labtax,lumptransfer,lumptransferpc,ssdebttogdp,corptax,assetdrop_A,assetdrop_B,caputil,govshock,transfershock,finwedge,labwedge,pricelev,prodgridscale,prodmarkovscale
integer 	:: neqmiter,nrhoiter
logical		:: converged,initialSS

!STATISTICS GLOBALS
real(8) 	:: Ea,Eb,Ec,Elabor,Elabor_N,Elabor_Y,Ed,Ewage,Enetlabinc,Egrosslabinc,Enetprofinc,Egrossprofinc,Einc,Ehours,Ehours_N,Ehours_Y,Enw,EbN,EbP,Eadjcost
real(8) 	:: FRACa0,FRACa0close,FRACb0,FRACb0close,FRACb0a0,FRACb0aP,FRACbN,FRACnw0,FRACnw0close,FRACb0a0close
real(8) 	:: PERCa(11),PERCb(11),PERCnw(11),PERCc(11),PERCinc(11)
real(8) 	:: GINIa,GINIb,GINInw,GINIc,GINIinc
real(8) 	:: Ea_nwQ(4),Eb_nwQ(4),Ec_nwQ(4),Einc_nwQ(4),Ea_incQ(4),Eb_incQ(4),Ec_incQ(4),Einc_incQ(4),Ec_nwQ_add(12)
real(8)		:: Ec_bN,Ec_b0close,Ec_b0far

!CALIBRATION GLOBALS
integer					:: nparam,nmoments,objeval,ndfls
real(8), allocatable	:: paramguess(:),paramout(:),paramscale(:),paramlb(:),paramub(:),diagweight(:),nobsvec(:)
real(8)		:: targetMeanIll,targetMeanLiq,targetMedianIll,targetP75Ill,targetMedianLiq,targetFracIll0,targetFracLiq0,targetFracIll0Liq0,targetFracLiqNEG,target_K_totoutput_ratio
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
	real(8), dimension(:,:,:), allocatable	:: V,c,s,h,d,u,gjoint,bdot,mpc,subeff1ass,subeff2ass,wealtheff1ass,wealtheff2ass
	real(8), dimension(:), allocatable 		:: gvec	
	real(8), dimension(:,:), allocatable 	:: gamarg,gbmarg,gmat
	type(CSR), dimension(:), allocatable 	:: A,B,AU
end type

type(SolutionType)		:: solnINITSS,solnFINALSS,solnTRANS(Ttransition)


!EQUILIBRIUM TYPE
type EquilibriumType
! 	real(8)		:: ra,rborr,rcapital,wage,netwage,KYratio,KNratio,mc,rb,tfp,pi,rnom,gap,bond,capital,labor,output,investment,govexp,taxrev,govbond,worldbond,labtax,borrwedge,rho,kappa0_w,kappa1_w,mpshock,prefshock,&
! 					priceadjust,fundlev,elast,gam,fundbond,profit,dividend, lumptransfer,equity,caputil,deprec,tfpadj,illassetdrop,govshock,transfershock,finwedge,labwedge,pricelev,prodgridscale,prodmarkovscale,yprodgrid(ngpy)
	real(8)     :: ra,rborr,rcapital,rb,pi,rnom,gap,bond,investment,govexp,taxrev,govbond,worldbond,profit,priceadjust,totoutput,varieties,output,&
					capital,K_totoutput_ratio,equity_A,equity_B,dividend_A,dividend_B,capital_Y,labor_Y,wage_Y,mc_Y,tfp_Y,capital_N,labor_N,wage_N,mc_N,tfp_N,price_W, &
					grossprofit_W,netprofit_W,grossprofit_R,netprofit_R,labtax,lumptransfer,lumptransferpc,ssdebttogdp,corptax,assetdrop_A,assetdrop_B,caputil,&
					govshock,transfershock,finwedge,labwedge,pricelev,prodgridscale,prodmarkovscale,yprodgrid(ngpy),&
					borrwedge,rho,kappa0_w,kappa1_w,mpshock,prefshock,gam,elast

end type

type(EquilibriumType)	:: equmINITSS,equmFINALSS,equmTRANS(Ttransition)


!DISTRIBUTION STATISTICS TYPE
type DistributionStatsType
	real(8) 	:: Ea,Eb,Ec,Elabor,Elabor_N,Elabor_Y,Ed,Ewage,Enetlabinc,Egrosslabinc,Enetprofinc,Egrossprofinc,Einc,Ehours,Ehours_N,Ehours_Y,Enw, &
					FRACa0,FRACa0close,FRACb0,FRACb0close,FRACb0a0,FRACb0aP,FRACbN,FRACnw0,FRACnw0close,FRACb0a0close, &
					EbN,EbP,Eadjcost,PERCa(11),PERCb(11),PERCnw(11),PERCc(11),PERCinc(11),GINIa,GINIb,GINInw,GINIc,GINIinc, &
					Ea_nwQ(4),Eb_nwQ(4),Ec_nwQ(4),Einc_nwQ(4),Ea_incQ(4),Eb_incQ(4),Ec_incQ(4),Einc_incQ(4),Ec_bN,Ec_b0close,Ec_b0far,Ec_nwQ_add(12)
end type

type(DistributionStatsType)	::statsINITSS,statsFINALSS,statsTRANS(Ttransition)

!CUMULATIVE POLICY FUNCTION TYPE
type CumulativePolicyType
	real(8), dimension(:,:,:), allocatable		:: ccum1,ccum2,ccum4,dcum1,dcum2,dcum4
end type	

type(CumulativePolicyType)		:: cumINITSS


!IMPULSE RESPONSE FUNCTION TYPE
type ImpulseResponseType
	type(SolutionType)			:: solnFLEX(Ttransition),solnSTICKY(Ttransition),solnZLB(Ttransition)
	type(EquilibriumType)		:: equmFLEX(Ttransition),equmSTICKY(Ttransition),equmZLB(Ttransition)
	type(DistributionStatsType)	:: statsFLEX(Ttransition),statsSTICKY(Ttransition),statsZLB(Ttransition)
	type(CumulativePolicyType)	:: cumFLEX,cumSTICKY,cumZLB
end type	

! type(ImpulseResponseType), dimension(0:0), target 	:: irfstruct
type(ImpulseResponseType), target 	:: irfstruct
type(ImpulseResponseType), target					:: irfpriceexp
type(ImpulseResponseType), pointer 					:: irfsave,irfpointer,irfpointer_fs 


!THREADPRIVATE GLOBALS
real(8)							:: gbdrift,gnetwage,gill,gVb,gidioprod,gidioscale,gidioprodss
!$OMP THREADPRIVATE(gbdrift,gnetwage,gill,gVb,gidioprod,gidioscale,gidioprodss)

END MODULE Globals
