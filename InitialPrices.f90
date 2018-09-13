SUBROUTINE InitialPrices

USE Parameters
USE Globals
USE Procedures

IMPLICIT NONE
REAL(8)	:: lmeanwage

!normalize steady state Y=N=P=P_R=1
output = 1.0
varieties = 1.0
totoutput = output*varieties

!steady-state prices and profits determined by elasticity of sub
price_W = (1.0-1.0/elast)
grossprofit_W = price_W*output*(1.0-drs_Y)
netprofit_W = varieties*grossprofit_W
grossprofit_R = (1.0-price_W)*output
netprofit_R = varieties*(1.0-drs_N)*grossprofit_R
profit = netprofit_R + netprofit_W

!compute steady-state target capital-output ratio
CALL KNYfromANY(targetMeanIll,target_K_totoutput_ratio,price_W)
K_totoutput_ratio = target_K_totoutput_ratio


!if solving for equilibrium, these are guesses
labor_Y = hourtarget*meanlabeff*DOT_PRODUCT(1.0-occgrid,occdist)/varieties
labor_N = hourtarget*meanlabeff*DOT_PRODUCT(occgrid,occdist)

capital = totoutput * K_totoutput_ratio
rcapital = (price_W*alpha_Y*drs_Y + (1.0-price_W)*alpha_N*drs_N)/K_totoutput_ratio

capital_Y = price_W*alpha_Y*drs_Y*output/rcapital
tfp_Y = output / (((capital_Y**alpha_Y)*(labor_Y**(1.0-alpha_Y))) ** drs_Y)
wage_Y = price_W*(1.0-alpha_Y)*drs_Y*output/labor_Y
mc_Y = (1.0/tfp_Y)*((rcapital/alpha_Y)**alpha_Y) *((wage_Y/(1.0-alpha_Y))**(1.0-alpha_Y))

capital_N = grossprofit_R*alpha_N*drs_N*varieties/rcapital
tfp_N = varieties / (((capital_N**alpha_N)*(labor_N**(1.0-alpha_N))) ** drs_N)
wage_N = grossprofit_R*(1.0-alpha_N)*drs_N*varieties/labor_N
mc_N = (1.0/tfp_N)*((rcapital/alpha_N)**alpha_N) *((wage_N/(1.0-alpha_N))**(1.0-alpha_N))

ra = rcapital - deprec
dividend_A = profdistfracA*profit*(1.0-corptax)
equity_A = dividend_A/ra
dividend_B = profdistfracB*profit*(1.0-corptax)
equity_B = dividend_B/rb

netwagegrid = (1.0-labtax) * yprodgrid * ( wage_N*yoccgrid + wage_Y*(1.0-yoccgrid) )

!guess labor disutility so that at average wages and average consumption hours =1/3 (sets C/Y = 0.70) for hand-to-mouth. true hours will be lower
lmeanwage = DOT_PRODUCT(1.0-occgrid,occdist)*wage_Y + DOT_PRODUCT(occgrid,occdist)*wage_N
IF (NoLaborSupply==1)	chi	= 0.0 
IF (LaborSupplySep==1)	chi	= lmeanwage*meanlabeff / (0.70 **(-gam) * hourtarget**(1.0/frisch))
IF (LaborSupplyGHH==1)	chi = lmeanwage*meanlabeff / (hourtarget**(1.0/frisch)) 


! IF(OneAssetNoCapital==1) THEN
! 	alpha = 0.0
! 	profdistfrac = 0.0
! 	rcapital = 0.0
! 	wage = mc*(1.0-alpha)*tfp
! 	netwage = (1.0-labtax)*wage
! 	ra = 0.0
! 	meanlabeff = 3.0
! 	kappafc_w = 1.0e8
! 	kappafc_d = 1.0e8
! 	kappa0_w = 100.0
! 	kappa0_d = 100.0
!
! END IF



END SUBROUTINE