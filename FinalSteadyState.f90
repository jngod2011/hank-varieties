SUBROUTINE FinalSteadyState

USE Parameters
USE Globals
USE Procedures

IMPLICIT NONE

INTEGER		:: iy

IF (Display>=1) write(*,*) "Solving for final steady state"

IF(PermanentShock==0) THEN
	DO iy = 1,ngpy
		ALLOCATE(solnFINALSS%A(iy)%val(solnINITSS%A(iy)%nz),solnFINALSS%A(iy)%row(solnINITSS%A(iy)%n+1),solnFINALSS%A(iy)%col(solnINITSS%A(iy)%nz))
		ALLOCATE(solnFINALSS%B(iy)%val(solnINITSS%B(iy)%nz), solnFINALSS%B(iy)%col(solnINITSS%B(iy)%nz), solnFINALSS%B(iy)%row(solnINITSS%B(iy)%n+1))
	END DO
	
	solnFINALSS = solnINITSS
	equmFINALSS = equmINITSS
	statsFINALSS = statsINITSS

ELSE IF(PermanentShock==1) THEN

	write(*,*) 'WARNING: THIS PART OF THE CODE MAY NOT WORK, NEEDS TO BE CHECKED'
	
	CALL SolveSteadyStateEqum
	
	DO iy = 1,ngpy
		ALLOCATE(solnFINALSS%A(iy)%val(ACSR(iy)%nz))
		ALLOCATE(solnFINALSS%A(iy)%row(ACSR(iy)%n+1))
		ALLOCATE(solnFINALSS%A(iy)%col(ACSR(iy)%nz))
		ALLOCATE(solnFINALSS%B(iy)%val(BCSR(iy)%nz))
		ALLOCATE(solnFINALSS%B(iy)%row(BCSR(iy)%n+1))
		ALLOCATE(solnFINALSS%B(iy)%col(BCSR(iy)%nz))
	END DO
	
	solnFINALSS%A = ACSR
	solnFINALSS%B = BCSR
	solnFINALSS%V = V
	solnFINALSS%u = u  
	solnFINALSS%c = c
	solnFINALSS%s = s
	solnFINALSS%h = h
	solnFINALSS%d = d  
	solnFINALSS%bdot = bdot

	solnFINALSS%gjoint = gjoint
	solnFINALSS%gamarg = gamarg
	solnFINALSS%gbmarg = gbmarg
	solnFINALSS%gvec = gvec
	solnFINALSS%gmat = gmat

	CALL DistributionStatistics

equmINITSS = EquilibriumType(	ra,rborr,rcapital,rb,pi,rnom,gap,bond,investment,govexp,taxrev,govbond,worldbond,profit,priceadjust,totoutput,varieties,output,&
					capital,K_totoutput_ratio,equity_A,equity_B,dividend_A,dividend_B,capital_Y,labor_Y,wage_Y,mc_Y,tfp_Y,capital_N,labor_N,wage_N,mc_N,tfp_N,price_W, &
					grossprofit_W,netprofit_W,grossprofit_R,netprofit_R,labtax,lumptransfer,lumptransferpc,ssdebttogdp,corptax,illassetdrop,caputil,&
					govshock,transfershock,finwedge,labwedge,pricelev,prodgridscale,prodmarkovscale,yprodgrid,&
					borrwedge,rho,kappa0_w,kappa1_w,mpshock,prefshock,gam,elast)

! 	statsFINALSS = DistributionStatsType(Ea,Eb,Ec,Elabor,Ed,Ewage,Enetlabinc,Egrosslabinc,Enetprofinc,Egrossprofinc,Einc,Ehours,Enw,FRACa0,FRACa0close,FRACb0,FRACb0close,FRACb0a0,FRACb0aP,FRACbN,FRACnw0,FRACnw0close,FRACb0a0close, &
! 										EbN,EbP,Eadjcost,PERCa,PERCb,PERCnw,PERCc,PERCinc,GINIa,GINIb,GINInw,GINIc,GINIinc, &
! 										Ea_nwQ,Eb_nwQ,Ec_nwQ,Einc_nwQ,Ea_incQ,Eb_incQ,Ec_incQ,Einc_incQ,Ec_bN,Ec_b0close,Ec_b0far,Ec_nwQ_add)

END IF



END SUBROUTINE FinalSteadyState