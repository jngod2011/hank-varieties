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
	solnFINALSS%p = p
	solnFINALSS%d = d  
	solnFINALSS%bdot = bdot

	solnFINALSS%gjoint = gjoint
	solnFINALSS%gamarg = gamarg
	solnFINALSS%gbmarg = gbmarg
	solnFINALSS%gvec = gvec
	solnFINALSS%gmat = gmat

	CALL DistributionStatistics
	equmFINALSS = EquilibriumType(	ra,rborr,rcapital,wage,netwage,KYratio,KNratio,mc,rb,tfp,pi,rnom,gap,bond,capital,labor,output,investment,govexp,taxrev,govbond,worldbond,labtax,&
									borrwedge,rho,kappafc_w,mpshock,ysig,prefshock,priceadjust,fundlev,elast,gam,fundbond,profit,dividend,divrate,intfirmbond,lumptransfer,equity,caputil,deprec,tfpadj)
	statsFINALSS = DistributionStatsType(Ea,Eb,Ec,Erent,Ed,Ewage,Enetlabinc,Egrosslabinc,Einc,Ehours,Enw,FRACa0,FRACa0close,FRACb0,FRACb0close,FRACb0a0,FRACb0aP,FRACbN,FRACnw0,FRACnw0close,FRACb0a0close, &
										EbN,EbP,Eadjcost,Efsp,PERCa,PERCb,PERCnw,PERCc,PERCinc,GINIa,GINIb,GINInw,GINIc,GINIinc, &
										Ea_nwQ,Eb_nwQ,Ec_nwQ,Einc_nwQ,Ea_incQ,Eb_incQ,Ec_incQ,Einc_incQ,Ec_bN,Ec_b0close,Ec_b0far)

END IF



END SUBROUTINE FinalSteadyState