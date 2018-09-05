SUBROUTINE SaveIRFOutput

USE Parameters
USE Globals
USE Procedures

IMPLICIT NONE
INTEGER		:: iy,it,ip,iq
CHARACTER	:: lstring*80,lt*80
REAL(8)		:: lmat(Ttransition,11),lmat4(Ttransition,4),lmat12(Ttransition,12)

!flexible price transition
IF(SolveFlexPriceTransition==1) THEN
	CALL system ("mkdir -p " // trim(OutputDirIRF) // "FLEX")	

	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'ra.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmFLEX%ra)
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'rcapital.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmFLEX%rcapital)
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'wage.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmFLEX%wage)
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'netwage.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmFLEX%netwage)
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'KYratio.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmFLEX%KYratio)
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'KNratio.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmFLEX%KNratio)
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'mc.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmFLEX%mc)
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'rb.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmFLEX%rb)
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'rborr.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmFLEX%rborr)
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'tfp.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmFLEX%tfp)
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'pi.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmFLEX%pi)
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'rnom.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmFLEX%rnom)
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'gap.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmFLEX%gap)
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'capital.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmFLEX%capital)
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'bond.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmFLEX%bond)
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'labor.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmFLEX%labor)
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'output.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmFLEX%output)
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'investment.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmFLEX%investment)
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'taxrev.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmFLEX%taxrev)
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'govexp.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmFLEX%govexp)
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'govbond.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmFLEX%govbond)
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'worldbond.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmFLEX%worldbond)
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'borrwedge.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmFLEX%borrwedge)
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'rho.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmFLEX%rho)
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'kappa0_w.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmFLEX%kappa0_w)
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'kappa1_w.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmFLEX%kappa1_w)
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'mpshock.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmFLEX%mpshock)
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'prefshock.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmFLEX%prefshock)
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'elast.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmFLEX%elast)
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'fundlev.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmFLEX%fundlev)
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'gam.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmFLEX%gam)
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'fundbond.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmFLEX%fundbond)
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'priceadjust.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmFLEX%priceadjust)
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'profit.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmFLEX%profit)
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'dividend.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmFLEX%dividend)
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'divrate.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmFLEX%divrate)
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'lumptransfer.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmFLEX%lumptransfer)
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'equity.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmFLEX%equity)
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'caputil.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmFLEX%caputil)
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'deprec.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmFLEX%deprec)
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'tfpadj.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmFLEX%tfpadj)
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'illassetdrop.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmFLEX%illassetdrop)
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'govshock.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmFLEX%govshock)
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'transfershock.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmFLEX%transfershock)
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'finwedge.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmFLEX%finwedge)
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'labwedge.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmFLEX%labwedge)
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'pricelev.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmFLEX%pricelev)
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'prodgridscale.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmFLEX%prodgridscale)
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'prodmarkovscale.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmFLEX%prodmarkovscale)

	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'Ea.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsFLEX%Ea)
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'Eb.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsFLEX%Eb)
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'Ec.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsFLEX%Ec)
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'Elabor.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsFLEX%Elabor)
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'Ed.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsFLEX%Ed)
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'Ewage.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsFLEX%Ewage)
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'Enetlabinc.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsFLEX%Enetlabinc)
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'Egrosslabinc.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsFLEX%Egrosslabinc)
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'Einc.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsFLEX%Einc)
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'Ehours.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsFLEX%Ehours)
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'Enw.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsFLEX%Enw)

	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'FRACa0.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsFLEX%FRACa0)
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'FRACa0close.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsFLEX%FRACa0close)
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'FRACb0.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsFLEX%FRACb0)
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'FRACb0close.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsFLEX%FRACb0close)
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'FRACb0a0.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsFLEX%FRACb0a0)
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'FRACb0aP.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsFLEX%FRACb0aP)
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'FRACbN.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsFLEX%FRACbN)
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'FRACnw0.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsFLEX%FRACnw0)
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'FRACnw0close.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsFLEX%FRACnw0close)

	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'EbN.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsFLEX%EbN)
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'EbP.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsFLEX%EbP)
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'Eadjcost.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsFLEX%Eadjcost)

	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'GINIa.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsFLEX%GINIa)
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'GINIb.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsFLEX%GINIb)
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'GINInw.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsFLEX%GINInw)
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'GINIc.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsFLEX%GINIc)
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'GINIinc.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsFLEX%GINIinc)

	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'Ec_bN.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsFLEX%Ec_bN)
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'Ec_b0close.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsFLEX%Ec_b0close)
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'Ec_b0far.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsFLEX%Ec_b0far)
	
	
	DO ip = 1,11;lmat(:,ip) = irfsave%statsFLEX%PERCa(ip); END DO	
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'PERCa.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,11,lmat)

	DO ip = 1,11;lmat(:,ip) = irfsave%statsFLEX%PERCb(ip); END DO	
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'PERCb.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,11,lmat)

	DO ip = 1,11;lmat(:,ip) = irfsave%statsFLEX%PERCnw(ip); END DO	
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'PERCnw.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,11,lmat)

	DO ip = 1,11;lmat(:,ip) = irfsave%statsFLEX%PERCc(ip); END DO	
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'PERCc.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,11,lmat)

	DO ip = 1,11;lmat(:,ip) = irfsave%statsFLEX%PERCinc(ip); END DO	
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'PERCinc.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,11,lmat)

	DO iq = 1,4;lmat4(:,iq) = irfsave%statsFLEX%Ea_nwQ(iq); END DO	
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'Ea_nwQ.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,4,lmat4)

	DO iq = 1,4;lmat4(:,iq) = irfsave%statsFLEX%Eb_nwQ(iq); END DO	
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'Eb_nwQ.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,4,lmat4)

	DO iq = 1,4;lmat4(:,iq) = irfsave%statsFLEX%Ec_nwQ(iq); END DO	
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'Ec_nwQ.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,4,lmat4)

	DO iq = 1,4;lmat4(:,iq) = irfsave%statsFLEX%Einc_nwQ(iq); END DO	
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'Einc_nwQ.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,4,lmat4)

	DO iq = 1,4;lmat4(:,iq) = irfsave%statsFLEX%Ea_incQ(iq); END DO	
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'Ea_incQ.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,4,lmat4)

	DO iq = 1,4;lmat4(:,iq) = irfsave%statsFLEX%Eb_incQ(iq); END DO	
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'Eb_incQ.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,4,lmat4)

	DO iq = 1,4;lmat4(:,iq) = irfsave%statsFLEX%Ec_incQ(iq); END DO	
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'Ec_incQ.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,4,lmat4)

	DO iq = 1,4;lmat4(:,iq) = irfsave%statsFLEX%Einc_incQ(iq); END DO	
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'Einc_incQ.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,4,lmat4)

	DO iq = 1,12;lmat12(:,iq) = irfsave%statsFLEX%Ec_nwQ_add(iq); END DO	
	OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'Ec_nwQ_add.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,12,lmat12)


	IF(SaveTime1PolicyFns==1) THEN
		DO iy = 1,ngpy
			IF (iy<10) WRITE(UNIT=lstring, FMT='(I1)') iy
			IF (iy >= 10) WRITE(UNIT=lstring, FMT='(I2)') iy
			
			OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'V_T1_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,irfsave%solnFLEX(1)%V(:,:,iy))
			OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/"  // 'dep_T1_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,irfsave%solnFLEX(1)%d(:,:,iy))
			OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/"  // 'con_T1_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,irfsave%solnFLEX(1)%c(:,:,iy))
			OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/"  // 'sav_T1_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,irfsave%solnFLEX(1)%s(:,:,iy))
			OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/"  // 'hour_T1_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,irfsave%solnFLEX(1)%h(:,:,iy))
			OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/"  // 'bdot_T1_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,irfsave%solnFLEX(1)%bdot(:,:,iy))
			OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/"  // 'gjoint_T1_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,irfsave%solnFLEX(1)%gjoint(:,:,iy))
			IF(ComputeDiscountedMPC==1) THEN
				OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/"  // 'mpc_T1_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,irfsave%solnFLEX(1)%mpc(:,:,iy))
				OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/"  // 'subeff1ass_T1_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,irfsave%solnFLEX(1)%subeff1ass(:,:,iy))
				OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/"  // 'subeff2ass_T1_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,irfsave%solnFLEX(1)%subeff2ass(:,:,iy))
				OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/"  // 'wealtheff1ass_T1_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,irfsave%solnFLEX(1)%wealtheff1ass(:,:,iy))
				OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/"  // 'wealtheff2ass_T1_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,irfsave%solnFLEX(1)%wealtheff2ass(:,:,iy))
			END IF
		
		
		END DO
	END IF
	
	IF (SaveCumPolicyFnsIRF==1) THEN
		DO iy = 1,ngpy
			IF (iy<10) WRITE(UNIT=lstring, FMT='(I1)') iy
			IF (iy >= 10) WRITE(UNIT=lstring, FMT='(I2)') iy
			
			OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'ccum1_T1_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,irfsave%cumFLEX%ccum1(:,:,iy))
			OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'ccum2_T1_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,irfsave%cumFLEX%ccum2(:,:,iy))
			OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'ccum4_T1_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,irfsave%cumFLEX%ccum4(:,:,iy))
			OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'dcum1_T1_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,irfsave%cumFLEX%dcum1(:,:,iy))
			OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'dcum2_T1_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,irfsave%cumFLEX%dcum2(:,:,iy))
			OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/" // 'dcum4_T1_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,irfsave%cumFLEX%dcum4(:,:,iy))
		END DO
		
	END IF	
	
	!transitions steady state distributions and policy functions
	IF(SaveTransitionPolicyFns==1) THEN
		CALL system ("mkdir -p " // trim(OutputDirIRF) // "FLEX/FuncDist/")	
	
		DO it = 1,Ttransition
		IF (it<10) WRITE(UNIT=lt, FMT='(I1)') it
		IF (it>=10 .and. it<100) WRITE(UNIT=lt, FMT='(I2)') it
		IF (it>=100) WRITE(UNIT=lt, FMT='(I3)') it

		DO iy = 1,ngpy
			IF (iy<10) WRITE(UNIT=lstring, FMT='(I1)') iy
			IF (iy >= 10) WRITE(UNIT=lstring, FMT='(I2)') iy

 			OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/FuncDist/" // 'V' // trim(lt) // '_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,irfsave%solnFLEX(it)%V(:,:,iy))
 			OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/FuncDist/" // 'dep' // trim(lt) // '_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,irfsave%solnFLEX(it)%d(:,:,iy))
 			OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/FuncDist/" // 'con' // trim(lt) // '_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,irfsave%solnFLEX(it)%c(:,:,iy))
 			OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/FuncDist/" // 'bdot' // trim(lt) // '_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,irfsave%solnFLEX(it)%bdot(:,:,iy))
			OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/FuncDist/" // 'gjoint' // trim(lt) // '_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,irfsave%solnFLEX(it)%gjoint(:,:,iy))
 			OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/FuncDist/" // 'B' // trim(lt) // '_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,irfsave%solnFLEX(it)%B(iy)%val)

		END DO

		OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/FuncDist/" // 'gamarg' // trim(lt) // '.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpy,irfsave%solnFLEX(it)%gamarg)
		OPEN(3, FILE = trim(OutputDirIRF) // "FLEX/FuncDist/" // 'gbmarg' // trim(lt) // '.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpb,ngpy,irfsave%solnFLEX(it)%gbmarg)

	END DO
	END IF
END IF

!sticky price transition
IF(SolveStickyPriceTransition==1) THEN
	CALL system ("mkdir -p " // trim(OutputDirIRF) // "STICKY")	

	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'ra.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%ra)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'rcapital.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%rcapital)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'wage.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%wage)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'netwage.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%netwage)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'KYratio.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%KYratio)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'KNratio.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%KNratio)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'mc.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%mc)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'rb.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%rb)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'rborr.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%rborr)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'tfp.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%tfp)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'pi.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%pi)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'rnom.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%rnom)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'gap.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%gap)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'capital.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%capital)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'bond.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%bond)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'labor.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%labor)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'output.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%output)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'investment.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%investment)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'taxrev.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%taxrev)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'govexp.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%govexp)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'govbond.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%govbond)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'worldbond.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%worldbond)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'borrwedge.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%borrwedge)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'rho.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%rho)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'kappa0_w.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%kappa0_w)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'kappa1_w.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%kappa1_w)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'mpshock.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%mpshock)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'prefshock.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%prefshock)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'elast.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%elast)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'fundlev.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%fundlev)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'gam.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%gam)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'fundbond.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%fundbond)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'priceadjust.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%priceadjust)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'profit.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%profit)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'dividend.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%dividend)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'divrate.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%divrate)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'lumptransfer.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%lumptransfer)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'equity.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%equity)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'caputil.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%caputil)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'deprec.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%deprec)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'tfpadj.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%tfpadj)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'illassetdrop.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%illassetdrop)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'govshock.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%govshock)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'transfershock.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%transfershock)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'finwedge.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%finwedge)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'labwedge.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%labwedge)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'pricelev.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%pricelev)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'prodgridscale.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%prodgridscale)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'prodmarkovscale.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%prodmarkovscale)

	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'Ea.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsSTICKY%Ea)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'Eb.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsSTICKY%Eb)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'Ec.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsSTICKY%Ec)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'Elabor.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsSTICKY%Elabor)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'Ed.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsSTICKY%Ed)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'Ewage.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsSTICKY%Ewage)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'Enetlabinc.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsSTICKY%Enetlabinc)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'Egrosslabinc.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsSTICKY%Egrosslabinc)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'Einc.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsSTICKY%Einc)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'Ehours.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsSTICKY%Ehours)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'Enw.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsSTICKY%Enw)

	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'FRACa0.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsSTICKY%FRACa0)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'FRACa0close.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsSTICKY%FRACa0close)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'FRACb0.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsSTICKY%FRACb0)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'FRACb0close.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsSTICKY%FRACb0close)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'FRACb0a0.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsSTICKY%FRACb0a0)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'FRACb0aP.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsSTICKY%FRACb0aP)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'FRACbN.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsSTICKY%FRACbN)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'FRACnw0.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsSTICKY%FRACnw0)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'FRACnw0close.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsSTICKY%FRACnw0close)

	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'EbN.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsSTICKY%EbN)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'EbP.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsSTICKY%EbP)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'Eadjcost.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsSTICKY%Eadjcost)

	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'GINIa.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsSTICKY%GINIa)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'GINIb.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsSTICKY%GINIb)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'GINInw.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsSTICKY%GINInw)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'GINIc.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsSTICKY%GINIc)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'GINIinc.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsSTICKY%GINIinc)

	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'Ec_bN.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsSTICKY%Ec_bN)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'Ec_b0close.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsSTICKY%Ec_b0close)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'Ec_b0far.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsSTICKY%Ec_b0far)

	DO ip = 1,11;lmat(:,ip) = irfsave%statsSTICKY%PERCa(ip); END DO	
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'PERCa.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,11,lmat)

	DO ip = 1,11;lmat(:,ip) = irfsave%statsSTICKY%PERCb(ip); END DO	
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'PERCb.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,11,lmat)

	DO ip = 1,11;lmat(:,ip) = irfsave%statsSTICKY%PERCnw(ip); END DO	
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'PERCnw.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,11,lmat)

	DO ip = 1,11;lmat(:,ip) = irfsave%statsSTICKY%PERCc(ip); END DO	
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'PERCc.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,11,lmat)

	DO ip = 1,11;lmat(:,ip) = irfsave%statsSTICKY%PERCinc(ip); END DO	
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'PERCinc.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,11,lmat)

	DO iq = 1,4;lmat4(:,iq) = irfsave%statsSTICKY%Ea_nwQ(iq); END DO	
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'Ea_nwQ.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,4,lmat4)

	DO iq = 1,4;lmat4(:,iq) = irfsave%statsSTICKY%Eb_nwQ(iq); END DO	
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'Eb_nwQ.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,4,lmat4)

	DO iq = 1,4;lmat4(:,iq) = irfsave%statsSTICKY%Ec_nwQ(iq); END DO	
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'Ec_nwQ.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,4,lmat4)

	DO iq = 1,4;lmat4(:,iq) = irfsave%statsSTICKY%Einc_nwQ(iq); END DO	
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'Einc_nwQ.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,4,lmat4)

	DO iq = 1,4;lmat4(:,iq) = irfsave%statsSTICKY%Ea_incQ(iq); END DO	
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'Ea_incQ.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,4,lmat4)

	DO iq = 1,4;lmat4(:,iq) = irfsave%statsSTICKY%Eb_incQ(iq); END DO	
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'Eb_incQ.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,4,lmat4)

	DO iq = 1,4;lmat4(:,iq) = irfsave%statsSTICKY%Ec_incQ(iq); END DO	
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'Ec_incQ.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,4,lmat4)

	DO iq = 1,4;lmat4(:,iq) = irfsave%statsSTICKY%Einc_incQ(iq); END DO	
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'Einc_incQ.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,4,lmat4)

	DO iq = 1,12;lmat12(:,iq) = irfsave%statsSTICKY%Ec_nwQ_add(iq); END DO	
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'Ec_nwQ_add.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,12,lmat12)


	IF(SaveTime1PolicyFns==1) THEN
		DO iy = 1,ngpy
			IF (iy<10) WRITE(UNIT=lstring, FMT='(I1)') iy
			IF (iy >= 10) WRITE(UNIT=lstring, FMT='(I2)') iy
			
			OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'V_T1_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,irfsave%solnSTICKY(1)%V(:,:,iy))
			OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/"  // 'dep_T1_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,irfsave%solnSTICKY(1)%d(:,:,iy))
			OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/"  // 'con_T1_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,irfsave%solnSTICKY(1)%c(:,:,iy))
			OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/"  // 'sav_T1_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,irfsave%solnSTICKY(1)%s(:,:,iy))
			OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/"  // 'hour_T1_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,irfsave%solnSTICKY(1)%h(:,:,iy))
			OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/"  // 'bdot_T1_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,irfsave%solnSTICKY(1)%bdot(:,:,iy))
			OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/"  // 'gjoint_T1_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,irfsave%solnSTICKY(1)%gjoint(:,:,iy))
			IF(ComputeDiscountedMPC==1) THEN
				OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/"  // 'mpc_T1_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,irfsave%solnSTICKY(1)%mpc(:,:,iy))
				OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/"  // 'subeff1ass_T1_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,irfsave%solnSTICKY(1)%subeff1ass(:,:,iy))
				OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/"  // 'subeff2ass_T1_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,irfsave%solnSTICKY(1)%subeff2ass(:,:,iy))
				OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/"  // 'wealtheff1ass_T1_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,irfsave%solnSTICKY(1)%wealtheff1ass(:,:,iy))
				OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/"  // 'wealtheff2ass_T1_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,irfsave%solnSTICKY(1)%wealtheff2ass(:,:,iy))
			END IF
		END DO
	END IF
	
	IF (SaveCumPolicyFnsIRF==1) THEN
		DO iy = 1,ngpy
			IF (iy<10) WRITE(UNIT=lstring, FMT='(I1)') iy
			IF (iy >= 10) WRITE(UNIT=lstring, FMT='(I2)') iy
			
			OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'ccum1_T1_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,irfsave%cumSTICKY%ccum1(:,:,iy))
			OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'ccum2_T1_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,irfsave%cumSTICKY%ccum2(:,:,iy))
			OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'ccum4_T1_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,irfsave%cumSTICKY%ccum4(:,:,iy))
			OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'dcum1_T1_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,irfsave%cumSTICKY%dcum1(:,:,iy))
			OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'dcum2_T1_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,irfsave%cumSTICKY%dcum2(:,:,iy))
			OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'dcum4_T1_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,irfsave%cumSTICKY%dcum4(:,:,iy))
		END DO
		
	END IF	
	
	!transitions steady state distributions and policy functions
	IF(SaveTransitionPolicyFns==1) THEN
		CALL system ("mkdir -p " // trim(OutputDirIRF) // "STICKY/FuncDist/")	
	
		DO it = 1,Ttransition
		IF (it<10) WRITE(UNIT=lt, FMT='(I1)') it
		IF (it>=10 .and. it<100) WRITE(UNIT=lt, FMT='(I2)') it
		IF (it>=100) WRITE(UNIT=lt, FMT='(I3)') it

		DO iy = 1,ngpy
			IF (iy<10) WRITE(UNIT=lstring, FMT='(I1)') iy
			IF (iy >= 10) WRITE(UNIT=lstring, FMT='(I2)') iy

			OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/FuncDist/" // 'V' // trim(lt) // '_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,irfsave%solnSTICKY(it)%V(:,:,iy))
			OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/FuncDist/" // 'dep' // trim(lt) // '_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,irfsave%solnSTICKY(it)%d(:,:,iy))
			OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/FuncDist/" // 'con' // trim(lt) // '_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,irfsave%solnSTICKY(it)%c(:,:,iy))
			OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/FuncDist/" // 'bdot' // trim(lt) // '_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,irfsave%solnSTICKY(it)%bdot(:,:,iy))
			OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/FuncDist/" // 'gjoint' // trim(lt) // '_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,irfsave%solnSTICKY(it)%gjoint(:,:,iy))

		END DO

		OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/FuncDist/" // 'gamarg' // trim(lt) // '.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpy,irfsave%solnSTICKY(it)%gamarg)
		OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/FuncDist/" // 'gbmarg' // trim(lt) // '.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpb,ngpy,irfsave%solnSTICKY(it)%gbmarg)

	END DO
	END IF
END IF

!sticky price transition (with ZLB)
IF(SolveZLBTransition==1) THEN
	CALL system ("mkdir -p " // trim(OutputDirIRF) // "ZLB")	

	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'ra.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmZLB%ra)
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'rcapital.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmZLB%rcapital)
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'wage.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmZLB%wage)
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'netwage.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmZLB%netwage)
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'KYratio.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmZLB%KYratio)
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'KNratio.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmZLB%KNratio)
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'mc.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmZLB%mc)
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'rb.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmZLB%rb)
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'rborr.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmZLB%rborr)
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'tfp.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmZLB%tfp)
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'pi.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmZLB%pi)
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'rnom.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmZLB%rnom)
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'gap.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmZLB%gap)
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'capital.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmZLB%capital)
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'bond.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmZLB%bond)
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'labor.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmZLB%labor)
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'output.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmZLB%output)
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'investment.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmZLB%investment)
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'taxrev.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmZLB%taxrev)
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'govexp.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmZLB%govexp)
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'govbond.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmZLB%govbond)
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'worldbond.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmZLB%worldbond)
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'borrwedge.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmZLB%borrwedge)
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'rho.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmZLB%rho)
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'kappa0_w.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmZLB%kappa0_w)
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'kappa1_w.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmZLB%kappa1_w)
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'mpshock.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmZLB%mpshock)
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'prefshock.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmZLB%prefshock)
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'elast.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmZLB%elast)
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'fundlev.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmZLB%fundlev)
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'gam.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmZLB%gam)
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'fundbond.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmZLB%fundbond)
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'priceadjust.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmZLB%priceadjust)
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'profit.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmZLB%profit)
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'dividend.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmZLB%dividend)
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'divrate.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmZLB%divrate)
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'lumptransfer.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmZLB%lumptransfer)
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'equity.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmZLB%equity)
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'caputil.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmZLB%caputil)
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'deprec.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmZLB%deprec)
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'tfpadj.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmZLB%tfpadj)
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'illassetdrop.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmZLB%illassetdrop)
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'govshock.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmZLB%govshock)
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'transfershock.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmZLB%transfershock)
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'finwedge.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmZLB%finwedge)
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'labwedge.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmZLB%labwedge)
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'pricelev.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmZLB%pricelev)
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'prodgridscale.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmZLB%prodgridscale)
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'prodmarkovscale.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmZLB%prodmarkovscale)

	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'Ea.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsZLB%Ea)
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'Eb.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsZLB%Eb)
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'Ec.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsZLB%Ec)
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'Elabor.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsZLB%Elabor)
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'Ed.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsZLB%Ed)
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'Ewage.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsZLB%Ewage)
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'Enetlabinc.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsZLB%Enetlabinc)
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'Egrosslabinc.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsZLB%Egrosslabinc)
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'Einc.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsZLB%Einc)
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'Ehours.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsZLB%Ehours)
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'Enw.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsZLB%Enw)

	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'FRACa0.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsZLB%FRACa0)
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'FRACa0close.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsZLB%FRACa0close)
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'FRACb0.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsZLB%FRACb0)
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'FRACb0close.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsZLB%FRACb0close)
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'FRACb0a0.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsZLB%FRACb0a0)
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'FRACb0aP.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsZLB%FRACb0aP)
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'FRACbN.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsZLB%FRACbN)
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'FRACnw0.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsZLB%FRACnw0)
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'FRACnw0close.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsZLB%FRACnw0close)

	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'EbN.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsZLB%EbN)
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'EbP.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsZLB%EbP)
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'Eadjcost.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsZLB%Eadjcost)

	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'GINIa.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsZLB%GINIa)
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'GINIb.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsZLB%GINIb)
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'GINInw.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsZLB%GINInw)
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'GINIc.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsZLB%GINIc)
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'GINIinc.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsZLB%GINIinc)

	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'Ec_bN.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsZLB%Ec_bN)
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'Ec_b0close.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsZLB%Ec_b0close)
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'Ec_b0far.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsZLB%Ec_b0far)

	DO ip = 1,11;lmat(:,ip) = irfsave%statsZLB%PERCa(ip); END DO	
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'PERCa.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,11,lmat)

	DO ip = 1,11;lmat(:,ip) = irfsave%statsZLB%PERCb(ip); END DO	
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'PERCb.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,11,lmat)

	DO ip = 1,11;lmat(:,ip) = irfsave%statsZLB%PERCnw(ip); END DO	
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'PERCnw.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,11,lmat)

	DO ip = 1,11;lmat(:,ip) = irfsave%statsZLB%PERCc(ip); END DO	
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'PERCc.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,11,lmat)

	DO ip = 1,11;lmat(:,ip) = irfsave%statsZLB%PERCinc(ip); END DO	
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'PERCinc.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,11,lmat)

	DO iq = 1,4;lmat4(:,iq) = irfsave%statsZLB%Ea_nwQ(iq); END DO	
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'Ea_nwQ.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,4,lmat4)

	DO iq = 1,4;lmat4(:,iq) = irfsave%statsZLB%Eb_nwQ(iq); END DO	
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'Eb_nwQ.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,4,lmat4)

	DO iq = 1,4;lmat4(:,iq) = irfsave%statsZLB%Ec_nwQ(iq); END DO	
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'Ec_nwQ.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,4,lmat4)

	DO iq = 1,4;lmat4(:,iq) = irfsave%statsZLB%Einc_nwQ(iq); END DO	
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'Einc_nwQ.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,4,lmat4)

	DO iq = 1,4;lmat4(:,iq) = irfsave%statsZLB%Ea_incQ(iq); END DO	
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'Ea_incQ.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,4,lmat4)

	DO iq = 1,4;lmat4(:,iq) = irfsave%statsZLB%Eb_incQ(iq); END DO	
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'Eb_incQ.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,4,lmat4)

	DO iq = 1,4;lmat4(:,iq) = irfsave%statsZLB%Ec_incQ(iq); END DO	
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'Ec_incQ.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,4,lmat4)

	DO iq = 1,4;lmat4(:,iq) = irfsave%statsZLB%Einc_incQ(iq); END DO	
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'Einc_incQ.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,4,lmat4)

	DO iq = 1,12;lmat12(:,iq) = irfsave%statsZLB%Ec_nwQ_add(iq); END DO	
	OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'Ec_nwQ_add.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,12,lmat12)

	IF(SaveTime1PolicyFns==1) THEN
		DO iy = 1,ngpy
			IF (iy<10) WRITE(UNIT=lstring, FMT='(I1)') iy
			IF (iy >= 10) WRITE(UNIT=lstring, FMT='(I2)') iy
			
			OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'V_T1_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,irfsave%solnZLB(1)%V(:,:,iy))
			OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/"  // 'dep_T1_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,irfsave%solnZLB(1)%d(:,:,iy))
			OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/"  // 'con_T1_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,irfsave%solnZLB(1)%c(:,:,iy))
			OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/"  // 'sav_T1_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,irfsave%solnZLB(1)%s(:,:,iy))
			OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/"  // 'hour_T1_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,irfsave%solnZLB(1)%h(:,:,iy))
			OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/"  // 'bdot_T1_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,irfsave%solnZLB(1)%bdot(:,:,iy))
			OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/"  // 'gjoint_T1_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,irfsave%solnZLB(1)%gjoint(:,:,iy))
		END DO
	END IF
	
	IF (SaveCumPolicyFnsIRF==1) THEN
		DO iy = 1,ngpy
			IF (iy<10) WRITE(UNIT=lstring, FMT='(I1)') iy
			IF (iy >= 10) WRITE(UNIT=lstring, FMT='(I2)') iy
			
			OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'ccum1_T1_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,irfsave%cumZLB%ccum1(:,:,iy))
			OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'ccum2_T1_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,irfsave%cumZLB%ccum2(:,:,iy))
			OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'ccum4_T1_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,irfsave%cumZLB%ccum4(:,:,iy))
			OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'dcum1_T1_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,irfsave%cumZLB%dcum1(:,:,iy))
			OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'dcum2_T1_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,irfsave%cumZLB%dcum2(:,:,iy))
			OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/" // 'dcum4_T1_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,irfsave%cumZLB%dcum4(:,:,iy))
		END DO
		
	END IF	
	
	!transitions steady state distributions and policy functions
	IF(SaveTransitionPolicyFns==1) THEN
		CALL system ("mkdir -p " // trim(OutputDirIRF) // "ZLB/FuncDist/")	
	
		DO it = 1,Ttransition
			IF (it<10) WRITE(UNIT=lt, FMT='(I1)') it
			IF (it>=10 .and. it<100) WRITE(UNIT=lt, FMT='(I2)') it
			IF (it>=100) WRITE(UNIT=lt, FMT='(I3)') it

			DO iy = 1,ngpy
				IF (iy<10) WRITE(UNIT=lstring, FMT='(I1)') iy
				IF (iy >= 10) WRITE(UNIT=lstring, FMT='(I2)') iy

				OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/FuncDist/" // 'V' // trim(lt) // '_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,irfsave%solnZLB(it)%V(:,:,iy))
				OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/FuncDist/" // 'dep' // trim(lt) // '_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,irfsave%solnZLB(it)%d(:,:,iy))
				OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/FuncDist/" // 'con' // trim(lt) // '_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,irfsave%solnZLB(it)%c(:,:,iy))
				OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/FuncDist/" // 'bdot' // trim(lt) // '_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,irfsave%solnZLB(it)%bdot(:,:,iy))
				OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/FuncDist/" // 'gjoint' // trim(lt) // '_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,irfsave%solnZLB(it)%gjoint(:,:,iy))

			END DO

			OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/FuncDist/" // 'gamarg' // trim(lt) // '.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpy,irfsave%solnZLB(it)%gamarg)
			OPEN(3, FILE = trim(OutputDirIRF) // "ZLB/FuncDist/" // 'gbmarg' // trim(lt) // '.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpb,ngpy,irfsave%solnZLB(it)%gbmarg)

		END DO
	END IF
END IF


END SUBROUTINE