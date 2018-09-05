SUBROUTINE FeedInPrices

USE Parameters
USE Globals
USE Procedures

IMPLICIT NONE
INTEGER		:: it,ipf,il,iv
CHARACTER(len=300)	:: skip
REAL(8), DIMENSION(Ttransition) 	:: wage_logdev, lumptransfer_logdev, rb_dev, ra_dev, profit_logdev, pref_dev
REAL(8) 	:: lillassetdrop


!allocate arrays
CALL AllocateArrays


!set up deltatransvec
IF (TransitionTimeStepType==1) THEN
	CALL PowerSpacedGrid (Ttransition,deltatransparam,deltatransmin,deltatransmax,deltatransvec)

	cumdeltatrans(1) = deltatransvec(1)
	DO it = 2,Ttransition
		cumdeltatrans(it) = cumdeltatrans(it-1) + deltatransvec(it)
	END DO
ELSEIF (TransitionTimeStepType==2) THEN
	cumdeltatrans(1) = 1.0/real(Ttransition+1) 
	DO it = 2,Ttransition
		cumdeltatrans(it) = cumdeltatrans(it-1) + 1.0/real(Ttransition+1) 
	END DO
	
	cumdeltatrans = (cumdeltatrans/(1.0-cumdeltatrans)) / deltatransnu
	deltatransvec(1) = cumdeltatrans(1)
	deltatransvec(2:Ttransition) = cumdeltatrans(2:Ttransition) - cumdeltatrans(1:Ttransition-1)

END IF
OPEN(3, FILE = trim(OutputDir) // 'deltatransvec.txt', STATUS = 'replace'); CALL WriteMatrix(3,Ttransition,1,deltatransvec)

!read in price sequeneces: must be correct length


DO ipf = 1,nfeed
	
	! get directories
	OPEN(1, FILE = trim(FeedPriceFile))
	IF (ipf>1) THEN
		DO il=1, ipf-1
		   READ (1,'(A)')  skip
		END DO
	END IF
	READ(1,'(A)') FeedPriceDir
	CLOSE(1)
	FeedPriceDir = trim(FeedPriceDir)
	WRITE(*,*) "Solving for partial equm response to fed in prices, directory: "
	WRITE(*,*) "  :",FeedPriceDir
	WRITE(*,*) " "
	
	
	
	!wage: log deviation
	OPEN(1, FILE = trim(FeedPriceDir) // '/wage_logdev.txt');READ(1,*) wage_logdev;CLOSE(1)

	!lumptransfer: log deviation
	OPEN(1, FILE = trim(FeedPriceDir) // '/lumptransfer_logdev.txt');READ(1,*) lumptransfer_logdev;CLOSE(1)

	!rb: deviation
	OPEN(1, FILE = trim(FeedPriceDir) // '/rb_dev.txt');READ(1,*) rb_dev;CLOSE(1)

	!ra: deviation
	OPEN(1, FILE = trim(FeedPriceDir) // '/ra_dev.txt');READ(1,*) ra_dev;CLOSE(1)

	!profit: log deviation
	OPEN(1, FILE = trim(FeedPriceDir) // '/profit_logdev.txt');READ(1,*) profit_logdev;CLOSE(1)

	!illiquid asset drop: fraction
	OPEN(1, FILE = trim(FeedPriceDir) // '/illassetdrop.txt');READ(1,*) lillassetdrop;CLOSE(1)

	!preference shock: fraction
	OPEN(1, FILE = trim(FeedPriceDir) // '/pref_dev.txt');READ(1,*) pref_dev;CLOSE(1)

	DO iv = 1,4
		equmTRANS(:) = equmINITSS	
		
		SELECT CASE(iv)

			CASE(1) !change all
				equmTRANS%wage = equmINITSS%wage*exp(wage_logdev)
				equmTRANS%lumptransfer = equmINITSS%lumptransfer*exp(lumptransfer_logdev)
				equmTRANS%rb = equmINITSS%rb + rb_dev
				equmTRANS%ra = equmINITSS%ra + ra_dev
				equmTRANS%profit = equmINITSS%profit*exp(profit_logdev)
				equmTRANS%illassetdrop  = lillassetdrop
				equmTRANS%prefshock  = pref_dev
				equmTRANS%netwage = (1.0-equmTRANS%labtax)*equmTRANS%wage
		
				OutputDirIRF = trim(FeedPriceDir) // '/all'
				CALL system ("mkdir -p " // trim(OutputDirIRF))
				
			CASE(2) !change rb
				equmTRANS%rb = equmINITSS%rb + rb_dev

				OutputDirIRF = trim(FeedPriceDir) // '/rb'
				CALL system ("mkdir -p " // trim(OutputDirIRF))

			CASE(3) !change wage, profit
				equmTRANS%wage = equmINITSS%wage*exp(wage_logdev)
				equmTRANS%profit = equmINITSS%profit*exp(profit_logdev)
				equmTRANS%netwage = (1.0-equmTRANS%labtax)*equmTRANS%wage
				
				OutputDirIRF = trim(FeedPriceDir) // '/wageprof'
				CALL system ("mkdir -p " // trim(OutputDirIRF))

			CASE(4) !change wage, profit, transfer
				equmTRANS%wage = equmINITSS%wage*exp(wage_logdev)
				equmTRANS%lumptransfer = equmINITSS%lumptransfer*exp(lumptransfer_logdev)				
				equmTRANS%profit = equmINITSS%profit*exp(profit_logdev)
				equmTRANS%netwage = (1.0-equmTRANS%labtax)*equmTRANS%wage
			
				OutputDirIRF = trim(FeedPriceDir) // '/wagetransfer'
				CALL system ("mkdir -p " // trim(OutputDirIRF))
				

		END SELECT

		CALL Transition


		! save output
		OPEN(3, FILE = trim(OutputDirIRF) // '/wage.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,equmTRANS%wage)
		OPEN(3, FILE = trim(OutputDirIRF) // '/lumptransfer.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,equmTRANS%lumptransfer)
		OPEN(3, FILE = trim(OutputDirIRF) // '/rb.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,equmTRANS%rb)
		OPEN(3, FILE = trim(OutputDirIRF) // '/ra.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,equmTRANS%ra)
		OPEN(3, FILE = trim(OutputDirIRF) // '/profit.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,equmTRANS%profit)
		OPEN(3, FILE = trim(OutputDirIRF) // '/prefshock.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,equmTRANS%prefshock)

		OPEN(3, FILE = trim(OutputDirIRF) // '/Ea.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,statsTRANS%Ea)
		OPEN(3, FILE = trim(OutputDirIRF) // '/Eb.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,statsTRANS%Eb)
		OPEN(3, FILE = trim(OutputDirIRF) // '/Ec.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,statsTRANS%Ec)
		OPEN(3, FILE = trim(OutputDirIRF) // '/Elabor.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,statsTRANS%Elabor)
		OPEN(3, FILE = trim(OutputDirIRF) // '/Ed.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,statsTRANS%Ed)
		OPEN(3, FILE = trim(OutputDirIRF) // '/Ewage.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,statsTRANS%Ewage)
		OPEN(3, FILE = trim(OutputDirIRF) // '/Enetlabinc.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,statsTRANS%Enetlabinc)
		OPEN(3, FILE = trim(OutputDirIRF) // '/Egrosslabinc.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,statsTRANS%Egrosslabinc)
		OPEN(3, FILE = trim(OutputDirIRF) // '/Einc.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,statsTRANS%Einc)
		OPEN(3, FILE = trim(OutputDirIRF) // '/Ehours.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,statsTRANS%Ehours)
		OPEN(3, FILE = trim(OutputDirIRF) // '/Enw.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,statsTRANS%Enw)
	END DO
END DO		

END SUBROUTINE FeedInPrices