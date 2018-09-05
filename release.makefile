FC = ifort
FCFLAGS = -m64 -traceback -O3 -openmp -implicitnone  -Wl,-stack_size,0x100000000 -L/usr/local/lib -lumfpack -lamd -lcholmod -lcolamd -lsuitesparseconfig -lblas 
LDFLAFS = -m64 -traceback -O3 -openmp -implicitnone  -Wl,-stack_size,0x100000000 -L/usr/local/lib -lumfpack -lamd -lcholmod -lcolamd -lsuitesparseconfig -lblas 


PROG = $(OUT)

MOD = Parameters.o Globals.o umfpack.o Procedures.o 

SUBR = 	AllocateArrays.o SetParameters.o Grids.o IterateBellman.o HJBUpdate.o cumnor.o rtsec.o StationaryDistribution.o SaveSteadyStateOutput.o DistributionStatistics.o rtbis.o rtflsp.o Transition.o InitialSteadyState.o FinalSteadyState.o IterateTransitionFlex.o IterateTransitionStickyB.o IterateTransitionStickyPi.o IterateTransitionStickyY.o IterateTransitionStickyRb.o SolveSteadyStateEqum.o Calibration.o MomentConditions.o dfovec.o newuoa-h.o newuob-h.o update.o trsapp-h.o biglag.o bigden.o mnbrak.o golden.o sort2.o FnGridAR1.o ImpulseResponses.o SaveIRFOutput.o SetupFiscalStimulus.o CumulativeConsumption.o IRFSequence.o FnDiscountRate.o Exploration.o sobol.o OptimalConsumption.o FnOptimalCon0Rent.o


OBJ = $(MOD) $(SUBR)

$(PROG).out: $(OBJ) Main.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)
Main.o: $(MOD)

%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.f90
	$(FC) $(FCFLAGS) -c $<
