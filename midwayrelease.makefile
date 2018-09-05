FC = ifort

FCFLAGS = -m64 -traceback -O2 -qopenmp -implicitnone -xSSE4.2 -axAVX -lumfpack -lamd -lcholmod -lcolamd -lsuitesparseconfig -lblas
LDFLAFS = -m64 -traceback -O2 -qopenmp -implicitnone -xSSE4.2 -axAVX -lumfpack -lamd -lcholmod -lcolamd -lsuitesparseconfig -lblas


PROG = $(OUT)

MOD = Parameters.o Globals.o umfpack.o Procedures.o

SUBR =  AllocateArrays.o SetParameters.o Grids.o IterateBellman.o HJBUpdate.o cumnor.o rtsec.o StationaryDistribution.o SaveSteadyStateOutput.o DistributionStatistics.o rtbis.o rtflsp.o InitialSteadyState.o FinalSteadyState.o SolveSteadyStateEqum.o Calibration.o MomentConditions.o dfovec.o newuoa-h.o newuob-h.o update.o trsapp-h.o biglag.o bigden.o mnbrak.o golden.o sort2.o  CumulativeConsumption.o  FnDiscountRate.o  OptimalConsumption.o FnHoursBC.o  ImpulseResponses.o IRFSequence.o Transition.o  SaveIRFOutput.o IterateTransOneAssetStickyRb.o IterateTransition.o Exploration.o sobol.o FnCapitalEquity.o CumulativeConsTransition.o DiscountedMPC.o DiscountedMPCTransition.o FeedInPrices.o




OBJ = $(MOD) $(SUBR)

$(PROG).out: $(OBJ) Main.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)
Main.o: $(MOD)

%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.f90
	$(FC) $(FCFLAGS) -c $<
