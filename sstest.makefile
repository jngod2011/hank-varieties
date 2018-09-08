FC = ifort
 FCFLAGS = -m64 -g -debug all -implicitnone -Wl,-stack_size,0x100000000 -save-temps -warn all -fp-stack-check -ftrapuv -traceback -L/Volumes/FILES/Projects/Fortran/SuiteSparse/lib -lumfpack -lamd -lcholmod -lcolamd -lsuitesparseconfig -lblas -check all
 LDFLAFS = -m64 -g -Wl -debug all,-stack_size,0x100000000  -L/Volumes/FILES/Projects/Fortran/SuiteSparse/lib -lumfpack -lamd -lcholmod -lcolamd -lsuitesparseconfig -lblas -check all

#FCFLAGS = -m64 -traceback -O3 -qopenmp -implicitnone  -Wl,-stack_size,0x100000000 -L/Volumes/FILES/Projects/Fortran/SuiteSparse/lib -lumfpack -lamd -lcholmod -lcolamd -lsuitesparseconfig -lblas
#LDFLAFS = -m64 -traceback -O3 -qopenmp -implicitnone  -Wl,-stack_size,0x100000000 -L/Volumes/FILES/Projects/Fortran/SuiteSparse/lib -lumfpack -lamd -lcholmod -lcolamd -lsuitesparseconfig -lblas

#-check noarg_temp_created -check arg_temp_created -gen-interfaces -warn interfaces

PROG = $(OUT)

MOD = Parameters.o Globals.o umfpack.o Procedures.o 

SUBR = 	AllocateArrays.o SetParameters.o Grids.o InitialPrices.o IterateBellman.o HJBUpdate.o  OptimalConsumption.o FnHoursBC.o cumnor.o rtsec.o sort2.o StationaryDistribution.o SaveSteadyStateOutput.o DistributionStatistics.o rtbis.o rtflsp.o InitialSteadyState.o FinalSteadyState.o SolveSteadyStateEqum.o CumulativeConsumption.o DiscountedMPC.o FnDiscountRate.o ImpulseResponses.o IRFSequence.o IterateTransition.o Transition.o SaveIRFOutput.o 

#Calibration.o MomentConditions.o dfovec.o newuoa-h.o newuob-h.o update.o trsapp-h.o biglag.o bigden.o


OBJ = $(MOD) $(SUBR)

$(PROG).out: $(OBJ) Main.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)
Main.o: $(MOD)

%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.f90
	$(FC) $(FCFLAGS) -c $<

