#source /opt/intel/bin/compilervars.sh intel64

release:
	make -f release.makefile

debug:
	make -f debug.makefile
	dsymutil $(OUT).out

sstest:
	make -f sstest.makefile
	# dsymutil $(OUT).out

hpcrelease:
	make -f hpcrelease.makefile

midwayrelease:
	make -f midwayrelease.makefile

pbs:
	make -f pbs.makefile

compilesubmit:
	make -f midwayrelease.makefile
	make -f slurm.makefile
	sbatch $(OUT).sbatch

compilesubmitNYU:
	make -f hpcrelease.makefile
	make -f pbs.makefile
	qsub $(OUT).pbs

.PHONY: clean

clean:
	rm -f *.o *.mod *.MOD *genmod* *~ 
	rm -fr *.dSYM
