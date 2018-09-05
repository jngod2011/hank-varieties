NAME = $(OUT)
#SRCDIR="${HOME}/git/hank-main"
SRCDIR="${HOME}/hank-main"
RUNDIR="${SCRATCH}/hank_runs/$(NAME)"
WORK="/project/gkaplan/hank_output"

.PHONY: $(NAME).pbs

$(NAME).sbatch:

	#create run directory on SCRATCH	
	rm -f $(NAME).sbatch
	mkdir -p $(RUNDIR)
	cp $(SRCDIR)/$(NAME).out $(RUNDIR)/
	cp $(SRCDIR)/SetParameters.f90 $(RUNDIR)/
	cp $(SRCDIR)/Parameters.f90 $(RUNDIR)/

	#make sbatch file
	echo "#! /bin/bash" >> $(NAME).sbatch

	# set up slurm environement
	echo "#SBATCH --nodes=1" >> $(NAME).sbatch
	echo "#SBATCH --partition=broadwl" >> $(NAME).sbatch
	echo "#SBATCH --exclusive" >> $(NAME).sbatch
	echo "#SBATCH --time=12:00:00"  >> $(NAME).sbatch
	echo "#SBATCH --mail-type BEGIN,END" >> $(NAME).sbatch

	#run model and zip results
	echo "cd $(RUNDIR)" >> $(NAME).sbatch
	echo "ulimit -s unlimited" >> $(NAME).sbatch
	echo "./$(NAME).out > output_$(NAME).txt" >> $(NAME).sbatch
	echo "cp -R Output/* ./" >> $(NAME).sbatch
	echo "rm -r Output" >> $(NAME).sbatch

	#copy output to directory on WORK
	echo "cd .." >> $(NAME).sbatch
	echo "zip -r results_$(NAME).zip $(NAME)" >> $(NAME).sbatch
	echo "mkdir -p $(WORK)" >> $(NAME).sbatch
	echo "cp results_$(NAME).zip $(WORK)/results_$(NAME).zip" >> $(NAME).sbatch
