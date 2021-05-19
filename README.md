* Compiling DQMC (determinantal QMC, essentially the same as AFQMC):

Get the code using 'git clone https://github.com/sanshar/VMC.git'. Checkout the dqmc branch using 'git checkout dqmc'. 
The Makefile on github uses some dependencies required by other parts of the code that are not required for DQMC, so you can use the Makefile in this directory (/projects/anma2640/VMC/test/VMC) instead. Provide paths for boost, eigen, and hdf5 in the Makefile on the following lines:

EIGEN=/projects/ilsa8974/apps/eigen/
BOOST=/projects/anma2640/boost_1_66_0/
HDF5=/curc/sw/hdf5/1.10.1/impi/17.3/intel/17.4/

Eigen is header-only so the above path should work since we all have access to Iliya's directory. The hdf5 path should also work as is after loading the hdf5 module as 'module load hdf5' (add to your ~/.bashrc, because it will be needed runtime) on the cluster. For boost you can copy the path from your Dice Makefile. To compile use 'make bin/DQMC -j'. 


* Running examples:

There are two small examples, n2 and h2o, in the /projects/anma2640/VMC/test/VMC/examples/DQMC directory. All required input files are already in the directories, so you can just run the code to check if it works right away. To run the code:

mpirun -np 4 your_path_to_the_code/bin/DQMC dqmc.json > dqmc.out


* DQMC input details:

The code requires two things from external sources: system Hamiltonian information and trial wave function information. Following are details about how different things are specified in the input file (dqmc.json) and how the corresponding data files are generated from other programs.

AFQMC uses integrals in the Cholesky format and the code reads them from the 'FCIDUMP_chol' file specified in the "system" block.  

Details of the trial wave functions are provided in the "wavefunction" block. For most cases we use "right": "ccsd" and "left": "multislater". Information about the ccsd wave function is read from the 'ccsd.h5' file and that about the multi-Slater wave function is read from the 'dets.bin' file specified in the input.

The "sampling" block contains details of the AFQMC sampling:

* "seed" is the number provided to the random number generator, so if you use the same seed the results should be exactly the same. 
* "dt" is the propagation time-step.
* "nsteps" is the number of steps. So nsteps * dt is the total propagation time. 
* "eneSteps" specifies the steps at which energy is measured. 
* "errorTargets" specifies the stochastic error threshold, so the program stops measuring the energy at the corresponding eneSteps when the noise goes below the threshold or if the number of iterations exceeds "stochasticIter". 
* "ene0Guess" is the guess for the ground state energy, one can use the CCSD(T) energy for this.

Other options have to do with small details.


* Generating data files:

We use python scripts to generate most of the data files (n2.py and h2o.py). A collection of useful functions are in the scripts/prepVMC.py file. To generate the integrals in the Cholesky format, a small part of another AFQMC code, PAUXY, is used. You can get it using 'git clone https://github.com/pauxy-qmc/pauxy.git'. Add both these to your python path using 'export PYTHONPATH=/projects/anma2640/VMC/master/VMC/scripts/:/projects/anma2640/pauxy/:$PYTHONPATH' (change to your paths if you
want and add to bashrc). By doing this these can be imported to your scripts (like 'import prepVMC'). 

If you look at the n2.py and h2o.py scripts they generate and write to disk:

* the cholesky integrals ('FCIDUMP_chol' in hdf5 format)
* rhf mo's ('rhf.txt' as a text file)
* ccsd amplitudes ('ccsd.h5' in hdf5 format)
* conventional eri integrals for canonical orbitals ('FCIDUMP_can' as a text file used by Dice, not used drectly by the afqmc program)

To generate the 'dets.bin' file that contains the multi-Slater / HCI / selected-CI wave function we do a Dice calculation. I added a small piece of code to Dice to dump the determinants in a binary format. I have not pushed it to the Dice repo, so you will have to copy it from /projects/anma2640/newDice/Dice. Compiling it is identical to the Dice you use. The 'writebestdeterminants #' line in the 'dice.dat' input file specifies how many determinants you want it to dump. 
