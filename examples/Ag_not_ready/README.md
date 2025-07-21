1. There are some basic input files one has to prepare when launching Quantum Espresso EPW calculations.
The details about the input parameters may be found e.g. at

https://www.quantum-espresso.org/Doc/INPUT_PW.html,
https://www.quantum-espresso.org/Doc/INPUT_PH.html,
https://docs.epw-code.org/doc/Inputs.html.

2. Relaxation step

The relax.sh script is to be configured by a user according to
the specific system parameters and the path of the parallel version of Quantum Espresso binaries. Here
only an example of batch stript is given:

```bash
sbatch relax.sh
```

Run the relax.sh script to obtain the relaxed fcc Ag lattice parameter to be used as the celldm(1) parameter
in the input files. The corresponding unit cell volume may be found in the relax.out file:

```bash
grep 'new unit-cell volume' relax.out
```

One has to choose the last (converged) value. The relax.sh script is to be configured by a user according to
the specific system parameters and the path of the parallel version of Quantum Espresso binaries.

3. Quantum Espresso calculations

Next, see the provided example of run_qe.sh file. It consists of 5 consecutive steps as described in the

https://docs.epw-code.org/doc/School2022.html (see exercise 2 of Wednesday Hands-On tutorial 1 prepared by S. Ponce).

The four first steps strongly follow the given manual, and the final step devoted to EPW calculations is a bit different.
The input file epw.in consists of minimal amount of lines which are enough to obtain all the necessary output files
used at the next stage in SKiES calculations. Please note special attention to the parameters connected with the
Wannier interpolation. The details might be found at 

https://docs.epw-code.org/doc/Inputs.html,
https://wannier.org/support/

After launching the following command

```bash
sbatch run_qe.sh
```

The wait until the end of all the stages. The output file epw.out should interrupt with the message

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Error in routine loadqmesh_serial (1):
     Cannot load fine q points
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     stopping ...

which is ok. After that one may continue with the SKiES part of calculations. See the example provided at the folder Ag_ready in the examples.
