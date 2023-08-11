# 2FE
Computes energy and force for gas-phase 2-fluoroethanol
=======================================================
This code reads Cartesian geometries from an input file and extracts energies and forces for each geometry. The input geometries should follow the atom number convention as found in the "2fe_ref.zmat" file, with each geometry written on a single line.

Compilation
-----------
Before running the code, compile it using the provided compile.sh script. The code is compatible with any Fortran compiler, and the default compiler used is gfortran. Make sure the PotFitParams2FE folder, containing the fitted parameter files, is present in the working directory.

Usage
-----
Run the compiled executable *getfe.exe*. Input and output file names should be specified in *input.param* file. The code will process the geometries, calculate energies and forces, and output the results to the output file. Each line in the output file will contain the energy followed by the forces for the corresponding geometry from the input file.

*input.param* format
--------------------
First and second line of *input.param* contains the input file and output file names, respectively.

Input Geometry Format
---------------------
Geometries in the input file should be provided in Angstrom units and adhere to the atom number convention as specified in the "2fe_ref.zmat" file.


Output Format
-------------
The output energies will be in units of $cm^{-1}$, and the forces will be in Hartree/Angstrom. Each line in the output file will follow this format:
> Energy, $F_x$, $F_y$, $F_z$ of atom 1, $F_x$, $F_y$, $F_z$ of atom 2 ....

Feel free to reach out if you have any questions or need further assistance.
