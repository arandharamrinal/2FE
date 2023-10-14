# 2FE
Computes energy and force for gas-phase 2-fluoroethanol
=======================================================
Usage in an MD code
-------------------
To utilize the PES in an MD code, follow these steps:
  1. Compile *UnitConversion.f90*, *sysvariables.f90*, *potvariables.f90*, *getEnForce2FE.f90*, *init2fe.f90* routines in this particular order.
  2.  Next, initialize the modules. Do this by invoking the *init2fe()* subroutine.
  3. Integrate the modules into the subroutine responsible for calculating Forces and Energies. To achieve this, include the line *use FEforceEMod* before specifying the variables.
  4. To compute the force, use the *getForce(X, F, V)* subroutine. Where X represents the input Cartesian geometry in Angstrom (a 1x3N array), F denotes the force array with a dimension of 1x3N, and V (a scalar) stands for the potential. N is the total number of atoms (N=9 for 2FE).
  5. If your objective is to compute only the energy, utilize the *getPot(X, V, Vrs)* subroutine. Here, X and V are as in step 3. Additionally, Vrs denotes the potential at the reaction surface.

Cartesian coordinate format
---------------------------
The input geometries should follow the atom number convention as found in the *2fe_ref.zmat* file.

PATH to parameter files
-----------------------
*PotFitParams2FE* contains all the fitted parameter files. By default, it must be in the same directory where the 
code is executed. This can be changed manually by editing the *potvariables.f90* file.

Units
-----
All input geometries must be in Angstrom. Output energies and forces are in $cm^{-1}$ and Hartee/Angstrom, respectively.

Example
=======
Let's try computing energies and forces for already-known geometries to see how it works.
Goal: Read cartesian geometries from an input file and write the corresponding energies and forces to an outputfile 
Input and output file names should be specified in *input.param* file. 

*input.param* format
--------------------
First and second line of *input.param* contains the input file and output file names, respectively.

Compilation
-----------
Before running the code, compile it using the provided *compile.sh* script. The code is compatible with any Fortran compiler, and the default compiler used is *gfortran*. Ensure the *PotFitParams2FE* folder containing the fitted parameter files is in the working directory.

Usage
-----
Run the compiled executable *getfe.exe*. The code will process the geometries, calculate energies and forces, and output the results to the output file. Each line in the output file will contain the energy followed by the forces for the corresponding geometry from the input file.

Input File Format
-----------------
The input file contains the Cartesian coordinates of each geometry written in a single line. 

Output Format
-------------
The output energies will be in units of $cm^{-1}$, and the forces will be in Hartree/Angstrom. Each line in the output file will follow this format:
>  Energy, F<sub>x</sub>, F<sub>y</sub>, F<sub>z</sub> of atom 1, F<sub>x</sub>, F<sub>y</sub>, F<sub>z</sub> of atom 2 ....

Feel free to reach out if you have any questions or need further assistance.
