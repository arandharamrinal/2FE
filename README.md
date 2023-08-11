# 2FE
Computes energy and force for gas-phase 2-fluoroethanol
=======================================================
Reads Cartesian geometries from a file and writes the energies and forces to an output file.
Each geometry should be in a single line.
Geometries must be in Angstrom and must follow the atom number convention as in "2fe_ref.zmat" file.
Output energies and forces are in $cm^{-1}$ and Hartee/\AA. 
Each line in the output file contains the energy followed by forces of corresponding geometries in the input file. 
> Input and output file names should be specified in "pot.param" file.
> "pot.param" has the following structure.
>> "inputfilename"
>> "outputfilename"
