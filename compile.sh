#!/bin/bash
#ifort -r8 -O4 UnitConversion.f90 sysvariables.f90 potvariables.f90 getEnForce2FE.f90 computeEnForce.f90 -o getfe.exe
gfortran -k8 -O4 -ffree-line-length-1024 UnitConversion.f90 sysvariables.f90 potvariables.f90 getEnForce2FE.f90 init2fe.f90 computeEnForce.f90 -o getfe.exe
