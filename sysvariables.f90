module molparams 
use constants
implicit none
integer::Natoms,Nbonds,Nang,Ndih,Ncarts,N_int
contains

subroutine initializeMolParams()
implicit none
integer::i,j,k
!Read param file 
Natoms = 9
Ncarts = 3*Natoms
N_int= Ncarts-6
Nbonds = Natoms-1
Nang = Natoms-2
Ndih = Natoms-3
end subroutine 

end module 
