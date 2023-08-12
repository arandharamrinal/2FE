program ComputeEnF
use FEforceEMod
implicit none 
character(len=200):: cartCoordFile,EnForceOutFile 
real(8)::Vpot 
real(8),allocatable :: Xcarts(:),Force(:) 
integer::i,j,iost,linectr

!Initialize 2fe modules
call init2fe()

!Allocate arrays
allocate(Xcarts(0:Ncarts-1),Force(0:Ncarts-1))

!Initialize variables to zero
Vpot = 0.d0;Xcarts=0.d0;Force=0.d0

!Read input and ouput file names
open(unit=10,file='input.param',status='old')
read(10,*)cartCoordFile
read(10,*)EnForceOutFile

!Input and output files 
open(unit=11,file=trim(adjustl(cartCoordFile)),status='old')
open(unit=12,file=trim(adjustl(EnForceOutFile)),status='unknown')

!Read input geometries from *cartCoordFile* -> compute energies and forces -> write to *EnForceOutFile* 
write(12,*)"En (cm^{-1}) and forces {hartee/Ang}"
linectr = 0
do
	read(11,*,iostat=iost)(Xcarts(i),i=0,Ncarts-1) 
	call getforce(Xcarts,Force,Vpot)
	linectr = linectr+1
    ! Check for end of file or insufficient characters
    if (iost /= 0) then
        if (iost == -1) then
            write(*, *) "End of file reached."
        else
            write(*, *) "Error reading line from file."
			write(*, *) "Check line ",linectr
        end if
        exit
    end if
	write(12,'(F10.2)',advance='no')Vpot
	do j = 0,Ncarts-1
		write(12,'(F16.8)',advance='no')Force(j)
	enddo
	write(12,*)
enddo
end program
