subroutine init2fe()
use molparams
use potvars
implicit none 
call initializeMolParams()
call initializePotvars()
endsubroutine 
