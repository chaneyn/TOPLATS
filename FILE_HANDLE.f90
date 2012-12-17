!##################################################################
!Subroutine to handle opening and closing the input files
!##################################################################

subroutine FILE_OPEN()

implicit none

character(len=200) :: filename
integer :: soilnvars,vegnvars,dvegnvars,nforcingvars,noutvars
integer :: nrow,ncol
soilnvars = 23
vegnvars = 20
noutvars = 1
dvegnvars = 2
nforcingvars = 7;
nrow = 60
ncol = 66

!Global Parameter File

filename = "/home/ice/nchaney/PROJECTS/TOPLATS_DEVELOPMENT/DATA/LittleRiver/GLOBAL_PARAMETER.txt"
open(1000,file=trim(filename))

!Soil Parameter File

filename = "/home/ice/nchaney/PROJECTS/TOPLATS_DEVELOPMENT/DATA/LittleRiver/GRADS/soil.bin"
open(1001,file=trim(filename),status='old',access='direct',form='unformatted',recl=ncol*nrow*soilnvars*4)

!Vegetation Static Parameter File

filename = "/home/ice/nchaney/PROJECTS/TOPLATS_DEVELOPMENT/DATA/LittleRiver/GRADS/veg.bin"
open(1002,file=trim(filename),status='old',access='direct',form='unformatted',recl=ncol*nrow*vegnvars*4)

!Vegetation Dynamic Parameter File

filename = "/home/ice/nchaney/PROJECTS/TOPLATS_DEVELOPMENT/DATA/LittleRiver/GRADS/dveg.bin"
open(1003,file=trim(filename),status='old',access='direct',form='unformatted',recl=ncol*nrow*dvegnvars*4)

!Forcing Data Set
filename = "/home/ice/nchaney/PROJECTS/TOPLATS_DEVELOPMENT/DATA/LittleRiver/GRADS/FORCING_01202005_10202005.bin"
open(1004,file=trim(filename),status='old',access='direct',form='unformatted',recl=ncol*nrow*nforcingvars*4)

!Output Data Set
filename = "/home/ice/nchaney/PROJECTS/TOPLATS_DEVELOPMENT/DATA/LittleRiver/GRADS/OUTPUT.bin"
open(1005,file=trim(filename),status='unknown',access='direct',form='unformatted',recl=ncol*nrow*noutvars*4)

!Regional Actual Energy Fluxes (VALIDATION FILE)
filename = "/home/ice/nchaney/PROJECTS/TOPLATS_DEVELOPMENT/DATA/LittleRiver/OLD/EF_fruit.txt"
open(2091,file=trim(filename))

end subroutine

subroutine FILE_CLOSE()

implicit none
integer :: FILEUNIT

!Global Parameter File
close(1000)

!Soil Parameter File
close(1001)

!Vegetation Static Parameter File
close(1002)

!Vegetation Dynamic Parameter File
close(1003)

!Forcing Data Set
close(1004)

!Regional Actual Energy Fluxes (VALIDATION FILE)
close(2091)

!Output Data Set
close(1005)

end subroutine
