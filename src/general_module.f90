module general
use datetime_module! , only : timedelta, datetime
IMPLICIT NONE
SAVE
type (datetime) :: simulation_start_date,file_start_date,current_date
type (timedelta) :: deltat_increment

integer, parameter :: max_num_levs=100
integer :: dayxfile
character(len=256) :: dirin ,dirout,filein,out_prefix
character(len=256) :: time_units ! read from wrf output
real   :: nhour 
logical :: rotate,fsstsk,fnrad ,fsgrad,flgrad,ftgrad,fhymet,fznt
logical :: falbedo,fsh,flh,fgh,fustar,flstar,fpbl,use_T2,use_W10 

logical :: utmcoordinate
logical :: flag_interp=.true.
logical :: if_bucket=.false. !set as true if attribute BUCKET_MM is found in input file in main.f90)
real :: xstart,ystart, dx,dy
integer :: nx,ny,nz,start_timestep
real, dimension(max_num_levs) :: zlev=-999



namelist/general_namelist/filein,out_prefix,dirin,dirout,   & 
start_timestep, rotate, fsstsk,fnrad ,fsgrad,     &
flgrad,ftgrad,falbedo,fsh,flh,fgh, fustar,flstar,fpbl,  &
fhymet,fznt,use_T2,use_W10

namelist/output_namelist/utmcoordinate,xstart,ystart,nx,ny,nz,dx,dy,zlev

contains

subroutine check(istatus)
use netcdf
IMPLICIT NONE
integer,intent (in) :: istatus
if (istatus /= nf90_noerr) then
   write(*,*) trim(adjustl(nf90_strerror(istatus)))
   write(*,*) 'PROGRAM STOPS'
  STOP
endif
end subroutine check


end module general

