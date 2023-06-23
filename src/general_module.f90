      module general
      use datetime_module, only : timedelta, datetime
      IMPLICIT NONE
      SAVE
      type (datetime) :: current_date, old_date
      type (timedelta) :: deltat

      integer, parameter :: max_num_levs=100
      character(len=252) :: dirin ,dirout,Time_startsim,filein,fileou
      integer:: dayxfile
      real   :: nhour 
      logical :: rotate,fsstsk,fnrad ,fsgrad,flgrad,ftgrad,fhymet,fznt
      logical :: falbedo,fsh,flh,fgh,fustar,flstar,fpbl,use_T2,use_W10 

      logical :: utmcoordinate
      logical :: flag_interp=.true.
      real :: xstart,ystart, dx,dy
      integer :: nx,ny,nz,start_timestep
      real, dimension(max_num_levs) :: zlev=-999

      

      namelist/general_namelist/filein,fileou,dirin,dirout,   & 
      Time_startsim,start_timestep,dayxfile,nhour,                   &
      rotate, fsstsk,fnrad ,fsgrad,flgrad,ftgrad,falbedo,fsh,flh,fgh,&
      fustar,flstar,fpbl,fhymet,fznt,use_T2,use_W10

      namelist/output_namelist/utmcoordinate,xstart,ystart,nx,ny,nz,dx,dy,zlev
      end module general


