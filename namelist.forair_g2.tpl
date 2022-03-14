&general_namelist
filein="@@filein@@"
fileou="@@fileou@@"
dirin="@@dirin@@"
dirout="@@dirout@@"
Time_startsim="@@timestartsim@@"
start_timestep=7
dayxfile=1
nhour=1
rotate = .true.   !wind rotation
fsstsk = .true.   ! skin temperature
fnrad = .true.    !net radiation
fsgrad = .true.    !showrtwave global radiation  (needed by surfpro)
flgrad = .true.    !longwave  global radiation
ftgrad = .true.    !short+long  global radiation
falbedo = .true.  !albedo
fsh = .true.      !sensible heat flux
flh = .true.      !latent heat flux
fgh = .true.      !ground heat flux
fustar = .true.   !frinction velocity
flstar = .true.   !Monin Obukov
fpbl = .true.   !PBL
fhymet = .true. ! hydrometeor
fznt = .true.    ! roughness length
use_T2 = .true.   !use T2 temperature for vertical temperature interpolation
use_W10 = .true.   !use W10 
/

&output_namelist
utmcoordinate = .true.
xstart = 266
ystart = 4016
nx=278
ny=308
nz=14
dx=4
dy=4
zlev= 20 65 125 210 325 480 690 975 1360 1880 2580 3525 4805 6290
/
