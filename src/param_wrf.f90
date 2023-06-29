      module param_wrf
      SAVE
!----- user must change ----------- BEGIN -------------------------
      integer, parameter :: DateStrLen=19
      character(len=DateStrLen)    :: Times
      real :: bucket_threshold_mm=0.
      integer :: bucket_count=0


      integer :: west_east 
      integer :: south_north
      integer :: bottom_top 
      integer :: bottom_top_stag
      integer :: soil_layers_stag
      integer :: west_east_stag 
      integer :: south_north_stag
      real, allocatable,dimension(:) :: XTIME
      real, allocatable,dimension(:,:) :: XLAT
      real, allocatable,dimension(:,:) :: XLONG
      real, allocatable,dimension(:,:) :: XLAT_U
      real, allocatable,dimension(:,:) :: XLONG_U
      real, allocatable,dimension(:,:) :: XLAT_V
      real, allocatable,dimension(:,:) :: XLONG_V
      real, allocatable,dimension(:,:,:) :: U
      real, allocatable,dimension(:,:,:) :: V
      real, allocatable,dimension(:,:,:) :: W
      real, allocatable,dimension(:,:,:) :: PH
      real, allocatable,dimension(:,:,:) :: PHB
      real, allocatable,dimension(:,:,:) :: T
      real, allocatable,dimension(:,:,:) :: PB
      real, allocatable,dimension(:,:,:) :: P
      real, allocatable,dimension(:,:,:) :: QVAPOR
      real, allocatable,dimension(:,:) :: RAINC
      real, allocatable,dimension(:,:) :: RAINNC
      integer, allocatable,dimension(:,:) :: I_RAINC
      integer, allocatable,dimension(:,:) :: I_RAINNC
      real, allocatable,dimension(:,:) :: RAINSH
      real, allocatable,dimension(:,:) :: SNOWNC
      real, allocatable,dimension(:,:) :: HAILNC
      real, allocatable,dimension(:,:) :: GRAUPELNC
      real, allocatable,dimension(:,:) :: SST
      real, allocatable,dimension(:,:) :: SSTSK
      real, allocatable,dimension(:,:) :: TSK
      real, allocatable,dimension(:,:) :: EMISS
      real, allocatable,dimension(:,:) :: HGT
      real, allocatable,dimension(:,:) :: COSALPHA
      real, allocatable,dimension(:,:) :: SINALPHA
      real, allocatable,dimension(:,:) :: T2
      real, allocatable,dimension(:,:) :: U10
      real, allocatable,dimension(:,:) :: V10
      real, allocatable,dimension(:,:) :: HFX
      real, allocatable,dimension(:,:) :: LH
      real, allocatable,dimension(:,:) :: GRDFLX
      real, allocatable,dimension(:,:) :: SWDOWN
      real, allocatable,dimension(:,:) :: GLW
      real, allocatable,dimension(:,:) :: SNOWH
      real, allocatable,dimension(:,:) :: ZNT

      real, allocatable,dimension(:,:,:) :: qcloud
      real, allocatable,dimension(:,:,:) :: qgraup
      real, allocatable,dimension(:,:,:) :: qice
      real, allocatable,dimension(:,:,:) :: qrain
      real, allocatable,dimension(:,:,:) :: qsnow

      real, allocatable,dimension(:,:) :: PBLH
      real, allocatable,dimension(:,:) :: PSFC
      real, allocatable,dimension(:,:) :: ALBEDO
      real, allocatable,dimension(:,:) :: UST
      real, allocatable,dimension(:,:) :: RMOL
!soubroutines output
      real, allocatable,dimension(:,:,:) :: insituT
      real, allocatable,dimension(:,:,:) :: press
      real, allocatable,dimension(:,:,:) :: rh
      real, allocatable,dimension(:,:,:) :: phi
      real, allocatable,dimension(:,:,:) :: geohgt_asl
      real, allocatable,dimension(:,:,:) :: geohgt_agd
      real, allocatable,dimension(:,:,:) :: zwrf_asl
      real, allocatable,dimension(:,:,:) :: zwrf_agd
      real, allocatable,dimension(:,:)   :: tcc
      real, allocatable,dimension(:,:)   :: totprec
      real, allocatable,dimension(:,:,:) :: uint
      real, allocatable,dimension(:,:,:) :: vint
      real, allocatable,dimension(:,:)   :: xwrf
      real, allocatable,dimension(:,:)   :: ywrf
      real, allocatable,dimension(:,:)   :: precflx
      real, allocatable,dimension(:,:)   :: nrad
      real, allocatable,dimension(:,:)   :: trad
      end module param_wrf
