       module param_minni
       real, allocatable, dimension(:,:)   :: xfarm,yfarm
       real, allocatable, dimension(:,:)   :: hgtfarm
       real, allocatable, dimension(:,:)   :: sstfarm
       real, allocatable, dimension(:,:)   :: sstskfarm
       real, allocatable, dimension(:,:)   :: precfarm
       real, allocatable, dimension(:,:)   :: tccfarm
       real, allocatable, dimension(:,:)   :: pblfarm
       real, allocatable, dimension(:,:,:) :: tfarm
       real, allocatable, dimension(:,:,:)   :: ufarm
       real, allocatable, dimension(:,:,:)   :: vfarm
       real, allocatable, dimension(:,:,:)   :: wfarm
       real, allocatable, dimension(:,:,:)   :: rhfarm
       real, allocatable, dimension(:,:,:)   :: phifarm
       real, allocatable, dimension(:,:,:)   :: pfarm
       real, allocatable, dimension(:,:)   :: netradfarm
       real, allocatable, dimension(:,:)   :: totradfarm
       real, allocatable, dimension(:,:)   :: totlradfarm
       real, allocatable, dimension(:,:)   :: totsradfarm
       real, allocatable, dimension(:,:)   :: shfarm
       real, allocatable, dimension(:,:)   :: lhfarm
       real, allocatable, dimension(:,:)   :: ghfarm
       real, allocatable, dimension(:,:)   :: albedofarm
       real, allocatable, dimension(:,:)   :: ustarfarm
       real, allocatable, dimension(:,:)   :: lstarfarm
       real, allocatable, dimension(:,:,:) :: qcloudfarm
       real, allocatable, dimension(:,:,:) :: qicefarm
       real, allocatable, dimension(:,:,:) :: qrainfarm
       real, allocatable, dimension(:,:,:) :: qsnowfarm
       real, allocatable, dimension(:,:,:) :: qgraupfarm
       real, allocatable, dimension(:,:) :: XA,XB,XC,XD,YA,YB,YC,YD,WA,WB,WC,WD
       real, allocatable, dimension(:,:)   :: snowfarm
       real, allocatable, dimension(:,:)   :: z0farm
       end module param_minni
