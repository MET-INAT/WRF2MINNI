
module load forair_lonlat

pathlib=pathlib=/gporq3/minni/FORAIR_LONLAT/software/libraries/LIB_INTEL19/
cat Makefile_CR6_ifort_tpl | sed "s!@@pathlib@@!$pathlib!g" > Makefile

make 
#make debug=1
