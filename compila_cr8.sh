module load forair_portici

#lib di FORAIR
pathlib=/lpor1/store_0/usr/forair/FORAIR/OPERATIONAL/software/libraries/LIB_INTEL_ONEAPI

#lib di CAMS
#pathlib=/lpor1/store_0/usr/forecast/CAMS/software/libraries/LIB_INTEL_ONEAPI

cat Makefile_CR8_ifort_tpl | sed "s!@@pathlib@@!$pathlib!g" > Makefile

make 
#make debug=1
