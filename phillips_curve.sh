#!/bin/sh

#make sure that MonetaryShockSize=0 in SetParameters.f90 before beginning (no spaces)
#make sure that 131 is indeed the relevant line number

name="$1"

pclist[1]=0.0125
pclist[2]=0.01
pclist[3]=0.0075
pclist[4]=0.005
pclist[5]=0.00375
pclist[6]=0.0025
pclist[7]=0.00125
pclist[8]=-0.00125
pclist[9]=-0.0025
pclist[10]=-0.00325
pclist[11]=-0.005
pclist[12]=-0.0075
pclist[13]=-0.01
pclist[14]=-0.0125


for ip in {1..14}
do
   echo MonetaryShockSize=${pclist[${ip}]}
   perl -pi -e "s/MonetaryShockSize=0/MonetaryShockSize=${pclist[${ip}]}/g" SetParameters.f90
   #sed -n -i.bak 's/MonetaryShockSize=0/MonetaryShockSize=${pclist[${ip}]}/g' SetParameters.f90
   make compilesubmit OUT=${name}_${ip}
   #echo "working on ${name}_${ip}"
   perl -pi -e "s/MonetaryShockSize=${pclist[${ip}]}/MonetaryShockSize=0/g" SetParameters.f90
   #sed -i '$131s/.*/MonetaryShockSize=0/' SetParameters.f90
done