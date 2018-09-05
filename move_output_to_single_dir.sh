#!/bin/sh

# first argumentis is name of stub for directory
STUB="$1"


#create enclosing directory
mkdir -p $STUB

shock[1]="TFP"
shock[2]="Pref"
shock[3]="Kappa0"
shock[4]="Markup"
shock[5]="FinWedge"
shock[6]="LabWedge"
shock[7]="GovExp"
shock[8]="BorrWedge"
shock[9]="ProdDisp"
shock[10]="ProdPers"
shock[11]="Monetary"
shock[12]="Kappa1"
shock[13]="News"
shock[14]="ForwardGuide"
shock[15]="Transfer"
shock[16]="TaylorPath"

# for remaining shocks, move only IRF director
FIRST=1
for is in {1..16}
do
	if [ -d "${STUB}_${shock[${is}]}" ]; then #checks if directory exists
		
		if [ "$FIRST" = "1" ] ; then
			# for first shock, move all files
			mv "${STUB}_${shock[${is}]}"/* "${STUB}"/
			rm -r "${STUB}_${shock[${is}]}"
			FIRST=0
		else
			mv "${STUB}_${shock[${is}]}"/IRF_"${shock[${is}]}" "${STUB}"/
			mv "${STUB}_${shock[${is}]}"/output*.txt "${STUB}"/
			rm -r "${STUB}_${shock[${is}]}"
		fi
	fi
done

# no need to keep executable(s)
rm "${STUB}"/*.out
