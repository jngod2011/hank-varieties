#!/bin/sh

# argument NAME is stub for directory
# make sure that shock indicators in SetParameters.f90 are set = 0 with no spaces

NAME="$1"

doTFP=1
doBorrWedge=1
doNews=1
doKappa0=1
doKappa1=1
doMonetary=1
doForwardGuide=0
doPref=1
doFundLev=0
doRiskAversion=0
doMarkup=1
doGovExp=1
doTransfer=0
doFinWedge=1
doLabWedge=1
doProdDisp=1
doProdPers=0
doTaylorPath=1

if [ "$doTFP" = "1" ] ; then
	echo "compiling, submiting ${NAME}_TFP" 
    perl -pi -e "s/IncludeTFPShock=0/IncludeTFPShock=1/g" SetParameters.f90
    make compilesubmit OUT=${NAME}_TFP
    perl -pi -e "s/IncludeTFPShock=1/IncludeTFPShock=0/g" SetParameters.f90
fi

if [ "$doBorrWedge" = "1" ] ; then
	echo "compiling, submiting ${NAME}_BorrWedge" 
    perl -pi -e "s/IncludeBorrWedgeShock=0/IncludeBorrWedgeShock=1/g" SetParameters.f90
    make compilesubmit OUT=${NAME}_BorrWedge
    perl -pi -e "s/IncludeBorrWedgeShock=1/IncludeBorrWedgeShock=0/g" SetParameters.f90
fi

if [ "$doNews" = "1" ] ; then
	echo "compiling, submiting ${NAME}_News" 
    perl -pi -e "s/IncludeNewsShock=0/IncludeNewsShock=1/g" SetParameters.f90
    make compilesubmit OUT=${NAME}_News
    perl -pi -e "s/IncludeNewsShock=1/IncludeNewsShock=0/g" SetParameters.f90
fi

if [ "$doKappa0" = "1" ] ; then
	echo "compiling, submiting ${NAME}_Kappa0" 
    perl -pi -e "s/IncludeKappa0Shock=0/IncludeKappa0Shock=1/g" SetParameters.f90
    make compilesubmit OUT=${NAME}_Kappa0
    perl -pi -e "s/IncludeKappa0Shock=1/IncludeKappa0Shock=0/g" SetParameters.f90
fi

if [ "$doKappa1" = "1" ] ; then
	echo "compiling, submiting ${NAME}_Kappa1" 
    perl -pi -e "s/IncludeKappa1Shock=0/IncludeKappa1Shock=1/g" SetParameters.f90
    make compilesubmit OUT=${NAME}_Kappa1
    perl -pi -e "s/IncludeKappa1Shock=1/IncludeKappa1Shock=0/g" SetParameters.f90
fi

if [ "$doMonetary" = "1" ] ; then
	echo "compiling, submiting ${NAME}_Monetary" 
    perl -pi -e "s/IncludeMonetaryShock=0/IncludeMonetaryShock=1/g" SetParameters.f90
    make compilesubmit OUT=${NAME}_Monetary
    perl -pi -e "s/IncludeMonetaryShock=1/IncludeMonetaryShock=0/g" SetParameters.f90
fi

if [ "$doForwardGuide" = "1" ] ; then
	echo "compiling, submiting ${NAME}_ForwardGuide" 
    perl -pi -e "s/IncludeForwardGuideShock=0/IncludeForwardGuideShock=1/g" SetParameters.f90
    make compilesubmit OUT=${NAME}_ForwardGuide
    perl -pi -e "s/IncludeForwardGuideShock=1/IncludeForwardGuideShock=0/g" SetParameters.f90
fi

if [ "$doPref" = "1" ] ; then
	echo "compiling, submiting ${NAME}_Pref" 
    perl -pi -e "s/IncludePrefShock=0/IncludePrefShock=1/g" SetParameters.f90
    make compilesubmit OUT=${NAME}_Pref
    perl -pi -e "s/IncludePrefShock=1/IncludePrefShock=0/g" SetParameters.f90
fi

if [ "$doFundLev" = "1" ] ; then
	echo "compiling, submiting ${NAME}_FundLev" 
    perl -pi -e "s/IncludeFundLevShock=0/IncludeFundLevShock=1/g" SetParameters.f90
    make compilesubmit OUT=${NAME}_FundLev
    perl -pi -e "s/IncludeFundLevShock=1/IncludeFundLevShock=0/g" SetParameters.f90
fi

if [ "$doRiskAversion" = "1" ] ; then
	echo "compiling, submiting ${NAME}_RiskAversion" 
    perl -pi -e "s/IncludeRiskAversionShock=0/IncludeRiskAversionShock=1/g" SetParameters.f90
    make compilesubmit OUT=${NAME}_RiskAversion
    perl -pi -e "s/IncludeRiskAversionShock=1/IncludeRiskAversionShock=0/g" SetParameters.f90
fi

if [ "$doMarkup" = "1" ] ; then
	echo "compiling, submiting ${NAME}_Markup" 
    perl -pi -e "s/IncludeMarkupShock=0/IncludeMarkupShock=1/g" SetParameters.f90
    make compilesubmit OUT=${NAME}_Markup
    perl -pi -e "s/IncludeMarkupShock=1/IncludeMarkupShock=0/g" SetParameters.f90
fi

if [ "$doGovExp" = "1" ] ; then
	echo "compiling, submiting ${NAME}_GovExp" 
    perl -pi -e "s/IncludeGovExpShock=0/IncludeGovExpShock=1/g" SetParameters.f90
    make compilesubmit OUT=${NAME}_GovExp
    perl -pi -e "s/IncludeGovExpShock=1/IncludeGovExpShock=0/g" SetParameters.f90
fi

if [ "$doTransfer" = "1" ] ; then
	echo "compiling, submiting ${NAME}_Transfer" 
    perl -pi -e "s/IncludeTransferShock=0/IncludeTransferShock=1/g" SetParameters.f90
    make compilesubmit OUT=${NAME}_Transfer
    perl -pi -e "s/IncludeTransferShock=1/IncludeTransferShock=0/g" SetParameters.f90
fi

if [ "$doFinWedge" = "1" ] ; then
	echo "compiling, submiting ${NAME}_FinWedge" 
    perl -pi -e "s/IncludeFinWedgeShock=0/IncludeFinWedgeShock=1/g" SetParameters.f90
    make compilesubmit OUT=${NAME}_FinWedge
    perl -pi -e "s/IncludeFinWedgeShock=1/IncludeFinWedgeShock=0/g" SetParameters.f90
fi

if [ "$doLabWedge" = "1" ] ; then
	echo "compiling, submiting ${NAME}_LabWedge" 
    perl -pi -e "s/IncludeLabWedgeShock=0/IncludeLabWedgeShock=1/g" SetParameters.f90
    make compilesubmit OUT=${NAME}_LabWedge
    perl -pi -e "s/IncludeLabWedgeShock=1/IncludeLabWedgeShock=0/g" SetParameters.f90
fi

if [ "$doProdDisp" = "1" ] ; then
	echo "compiling, submiting ${NAME}_ProdDisp" 
    perl -pi -e "s/IncludeProdDispShock=0/IncludeProdDispShock=1/g" SetParameters.f90
    make compilesubmit OUT=${NAME}_ProdDisp
    perl -pi -e "s/IncludeProdDispShock=1/IncludeProdDispShock=0/g" SetParameters.f90
fi

if [ "$doProdPers" = "1" ] ; then
	echo "compiling, submiting ${NAME}_ProdPers" 
    perl -pi -e "s/IncludeProdPersShock=0/IncludeProdPersShock=1/g" SetParameters.f90
    make compilesubmit OUT=${NAME}_ProdPers
    perl -pi -e "s/IncludeProdPersShock=1/IncludeProdPersShock=0/g" SetParameters.f90
fi

if [ "$doTaylorPath" = "1" ] ; then
	echo "compiling, submiting ${NAME}_TaylorPath" 
    perl -pi -e "s/IncludeTaylorPathShock=0/IncludeTaylorPathShock=1/g" SetParameters.f90
    make compilesubmit OUT=${NAME}_TaylorPath
    perl -pi -e "s/IncludeTaylorPathShock=1/IncludeTaylorPathShock=0/g" SetParameters.f90
fi