#!/bin/tcsh

set rootCmd = "root -b -l -q -x macros/fitWrapper.C"
set fitDir = fitFiles
set outDir = pdfs

set augerAndProton = 0
set galacticMass = 1

if ($galacticMass) then
  set masses = "12 16 20 24 28 32 36 40 44 48 52 56"
  foreach mass ($masses)
     sed -e 's/@MASS@/'$mass'/' $fitDir/GalacticMass.txt.in > $fitDir/GalacticMass$mass.txt
     setenv FITOPTION GalacticMass$mass
     $rootCmd |& tee $outDir/$FITOPTION.log
     mv $fitDir/$FITOPTION.txt $outDir
  end
endif

if ($augerAndProton) then
   setenv FITOPTION protonAGN
   $rootCmd |& tee $outDir/$FITOPTION.log
   cp $fitDir/$FITOPTION.txt $outDir

   setenv FITOPTION protonSFR2
   $rootCmd |& tee $outDir/$FITOPTION.log
   cp $fitDir/$FITOPTION.txt $outDir
endif
