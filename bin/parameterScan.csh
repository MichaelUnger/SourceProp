#!/bin/tcsh

set rootCmd = "root -b -l -q -x macros/fitWrapper.C"
set fitDir = fitFiles
set outDir = pdfs

set augerAndProton = 0
set systematics = 0
set evolution = 0
set galacticMass = 0
set outliers = 0
set standard = 0
set photonfield = 0
set escape = 0
set irb = 1


if ($photonfield) then
  set eps0 = "0.03"
  set beta = "2.0"
  set alpha = "32"

  set eps0s = "0.025 0.01 0.03 0.05 0.07 0.10 0.15"
  foreach eeps0 ($eps0s)
    setenv FITOPTION PhotonField${eeps0}_${beta}_${alpha}
    sed -e 's/@EPS@/'$eeps0'/g;s/@ALPHA@/'$alpha'/g;s/@BETA@/'$beta'/g' \
    $fitDir/PhotonField.txt.in > $fitDir/$FITOPTION.txt
    $rootCmd |& tee $outDir/$FITOPTION.log
    mv $fitDir/$FITOPTION.txt $outDir
  end

  set betas = "1.5 1.8 2.2 2.5"
  foreach bbeta ($betas)
    setenv FITOPTION PhotonField${eps0}_${bbeta}_${alpha}
    sed -e 's/@EPS@/'$eps0'/g;s/@ALPHA@/'$alpha'/g;s/@BETA@/'$bbeta'/g' \
    $fitDir/PhotonField.txt.in > $fitDir/$FITOPTION.txt
    $rootCmd |& tee $outDir/$FITOPTION.log
    mv $fitDir/$FITOPTION.txt $outDir
  end

  set alphas = "12 22 42 52"
  foreach aalpha ($alphas)
    setenv FITOPTION PhotonField${eps0}_${beta}_${aalpha}
    sed -e 's/@EPS@/'$eps0'/g;s/@ALPHA@/'$aalpha'/g;s/@BETA@/'$beta'/g' \
    $fitDir/PhotonField.txt.in > $fitDir/$FITOPTION.txt
    $rootCmd |& tee $outDir/$FITOPTION.log
    mv $fitDir/$FITOPTION.txt $outDir
  end
endif

if ($standard) then

  setenv FITOPTION StandardParScan
  $rootCmd |& tee $outDir/$FITOPTION.log
  cp $fitDir/$FITOPTION.txt $outDir

endif

if ($outliers) then

  setenv FITOPTION SpectrumOutliers
  $rootCmd |& tee $outDir/$FITOPTION.log
  cp $fitDir/$FITOPTION.txt $outDir

endif


if ($irb) then
  set irbs = "Kneiske04 Stecker05 Kneiske04G4"
  foreach irb ($irbs)
    setenv FITOPTION IRB$irb
    sed -e 's/@IRB@/'$irb'/' $fitDir/IRB.txt.in > \
      $fitDir/$FITOPTION.txt
    $rootCmd |& tee $outDir/$FITOPTION.log
    mv $fitDir/$FITOPTION.txt $outDir
  end
endif

if ($escape) then

  set deltas = "1.000 0.875 0.750 0.625 0.500 0.375 0.250"
  foreach delta ($deltas)
    echo $delta
    setenv FITOPTION Escape$delta
    sed -e 's/@ESC@/'$delta'/' $fitDir/Escape.txt.in > \
      $fitDir/$FITOPTION.txt
    $rootCmd |& tee $outDir/$FITOPTION.log
    mv $fitDir/$FITOPTION.txt $outDir
  end
endif


if ($evolution) then

  set m = 10
  while ($m < 55)
    echo $m
    setenv FITOPTION Evolution$m
    sed -e 's/@EVO@/'$m'/' $fitDir/Evolution.txt.in > \
      $fitDir/$FITOPTION.txt
    $rootCmd |& tee $outDir/$FITOPTION.log
    mv $fitDir/$FITOPTION.txt $outDir
    @ m += 5
  end
endif

if ($systematics) then
  setenv FITOPTION Systematics00
  sed -e 's/@ESHIFT@/0/;s/@XMAXSHIFT@/0/' $fitDir/Systematics.txt.in > \
    $fitDir/$FITOPTION.txt
  $rootCmd |& tee $outDir/$FITOPTION.log
  mv $fitDir/$FITOPTION.txt $outDir

  setenv FITOPTION SystematicsNN
  sed -e 's/@ESHIFT@/-1/;s/@XMAXSHIFT@/-1/' $fitDir/Systematics.txt.in > \
    $fitDir/$FITOPTION.txt
  $rootCmd |& tee $outDir/$FITOPTION.log
  mv $fitDir/$FITOPTION.txt $outDir

  setenv FITOPTION SystematicsPP
  sed -e 's/@ESHIFT@/1/;s/@XMAXSHIFT@/1/' $fitDir/Systematics.txt.in > \
    $fitDir/$FITOPTION.txt
  $rootCmd |& tee $outDir/$FITOPTION.log
  mv $fitDir/$FITOPTION.txt $outDir

  setenv FITOPTION SystematicsNP
  sed -e 's/@ESHIFT@/-1/;s/@XMAXSHIFT@/1/' $fitDir/Systematics.txt.in > \
    $fitDir/$FITOPTION.txt
  $rootCmd |& tee $outDir/$FITOPTION.log
  mv $fitDir/$FITOPTION.txt $outDir

  setenv FITOPTION SystematicsPN
  sed -e 's/@ESHIFT@/1/;s/@XMAXSHIFT@/-1/' $fitDir/Systematics.txt.in > \
    $fitDir/$FITOPTION.txt
  $rootCmd |& tee $outDir/$FITOPTION.log
  mv $fitDir/$FITOPTION.txt $outDir

endif

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
