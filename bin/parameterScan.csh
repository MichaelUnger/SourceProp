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
set irb = 0
set massScan = 0
set bigScan = 1

if ($bigScan) then
  set temperatures = "80 100 120 140 160 180 200 220 240 260 280 300 320 340 360 380 400"
  set temperatures = "50 100 150 200 250 300 350 400 450 500 550 600 650"
#  set temperatures = "300 350 400 450 500 550"
  set masses = "20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35"
#  set masses = "28"
  set gammas = "1 1.5 2"
  set deltas = "1 0.66667 0.333333"
  set gammas = "1"
  set deltas = "1"
  set nTot = 2448
  set nCurr = 0
  foreach temperature ($temperatures)
    foreach mass ($masses)
      foreach gamma ($gammas)
        foreach delta ($deltas)
           @ nCurr ++
           echo "\n----------- $nCurr of $nTot"
           setenv FITOPTION BigScan_${temperature}_${mass}_${gamma}_${delta}
           sed -e 's/@TEMP@/'$temperature'/g;s/@MASS@/'$mass'/g;s/@GAMMA@/'$gamma'/g;s/@DELTA@/'$delta'/g' \
           $fitDir/BigScan.txt.in > $fitDir/$FITOPTION.txt
           time $rootCmd |& tee $outDir/$FITOPTION.log
           mv $fitDir/$FITOPTION.txt $outDir
        end
      end
    end
  end
endif

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
  set irbs = "Kneiske04 Stecker05 Kneiske04G4 Kneiske10"
  set irbs = "Kneiske04G4"
  foreach irb ($irbs)
    setenv FITOPTION IRB$irb
    sed -e 's/@IRB@/'$irb'/' $fitDir/IRB.txt.in > \
      $fitDir/$FITOPTION.txt
    $rootCmd |& tee $outDir/$FITOPTION.log
    mv $fitDir/$FITOPTION.txt $outDir
  end
endif

if ($massScan) then
  set irb = "Kneiske04"
  set masses = "23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40"
  foreach mass ($masses)
    setenv FITOPTION MASS${irb}_${mass}_NP
    sed -e 's/@IRB@/'$irb'/;s/@MASS@/'$mass'/' $fitDir/MASS.txt.in > \
      $fitDir/$FITOPTION.txt
    $rootCmd |& tee $outDir/$FITOPTION.log
    mv $fitDir/$FITOPTION.txt $outDir
  end
  grep FCN pdfs/MASS${irb}_*_NP.log | awk -F FROM '{print $1}'
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

  setenv FITOPTION Systematics0N
  sed -e 's/@ESHIFT@/0/;s/@XMAXSHIFT@/-1/' $fitDir/Systematics.txt.in > \
    $fitDir/$FITOPTION.txt
  $rootCmd |& tee $outDir/$FITOPTION.log
  mv $fitDir/$FITOPTION.txt $outDir

  setenv FITOPTION Systematics0P
  sed -e 's/@ESHIFT@/0/;s/@XMAXSHIFT@/1/' $fitDir/Systematics.txt.in > \
    $fitDir/$FITOPTION.txt
  $rootCmd |& tee $outDir/$FITOPTION.log
  mv $fitDir/$FITOPTION.txt $outDir

  setenv FITOPTION SystematicsN0
  sed -e 's/@ESHIFT@/-1/;s/@XMAXSHIFT@/0/' $fitDir/Systematics.txt.in > \
    $fitDir/$FITOPTION.txt
  $rootCmd |& tee $outDir/$FITOPTION.log
  mv $fitDir/$FITOPTION.txt $outDir

  setenv FITOPTION SystematicsP0
  sed -e 's/@ESHIFT@/1/;s/@XMAXSHIFT@/0/' $fitDir/Systematics.txt.in > \
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
