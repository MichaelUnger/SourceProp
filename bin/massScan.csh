#!/bin/tcsh

#set masses = "22 23 24 25 26 27 28 29 30 31 32 33 34"
#set fractions = "0.00 0.05 0.1 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50"
#set fractions = "0.55 0.60 0.65 0.70 0.75 0.80 0.85 0.90 0.95 1.00"
set masses = "23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40"
set fractions = "0.00 0.05 0.1 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50"

set masses = "1 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45"
#set masses = "22 23 24 32 33 34"
#set masses = "45 44 43 42 41"
#set masses = "35 36 37 38 39 40"
set masses = "41"
set fractions = "0.00 0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50 0.55 0.60 0.65 0.70 0.75 0.80 0.85 0.90 0.95 1.00"
set fractions = "0.40 0.45 0.50 0.55 0.60 0.65 0.70 0.75 0.80 0.85 0.90 0.95 1.00"
set testMass = 4
setenv EXEDIR $PWD
foreach mass ($masses)
  foreach fraction ($fractions)
     set fitname = massScan_${mass}_${testMass}_${fraction}
     if (! -e pdfs/${fitname}.root) then
       setenv FITFILE ${fitname}.txt
       echo "OutDir pdfs" > $FITFILE
       echo "DataDir data" >> $FITFILE
       echo "evolution SFR2" >> $FITFILE
       echo "IRB Gilmore12" >> $FITFILE
       echo "par gammaInj -1.0 0.1 0 0 1" >> $FITFILE
       echo "par deltaEsc   -0.9 0.1 -1.01 -0.2 0" >> $FITFILE
       echo "mass $testMass $fraction 1 56 1 1" >> $FITFILE
       echo "mass $mass 0.1 1 56 1 1" >> $FITFILE
       echo "spectrumData Auger2013" >> $FITFILE
       echo "PhotonBPL 0.07 32 2.0" >> $FITFILE
       echo "energyBinShift +1" >> $FITFILE
       echo "xmaxSigmaShift -1" >> $FITFILE
       root.exe -b -l -q -x macros/fitWrapper.C
       mv $FITFILE pdfs
    else
       echo '****************************************'
       echo "pdfs/${fitname}.root exists"
       echo '****************************************'
    endif
  end
end
