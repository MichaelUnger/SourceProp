#!/bin/tcsh

set options = "SFR2 AAGHRW05 AGN eMm40z10 eMm40z20  eMm40z30  eMm40z40 eMm40z50 eMm20z10 eMm20z20 eMm20z30 eMm20z40 eMm20z50 eM00z10 eM00z20 eM00z30 eM00z40 eM00z50 eMp20z10 eMp20z20 eMp20z30 eMp20z40 eMp20z50 eMp40z10 eMp40z20 eMp40z30 eMp40z40  eMp40z50"

foreach option ($options)
   set filename = evo_${option}.txt
   echo "fit("\""Marco/$filename"\"", true, true)"
   echo "# generated at `date` by $0" > $filename
   echo "OutDir Marco/pdfs" >> $filename
   echo "DataDir Data" >> $filename
   echo "evolution $option">> $filename
   echo "IRB CRPropaG12" >> $filename
   echo "energyBinShift +1" >> $filename
   echo "xmaxSigmaShift -1" >> $filename
   echo "par gammaInj -1.0 0.1 0 0 0" >> $filename
   echo "par deltaEsc   -6.53346e-01 0.1 -1.01 -0.2 0" >> $filename
   echo "mass 28 0.1 1 56 0 1" >> $filename
   echo "spectrumData Auger2013" >> $filename
   echo "PhotonBPL 0.07 32 2.0" >> $filename
end   
