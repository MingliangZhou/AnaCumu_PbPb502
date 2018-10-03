#!/bin/sh
echo "RUNNING PROGRAM"

job=2
for ((I1=0;I1<110;I1+=1)) # 20*33
do
I2=$[$I1+1]
Filenm01=jobPhase1_$I1"_"$I2.job
Filenm02=jobPhase1_$I1"_"$I2.sh

cp condor_script $Filenm01
echo "Executable      = $Filenm02">>$Filenm01
echo "queue"                        >>$Filenm01

echo "#!/bin/sh    ">$Filenm02
echo 'echo $0'>>$Filenm02
echo "root -b -l <<EOF ">>$Filenm02
echo "gSystem->Load(\"../Tool_cxx.so\")">>$Filenm02
echo "gSystem->Load(\"../Event_cxx.so\")">>$Filenm02
echo "gSystem->Load(\"../MultiCorr_cxx.so\")">>$Filenm02
echo "gSystem->Load(\"../Phase1_cxx.so\")">>$Filenm02
echo "Phase1($I1)">>$Filenm02
echo ".q">>$Filenm02
echo "EOF">>$Filenm02
echo "rm $Filenm01">>$Filenm02
echo "rm $Filenm02">>$Filenm02

chmod 711 $Filenm02
#sh  $Filenm02
condor_submit $Filenm01 
done

