#! /bin/bash
# string to run this:
# chmod 744 ./SubmitKMC.sh;./SubmitKMC.sh
if [ $SGE_CLUSTER_NAME = farber ]; then
  JobScriptName="Farber_Submit.qs"
elif [ $SGE_CLUSTER_NAME = squidward ]; then
  JobScriptName="Squidward_Submit.ps"
fi

for i in $(echo ./*/);do
 submit="${i}Zacros.ps"
 cp ./$JobScriptName $submit
 cd $i;qsub Zacros.ps;cd ..
done
