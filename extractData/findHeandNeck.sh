rm -f /home/tbaudier/studies/cbctHeadNeck/code/listPatientHeadNeck.log
db="xvi_synergy_db_45"
cd /mnt/clb_storage/$db/

#for D in `find . -maxdepth 1 -type d -name "patient_*"`
patientFolder=`ls -rtd patient_*`
for D in $patientFolder
do
    iniFile=""
    tag=""
    echo $D
    if [ -d "/mnt/clb_storage/$db/$D/CT_SET/" ]; then
	cd /mnt/clb_storage/$db/$D/CT_SET/
    	iniFile=`find . -maxdepth 2 -type f -name "1.*.INI"`
	if [ -n "$iniFile" ]
    	then
      		tag=`cat $iniFile | grep -i arynx`
      		if [ -n "$tag" ]
      		then
			treatmentName=`cat $iniFile | grep -i "TreatmentID="`
        		echo "$D $treatmentName" >> /home/tbaudier/studies/cbctHeadNeck/code/listPatientHeadNeck.log
      		fi
    	fi
    	cd /mnt/clb_storage/$db/
    fi
done


