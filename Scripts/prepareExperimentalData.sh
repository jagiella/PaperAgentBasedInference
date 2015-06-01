# Growth curves
C=(I II III IV);
O=(0.28 0.28 0.28 0.07);
G=(1 5 25 25);
for i in 0 1 2 3; do
	cp /Users/jagiella/Dropbox/Work/PUBLICATIONS/Paper-wp2-lungsys/data/SK-MES-1/SK-MES_Radius_O${O[$i]}_G${G[$i]}.dat SK-MES1_${C[$i]}_GC.dat
done

# Profiles
for c in II III; do
	for t in 3 4; do
		cat ~/Dropbox/Work/PUBLICATIONS/Paper-wp2-lungsys/data/SK-MES-1-CryoSections-MACBOOK/Condition${c}/T${t}/Ki67_MEDIAN_histDiv_StandardDerivation.dat | awk '{print $1, $4, $12}' > SK-MES1_${c}_T${t}_Ki67.dat
		cat ~/Dropbox/Work/PUBLICATIONS/Paper-wp2-lungsys/data/SK-MES-1-CryoSections-MACBOOK/Condition${c}/T${t}/TUNEL_MEDIAN_histDiv_StandardDerivation.dat | awk '{print $1, $4, $12}' > SK-MES1_${c}_T${t}_TUNEL.dat
		cat ~/Dropbox/Work/PUBLICATIONS/Paper-wp2-lungsys/data/SK-MES-1-CryoSections-MACBOOK/Condition${c}/T${t}/ColIV_ECM_histDiv_StandardDerivation.dat | awk '{print $1, $2/256., $5/256.}' > SK-MES1_${c}_T${t}_ECM.dat
	done
done