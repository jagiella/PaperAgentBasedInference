# sigma = mean(f)
#for i in SK-MES1_I*.dat; do echo $i; mean=`cat $i | awk 'BEGIN{x=0; n=0}{if($1!=""&&$2!=0&&$3!=0) {x=x+$2;n++;}}END{print x/n}'`; cat $i | awk '{if($1!=""&&$2!=0&&$3!=0) print $1, $2, '${mean}'}' > $i.FixErr ; done

# sigma = max(f) - min(f)
#for i in SK-MES1_I*.dat; do echo $i; min=`cat $i | awk 'BEGIN{x=10000; n=0}{if($1!=""&&$2!=0&&$3!=0) {x=(x<$2?x:$2);n++;}}END{print x}'`; max=`cat $i | awk 'BEGIN{x=0; n=0}{if($1!=""&&$2!=0&&$3!=0) {x=(x>$2?x:$2);n++;}}END{print x}'`; cat $i | awk '{if($1!=""&&$2!=0&&$3!=0) print $1, $2, '${max}-${min}'}' > $i.FixErr ; done

# sigma = std
for i in SK-MES1_I*.dat; do echo $i; std=`cat $i | awk 'BEGIN{mean=0; mean2=0; n=0}{if($1!=""&&$2!=0&&$3!=0) {mean=mean+$2; mean2=mean2+$2*$2; n++;}}END{print sqrt(mean2/n - mean*mean/(n*n))}'`; cat $i | awk '{if($1!=""&&$2!=0&&$3!=0) print $1, $2, '${std}'}'  > $i.Std ; done
