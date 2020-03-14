#!/usr/bin/env bash

END=7;
filestr="";
searchstr="";
acid="E";
for i in $(seq 1 $END); 
do 

added=$(( $i + $END ))

if [ $added -lt 10 ]
then
filestr=$acid"0";
else
filestr="$acid";
fi 

if [ $i -lt 10 ]
then
searchstr=$acid"0";
else
searchstr="$acid";
fi

echo "$searchstr";
echo "$filestr";

val="$i";

regexSearch="$searchstr$val";



regexReplace="$filestr$added";

filename="$regexReplace.off";

echo "added: $added"; 
echo "regexSearch: $regexSearch"; 
echo "regexReplace: $regexReplace"; 
echo "fileName: $filename";
sed -i "s/$regexSearch/$regexReplace/g" "$filename";


done


