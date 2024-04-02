#!/bin/bash -ue
touch genomescope.flag.txt
xvfb-run genomescope.R -i Maut_F-HiCphasing_SimplePolish_Yahs-b.histo -o Maut_F-HiCphasing_SimplePolish_Yahs-b -k 21 -p 2 --fitted_hist
genomescope.R --version > version.txt
awk '/kmercov [0-9]/ { print $2 }' Maut_F-HiCphasing_SimplePolish_Yahs-b/model.txt >> kcov.txt
echo "finished genomescope"
sleep 120;
exit 0;
