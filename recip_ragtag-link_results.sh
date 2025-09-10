#!/bin/bash

# recip_ragtag-link_results.sh collects final output from 
# makes RT_final/ output directory and writes s links to final output assemblies

mkdir RT_final
cd RT_final
ln -s ../{input.FpA} Fp_RTfinal.fa
ln -s ../{input.Fh1} Fh1_RTfinal.fa
ln -s ../{input.Fh2} Fh2_RTfinal.fa
ln -s ../{input.MpA} Mp_RTfinal.fa
ln -s ../{input.Mh1} Mh1_RTfinal.fa
ln -s ../{input.Mh2} Mh2_RTfinal.fa
ln -s ../{input.FpB} ALT-Fp_RTfinal.fa
ln -s ../{input.MpB} ALT-Mp_RTfinal.fa
cd ..
