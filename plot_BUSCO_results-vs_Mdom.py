#!/bin/python3

# plot_BUSCO_results-vs_Mdom.py uses the buscoplotpy package 
# to plot Fcan V2 (female primary) synteny pair vs M domesticus M3
# For establishing numbering and orientation of autosomes.

# Takes busco output full_table.tsv (all against same lineage)
# and karyotype.tsv file with "chr" "sequence" and "end" headers, end=length
# Row per sequence to plot (eg with busco hits), can make from .fai

# see repository at github.com/lorenzo-arcioni/BUSCO-Plot-Py
# with example scripts in test.ipynb

# first need to activate the environment (in salloc):
# $ module load python_3
# $ source /project/vpgru/software/python_envs/buscoplotpy/bin/activate

# then either open python interactive session and paste commands or run whole script via:
# $ python /path/to/plot_BUSCO_results.py


import pandas as pd # type: ignore
from buscoplotpy.utils.load_busco_fulltable import load_busco_fulltable # type: ignore
from buscoplotpy.graphics.karyoplot import karyoplot # type: ignore
from buscoplotpy.graphics.synteny import horizontal_synteny_plot, vertical_synteny_plot # type: ignore

# Mdom
FT_MdomM3 = load_busco_fulltable('/90daydata/vpgru/DavidLuecke/MusD/M3/OUT/busco/MdomM3-pctg_v3_diptera/run_diptera_odb10/full_table.tsv')
KT_MdomM3 = pd.read_csv('/90daydata/vpgru/DavidLuecke/MusD/M3/OUT/busco/MdomM3-pctg_v3-karyotype.tsv', sep='\t')
KT_MdomM3['organism'] = 'M. domestica M3 male'
KT_MdomM3['color'] = "#2C2C2C"

# Fcan
FT_FcanF = load_busco_fulltable('/90daydata/vpgru/DavidLuecke/Fcan/V2-frozen/female/pri/FcanF_V2-pri_diptera/run_diptera_odb10/full_table.tsv')
KT_FcanF = pd.read_csv('/90daydata/vpgru/DavidLuecke/Fcan/V2-frozen/female/pri/FcanF_V2-pri-karyotype.tsv', sep='\t')
KT_FcanF['organism'] = 'F. canicularis female'
KT_FcanF['color'] = "#036314"

# Dmel
FT_Dmel = load_busco_fulltable('/90daydata/vpgru/DavidLuecke/Dmel/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic_diptera/run_diptera_odb10/full_table.tsv')
KT_Dmel = pd.read_csv('/90daydata/vpgru/DavidLuecke/Dmel/Dmel-karyotype.tsv', sep='\t')
KT_Dmel['organism'] = 'D. melanogaster'
KT_Dmel['color'] = "#020077"

# M dom chromosome colors
link_colors = {
    'Chromosome1': "#772877",
    'Chromosome2': "#6b77bd",
    'Chromosome3': "#8a1515",
    'Chromosome4': '#ff9900',
    'Chromosome5': "#58A866",
}

horizontal_synteny_plot(ft_1 = FT_MdomM3, ft_2 = FT_Dmel, karyotype_1 = KT_MdomM3, karyotype_2 = KT_Dmel, link_colors=link_colors, output_path = 'synteny_MdomM3-Dmel.png')
horizontal_synteny_plot(ft_1 = FT_MdomM3, ft_2 = FT_FcanF, karyotype_1 = KT_MdomM3, karyotype_2 = KT_FcanF, link_colors=link_colors, output_path = 'synteny_MdomM3-FcanFV2.png')
