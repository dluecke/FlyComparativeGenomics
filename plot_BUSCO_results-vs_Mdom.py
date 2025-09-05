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

# Fcan
FT_FcanF = load_busco_fulltable('/90daydata/vpgru/DavidLuecke/Fcan/V2-frozen/female/pri/FcanF_V2-pri_diptera/run_diptera_odb10/full_table.tsv')
KT_FcanF = pd.read_csv('/90daydata/vpgru/DavidLuecke/Fcan/V2-frozen/female/pri/FcanF_V2-pri-karyotype.tsv', sep='\t')
KT_FcanF['organism'] = 'F. canicularis female'
KT_FcanF['color'] = "#036314"

# Mdom
FT_MdomM3 = load_busco_fulltable('/90daydata/vpgru/DavidLuecke/MusD/M3/OUT/busco/MdomM3-pctg_v3_diptera/run_diptera_odb10/full_table.tsv')
KT_MdomM3 = pd.read_csv('/90daydata/vpgru/DavidLuecke/MusD/M3/OUT/busco/MdomM3-pctg_v3-karyotype.tsv', sep='\t')
KT_MdomM3['organism'] = 'M. domestica M3 male'
KT_MdomM3['color'] = "#2C2C2C"

# Dmel
FT_Dmel = load_busco_fulltable('/90daydata/vpgru/DavidLuecke/Dmel/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic_diptera/run_diptera_odb10/full_table.tsv')
KT_Dmel = pd.read_csv('/90daydata/vpgru/DavidLuecke/Dmel/Dmel-karyotype.tsv', sep='\t')
KT_Dmel['organism'] = 'D. melanogaster'
KT_Dmel['color'] = "#020077"

# confirm Muller elements on M. dom
link_colors = {
    'NC_004354.4': "#772877",
    'NT_033779.5': "#6b77bd",
    'NT_033778.4': "#8a1515",
    'NT_037436.4': '#ff9900',
    'NT_033777.3': "#58A866",
    'NC_004353.4': "#E98484",
}
horizontal_synteny_plot(ft_1 = FT_Dmel, ft_2 = FT_MdomM3, karyotype_1 = KT_Dmel, karyotype_2 = KT_MdomM3, link_colors=link_colors, output_path = 'synteny_MdomM3-Dmel.png')

# synteny between Mdom and Fcan
horizontal_synteny_plot(ft_1 = FT_FcanF, ft_2 = FT_MdomM3, karyotype_1 = KT_FcanF, karyotype_2 = KT_MdomM3, output_path = 'synteny_MdomM3-FcanF.png')
