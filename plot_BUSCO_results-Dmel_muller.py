#!/bin/python3

# plot_BUSCO_results-Dmel_muller.py uses the buscoplotpy package 
# to plot Muscid busco results as synteny pairs against D. melanogaster to track Muller elements

# first need to activate the environment (in salloc):
# $ module load python_3
# $ source /project/vpgru/software/python_envs/buscoplotpy/bin/activate

# then either open python interactive session and paste commands or run whole script via:
# $ python /path/to/plot_BUSCO_results-Dmel_muller.py


import pandas as pd # type: ignore
from buscoplotpy.utils.load_busco_fulltable import load_busco_fulltable # type: ignore
from buscoplotpy.graphics.karyoplot import karyoplot # type: ignore
from buscoplotpy.graphics.synteny import horizontal_synteny_plot, vertical_synteny_plot # type: ignore

# Dmel
FT_Dmel = load_busco_fulltable('/90daydata/vpgru/DavidLuecke/Dmel/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic_diptera/run_diptera_odb10/full_table.tsv')
KT_Dmel = pd.read_csv('/90daydata/vpgru/DavidLuecke/Dmel/Dmel-karyotype.tsv', sep='\t')
KT_Dmel['organism'] = 'D. melanogaster'
KT_Dmel['color'] = '#808080'

# Fcan
FT_FcanF = load_busco_fulltable('/90daydata/vpgru/DavidLuecke/Finished/Fcan/female/Fcan_F-v1_diptera/run_diptera_odb10/full_table.tsv')
KT_FcanF = pd.read_csv('/90daydata/vpgru/DavidLuecke/Finished/Fcan/female/FcanF_karyotype.tsv', sep='\t')
KT_FcanF['organism'] = 'F. canicularis female'
KT_FcanF['color'] = '#52b765'

# Haen
FT_HaenF = load_busco_fulltable('/90daydata/vpgru/DavidLuecke/Finished/Haen/female/Haen_F-v1_diptera/run_diptera_odb10/full_table.tsv')
KT_HaenF = pd.read_csv('/90daydata/vpgru/DavidLuecke/Finished/Haen/female/HaenF_karyotype.tsv', sep='\t')
KT_HaenF['organism'] = 'H. aenescens female'
KT_HaenF['color'] = '#fd9020'

# Hirr
FT_HirrF = load_busco_fulltable('/90daydata/vpgru/DavidLuecke/Finished/Hirr/female/Hirr_F-v2.masked_diptera/run_diptera_odb10/full_table.tsv')
KT_HirrF = pd.read_csv('/90daydata/vpgru/DavidLuecke/Finished/Hirr/female/HirrF_karyotype.tsv', sep='\t')
KT_HirrF['organism'] = 'H. irritans female'
KT_HirrF['color'] = '#dccd45'

# Maut
FT_MautF = load_busco_fulltable('/90daydata/vpgru/DavidLuecke/Finished/Maut/female/Maut_F-v2.masked_diptera/run_diptera_odb10/full_table.tsv')
KT_MautF = pd.read_csv('/90daydata/vpgru/DavidLuecke/Finished/Maut/female/MautF_karyotype.tsv', sep='\t')
KT_MautF['organism'] = 'M autumnalis female'
KT_MautF['color'] = '#006994'

# Scal
FT_ScalF = load_busco_fulltable('/90daydata/vpgru/DavidLuecke/Finished/Scal/female/Scal_F-v1_diptera/run_diptera_odb10/full_table.tsv')
KT_ScalF = pd.read_csv('/90daydata/vpgru/DavidLuecke/Finished/Scal/female/ScalF_karyotype.tsv', sep='\t')
KT_ScalF['organism'] = 'S. calcitrans female'
KT_ScalF['color'] = '#d43a3f'

# all Muscids vs Dmel

link_colors = {
    'NC_004354.4': '#800080',
    'NT_033779.5': '#3f4c99',
    'NT_033778.4': '#851d1d',
    'NT_037436.4': '#ff9900',
    'NT_033777.3': '#368244',
    'NC_004353.4': '#000000',
}

horizontal_synteny_plot(ft_1 = FT_Dmel, ft_2 = FT_FcanF, karyotype_1 = KT_Dmel, karyotype_2 = KT_FcanF, link_colors=link_colors, output_path = 'synteny_Dmel-FcanF.png')

horizontal_synteny_plot(ft_1 = FT_Dmel, ft_2 = FT_HaenF, karyotype_1 = KT_Dmel, karyotype_2 = KT_HaenF, link_colors=link_colors, output_path = 'synteny_Dmel-HaenF.png')

horizontal_synteny_plot(ft_1 = FT_Dmel, ft_2 = FT_HirrF, karyotype_1 = KT_Dmel, karyotype_2 = KT_HirrF, link_colors=link_colors, output_path = 'synteny_Dmel-HirrF.png')

horizontal_synteny_plot(ft_1 = FT_Dmel, ft_2 = FT_MautF, karyotype_1 = KT_Dmel, karyotype_2 = KT_MautF, link_colors=link_colors, output_path = 'synteny_Dmel-MautF.png')

horizontal_synteny_plot(ft_1 = FT_Dmel, ft_2 = FT_ScalF, karyotype_1 = KT_Dmel, karyotype_2 = KT_ScalF, link_colors=link_colors, output_path = 'synteny_Dmel-ScalF.png')

