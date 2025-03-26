#!/bin/python3

# plot_BUSCO_results-2025.py uses the buscoplotpy package 
# to plot Muscid busco results as karyotype and synteny pairs

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
FT_FcanF = load_busco_fulltable('/90daydata/vpgru/DavidLuecke/Fcan/female/pctg/yahs/Fcan_01b-pctg.review.FINAL_diptera/run_diptera_odb10/full_table.tsv')
KT_FcanF = pd.read_csv('/90daydata/vpgru/DavidLuecke/Fcan/female/pctg/yahs/Fcan_01b-pctg.review.FINAL_diptera/FcanF2025_karyotype.tsv', sep='\t')
KT_FcanF['organism'] = 'F. canicularis female'
KT_FcanF['color'] = '#52b765'

karyoplot(karyotype = KT_FcanF, fulltable = FT_FcanF, output_file = "karyoplot_FcanF2025.png")

FT_FcanM = load_busco_fulltable('/90daydata/vpgru/DavidLuecke/Fcan/male/pctg/yahs/FcanM_01b-pctg.review.FINAL_diptera/run_diptera_odb10/full_table.tsv')
KT_FcanM = pd.read_csv('/90daydata/vpgru/DavidLuecke/Fcan/male/pctg/yahs/FcanM_01b-pctg.review.FINAL_diptera/FcanM2025_karyotype.tsv', sep='\t')
KT_FcanM['organism'] = 'F. canicularis male'
KT_FcanM['color'] = '#368244'

karyoplot(karyotype = KT_FcanM, fulltable = FT_FcanM, output_file = "karyoplot_FcanM2025.png")

horizontal_synteny_plot(ft_1 = FT_FcanF, ft_2 = FT_FcanM, karyotype_1 = KT_FcanF, karyotype_2 = KT_FcanM, output_path = 'synteny_Fcan2025.png')

# Haen
FT_HaenF = load_busco_fulltable('/90daydata/vpgru/DavidLuecke/Haen/female/pctg/yahs/HaenF_01-pctg.review.FINAL_diptera/run_diptera_odb10/full_table.tsv')
KT_HaenF = pd.read_csv('/90daydata/vpgru/DavidLuecke/Haen/female/pctg/yahs/HaenF_01-pctg.review.FINAL_diptera/HaenF2025_karyotype.tsv', sep='\t')
KT_HaenF['organism'] = 'H. aenescens female'
KT_HaenF['color'] = '#fd9020'

karyoplot(karyotype = KT_HaenF, fulltable = FT_HaenF, output_file = "karyoplot_HaenF2025.png")

FT_HaenM = load_busco_fulltable('/90daydata/vpgru/DavidLuecke/Haen/male/pctg/yahs/HaenM_01-pctg.review.FINAL_diptera/run_diptera_odb10/full_table.tsv')
KT_HaenM = pd.read_csv('/90daydata/vpgru/DavidLuecke/Haen/male/pctg/yahs/HaenM_01-pctg.review.FINAL_diptera/HaenM2025_karyotype.tsv', sep='\t')
KT_HaenM['organism'] = 'H. aenescens male'
KT_HaenM['color'] = '#ca6702'

karyoplot(karyotype = KT_HaenM, fulltable = FT_HaenM, output_file = "karyoplot_HaenM2025.png")

horizontal_synteny_plot(ft_1 = FT_HaenF, ft_2 = FT_HaenM, karyotype_1 = KT_HaenF, karyotype_2 = KT_HaenM, output_path = 'synteny_Haen2025.png')

# Hirr
FT_HirrF = load_busco_fulltable('/90daydata/vpgru/DavidLuecke/Hirr/female/pctg/yahs/HirrF_01-pctg.review.FINAL_diptera/run_diptera_odb10/full_table.tsv')
KT_HirrF = pd.read_csv('/90daydata/vpgru/DavidLuecke/Hirr/female/pctg/yahs/HirrF_01-pctg.review.FINAL_diptera/HirrF2025_karyotype.tsv', sep='\t')
KT_HirrF['organism'] = 'H. irritans female'
KT_HirrF['color'] = '#dccd45'

karyoplot(karyotype = KT_HirrF, fulltable = FT_HirrF, output_file = "karyoplot_HirrF2025.png")

FT_HirrM = load_busco_fulltable('/90daydata/vpgru/DavidLuecke/Hirr/male/pctg/yahs/busco/HirrM_01-pctg.review.FINAL_diptera/run_diptera_odb10/full_table.tsv')
KT_HirrM = pd.read_csv('/90daydata/vpgru/DavidLuecke/Hirr/male/pctg/yahs/busco/HirrM_01-pctg.review.FINAL_diptera/HirrM2025_karyotype.tsv', sep='\t')
KT_HirrM['organism'] = 'H. irritans male'
KT_HirrM['color'] = '#afa021'

karyoplot(karyotype = KT_HirrM, fulltable = FT_HirrM, output_file = "karyoplot_HirrM2025.png")

horizontal_synteny_plot(ft_1 = FT_HirrF, ft_2 = FT_HirrM, karyotype_1 = KT_HirrF, karyotype_2 = KT_HirrM, output_path = 'synteny_Hirr2025.png')

# Maut
FT_MautF = load_busco_fulltable('/90daydata/vpgru/DavidLuecke/Maut/female/pctg/yahs/MautF_01-pctg.review.FINAL_diptera/run_diptera_odb10/full_table.tsv')
KT_MautF = pd.read_csv('/90daydata/vpgru/DavidLuecke/Maut/female/pctg/yahs/MautF_01-pctg.review.FINAL_diptera/MautF2025_karyotype.tsv', sep='\t')
KT_MautF['organism'] = 'M autumnalis female'
KT_MautF['color'] = '#006994'

karyoplot(karyotype = KT_MautF, fulltable = FT_MautF, output_file = "karyoplot_MautF2025.png")

FT_MautM = load_busco_fulltable('/90daydata/vpgru/DavidLuecke/Maut/male/pctg/yahs/MautM_01-pctg.review.FINAL_diptera/run_diptera_odb10/full_table.tsv')
KT_MautM = pd.read_csv('/90daydata/vpgru/DavidLuecke/Maut/male/pctg/yahs/MautM_01-pctg.review.FINAL_diptera/MautM2025_karyotype.tsv', sep='\t')
KT_MautM['organism'] = 'M autumnalis male'
KT_MautM['color'] = '#002f42'

karyoplot(karyotype = KT_MautM, fulltable = FT_MautM, output_file = "karyoplot_MautM2025.png")

horizontal_synteny_plot(ft_1 = FT_MautF, ft_2 = FT_MautM, karyotype_1 = KT_MautF, karyotype_2 = KT_MautM, output_path = 'synteny_Maut2025.png')

# Scal
FT_ScalF = load_busco_fulltable('/90daydata/vpgru/DavidLuecke/Scal/female/pctg/yahs/ScalF_01-pctg.review.FINAL_diptera/run_diptera_odb10/full_table.tsv')
KT_ScalF = pd.read_csv('/90daydata/vpgru/DavidLuecke/Scal/female/pctg/yahs/ScalF_01-pctg.review.FINAL_diptera/ScalF2025_karyotype.tsv', sep='\t')
KT_ScalF['organism'] = 'S. calcitrans female'
KT_ScalF['color'] = '#d43a3f'

karyoplot(karyotype = KT_ScalF, fulltable = FT_ScalF, output_file = "karyoplot_ScalF2025.png")

FT_ScalM = load_busco_fulltable('/90daydata/vpgru/DavidLuecke/Scal/male/pctg/yahs/ScalM_01-pctg.review.FINAL_diptera/run_diptera_odb10/full_table.tsv')
KT_ScalM = pd.read_csv('/90daydata/vpgru/DavidLuecke/Scal/male/pctg/yahs/ScalM_01-pctg.review.FINAL_diptera/ScalM2025_karyotype.tsv', sep='\t')
KT_ScalM['organism'] = 'S. calcitrans male'
KT_ScalM['color'] = '#9b2226'

karyoplot(karyotype = KT_ScalM, fulltable = FT_ScalM, output_file = "karyoplot_ScalM2025.png")

horizontal_synteny_plot(ft_1 = FT_ScalF, ft_2 = FT_ScalM, karyotype_1 = KT_ScalF, karyotype_2 = KT_ScalM, output_path = 'synteny_Scal2025.png')

# Phylogeny via Li et al 2023: (Fcan(Haen((Maut,Scal)Hirr)))
horizontal_synteny_plot(ft_1 = FT_FcanF, ft_2 = FT_HaenF, karyotype_1 = KT_FcanF, karyotype_2 = KT_HaenF, output_path = 'synteny_FcanF-HaenF2025.png')
horizontal_synteny_plot(ft_1 = FT_FcanM, ft_2 = FT_HaenM, karyotype_1 = KT_FcanM, karyotype_2 = KT_HaenM, output_path = 'synteny_FcanM-HaenM2025.png')

horizontal_synteny_plot(ft_1 = FT_FcanF, ft_2 = FT_MautF, karyotype_1 = KT_FcanF, karyotype_2 = KT_MautF, output_path = 'synteny_FcanF-MautF2025.png')
horizontal_synteny_plot(ft_1 = FT_FcanM, ft_2 = FT_MautM, karyotype_1 = KT_FcanM, karyotype_2 = KT_MautM, output_path = 'synteny_FcanM-MautM2025.png')

horizontal_synteny_plot(ft_1 = FT_FcanF, ft_2 = FT_ScalF, karyotype_1 = KT_FcanF, karyotype_2 = KT_ScalF, output_path = 'synteny_FcanF-ScalF2025.png')
horizontal_synteny_plot(ft_1 = FT_FcanM, ft_2 = FT_ScalM, karyotype_1 = KT_FcanM, karyotype_2 = KT_ScalM, output_path = 'synteny_FcanM-ScalM2025.png')

horizontal_synteny_plot(ft_1 = FT_HaenF, ft_2 = FT_MautF, karyotype_1 = KT_HaenF, karyotype_2 = KT_MautF, output_path = 'synteny_HaenF-MautF2025.png')
horizontal_synteny_plot(ft_1 = FT_HaenM, ft_2 = FT_MautM, karyotype_1 = KT_HaenM, karyotype_2 = KT_MautM, output_path = 'synteny_HaenM-MautM2025.png')

horizontal_synteny_plot(ft_1 = FT_HaenF, ft_2 = FT_HirrF, karyotype_1 = KT_HaenF, karyotype_2 = KT_HirrF, output_path = 'synteny_HaenF-HirrF2025.png')
horizontal_synteny_plot(ft_1 = FT_HaenM, ft_2 = FT_HirrM, karyotype_1 = KT_HaenM, karyotype_2 = KT_HirrM, output_path = 'synteny_HaenM-HirrM2025.png')

horizontal_synteny_plot(ft_1 = FT_MautF, ft_2 = FT_ScalF, karyotype_1 = KT_MautF, karyotype_2 = KT_ScalF, output_path = 'synteny_MautF-ScalF2025.png')
horizontal_synteny_plot(ft_1 = FT_MautM, ft_2 = FT_ScalM, karyotype_1 = KT_MautM, karyotype_2 = KT_ScalM, output_path = 'synteny_MautM-ScalM2025.png')

horizontal_synteny_plot(ft_1 = FT_ScalF, ft_2 = FT_HirrF, karyotype_1 = KT_ScalF, karyotype_2 = KT_HirrF, output_path = 'synteny_ScalF-HirrF2025.png')
horizontal_synteny_plot(ft_1 = FT_ScalM, ft_2 = FT_HirrM, karyotype_1 = KT_ScalM, karyotype_2 = KT_HirrM, output_path = 'synteny_ScalM-HirrM2025.png')

# Get Muller elements
# Dmel
FT_Dmel = load_busco_fulltable('/90daydata/vpgru/DavidLuecke/Dmel/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic_diptera/run_diptera_odb10/full_table.tsv')
KT_Dmel = pd.read_csv('/90daydata/vpgru/DavidLuecke/Dmel/Dmel-karyotype.tsv', sep='\t')
KT_Dmel['organism'] = 'D. melanogaster'
KT_Dmel['color'] = '#808080'

# all Muscids vs Dmel
link_colors = {
    'NC_004354.4': '#800080',
    'NT_033779.5': '#3f4c99',
    'NT_033778.4': '#851d1d',
    'NT_037436.4': '#ff9900',
    'NT_033777.3': '#368244',
    'NC_004353.4': '#000000',
}

horizontal_synteny_plot(ft_1 = FT_Dmel, ft_2 = FT_FcanF, karyotype_1 = KT_Dmel, karyotype_2 = KT_FcanF, link_colors=link_colors, output_path = 'synteny_Dmel-FcanF2025.png')

horizontal_synteny_plot(ft_1 = FT_Dmel, ft_2 = FT_HaenF, karyotype_1 = KT_Dmel, karyotype_2 = KT_HaenF, link_colors=link_colors, output_path = 'synteny_Dmel-HaenF2025.png')

horizontal_synteny_plot(ft_1 = FT_Dmel, ft_2 = FT_HirrF, karyotype_1 = KT_Dmel, karyotype_2 = KT_HirrF, link_colors=link_colors, output_path = 'synteny_Dmel-HirrF2025.png')

horizontal_synteny_plot(ft_1 = FT_Dmel, ft_2 = FT_MautF, karyotype_1 = KT_Dmel, karyotype_2 = KT_MautF, link_colors=link_colors, output_path = 'synteny_Dmel-MautF2025.png')

horizontal_synteny_plot(ft_1 = FT_Dmel, ft_2 = FT_ScalF, karyotype_1 = KT_Dmel, karyotype_2 = KT_ScalF, link_colors=link_colors, output_path = 'synteny_Dmel-ScalF2025.png')
