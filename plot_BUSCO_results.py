#!/bin/python3

# plot_BUSCO_results.py uses the buscoplotpy package 
# to plot Muscid busco results as karyotype and synteny pairs
# first need to activate the environment (in salloc):
# $ module load python_3
# $ source /project/vpgru/software/python_envs/buscoplotpy/bin/activate

import pandas as pd
from buscoplotpy.utils.load_busco_fulltable import load_busco_fulltable
from buscoplotpy.graphics.karyoplot import karyoplot
from buscoplotpy.graphics.synteny import horizontal_synteny_plot, vertical_synteny_plot

# Fcan
FT_FcanF = load_busco_fulltable('/90daydata/vpgru/DavidLuecke/Finished/Fcan/female/Fcan_F-v1_diptera/run_diptera_odb10/full_table.tsv')
KT_FcanF = pd.read_csv('/90daydata/vpgru/DavidLuecke/Finished/Fcan/female/FcanF_karyotype.tsv', sep='\t')
KT_FcanF['organism'] = 'F. canicularis female'
KT_FcanF['color'] = '#52b765'

karyoplot(karyotype = KT_FcanF, fulltable = FT_FcanF, output_file = "karyoplot_FcanF.png")

FT_FcanM = load_busco_fulltable('/90daydata/vpgru/DavidLuecke/Finished/Fcan/male/Fcan_M-v1_diptera/run_diptera_odb10/full_table.tsv')
KT_FcanM = pd.read_csv('/90daydata/vpgru/DavidLuecke/Finished/Fcan/male/FcanM_karyotype.tsv', sep='\t')
KT_FcanM['organism'] = 'F. canicularis male'
KT_FcanM['color'] = '#368244'

karyoplot(karyotype = KT_FcanM, fulltable = FT_FcanM, output_file = "karyoplot_FcanM.png")

horizontal_synteny_plot(ft_1 = FT_FcanF, ft_2 = FT_FcanM, karyotype_1 = KT_FcanF, karyotype_2 = KT_FcanM, output_path = 'synteny_Fcan.png')

# Haen
FT_HaenF = load_busco_fulltable('/90daydata/vpgru/DavidLuecke/Finished/Haen/female/Haen_F-v1_diptera/run_diptera_odb10/full_table.tsv')
KT_HaenF = pd.read_csv('/90daydata/vpgru/DavidLuecke/Finished/Haen/female/HaenF_karyotype.tsv', sep='\t')
KT_HaenF['organism'] = 'H. aenescens female'
KT_HaenF['color'] = '#fd9020'

karyoplot(karyotype = KT_HaenF, fulltable = FT_HaenF, output_file = "karyoplot_HaenF.png")

FT_HaenM = load_busco_fulltable('/90daydata/vpgru/DavidLuecke/Finished/Haen/male/Haen_M-v1_diptera/run_diptera_odb10/full_table.tsv')
KT_HaenM = pd.read_csv('/90daydata/vpgru/DavidLuecke/Finished/Haen/male/HaenM_karyotype.tsv', sep='\t')
KT_HaenM['organism'] = 'H. aenescens male'
KT_HaenM['color'] = '#ca6702'

karyoplot(karyotype = KT_HaenM, fulltable = FT_HaenM, output_file = "karyoplot_HaenM.png")

horizontal_synteny_plot(ft_1 = FT_HaenF, ft_2 = FT_HaenM, karyotype_1 = KT_HaenF, karyotype_2 = KT_HaenM, output_path = 'synteny_Haen.png')

# Hirr
FT_HirrF = load_busco_fulltable('/90daydata/vpgru/DavidLuecke/Finished/Hirr/female/Hirr_F-v2.masked_diptera/run_diptera_odb10/full_table.tsv')
KT_HirrF = pd.read_csv('/90daydata/vpgru/DavidLuecke/Finished/Hirr/female/HirrF_karyotype.tsv', sep='\t')
KT_HirrF['organism'] = 'H. irritans female'
KT_HirrF['color'] = '#dccd45'

karyoplot(karyotype = KT_HirrF, fulltable = FT_HirrF, output_file = "karyoplot_HirrF.png")

# FT_MautM = load_busco_fulltable('/90daydata/vpgru/DavidLuecke/Finished/Maut/male/Maut_M-v3.masked_diptera/run_diptera_odb10/full_table.tsv')
# KT_MautM = pd.read_csv('/90daydata/vpgru/DavidLuecke/Finished/Maut/male/MautM_karyotype.tsv', sep='\t')
# KT_MautM['organism'] = 'M autumnalis male'
# KT_HirrM['color'] = '#afa021'

# karyoplot(karyotype = KT_MautM, fulltable = FT_MautM, output_file = "karyoplot_MautM.png")

# horizontal_synteny_plot(ft_1 = FT_MautF, ft_2 = FT_MautM, karyotype_1 = KT_MautF, karyotype_2 = KT_MautM, output_path = 'synteny_Maut.png')

# Maut
FT_MautF = load_busco_fulltable('/90daydata/vpgru/DavidLuecke/Finished/Maut/female/Maut_F-v2.masked_diptera/run_diptera_odb10/full_table.tsv')
KT_MautF = pd.read_csv('/90daydata/vpgru/DavidLuecke/Finished/Maut/female/MautF_karyotype.tsv', sep='\t')
KT_MautF['organism'] = 'M autumnalis female'
KT_MautF['color'] = '#006994'

karyoplot(karyotype = KT_MautF, fulltable = FT_MautF, output_file = "karyoplot_MautF.png")

FT_MautM = load_busco_fulltable('/90daydata/vpgru/DavidLuecke/Finished/Maut/male/Maut_M-v3.masked_diptera/run_diptera_odb10/full_table.tsv')
KT_MautM = pd.read_csv('/90daydata/vpgru/DavidLuecke/Finished/Maut/male/MautM_karyotype.tsv', sep='\t')
KT_MautM['organism'] = 'M autumnalis male'
KT_MautM['color'] = '#002f42'

karyoplot(karyotype = KT_MautM, fulltable = FT_MautM, output_file = "karyoplot_MautM.png")

horizontal_synteny_plot(ft_1 = FT_MautF, ft_2 = FT_MautM, karyotype_1 = KT_MautF, karyotype_2 = KT_MautM, output_path = 'synteny_Maut.png')

# Scal
FT_ScalF = load_busco_fulltable('/90daydata/vpgru/DavidLuecke/Finished/Scal/female/Scal_F-v1_diptera/run_diptera_odb10/full_table.tsv')
KT_ScalF = pd.read_csv('/90daydata/vpgru/DavidLuecke/Finished/Scal/female/ScalF_karyotype.tsv', sep='\t')
KT_ScalF['organism'] = 'S. calcitrans female'
KT_ScalF['color'] = '#d43a3f'

karyoplot(karyotype = KT_ScalF, fulltable = FT_ScalF, output_file = "karyoplot_ScalF.png")

FT_ScalM = load_busco_fulltable('/90daydata/vpgru/DavidLuecke/Finished/Scal/male/Scal_M-v1_diptera/run_diptera_odb10/full_table.tsv')
KT_ScalM = pd.read_csv('/90daydata/vpgru/DavidLuecke/Finished/Scal/male/ScalM_karyotype.tsv', sep='\t')
KT_ScalM['organism'] = 'S. calcitrans male'
KT_ScalM['color'] = '#9b2226'

karyoplot(karyotype = KT_ScalM, fulltable = FT_ScalM, output_file = "karyoplot_ScalM.png")

horizontal_synteny_plot(ft_1 = FT_ScalF, ft_2 = FT_ScalM, karyotype_1 = KT_ScalF, karyotype_2 = KT_ScalM, output_path = 'synteny_Scal.png')

# Phylogeny via Li et al 2023: (Fcan(Haen((Maut,Scal)Hirr)))
horizontal_synteny_plot(ft_1 = FT_FcanF, ft_2 = FT_HaenF, karyotype_1 = KT_FcanF, karyotype_2 = KT_HaenF, output_path = 'synteny_FcanF-HaenF.png')
horizontal_synteny_plot(ft_1 = FT_FcanM, ft_2 = FT_HaenM, karyotype_1 = KT_FcanM, karyotype_2 = KT_HaenM, output_path = 'synteny_FcanM-HaenM.png')

horizontal_synteny_plot(ft_1 = FT_HaenF, ft_2 = FT_MautF, karyotype_1 = KT_HaenF, karyotype_2 = KT_MautF, output_path = 'synteny_HaenF-MautF.png')
horizontal_synteny_plot(ft_1 = FT_HaenM, ft_2 = FT_MautM, karyotype_1 = KT_HaenM, karyotype_2 = KT_MautM, output_path = 'synteny_HaenM-MautM.png')

horizontal_synteny_plot(ft_1 = FT_MautF, ft_2 = FT_ScalF, karyotype_1 = KT_MautF, karyotype_2 = KT_ScalF, output_path = 'synteny_MautF-ScalF.png')
horizontal_synteny_plot(ft_1 = FT_MautM, ft_2 = FT_ScalM, karyotype_1 = KT_MautM, karyotype_2 = KT_ScalM, output_path = 'synteny_MautM-ScalM.png')

horizontal_synteny_plot(ft_1 = FT_ScalF, ft_2 = FT_HirrF, karyotype_1 = KT_ScalF, karyotype_2 = KT_HirrF, output_path = 'synteny_ScalF-HirrF.png')


