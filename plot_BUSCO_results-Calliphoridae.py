#!/bin/python3

# plot_BUSCO_results-Calliphoridae.py uses the buscoplotpy package 
# to plot Muscid busco results as karyotype and synteny pairs
# first need to activate the environment (in salloc):
# $ module load python_3
# $ source /project/vpgru/software/python_envs/buscoplotpy/bin/activate

import pandas as pd
from buscoplotpy.utils.load_busco_fulltable import load_busco_fulltable
from buscoplotpy.graphics.karyoplot import karyoplot
from buscoplotpy.graphics.synteny import horizontal_synteny_plot, vertical_synteny_plot

# Dmel
FT_Dmel = load_busco_fulltable('/90daydata/vpgru/DavidLuecke/Dmel/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic_diptera/run_diptera_odb10/full_table.tsv')
KT_Dmel = pd.read_csv('/90daydata/vpgru/DavidLuecke/Dmel/Dmel-karyotype.tsv', sep='\t')
KT_Dmel['organism'] = 'D. melanogaster'
KT_Dmel['color'] = '#808080'

karyoplot(karyotype = KT_Dmel, fulltable = FT_Dmel, output_file = "karyoplot_Dmel.png")


# Bpan
FT_Bpan = load_busco_fulltable('/90daydata/vpgru/DavidLuecke/Bpan/Bpan_diptera/run_diptera_odb10/full_table.tsv')
KT_Bpan = pd.read_csv('/90daydata/vpgru/DavidLuecke/Bpan/Bpan-karyotype.tsv', sep='\t')
KT_Bpan['organism'] = 'B. pandia'
KT_Bpan['color'] = '#800000'

karyoplot(karyotype = KT_Bpan, fulltable = FT_Bpan, output_file = "karyoplot_Bpan.png")


# Chom_old
FT_Chom_old = load_busco_fulltable('/90daydata/vpgru/DavidLuecke/Chom/Chom_old_diptera/run_diptera_odb10/full_table.tsv')
KT_Chom_old = pd.read_csv('/90daydata/vpgru/DavidLuecke/Chom/Chom_old-karyotype.tsv', sep='\t')
KT_Chom_old['organism'] = 'C. hominivorax - old assembly'
KT_Chom_old['color'] = '#808000'

karyoplot(karyotype = KT_Chom_old, fulltable = FT_Chom_old, output_file = "karyoplot_Chom_old.png")

# Chom_v1
FT_Chom_v1 = load_busco_fulltable('/90daydata/vpgru/DavidLuecke/Chom/Chom_v1_diptera/run_diptera_odb10/full_table.tsv')
KT_Chom_v1 = pd.read_csv('/90daydata/vpgru/DavidLuecke/Chom/Chom_v1-karyotype.tsv', sep='\t')
KT_Chom_v1['organism'] = 'C. hominivorax - v1'
KT_Chom_v1['color'] = '#808000'

karyoplot(karyotype = KT_Chom_v1, fulltable = FT_Chom_v1, output_file = "karyoplot_Chom_v1.png")


# Cmac
FT_Cmac = load_busco_fulltable('/90daydata/vpgru/DavidLuecke/Cmac/Cmac_v2_diptera/run_diptera_odb10/full_table.tsv')
KT_Cmac = pd.read_csv('/90daydata/vpgru/DavidLuecke/Cmac/Cmac-karyotype.tsv', sep='\t')
KT_Cmac['organism'] = 'C. macellaria'
KT_Cmac['color'] = '#008000'

karyoplot(karyotype = KT_Cmac, fulltable = FT_Cmac, output_file = "karyoplot_Cmac.png")


# Cvic
FT_Cvic = load_busco_fulltable('/90daydata/vpgru/DavidLuecke/Cvic/Cvic_diptera/run_diptera_odb10/full_table.tsv')
KT_Cvic = pd.read_csv('/90daydata/vpgru/DavidLuecke/Cvic/Cvic-karyotype.tsv', sep='\t')
KT_Cvic['organism'] = 'Ca. vicina'
KT_Cvic['color'] = '#800080'

karyoplot(karyotype = KT_Cvic, fulltable = FT_Cvic, output_file = "karyoplot_Cvic.png")


# Lcup
FT_Lcup = load_busco_fulltable('/90daydata/vpgru/DavidLuecke/Lcup/Lcup_diptera/run_diptera_odb10/full_table.tsv')
KT_Lcup = pd.read_csv('/90daydata/vpgru/DavidLuecke/Lcup/Lcup-karyotype.tsv', sep='\t')
KT_Lcup['organism'] = 'L. cuprina'
KT_Lcup['color'] = '#008080'

karyoplot(karyotype = KT_Lcup, fulltable = FT_Lcup, output_file = "karyoplot_Lcup.png")


# Pter
FT_Pter = load_busco_fulltable('/90daydata/vpgru/DavidLuecke/Pter/Pter_diptera/run_diptera_odb10/full_table.tsv')
KT_Pter = pd.read_csv('/90daydata/vpgru/DavidLuecke/Pter/Pter-karyotype.tsv', sep='\t')
KT_Pter['organism'] = 'P. terraenovae'
KT_Pter['color'] = '#000080'

karyoplot(karyotype = KT_Pter, fulltable = FT_Pter, output_file = "karyoplot_Pter.png")


# cross-species synteny, branch order: (Dmel((Lcup(Bpan,Cvic))(Pter(Cmac,Chom))))

horizontal_synteny_plot(ft_1 = FT_Dmel, ft_2 = FT_Lcup, karyotype_1 = KT_Dmel, karyotype_2 = KT_Lcup, output_path = 'synteny_Dmel-Lcup.png')
horizontal_synteny_plot(ft_1 = FT_Lcup, ft_2 = FT_Bpan, karyotype_1 = KT_Lcup, karyotype_2 = KT_Bpan, output_path = 'synteny_Lcup-Bpan.png')
horizontal_synteny_plot(ft_1 = FT_Bpan, ft_2 = FT_Cvic, karyotype_1 = KT_Bpan, karyotype_2 = KT_Cvic, output_path = 'synteny_Bpan-Cvic.png')
horizontal_synteny_plot(ft_1 = FT_Cvic, ft_2 = FT_Pter, karyotype_1 = KT_Cvic, karyotype_2 = KT_Pter, output_path = 'synteny_Cvic-Pter.png')
horizontal_synteny_plot(ft_1 = FT_Pter, ft_2 = FT_Cmac, karyotype_1 = KT_Pter, karyotype_2 = KT_Cmac, output_path = 'synteny_Pter-Cmac.png')
horizontal_synteny_plot(ft_1 = FT_Cmac, ft_2 = FT_Chom_old, karyotype_1 = KT_Cmac, karyotype_2 = KT_Chom_old, output_path = 'synteny_Cmac-Chom_old.png')
horizontal_synteny_plot(ft_1 = FT_Cmac, ft_2 = FT_Chom_v1, karyotype_1 = KT_Cmac, karyotype_2 = KT_Chom_v1, output_path = 'synteny_Cmac-Chom_v1.png')
