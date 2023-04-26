import pandas as pd
import numpy as np
import sys, os

if not os.path.exists("OSAT") and os.path.exists(
    "../OSAT"
):  # hack to allow scripts to be placed in subdirectories next to AlOlivThermo:
    sys.path.insert(1, os.path.abspath(".."))

from OSAT import thermometers
from OSAT.tools import vec_mc


df_ol = pd.read_excel("input_ol.xlsx")
df_sp = pd.read_excel("input_sp.xlsx")
df_ol_err =  pd.read_excel("input_ol_err.xlsx")
df_sp_err = pd.read_excel("input_sp_err.xlsx")

mc_col_sp = df_sp_err.columns[1::]

mc_col_ol = df_ol_err.columns[1:-1]


mc = 1000

sp_mc_collect = pd.DataFrame()
for idx in range(len(df_sp)):
    df_iter = vec_mc(df_sp[mc_col_sp], df_sp_err[mc_col_sp], idx, mc)
    sp_mc_collect = sp_mc_collect.append(df_iter)

ol_mc_collect = pd.DataFrame()
for idx in range(len(df_ol)):
    df_iter = vec_mc(df_ol[mc_col_ol], df_ol_err[mc_col_ol], idx, mc)
    ol_mc_collect = ol_mc_collect.append(df_iter)

sp_mc_collect.columns = mc_col_sp
sp_mc_collect[sp_mc_collect < 0] = 0

ol_mc_collect.columns = mc_col_ol
ol_mc_collect[ol_mc_collect < 0] = 0

sp_mc_collect['sample_group'] = df_sp.loc[df_sp.index.repeat(mc)]['sample_no']  # avoid merge error in the script of alolivthermo
ol_mc_collect['sample_group'] = df_ol.loc[df_ol.index.repeat(mc)]['sample_no']

ol_mc_collect.reset_index(drop=True, inplace=True)
sp_mc_collect.reset_index(drop=True, inplace=True)

ol_mc_collect['sample_no'] = ol_mc_collect.index  # sampel_no match use
sp_mc_collect['sample_no'] = sp_mc_collect.index


reg = thermometers.models(ol_mc_collect, sp_mc_collect).compute()

print(
reg.groupby('sample_group')[['t_thermo', 't_kdcr', 't_kdal', 't_coogan','t_zThermoAlCr', 't_zThermoAl']].agg(['mean', 'std'])
)
